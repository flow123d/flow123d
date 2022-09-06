/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    transport.cc
 * @ingroup transport
 * @brief   Transport
 */

#include <memory>

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/index_types.hh"

#include "mesh/mesh.h"
#include "mesh/partitioning.hh"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "mesh/neighbours.h"
#include "transport/transport.h"
#include "transport/assembly_convection.hh"

#include "la/distribution.hh"

#include "la/sparse_graph.hh"
#include <iostream>
#include <iomanip>
#include <string>

#include "io/output_time.hh"
#include "tools/time_governor.hh"
#include "tools/mixed.hh"
#include "coupling/balance.hh"
#include "input/accessors.hh"
#include "input/input_type.hh"

#include "fields/field_algo_base.hh"
#include "fields/field_values.hh"
#include "fields/field_fe.hh"
#include "fields/fe_value_handler.hh"
#include "fields/generic_field.hh"

#include "reaction/isotherm.hh" // SorptionType enum

#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "quadrature/quadrature_lib.hh"


FLOW123D_FORCE_LINK_IN_CHILD(convectionTransport)


namespace IT = Input::Type;

/********************************************************************************
 * Static methods and data members
 */
const string _equation_name = "Solute_Advection_FV";

const int ConvectionTransport::registrar =
		Input::register_class< ConvectionTransport, Mesh &, const Input::Record >(_equation_name) +
		ConvectionTransport::get_input_type().size();

const IT::Record &ConvectionTransport::get_input_type()
{
	return IT::Record(_equation_name, "Finite volume method, explicit in time, for advection only solute transport.")
			.derive_from(ConcentrationTransportBase::get_input_type())
			.copy_keys(EquationBase::user_fields_template(_equation_name))
			.declare_key("input_fields", IT::Array(
			        EqFields().make_field_descriptor_type(_equation_name)),
			        IT::Default::obligatory(),
			        "")
            .declare_key("output",
                    EqFields().output_fields.make_output_type(_equation_name, ""),
                    IT::Default("{ \"fields\": [ \"conc\" ] }"),
                    "Specification of output fields and output times.")
			.close();
}


ConvectionTransport::EqFields::EqFields() : TransportEqFields()
{
    *this += bc_conc.name("bc_conc")
            .description("Boundary condition for concentration of substances.")
            .input_default("0.0")
            .units( UnitSI().kg().m(-3) );

    *this += init_conc.name("init_conc")
            .description("Initial values for concentration of substances.")
            .input_default("0.0")
            .units( UnitSI().kg().m(-3) );

    output_fields += *this;
    output_fields += conc_mobile.name("conc")
            .units( UnitSI().kg().m(-3) )
            .flags(FieldFlag::equation_result)
            .description("Concentration solution.");
	output_fields += region_id.name("region_id")
	        .units( UnitSI::dimensionless())
	        .flags(FieldFlag::equation_external_output)
            .description("Region ids.");
    output_fields += subdomain.name("subdomain")
            .units( UnitSI::dimensionless() )
            .flags(FieldFlag::equation_external_output)
            .description("Subdomain ids of the domain decomposition.");

    this->add_coords_field();
    this->set_default_fieldset();
}


/********************************************************************************
 * ConvectionTransport
 */
ConvectionTransport::ConvectionTransport(Mesh &init_mesh, const Input::Record in_rec)
: ConcentrationTransportBase(init_mesh, in_rec),
  input_rec(in_rec)
{
	START_TIMER("ConvectionTransport");
    eq_data_ = make_shared<EqData>();
    eq_fields_ = make_shared<EqFields>();
	this->eq_fieldset_ = eq_fields_.get();

	eq_data_->transport_matrix_time = -1.0; // or -infty
    eq_data_->transport_bc_time = -1.0;
    eq_data_->is_convection_matrix_scaled = false;
    is_src_term_scaled = false;
    is_bc_term_scaled = false;

    //initialization of DOF handler
    MixedPtr<FE_P_disc> fe(0);
    eq_data_->dh_ = make_shared<DOFHandlerMultiDim>(init_mesh);
    shared_ptr<DiscreteSpace> ds = make_shared<EqualOrderDiscreteSpace>( &init_mesh, fe);
    eq_data_->dh_->distribute_dofs(ds);
    vcumulative_corr = nullptr;

}

void ConvectionTransport::initialize()
{
    target_mark_type = time_->equation_fixed_mark_type();

    cfl_max_step = time_->end_time();

    eq_fields_->set_components(eq_data_->substances_.names());
    eq_fields_->set_input_list( input_rec.val<Input::Array>("input_fields"), *time_ );
    eq_fields_->set_mesh(*mesh_);

    alloc_transport_vectors();
    alloc_transport_structs_mpi();

    Input::Array user_fields_arr;
    if (input_rec.opt_val("user_fields", user_fields_arr)) {
       	this->init_user_fields(user_fields_arr, time().step().end(), eq_fields_->output_fields);
    }

	// register output vectors
    eq_fields_->output_fields.set_components(eq_data_->substances_.names());
    eq_fields_->output_fields.set_mesh(*mesh_);
    eq_fields_->output_fields.output_type(OutputTime::ELEM_DATA);
    eq_fields_->conc_mobile.setup_components();
    eq_fields_->region_id = GenericField<3>::region_id(*mesh_);
    eq_fields_->subdomain = GenericField<3>::subdomain(*mesh_);
	for (unsigned int sbi=0; sbi<n_substances(); sbi++)
	{
		// create shared pointer to a FieldFE and push this Field to output_field on all regions
	    eq_fields_->conc_mobile_fe[sbi] = create_field_fe< 3, FieldValue<3>::Scalar >(eq_data_->dh_);
	    eq_fields_->conc_mobile[sbi].set(eq_fields_->conc_mobile_fe[sbi], 0);
	}
	//output_stream_->add_admissible_field_names(input_rec.val<Input::Array>("output_fields"));
    //output_stream_->mark_output_times(*time_);

	eq_fields_->output_fields.initialize(output_stream_, mesh_, input_rec.val<Input::Record>("output"), time() );
	//cout << "Transport." << endl;
	//cout << time().marks();

    balance_->allocate(mesh_->get_el_ds()->lsize(), 1);
    eq_data_->balance_ = this->balance();
    eq_data_->set_time_governor(this->time_);
    eq_data_->max_edg_sides = max(this->mesh_->max_edge_sides(1), max(this->mesh_->max_edge_sides(2), this->mesh_->max_edge_sides(3)));

    mass_assembly_ = new GenericAssembly< MassAssemblyConvection >(eq_fields_.get(), eq_data_.get());
    init_cond_assembly_ = new GenericAssembly< InitCondAssemblyConvection >(eq_fields_.get(), eq_data_.get());
    conc_sources_bdr_assembly_ = new GenericAssembly< ConcSourcesBdrAssemblyConvection >(eq_fields_.get(), eq_data_.get());
    matrix_mpi_assembly_ = new GenericAssembly< MatrixMpiAssemblyConvection >(eq_fields_.get(), eq_data_.get());
    matrix_mpi_assembly_->set_min_edge_sides(1);
}


Vec ConvectionTransport::get_component_vec(unsigned int sbi)
{
    return eq_fields_->conc_mobile_fe[sbi]->vec().petsc_vec();
}



ConvectionTransport::~ConvectionTransport()
{
    unsigned int sbi;

    if (vcumulative_corr) {
        //Destroy mpi vectors at first
        chkerr(MatDestroy(&eq_data_->tm));
        chkerr(VecDestroy(&eq_data_->mass_diag));
        chkerr(VecDestroy(&vpmass_diag));

        for (sbi = 0; sbi < n_substances(); sbi++) {
            // mpi vectors
        	chkerr(VecDestroy(&vpconc[sbi]));
        	chkerr(VecDestroy(&eq_data_->bcvcorr[sbi]));
        	chkerr(VecDestroy(&vcumulative_corr[sbi]));
        }

        // arrays of mpi vectors
        delete vpconc;
        delete eq_data_->bcvcorr;
        delete vcumulative_corr;
        
        // assembly objects
        delete mass_assembly_;
        delete init_cond_assembly_;
        delete conc_sources_bdr_assembly_;
        delete matrix_mpi_assembly_;
    }
}





//=============================================================================
//	ALLOCATE OF TRANSPORT VARIABLES (ELEMENT & NODES)
//=============================================================================
void ConvectionTransport::alloc_transport_vectors() {

    eq_fields_->conc_mobile_fe.resize(this->n_substances());
}

//=============================================================================
//	ALLOCATION OF TRANSPORT VECTORS (MPI)
//=============================================================================
void ConvectionTransport::alloc_transport_structs_mpi() {

    int sbi, n_subst, lsize;
    n_subst = n_substances();
    lsize = mesh_->get_el_ds()->lsize();

    MPI_Barrier(PETSC_COMM_WORLD);

    vpconc = new Vec[n_subst];
    eq_data_->bcvcorr = new Vec[n_subst];
    vcumulative_corr = new Vec[n_subst];
    eq_data_->tm_diag.reserve(eq_data_->n_substances());
    eq_data_->corr_vec.reserve(eq_data_->n_substances());
    

    for (sbi = 0; sbi < n_subst; sbi++) {
        VecCreateMPI(PETSC_COMM_WORLD, lsize, mesh_->n_elements(), &eq_data_->bcvcorr[sbi]);
        VecZeroEntries(eq_data_->bcvcorr[sbi]);

        VecCreateMPI(PETSC_COMM_WORLD, lsize, mesh_->n_elements(), &vpconc[sbi]);
        VecZeroEntries(vpconc[sbi]);

        // SOURCES
        VecCreateMPI(PETSC_COMM_WORLD, lsize, mesh_->n_elements(), &vcumulative_corr[sbi]);
        
        eq_data_->corr_vec.emplace_back(lsize, PETSC_COMM_WORLD);
        
        eq_data_->tm_diag.emplace_back(lsize, PETSC_COMM_WORLD);

        VecZeroEntries(vcumulative_corr[sbi]);
    }


    MatCreateAIJ(PETSC_COMM_WORLD, lsize, lsize, mesh_->n_elements(),
            mesh_->n_elements(), 16, PETSC_NULL, 4, PETSC_NULL, &eq_data_->tm);
    
    VecCreateMPI(PETSC_COMM_WORLD, lsize, mesh_->n_elements(), &eq_data_->mass_diag);
    VecCreateMPI(PETSC_COMM_WORLD, lsize, mesh_->n_elements(), &vpmass_diag);

    eq_data_->alloc_transport_structs_mpi(lsize);
}


void ConvectionTransport::zero_time_step()
{
	ASSERT_EQ(time_->tlevel(), 0);

	eq_fields_->mark_input_times(*time_);
	eq_fields_->set_time(time_->step(), LimitSide::right);
	std::stringstream ss; // print warning message with table of uninitialized fields
	if ( FieldCommon::print_message_table(ss, "convection transport") ) {
		WarningOut() << ss.str();
	}

    //set_initial_condition();
    init_cond_assembly_->assemble(eq_data_->dh_);
    //create_mass_matrix();
    mass_assembly_->assemble(eq_data_->dh_);
    
    START_TIMER("Convection balance zero time step");

    START_TIMER("convection_matrix_assembly");
    //create_transport_matrix_mpi();
    matrix_mpi_assembly_->assemble(eq_data_->dh_);
    END_TIMER("convection_matrix_assembly");
	START_TIMER("sources_reinit_set_bc");
    //conc_sources_bdr_conditions();
    conc_sources_bdr_assembly_->assemble(eq_data_->dh_);
	END_TIMER("sources_reinit_set_bc");

    // write initial condition
	output_data();
}


bool ConvectionTransport::evaluate_time_constraint(double& time_constraint)
{
	ASSERT_PTR(eq_data_->dh_).error( "Null DOF handler object.\n" );
    // read changed status before setting time
    bool changed_flux = eq_fields_->flow_flux.changed();
    eq_fields_->set_time(time_->step(), LimitSide::right); // set to the last computed time

    START_TIMER("data reinit");
    
    bool cfl_changed = false;
    
    // if FLOW or DATA changed ---------------------> recompute transport matrix
    if (changed_flux)
    {
        START_TIMER("convection_matrix_assembly");
        //create_transport_matrix_mpi();
        matrix_mpi_assembly_->assemble(eq_data_->dh_);
        END_TIMER("convection_matrix_assembly");
        eq_data_->is_convection_matrix_scaled=false;
        cfl_changed = true;
        DebugOut() << "CFL changed - flow.\n";
    }
    
    if (eq_data_->is_mass_diag_changed)
    {
        //create_mass_matrix();
        mass_assembly_->assemble(eq_data_->dh_);
        cfl_changed = true;
        DebugOut() << "CFL changed - mass matrix.\n";
    }
    
    // if DATA changed ---------------------> recompute concentration sources (rhs and matrix diagonal)
    if( eq_fields_->sources_density.changed() || eq_fields_->sources_conc.changed() || eq_fields_->sources_sigma.changed()
       || eq_fields_->cross_section.changed() || eq_fields_->flow_flux.changed() || eq_fields_->porosity.changed()
       || eq_fields_->water_content.changed() || eq_fields_->bc_conc.changed() )
    {
    	START_TIMER("sources_reinit_set_bc");
        //conc_sources_bdr_conditions();
        conc_sources_bdr_assembly_->assemble(eq_data_->dh_);
    	END_TIMER("sources_reinit_set_bc");
        if( eq_data_->sources_changed_ ) {
            is_src_term_scaled = false;
            cfl_changed = true;
        }
        if (eq_fields_->flow_flux.changed() || eq_fields_->porosity.changed()
        	       || eq_fields_->water_content.changed() || eq_fields_->bc_conc.changed() ) is_bc_term_scaled = false;
        DebugOut() << "CFL changed - source.\n";
    }
    
    // now resolve the CFL condition
    if(cfl_changed)
    {
        // find maximum of sum of contribution from flow and sources: MAX(vcfl_flow_ + vcfl_source_)
        Vec cfl;
        VecCreateMPI(PETSC_COMM_WORLD, mesh_->get_el_ds()->lsize(), PETSC_DETERMINE, &cfl);
        VecWAXPY(cfl, 1.0, eq_data_->cfl_flow_.petsc_vec(), eq_data_->cfl_source_.petsc_vec());
        VecMaxPointwiseDivide(cfl, eq_data_->mass_diag, &cfl_max_step);
        // get a reciprocal value as a time constraint
        cfl_max_step = 1 / cfl_max_step;
        DebugOut().fmt("CFL constraint (transport): {}\n", cfl_max_step);
    }
    
    END_TIMER("data reinit");
    
    // return time constraint
    time_constraint = cfl_max_step;
    return cfl_changed;
}

void ConvectionTransport::update_solution() {

    START_TIMER("convection-one step");
    
    // proceed to next time - which we are about to compute
    // explicit scheme looks one step back and uses data from previous time
    // (data time set previously in assess_time_constraint())
    time_->next_time();
    
    double dt_new = time_->dt(),                    // current time step we are about to compute
           dt_scaled = dt_new / time_->last_dt();   // scaling ratio to previous time step
    
    START_TIMER("time step rescaling");
    
    // if FLOW or DATA or BC or DT changed ---------------------> rescale boundary condition
    if( ! is_bc_term_scaled || time_->is_changed_dt() )
    {
    	DebugOut() << "BC - rescale dt.\n";
        //choose between fresh scaling with new dt or rescaling to a new dt
        double dt = (!is_bc_term_scaled) ? dt_new : dt_scaled;
        for (unsigned int  sbi=0; sbi<n_substances(); sbi++)
            VecScale(eq_data_->bcvcorr[sbi], dt);
        is_bc_term_scaled = true;
    }
    

    // if DATA or TIME STEP changed -----------------------> rescale source term
    if( !is_src_term_scaled || time_->is_changed_dt()) {
    	DebugOut() << "SRC - rescale dt.\n";
        //choose between fresh scaling with new dt or rescaling to a new dt
        double dt = (!is_src_term_scaled) ? dt_new : dt_scaled;
        for (unsigned int sbi=0; sbi<n_substances(); sbi++)
        {
            VecScale(eq_data_->corr_vec[sbi].petsc_vec(), dt);
            VecScale(eq_data_->tm_diag[sbi].petsc_vec(), dt);
        }
        is_src_term_scaled = true;
    }
    
    // if DATA or TIME STEP changed -----------------------> rescale transport matrix
    if ( !eq_data_->is_convection_matrix_scaled || time_->is_changed_dt()) {
    	DebugOut() << "TM - rescale dt.\n";
        //choose between fresh scaling with new dt or rescaling to a new dt
        double dt = (!eq_data_->is_convection_matrix_scaled) ? dt_new : dt_scaled;
        
        MatScale(eq_data_->tm, dt);
        eq_data_->is_convection_matrix_scaled = true;
    }
    
    END_TIMER("time step rescaling");
    
    
    eq_fields_->set_time(time_->step(), LimitSide::right); // set to the last computed time
    if (eq_fields_->cross_section.changed() || eq_fields_->water_content.changed() || eq_fields_->porosity.changed())
    {
        VecCopy(eq_data_->mass_diag, vpmass_diag);
        //create_mass_matrix();
        mass_assembly_->assemble(eq_data_->dh_);
    } else eq_data_->is_mass_diag_changed = false;
    

    // Compute new concentrations for every substance.
    
    for (unsigned int sbi = 0; sbi < n_substances(); sbi++) {
      // one step in MOBILE phase
      START_TIMER("mat mult");
      Vec vconc = eq_fields_->conc_mobile_fe[sbi]->vec().petsc_vec();
      
      // tm_diag is a diagonal part of transport matrix, which depends on substance data (sources_sigma)
      // Wwe need keep transport matrix independent of substance, therefore we keep this diagonal part
      // separately in a vector tm_diag.
      // Computation: first, we compute this diagonal addition D*pconc and save it temporaly into RHS
        
      // RHS = D*pconc, where D is diagonal matrix represented by a vector
      VecPointwiseMult(vcumulative_corr[sbi], eq_data_->tm_diag[sbi].petsc_vec(), vconc); //w = x.*y
      
      // Then we add boundary terms ans other source terms into RHS.
      // RHS = 1.0 * bcvcorr + 1.0 * corr_vec + 1.0 * rhs
      VecAXPBYPCZ(vcumulative_corr[sbi], 1.0, 1.0, 1.0, eq_data_->bcvcorr[sbi], eq_data_->corr_vec[sbi].petsc_vec());   //z = ax + by + cz
      
      // Then we set the new previous concentration.
      VecCopy(vconc, vpconc[sbi]); // pconc = conc
      // And finally proceed with transport matrix multiplication.
      if (eq_data_->is_mass_diag_changed) {
        VecPointwiseMult(vconc, vconc, vpmass_diag);         // vconc*=vpmass_diag
        MatMultAdd(eq_data_->tm, vpconc[sbi], vconc, vconc);           // vconc+=tm*vpconc
        VecAXPY(vconc, 1, vcumulative_corr[sbi]);            // vconc+=vcumulative_corr
        VecPointwiseDivide(vconc, vconc, eq_data_->mass_diag); // vconc/=mass_diag
      } else {
        MatMultAdd(eq_data_->tm, vpconc[sbi], vcumulative_corr[sbi], vconc);  // vconc =tm*vpconc+vcumulative_corr
        VecPointwiseDivide(vconc, vconc, eq_data_->mass_diag);        // vconc/=mass_diag
        VecAXPY(vconc, 1, vpconc[sbi]);                             // vconc+=vpconc
      }

      END_TIMER("mat mult");
    }
    
    for (unsigned int sbi=0; sbi<n_substances(); ++sbi)
      balance_->calculate_cumulative(sbi, vpconc[sbi]);
    
    END_TIMER("convection-one step");
}


void ConvectionTransport::set_target_time(double target_time)
{

    //sets target_mark_type (it is fixed) to be met in next_time()
    time_->marks().add(TimeMark(target_time, target_mark_type));

    // This is done every time TOS calls update_solution.
    // If CFL condition is changed, time fixation will change later from TOS.
    
    // Set the same constraint as was set last time.
    
    // TODO: fix this hack, remove this method completely, leaving just the first line at the calling point
    // in TransportOperatorSplitting::update_solution()
    // doing this directly leads to choose of large time step violationg CFL condition
    if (cfl_max_step > time_->dt()*1e-10)
    time_->set_upper_constraint(cfl_max_step, "CFL condition used from previous step.");

    // fixing convection time governor till next target_mark_type (got from TOS or other)
    // may have marks for data changes
    time_->fix_dt_until_mark();
}


ConvectionTransport::FieldFEScalarVec& ConvectionTransport::get_p0_interpolation() {
    return eq_fields_->conc_mobile_fe;
}


void ConvectionTransport::output_data() {

    eq_fields_->output_fields.set_time(time().step(), LimitSide::right);
    eq_fields_->output_fields.output(time().step());
    
    START_TIMER("TOS-balance");
    for (unsigned int sbi=0; sbi<n_substances(); ++sbi)
        balance_->calculate_instant(sbi, eq_fields_->conc_mobile_fe[sbi]->vec().petsc_vec());
    balance_->output();
    END_TIMER("TOS-balance");
}

void ConvectionTransport::set_balance_object(std::shared_ptr<Balance> balance)
{
	balance_ = balance;
	eq_data_->subst_idx = balance_->add_quantities(eq_data_->substances_.names());
    eq_data_->balance_ = this->balance();
}
