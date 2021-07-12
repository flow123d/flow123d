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
			.declare_key("input_fields", IT::Array(
			        EqData().make_field_descriptor_type(_equation_name)),
			        IT::Default::obligatory(),
			        "")
            .declare_key("output",
                    EqData().output_fields.make_output_type(_equation_name, ""),
                    IT::Default("{ \"fields\": [ \"conc\" ] }"),
                    "Specification of output fields and output times.")
			.close();
}


ConvectionTransport::EqData::EqData() : TransportEqFields()
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
}


/********************************************************************************
 * Helper class FETransportObjects
 */
FETransportObjects::FETransportObjects()
: q_(QGauss::make_array(0))
{
    fe0_ = new FE_P_disc<0>(0);
    fe1_ = new FE_P_disc<1>(0);
    fe2_ = new FE_P_disc<2>(0);
    fe3_ = new FE_P_disc<3>(0);


    fe_values_[0].initialize(q(0), *fe1_,
            update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
    fe_values_[1].initialize(q(1), *fe2_,
            update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
    fe_values_[2].initialize(q(2), *fe3_,
            update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
}


FETransportObjects::~FETransportObjects()
{
    delete fe0_;
    delete fe1_;
    delete fe2_;
    delete fe3_;
}

template<> FiniteElement<0> *FETransportObjects::fe<0>() { return fe0_; }
template<> FiniteElement<1> *FETransportObjects::fe<1>() { return fe1_; }
template<> FiniteElement<2> *FETransportObjects::fe<2>() { return fe2_; }
template<> FiniteElement<3> *FETransportObjects::fe<3>() { return fe3_; }

Quadrature &FETransportObjects::q(unsigned int dim) { return q_[dim]; }



/********************************************************************************
 * ConvectionTransport
 */
ConvectionTransport::ConvectionTransport(Mesh &init_mesh, const Input::Record in_rec)
: ConcentrationTransportBase(init_mesh, in_rec),
  is_mass_diag_changed(false),
  tm_diag(nullptr),
  input_rec(in_rec)
{
	START_TIMER("ConvectionTransport");
	this->eq_fieldset_ = &data_;

    transport_matrix_time = -1.0; // or -infty
    transport_bc_time = -1.0;
    is_convection_matrix_scaled = false;
    is_src_term_scaled = false;
    is_bc_term_scaled = false;

    //initialization of DOF handler
    MixedPtr<FE_P_disc> fe(0);
    dh_ = make_shared<DOFHandlerMultiDim>(init_mesh);
    shared_ptr<DiscreteSpace> ds = make_shared<EqualOrderDiscreteSpace>( &init_mesh, fe);
    dh_->distribute_dofs(ds);

}

void ConvectionTransport::initialize()
{
    target_mark_type = time_->equation_fixed_mark_type();

    cfl_max_step = time_->end_time();

    data_.set_components(substances_.names());
    data_.set_input_list( input_rec.val<Input::Array>("input_fields"), *time_ );
    data_.set_mesh(*mesh_);

    make_transport_partitioning();
    alloc_transport_vectors();
    alloc_transport_structs_mpi();

	// register output vectors
	data_.output_fields.set_components(substances_.names());
	data_.output_fields.set_mesh(*mesh_);
    auto output_types = OutputTime::empty_discrete_flags();
    output_types[OutputTime::ELEM_DATA] = true;
	data_.output_fields.output_type(output_types);
	data_.conc_mobile.setup_components();
	data_.region_id = GenericField<3>::region_id(*mesh_);
    data_.subdomain = GenericField<3>::subdomain(*mesh_);
	for (unsigned int sbi=0; sbi<n_substances(); sbi++)
	{
		// create shared pointer to a FieldFE and push this Field to output_field on all regions
		data_.conc_mobile_fe[sbi] = create_field_fe< 3, FieldValue<3>::Scalar >(dh_);
		data_.conc_mobile[sbi].set(data_.conc_mobile_fe[sbi], 0);
	}
	//output_stream_->add_admissible_field_names(input_rec.val<Input::Array>("output_fields"));
    //output_stream_->mark_output_times(*time_);

	data_.output_fields.initialize(output_stream_, mesh_, input_rec.val<Input::Record>("output"), time() );
	//cout << "Transport." << endl;
	//cout << time().marks();

    balance_->allocate(el_ds->lsize(), 1);

}


Vec ConvectionTransport::get_component_vec(unsigned int sbi)
{
    return data_.conc_mobile_fe[sbi]->vec().petsc_vec(); 
}


//=============================================================================
// MAKE TRANSPORT
//=============================================================================
void ConvectionTransport::make_transport_partitioning() {

//    int * id_4_old = new int[mesh_->n_elements()];
//    int i = 0;
//    for (auto ele : mesh_->elements_range()) id_4_old[i++] = ele.index();
//    mesh_->get_part()->id_maps(mesh_->n_elements(), id_4_old, el_ds, el_4_loc, row_4_el);
//    delete[] id_4_old;
	el_ds = mesh_->get_el_ds();
	el_4_loc = mesh_->get_el_4_loc();
	row_4_el = mesh_->get_row_4_el();

    // TODO: make output of partitioning is usefull but makes outputs different
    // on different number of processors, which breaks tests.
    //
    // Possible solution:
    // - have flag in ini file to turn this output ON
    // - possibility to have different ref_output for different num of proc.
    // - or do not test such kind of output
    //
    //for (auto ele : mesh_->elements_range()) {
    //    ele->pid()=el_ds->get_proc(row_4_el[ele.index()]);
    //}

}



ConvectionTransport::~ConvectionTransport()
{
    unsigned int sbi;

    if (tm_diag) {
        //Destroy mpi vectors at first
        chkerr(MatDestroy(&tm));
        chkerr(VecDestroy(&mass_diag));
        chkerr(VecDestroy(&vpmass_diag));
        chkerr(VecDestroy(&vcfl_flow_));
        chkerr(VecDestroy(&vcfl_source_));
        delete cfl_flow_;
        delete cfl_source_;

        for (sbi = 0; sbi < n_substances(); sbi++) {
            // mpi vectors
        	chkerr(VecDestroy(&vpconc[sbi]));
        	chkerr(VecDestroy(&bcvcorr[sbi]));
        	chkerr(VecDestroy(&vcumulative_corr[sbi]));
        	chkerr(VecDestroy(&v_tm_diag[sbi]));

            // arrays of arrays
            delete cumulative_corr[sbi];
            delete tm_diag[sbi];
        }

        // arrays of mpi vectors
        delete vpconc;
        delete bcvcorr;
        delete vcumulative_corr;
        delete v_tm_diag;
        
        // arrays of arrays
        delete cumulative_corr;
        delete tm_diag;
    }
}





void ConvectionTransport::set_initial_condition()
{
    DebugOut() << "ConvectionTransport::set_initial_condition()\n";
    // get vecs
    std::vector<VectorMPI> vecs;
    for (unsigned int sbi=0; sbi<n_substances(); sbi++)
        vecs.push_back(data_.conc_mobile_fe[sbi]->vec());

	for ( DHCellAccessor dh_cell : dh_->own_range() ) {
		LongIdx index = dh_cell.local_idx();
		ElementAccessor<3> ele_acc = mesh_->element_accessor( dh_cell.elm_idx() );
        
		for (unsigned int sbi=0; sbi<n_substances(); sbi++) { // Optimize: SWAP LOOPS
			vecs[sbi].set( index, data_.init_conc[sbi].value(ele_acc.centre(), ele_acc) );
        }
	}
}

//=============================================================================
//	ALLOCATE OF TRANSPORT VARIABLES (ELEMENT & NODES)
//=============================================================================
void ConvectionTransport::alloc_transport_vectors() {

    unsigned int sbi, n_subst;
    n_subst = n_substances();
    
    tm_diag = new double*[n_subst];
    cumulative_corr = new double*[n_subst];
    for (sbi = 0; sbi < n_subst; sbi++) {
      cumulative_corr[sbi] = new double[el_ds->lsize()];
      tm_diag[sbi] = new double[el_ds->lsize()];
    }

    data_.conc_mobile_fe.resize(n_subst);
    
    cfl_flow_ = new double[el_ds->lsize()];
    cfl_source_ = new double[el_ds->lsize()];
}

//=============================================================================
//	ALLOCATION OF TRANSPORT VECTORS (MPI)
//=============================================================================
void ConvectionTransport::alloc_transport_structs_mpi() {

    int sbi, n_subst;
    n_subst = n_substances();

    MPI_Barrier(PETSC_COMM_WORLD);

    vpconc = new Vec[n_subst];
    bcvcorr = new Vec[n_subst];
    vcumulative_corr = new Vec[n_subst];
    v_tm_diag = new Vec[n_subst];
    

    for (sbi = 0; sbi < n_subst; sbi++) {
        VecCreateMPI(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(), &bcvcorr[sbi]);
        VecZeroEntries(bcvcorr[sbi]);

        VecCreateMPI(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(), &vpconc[sbi]);
        VecZeroEntries(vpconc[sbi]);

        // SOURCES
        VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), mesh_->n_elements(),
        		cumulative_corr[sbi],&vcumulative_corr[sbi]);
        
        corr_vec.emplace_back(el_ds->lsize(), PETSC_COMM_WORLD);
        
        VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), mesh_->n_elements(),
                tm_diag[sbi],&v_tm_diag[sbi]);

        VecZeroEntries(vcumulative_corr[sbi]);
    }


    MatCreateAIJ(PETSC_COMM_WORLD, el_ds->lsize(), el_ds->lsize(), mesh_->n_elements(),
            mesh_->n_elements(), 16, PETSC_NULL, 4, PETSC_NULL, &tm);
    
    VecCreateMPI(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(), &mass_diag);
    VecCreateMPI(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(), &vpmass_diag);

    VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), mesh_->n_elements(),
            cfl_flow_,&vcfl_flow_);
    VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), mesh_->n_elements(),
            cfl_source_,&vcfl_source_);
}


//=============================================================================
// COMPUTE SOURCES, SET BOUNDARY CONDITIONS
//=============================================================================
void ConvectionTransport::conc_sources_bdr_conditions() {

    //temporary variables
    double csection, source, diag;
    unsigned int sbi;
    bool sources_changed = ( (data_.sources_density.changed() )
            || (data_.sources_conc.changed() )
            || (data_.sources_sigma.changed() )
            || (data_.cross_section.changed()));
    
    //TODO: would it be possible to check the change in data for chosen substance? (may be in multifields?)
  
	START_TIMER("sources_reinit_set_bc");
	if (sources_changed) balance_->start_source_assembly(subst_idx);
    // Assembly bcvcorr vector
    for(sbi=0; sbi < n_substances(); sbi++) VecZeroEntries(bcvcorr[sbi]);

   	balance_->start_flux_assembly(subst_idx);

	for ( DHCellAccessor dh_cell : dh_->own_range() )
	{
		ElementAccessor<3> elm = dh_cell.elm();
		// we have currently zero order P_Disc FE
		ASSERT_DBG(dh_cell.get_loc_dof_indices().size() == 1);
		IntIdx local_p0_dof = dh_cell.get_loc_dof_indices()[0];
        LongIdx glob_p0_dof = dh_->get_local_to_global_map()[local_p0_dof];

		arma::vec3 center = elm.centre();
		csection = data_.cross_section.value(center, elm);

		// SET SOURCES

		if (sources_changed) {
            // read for all substances
            double max_cfl=0;
            for (sbi = 0; sbi < n_substances(); sbi++)
            {
                double src_sigma = data_.sources_sigma[sbi].value(center, elm);

                source = csection * (data_.sources_density[sbi].value(center, elm)
                     + src_sigma * data_.sources_conc[sbi].value(center, elm));
                // addition to RHS
                corr_vec[sbi].set(local_p0_dof, source);
                // addition to diagonal of the transport matrix
                diag = src_sigma * csection;
                tm_diag[sbi][local_p0_dof] = - diag;

                // compute maximal cfl condition over all substances
                max_cfl = std::max(max_cfl, fabs(diag));

                balance_->add_source_values(sbi, elm.region().bulk_idx(), {local_p0_dof},
                                            {- src_sigma * elm.measure() * csection},
                                            {source * elm.measure()});
            }

            cfl_source_[local_p0_dof] = max_cfl;
		}

		// BOUNDARY CONDITIONS

        for(DHCellSide dh_side: dh_cell.side_range()) {
            if (dh_side.side().is_boundary()) {
                ElementAccessor<3> bc_elm = dh_side.cond().element_accessor();
                double flux = this->side_flux(dh_side);
                if (flux < 0.0) {
                    double aij = -(flux / elm.measure() );

                    for (sbi=0; sbi<n_substances(); sbi++)
                    {
                        double value = data_.bc_conc[sbi].value( bc_elm.centre(), bc_elm );

                        VecSetValue(bcvcorr[sbi], glob_p0_dof, value * aij, ADD_VALUES);

                        // CAUTION: It seems that PETSc possibly optimize allocated space during assembly.
                        // So we have to add also values that may be non-zero in future due to changing velocity field.
                        balance_->add_flux_values(subst_idx[sbi], dh_side,
                                                  {local_p0_dof}, {0.0}, flux*value);
                    }
                } else {
                    for (sbi=0; sbi<n_substances(); sbi++)
                        VecSetValue(bcvcorr[sbi], glob_p0_dof, 0, ADD_VALUES);

                    for (unsigned int sbi=0; sbi<n_substances(); sbi++)
                        balance_->add_flux_values(subst_idx[sbi], dh_side,
                                                  {local_p0_dof}, {flux}, 0.0);
                }
            }
        }
    }

    balance_->finish_flux_assembly(subst_idx);
    if (sources_changed) balance_->finish_source_assembly(subst_idx);

    for (sbi=0; sbi<n_substances(); sbi++)  	VecAssemblyBegin(bcvcorr[sbi]);
    for (sbi=0; sbi<n_substances(); sbi++)   	VecAssemblyEnd(bcvcorr[sbi]);

    // we are calling set_boundary_conditions() after next_time() and
    // we are using data from t() before, so we need to set corresponding bc time
    transport_bc_time = time_->last_t();
	END_TIMER("sources_reinit_set_bc");
}



void ConvectionTransport::zero_time_step()
{
	OLD_ASSERT_EQUAL(time_->tlevel(), 0);

	data_.mark_input_times(*time_);
	data_.set_time(time_->step(), LimitSide::right);
	std::stringstream ss; // print warning message with table of uninitialized fields
	if ( FieldCommon::print_message_table(ss, "convection transport") ) {
		WarningOut() << ss.str();
	}

    set_initial_condition();
    create_mass_matrix();
    
    START_TIMER("Convection balance zero time step");

    create_transport_matrix_mpi();
    conc_sources_bdr_conditions();

    // write initial condition
	output_data();
}


bool ConvectionTransport::evaluate_time_constraint(double& time_constraint)
{
	ASSERT_PTR(dh_).error( "Null DOF handler object.\n" );
    // read changed status before setting time
    bool changed_flux = data_.flow_flux.changed();
    data_.set_time(time_->step(), LimitSide::right); // set to the last computed time

    START_TIMER("data reinit");
    
    bool cfl_changed = false;
    
    // if FLOW or DATA changed ---------------------> recompute transport matrix
    if (changed_flux)
    {
        create_transport_matrix_mpi();
        is_convection_matrix_scaled=false;
        cfl_changed = true;
        DebugOut() << "CFL changed - flow.\n";
    }
    
    if (is_mass_diag_changed)
    {
        create_mass_matrix();
        cfl_changed = true;
        DebugOut() << "CFL changed - mass matrix.\n";
    }
    
    // if DATA changed ---------------------> recompute concentration sources (rhs and matrix diagonal)
    if( data_.sources_density.changed() || data_.sources_conc.changed() || data_.sources_sigma.changed()
       || data_.cross_section.changed() || data_.flow_flux.changed() || data_.porosity.changed()
       || data_.water_content.changed() || data_.bc_conc.changed() )
    {
        conc_sources_bdr_conditions();
        if( data_.sources_density.changed() || data_.sources_conc.changed() || data_.sources_sigma.changed()
               || data_.cross_section.changed() ) {
            is_src_term_scaled = false;
            cfl_changed = true;
        }
        if (data_.flow_flux.changed() || data_.porosity.changed()
        	       || data_.water_content.changed() || data_.bc_conc.changed() ) is_bc_term_scaled = false;
        DebugOut() << "CFL changed - source.\n";
    }
    
    // now resolve the CFL condition
    if(cfl_changed)
    {
        // find maximum of sum of contribution from flow and sources: MAX(vcfl_flow_ + vcfl_source_)
        Vec cfl;
        VecCreateMPI(PETSC_COMM_WORLD, el_ds->lsize(),PETSC_DETERMINE, &cfl);
        VecWAXPY(cfl, 1.0, vcfl_flow_, vcfl_source_);
        VecMaxPointwiseDivide(cfl,mass_diag, &cfl_max_step);
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
            VecScale(bcvcorr[sbi], dt);
        is_bc_term_scaled = true;
    }
    

    // if DATA or TIME STEP changed -----------------------> rescale source term
    if( !is_src_term_scaled || time_->is_changed_dt()) {
    	DebugOut() << "SRC - rescale dt.\n";
        //choose between fresh scaling with new dt or rescaling to a new dt
        double dt = (!is_src_term_scaled) ? dt_new : dt_scaled;
        for (unsigned int sbi=0; sbi<n_substances(); sbi++)
        {
            VecScale(corr_vec[sbi].petsc_vec(), dt);
            VecScale(v_tm_diag[sbi], dt);
        }
        is_src_term_scaled = true;
    }
    
    // if DATA or TIME STEP changed -----------------------> rescale transport matrix
    if ( !is_convection_matrix_scaled || time_->is_changed_dt()) {
    	DebugOut() << "TM - rescale dt.\n";
        //choose between fresh scaling with new dt or rescaling to a new dt
        double dt = (!is_convection_matrix_scaled) ? dt_new : dt_scaled;
        
        MatScale(tm, dt);
        is_convection_matrix_scaled = true;
    }
    
    END_TIMER("time step rescaling");
    
    
    data_.set_time(time_->step(), LimitSide::right); // set to the last computed time
    if (data_.cross_section.changed() || data_.water_content.changed() || data_.porosity.changed())
    {
        VecCopy(mass_diag, vpmass_diag);
        create_mass_matrix();
    } else is_mass_diag_changed = false;
    

    // Compute new concentrations for every substance.
    
    for (unsigned int sbi = 0; sbi < n_substances(); sbi++) {
      // one step in MOBILE phase
      START_TIMER("mat mult");
      Vec vconc = data_.conc_mobile_fe[sbi]->vec().petsc_vec();
      
      // tm_diag is a diagonal part of transport matrix, which depends on substance data (sources_sigma)
      // Wwe need keep transport matrix independent of substance, therefore we keep this diagonal part
      // separately in a vector tm_diag.
      // Computation: first, we compute this diagonal addition D*pconc and save it temporaly into RHS
        
      // RHS = D*pconc, where D is diagonal matrix represented by a vector
      VecPointwiseMult(vcumulative_corr[sbi], v_tm_diag[sbi], vconc); //w = x.*y
      
      // Then we add boundary terms ans other source terms into RHS.
      // RHS = 1.0 * bcvcorr + 1.0 * corr_vec + 1.0 * rhs
      VecAXPBYPCZ(vcumulative_corr[sbi], 1.0, 1.0, 1.0, bcvcorr[sbi], corr_vec[sbi].petsc_vec());   //z = ax + by + cz
      
      // Then we set the new previous concentration.
      VecCopy(vconc, vpconc[sbi]); // pconc = conc
      // And finally proceed with transport matrix multiplication.
      if (is_mass_diag_changed) {
        VecPointwiseMult(vconc, vconc, vpmass_diag); // vconc*=vpmass_diag
        MatMultAdd(tm, vpconc[sbi], vconc, vconc);   // vconc+=tm*vpconc
        VecAXPY(vconc, 1, vcumulative_corr[sbi]);    // vconc+=vcumulative_corr
        VecPointwiseDivide(vconc, vconc, mass_diag); // vconc/=mass_diag
      } else {
        MatMultAdd(tm, vpconc[sbi], vcumulative_corr[sbi], vconc);  // vconc =tm*vpconc+vcumulative_corr
        VecPointwiseDivide(vconc, vconc, mass_diag);                // vconc/=mass_diag
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


void ConvectionTransport::create_mass_matrix()
{
    ElementAccessor<3> elm;
    
    VecZeroEntries(mass_diag);
    
    balance_->start_mass_assembly(subst_idx);

    for ( DHCellAccessor dh_cell : dh_->own_range() ) {
        ElementAccessor<3> elm = dh_cell.elm();
        // we have currently zero order P_Disc FE
        ASSERT_DBG(dh_cell.get_loc_dof_indices().size() == 1);
        IntIdx local_p0_dof = dh_cell.get_loc_dof_indices()[0];

        double csection = data_.cross_section.value(elm.centre(), elm);
        //double por_m = data_.porosity.value(elm.centre(), elm->element_accessor());
        double por_m = data_.water_content.value(elm.centre(), elm);

        for (unsigned int sbi=0; sbi<n_substances(); ++sbi)
            balance_->add_mass_values(subst_idx[sbi], dh_cell, {local_p0_dof}, {csection*por_m*elm.measure()}, 0);
        
        VecSetValue(mass_diag, dh_->get_local_to_global_map()[local_p0_dof], csection*por_m, INSERT_VALUES);
    }
    
    balance_->finish_mass_assembly(subst_idx);
    
    VecAssemblyBegin(mass_diag);
    VecAssemblyEnd(mass_diag);
    
    is_mass_diag_changed = true;
}


//=============================================================================
// CREATE TRANSPORT MATRIX
//=============================================================================
void ConvectionTransport::create_transport_matrix_mpi() {

    START_TIMER("convection_matrix_assembly");

    ElementAccessor<3> el2;
    ElementAccessor<3> elm;
    int j;
    LongIdx new_j, new_i;
    double aij, aii;
        
    MatZeroEntries(tm);

    double flux, flux2, edg_flux;

    aii = 0.0;

    unsigned int loc_el = 0;
    for ( DHCellAccessor dh_cell : dh_->own_range() ) {
        new_i = row_4_el[ dh_cell.elm_idx() ];
        elm = dh_cell.elm();
        for( DHCellSide cell_side : dh_cell.side_range() ) {
            flux = this->side_flux(cell_side);
            if (! cell_side.side().is_boundary()) {
                edg_flux = 0;
                for( DHCellSide edge_side : cell_side.edge_sides() ) {
                    el2 = edge_side.element();
                    flux2 = this->side_flux(edge_side);
                    if ( flux2 > 0)  edg_flux+= flux2;
                }
                for( DHCellSide edge_side : cell_side.edge_sides() )
                    if (edge_side != cell_side) {
                        j = edge_side.element().idx();
                        new_j = row_4_el[j];

                        el2 = edge_side.element();
                        flux2 = this->side_flux(edge_side);
                        if ( flux2 > 0.0 && flux <0.0)
                            aij = -(flux * flux2 / ( edg_flux * dh_cell.elm().measure() ) );
                        else aij =0;
                        MatSetValue(tm, new_i, new_j, aij, INSERT_VALUES);
                    }
            }
            if (flux > 0.0)
              aii -= (flux / dh_cell.elm().measure() );
        }  // end same dim     //ELEMENT_SIDES

        for( DHCellSide neighb_side : dh_cell.neighb_sides() ) // dh_cell lower dim
        {
            ASSERT( neighb_side.elem_idx() != dh_cell.elm_idx() ).error("Elm. same\n");
            new_j = row_4_el[ neighb_side.elem_idx() ];
            el2 = neighb_side.element();
            flux = this->side_flux(neighb_side);

            // volume source - out-flow from higher dimension
            if (flux > 0.0)  aij = flux / dh_cell.elm().measure();
            else aij=0;
            MatSetValue(tm, new_i, new_j, aij, INSERT_VALUES);
            // out flow from higher dim. already accounted

            // volume drain - in-flow to higher dimension
            if (flux < 0.0) {
                aii -= (-flux) / dh_cell.elm().measure();                           // diagonal drain
                aij = (-flux) / neighb_side.element().measure();
            } else aij=0;
            MatSetValue(tm, new_j, new_i, aij, INSERT_VALUES);
        }

    MatSetValue(tm, new_i, new_i, aii, INSERT_VALUES);

    cfl_flow_[loc_el++] = fabs(aii);
    aii = 0.0;
    }

    MatAssemblyBegin(tm, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(tm, MAT_FINAL_ASSEMBLY);

    is_convection_matrix_scaled = false;
    END_TIMER("convection_matrix_assembly");

    transport_matrix_time = time_->t();
}


ConvectionTransport::FieldFEScalarVec& ConvectionTransport::get_p0_interpolation() {
    return data_.conc_mobile_fe;
}

void ConvectionTransport::get_par_info(LongIdx * &el_4_loc_out, Distribution * &el_distribution_out){
	el_4_loc_out = this->el_4_loc;
	el_distribution_out = this->el_ds;
	return;
}

//int *ConvectionTransport::get_el_4_loc(){
//	return el_4_loc;
//}

LongIdx *ConvectionTransport::get_row_4_el(){
	return row_4_el;
}



void ConvectionTransport::output_data() {

    data_.output_fields.set_time(time().step(), LimitSide::right);
	data_.output_fields.output(time().step());
    
    START_TIMER("TOS-balance");
    for (unsigned int sbi=0; sbi<n_substances(); ++sbi)
        balance_->calculate_instant(sbi, data_.conc_mobile_fe[sbi]->vec().petsc_vec());
    balance_->output();
    END_TIMER("TOS-balance");
}

void ConvectionTransport::set_balance_object(std::shared_ptr<Balance> balance)
{
	balance_ = balance;
    subst_idx = balance_->add_quantities(substances_.names());
}

double ConvectionTransport::side_flux(const DHCellSide &cell_side)
{
    if (cell_side.dim()==3) return calculate_side_flux<3>(cell_side);
    else if (cell_side.dim()==2) return calculate_side_flux<2>(cell_side);
    else return calculate_side_flux<1>(cell_side);
}

template<unsigned int dim>
double ConvectionTransport::calculate_side_flux(const DHCellSide &cell_side)
{
    ASSERT_EQ(cell_side.dim(), dim).error("Element dimension mismatch!");

    feo_.fe_values(dim).reinit(cell_side.side());
    auto vel = data_.flow_flux.value(cell_side.centre(), cell_side.element());
    double side_flux = arma::dot(vel, feo_.fe_values(dim).normal_vector(0)) * feo_.fe_values(dim).JxW(0);
    return side_flux;
}
