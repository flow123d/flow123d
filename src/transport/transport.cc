/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @ingroup transport
 * @brief  Transport
 *
 *
 */


#include <memory>

#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/partitioning.hh"
#include "transport/transport.h"

#include "la/distribution.hh"

#include "la/sparse_graph.hh"
#include <iostream>
#include <iomanip>
#include <string>

#include "io/output_time.hh"
#include "tools/time_governor.hh"
#include "flow/old_bcd.hh"
#include "coupling/balance.hh"
#include "input/accessors.hh"
#include "input/input_type.hh"

#include "fields/field_algo_base.hh"
#include "fields/field_values.hh"
#include "fields/field_elementwise.hh" 
#include "fields/generic_field.hh"

#include "reaction/isotherm.hh" // SorptionType enum

namespace IT = Input::Type;


const IT::Selection & ConvectionTransport::EqData::get_output_selection() {
	return EqData().output_fields
		.make_output_field_selection("ConvectionTransport_Output")
		.close();
}


ConvectionTransport::EqData::EqData() : TransportBase::TransportEqData()
{
	ADD_FIELD(bc_conc, "Boundary conditions for concentrations.", "0.0");
    	bc_conc.add_factory( OldBcdInput::instance()->trans_conc_factory );
    	bc_conc.units( UnitSI().kg().m(-3) );
	ADD_FIELD(init_conc, "Initial concentrations.", "0.0");
    	init_conc.units( UnitSI().kg().m(-3) );

    output_fields += *this;
    output_fields += conc_mobile.name("conc").units( UnitSI().kg().m(-3) );
	output_fields += region_id.name("region_id")
	        .units( UnitSI::dimensionless())
	        .flags(FieldFlag::equation_external_output);
}


ConvectionTransport::ConvectionTransport(Mesh &init_mesh, const Input::Record in_rec)
: TransportBase(init_mesh, in_rec), time_constraint_cfl("Convection CFL")
{
	START_TIMER("ConvectionTransport");
	this->eq_data_ = &data_;

    //mark type of the equation of convection transport (created in EquationBase constructor) and it is fixed
    time_ = new TimeGovernor(in_rec.val<Input::Record>("time"));
    target_mark_type = time_->equation_fixed_mark_type();
    time_->define_constraint(time_constraint_cfl, 
                             "Time step constrained due to CFL condition (including both flow and sources).");

    cfl_max_step = time_->end_time();

	// Initialize list of substances.
	substances_.initialize(in_rec.val<Input::Array>("substances"));
    n_subst_ = substances_.size();
    INPUT_CHECK(n_subst_ >= 1 ,"Number of substances must be positive.\n");

    data_.set_components(substances_.names());
    data_.set_mesh(init_mesh);
    data_.set_input_list( in_rec.val<Input::Array>("input_fields") );
    data_.set_limit_side(LimitSide::right);

    make_transport_partitioning();
    alloc_transport_vectors();
    alloc_transport_structs_mpi();
    transport_matrix_time = -1.0; // or -infty
    transport_bc_time = -1.0;
    is_convection_matrix_scaled = false;
    is_src_term_scaled = false;
    is_bc_term_scaled = false;

    // register output vectors
    output_rec = in_rec.val<Input::Record>("output_stream");
	data_.output_fields.set_components(substances_.names());
	data_.output_fields.set_mesh(*mesh_);
	data_.output_fields.set_limit_side(LimitSide::right);
	data_.output_fields.output_type(OutputTime::ELEM_DATA);
	data_.conc_mobile.set_up_components();
	data_.region_id = GenericField<3>::region_id(*mesh_);
	for (unsigned int sbi=0; sbi<n_subst_; sbi++)
	{
		// create shared pointer to a FieldElementwise and push this Field to output_field on all regions
		auto output_field_ptr = out_conc[sbi].create_field<3, FieldValue<3>::Scalar>(n_subst_);
		data_.conc_mobile[sbi].set_field(mesh_->region_db().get_region_set("ALL"), output_field_ptr, 0);
	}
	output_stream_ = OutputTime::create_output_stream(output_rec);
	output_stream_->add_admissible_field_names(in_rec.val<Input::Array>("output_fields"));
	output_stream_->mark_output_times(*time_);
}


//=============================================================================
// MAKE TRANSPORT
//=============================================================================
void ConvectionTransport::make_transport_partitioning() {

    int * id_4_old = new int[mesh_->n_elements()];
    int i = 0;
    FOR_ELEMENTS(mesh_, ele) id_4_old[i++] = ele.index();
    mesh_->get_part()->id_maps(mesh_->n_elements(), id_4_old, el_ds, el_4_loc, row_4_el);
    delete[] id_4_old;

    // TODO: make output of partitioning is usefull but makes outputs different
    // on different number of processors, which breaks tests.
    //
    // Possible solution:
    // - have flag in ini file to turn this output ON
    // - possibility to have different ref_output for different num of proc.
    // - or do not test such kind of output
    //
    //FOR_ELEMENTS(mesh_, ele) {
    //    ele->pid=el_ds->get_proc(row_4_el[ele.index()]);
    //}

}



ConvectionTransport::~ConvectionTransport()
{
    unsigned int sbi;

    //Destroy mpi vectors at first
    MatDestroy(&tm);
    VecDestroy(&vcfl_flow_);
    VecDestroy(&vcfl_source_);
    delete cfl_flow_;
    delete cfl_source_;

    for (sbi = 0; sbi < n_subst_; sbi++) {
        // mpi vectors
        VecDestroy(&(vconc[sbi]));
        VecDestroy(&(vpconc[sbi]));
        VecDestroy(&(bcvcorr[sbi]));
        VecDestroy(&(vcumulative_corr[sbi]));
        VecDestroy(&(v_tm_diag[sbi]));
        VecDestroy(&(v_sources_corr[sbi]));
        
        // arrays of arrays
        delete conc[sbi];
        delete cumulative_corr[sbi];
        delete tm_diag[sbi];
        delete sources_corr[sbi];
    }
    
    // arrays of mpi vectors
    delete vconc;
    delete vpconc;
    delete bcvcorr;
    delete vcumulative_corr;
    delete v_tm_diag;
    delete v_sources_corr;
    
    // arrays of arrays
    delete conc;
    delete cumulative_corr;
    delete tm_diag;
    delete sources_corr;
}





void ConvectionTransport::set_initial_condition()
{
    FOR_ELEMENTS(mesh_, elem)
    {
    	if (!el_ds->is_local(row_4_el[elem.index()])) continue;

    	unsigned int index = row_4_el[elem.index()] - el_ds->begin();
    	ElementAccessor<3> ele_acc = mesh_->element_accessor(elem.index());
		arma::vec value = data_.init_conc.value(elem->centre(), ele_acc);

		for (unsigned int sbi=0; sbi<n_subst_; sbi++)
			conc[sbi][index] = value(sbi);
    }

}

//=============================================================================
//	ALLOCATE OF TRANSPORT VARIABLES (ELEMENT & NODES)
//=============================================================================
void ConvectionTransport::alloc_transport_vectors() {

    unsigned int i, sbi;
    
    sources_corr = new double*[n_subst_];
    tm_diag = new double*[n_subst_];
    cumulative_corr = new double*[n_subst_];
    for (sbi = 0; sbi < n_subst_; sbi++) {
      cumulative_corr[sbi] = new double[el_ds->lsize()];
      sources_corr[sbi] = new double[el_ds->lsize()];
      tm_diag[sbi] = new double[el_ds->lsize()];
    }

    conc = new double*[n_subst_];
    out_conc.clear();
    out_conc.resize(n_subst_);
    for (sbi = 0; sbi < n_subst_; sbi++) {
        conc[sbi] = new double[el_ds->lsize()];
        out_conc[sbi].resize( el_ds->size() );
        for (i = 0; i < el_ds->lsize(); i++) {
            conc[sbi][i] = 0.0;
        }
    }
    
    cfl_flow_ = new double[el_ds->lsize()];
    cfl_source_ = new double[el_ds->lsize()];
}

//=============================================================================
//	ALLOCATION OF TRANSPORT VECTORS (MPI)
//=============================================================================
void ConvectionTransport::alloc_transport_structs_mpi() {

    unsigned int sbi;
    int rank, np;

    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);

    vconc = new Vec[n_subst_];
    vpconc = new Vec[n_subst_];
    bcvcorr = new Vec[n_subst_];
    vcumulative_corr = new Vec[n_subst_];
    v_tm_diag = new Vec[n_subst_];
    v_sources_corr = new Vec[n_subst_];
    

    for (sbi = 0; sbi < n_subst_; sbi++) {
        VecCreateMPI(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(), &bcvcorr[sbi]);
        VecZeroEntries(bcvcorr[sbi]);
        VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), mesh_->n_elements(), conc[sbi],
                &vconc[sbi]);

        VecCreateMPI(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(), &vpconc[sbi]);
        VecZeroEntries(vconc[sbi]);
        VecZeroEntries(vpconc[sbi]);

        // SOURCES
        VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), mesh_->n_elements(),
        		cumulative_corr[sbi],&vcumulative_corr[sbi]);
        
        VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), mesh_->n_elements(),
                sources_corr[sbi],&v_sources_corr[sbi]);
        
        VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), mesh_->n_elements(),
                tm_diag[sbi],&v_tm_diag[sbi]);

        VecZeroEntries(vcumulative_corr[sbi]);
        VecZeroEntries(out_conc[sbi].get_data_petsc());
    }


    MatCreateAIJ(PETSC_COMM_WORLD, el_ds->lsize(), el_ds->lsize(), mesh_->n_elements(),
            mesh_->n_elements(), 16, PETSC_NULL, 4, PETSC_NULL, &tm);

    VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), mesh_->n_elements(),
            cfl_flow_,&vcfl_flow_);
    VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), mesh_->n_elements(),
            cfl_source_,&vcfl_source_);
}


void ConvectionTransport::set_boundary_conditions()
{
    START_TIMER ("set_boundary_conditions");

    ElementFullIter elm = ELEMENT_FULL_ITER_NULL(mesh_);

    unsigned int sbi, loc_el, loc_b = 0;
    
    // Assembly bcvcorr vector
    for(sbi=0; sbi < n_subst_; sbi++) VecZeroEntries(bcvcorr[sbi]);

    if (balance_ != nullptr)
    	balance_->start_flux_assembly(subst_idx);

    for (loc_el = 0; loc_el < el_ds->lsize(); loc_el++) {
        elm = mesh_->element(el_4_loc[loc_el]);
        if (elm->boundary_idx_ != NULL) {
            unsigned int new_i = row_4_el[elm.index()];
            double csection = data_.cross_section.value(elm->centre(), elm->element_accessor());
            double por_m = data_.porosity.value(elm->centre(), elm->element_accessor());

            FOR_ELEMENT_SIDES(elm,si) {
                Boundary *b = elm->side(si)->cond();
                if (b != NULL) {
                    double flux = mh_dh->side_flux( *(elm->side(si)) );
                    if (flux < 0.0) {
                        double aij = -(flux / (elm->measure() * csection * por_m) );

                        arma::vec value = data_.bc_conc.value( b->element()->centre(), b->element_accessor() );
                        for (sbi=0; sbi<n_subst_; sbi++)
                            VecSetValue(bcvcorr[sbi], new_i, value[sbi] * aij, ADD_VALUES);

                        if (balance_ != nullptr)
                        {
                        	for (unsigned int sbi=0; sbi<n_substances(); sbi++)
                        	{
                        		balance_->add_flux_matrix_values(subst_idx[sbi], loc_b, {row_4_el[el_4_loc[loc_el]]}, {0.});
                        		balance_->add_flux_vec_value(subst_idx[sbi], loc_b, flux*value[sbi]);
                        	}
                        }
                    } else {
                    	if (balance_ != nullptr)
						{
							for (unsigned int sbi=0; sbi<n_substances(); sbi++)
							{
								balance_->add_flux_matrix_values(subst_idx[sbi], loc_b, {row_4_el[el_4_loc[loc_el]]}, {flux});
							}
						}
                    }
                    ++loc_b;
                }
            }

        }
    }

    if (balance_ != nullptr)
    	balance_->finish_flux_assembly(subst_idx);

    for (sbi=0; sbi<n_subst_; sbi++)  	VecAssemblyBegin(bcvcorr[sbi]);
    for (sbi=0; sbi<n_subst_; sbi++)   	VecAssemblyEnd(bcvcorr[sbi]);

    // we are calling set_boundary_conditions() after next_time() and
    // we are using data from t() before, so we need to set corresponding bc time
    transport_bc_time = time_->last_t();
}


//=============================================================================
// COMPUTE SOURCES
//=============================================================================
void ConvectionTransport::compute_concentration_sources() {

  //temporary variables
  unsigned int loc_el, sbi;
  double measure, por_m, source, diag;
  double max_cfl;
  Element *ele;
  ElementAccessor<3> ele_acc;
  arma::vec3 p;
  arma::vec src_density(n_subst_), src_conc(n_subst_), src_sigma(n_subst_);
    
  //TODO: would it be possible to check the change in data for chosen substance? (may be in multifields?)
  
  //checking if the data were changed
    if( (data_.sources_density.changed() )
          || (data_.sources_conc.changed() )
          || (data_.sources_sigma.changed() )
          || (data_.cross_section.changed())
          || (data_.porosity.changed() ))
    {
        START_TIMER("sources_reinit"); 
        if (balance_ != nullptr) balance_->start_source_assembly(subst_idx);
        
        for (loc_el = 0; loc_el < el_ds->lsize(); loc_el++) 
        {
            ele = mesh_->element(el_4_loc[loc_el]);
            ele_acc = ele->element_accessor();
            p = ele_acc.centre();
            por_m = data_.porosity.value(p, ele_acc);
            
            if (balance_ != nullptr) 
                measure = ele->measure() *
                          data_.cross_section.value(p, ele_acc);
            
            // read for all substances
            src_density = data_.sources_density.value(p, ele_acc);
            src_conc = data_.sources_conc.value(p, ele_acc);
            src_sigma = data_.sources_sigma.value(p, ele_acc);
                
            for (sbi = 0; sbi < n_subst_; sbi++) 
            {      
                source = src_density(sbi) + src_sigma(sbi) * src_conc(sbi);
                // addition to RHS
                sources_corr[sbi][loc_el] = source / por_m;
                // addition to diagonal of the transport matrix
                diag = src_sigma(sbi) / por_m;
                tm_diag[sbi][loc_el] = - diag;
                
                // compute maximal cfl condition over all substances
                max_cfl = std::max(max_cfl, diag);
                
                if (balance_ != nullptr)
                {
                    balance_->add_source_matrix_values(sbi, ele_acc.region().bulk_idx(), {row_4_el[el_4_loc[loc_el]]}, 
                                                       {- src_sigma(sbi) * measure});
                    balance_->add_source_rhs_values(sbi, ele_acc.region().bulk_idx(), {row_4_el[el_4_loc[loc_el]]}, 
                                                    {source * measure});
                }
            }
            
            cfl_source_[loc_el] = max_cfl;
            max_cfl = 0;
        }
        
        if (balance_ != nullptr) balance_->finish_source_assembly(subst_idx);
        
        END_TIMER("sources_reinit");
    }
}



void ConvectionTransport::zero_time_step()
{
	ASSERT_EQUAL(time_->tlevel(), 0);

	data_.mark_input_times(target_mark_type);
	data_.set_time(time_->step());

    set_initial_condition();

    if (balance_ != nullptr)
    {
    	START_TIMER("Convection balance zero time step");

    	create_transport_matrix_mpi();
        compute_concentration_sources();
    	set_boundary_conditions();

    	calculate_instant_balance();
    }

    // write initial condition
	output_data();
}


bool ConvectionTransport::evaluate_time_constraint(double& time_constraint)
{
    ASSERT(mh_dh, "Null MH object.\n" );
    data_.set_time(time_->step()); // set to the last computed time
    
    START_TIMER("data reinit");
    
    bool cfl_changed = false;
    
    // if FLOW or DATA changed ---------------------> recompute transport matrix
    if (mh_dh->time_changed() > transport_matrix_time  || data_.porosity.changed())
    {
        create_transport_matrix_mpi();
        is_convection_matrix_scaled=false;
        cfl_changed = true;
        DBGMSG("CFL changed - flow.\n");
    }
    
    // if DATA changed ---------------------> recompute concentration sources (rhs and matrix diagonal)
    if( data_.sources_density.changed() || data_.sources_conc.changed() || data_.sources_sigma.changed()
       || data_.cross_section.changed() || data_.porosity.changed() )
    {
        compute_concentration_sources();
        is_src_term_scaled = false;
        cfl_changed = true;
        DBGMSG("CFL changed - source.\n");
    }
    
    // now resolve the CFL condition
    if(cfl_changed)
    {
        // find maximum of sum of contribution from flow and sources: MAX(vcfl_flow_ + vcfl_source_)
        Vec cfl;
        VecCreateMPI(PETSC_COMM_WORLD, el_ds->lsize(),PETSC_DETERMINE, &cfl);
        VecWAXPY(cfl, 1.0, vcfl_flow_, vcfl_source_);
        VecMax(cfl,nullptr, &cfl_max_step);
        // get a reciprocal value as a time constraint
        cfl_max_step = 1 / cfl_max_step;
    }
    
    // although it does not influence CFL, compute BC so the full system is assembled
    if ( (mh_dh->time_changed() > transport_bc_time)
        || data_.porosity.changed()
        || data_.bc_conc.changed() )
    {
        set_boundary_conditions();
        is_bc_term_scaled = false;
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
    if ( ! is_bc_term_scaled || time_->is_changed_dt())
    {
        DBGMSG("BC - rescale dt.\n");
        //choose between fresh scaling with new dt or rescaling to a new dt
        double dt = (!is_bc_term_scaled) ? dt_new : dt_scaled;
        for (unsigned int  sbi=0; sbi<n_subst_; sbi++) 
            VecScale(bcvcorr[sbi], dt);
        is_bc_term_scaled = true;
    }

    // if DATA or TIME STEP changed -----------------------> rescale source term
    if( !is_src_term_scaled  || time_->is_changed_dt())
    {
        DBGMSG("SRC - rescale dt.\n");
        //choose between fresh scaling with new dt or rescaling to a new dt
        double dt = (!is_src_term_scaled) ? dt_new : dt_scaled;
        for (unsigned int sbi=0; sbi<n_subst_; sbi++)
        {
            VecScale(v_sources_corr[sbi], dt);
            VecScale(v_tm_diag[sbi], dt);
        }
        is_src_term_scaled = true;
    }
    
    // if DATA or TIME STEP changed -----------------------> rescale transport matrix
    if ( !is_convection_matrix_scaled || time_->is_changed_dt()) 
    {
        DBGMSG("TM - rescale dt.\n");
        //choose between fresh scaling with new dt or rescaling to a new dt
        double dt = (!is_convection_matrix_scaled) ? dt_new : dt_scaled;
        
        if(is_convection_matrix_scaled) MatShift(tm, -1.0);
        MatScale(tm, dt);
        MatShift(tm, 1.0);
        is_convection_matrix_scaled = true;
    }
    
    END_TIMER("time step rescaling");
    

    // Compute new concentrations for every substance.
    
    for (unsigned int sbi = 0; sbi < n_subst_; sbi++) {
      // one step in MOBILE phase
      START_TIMER("mat mult");
      
      // tm_diag is a diagonal part of transport matrix, which depends on substance data (sources_sigma)
      // Wwe need keep transport matrix independent of substance, therefore we keep this diagonal part
      // separately in a vector tm_diag.
      // Computation: first, we compute this diagonal addition D*pconc and save it temporaly into RHS
        
      // RHS = D*pconc, where D is diagonal matrix represented by a vector
      VecPointwiseMult(vcumulative_corr[sbi], v_tm_diag[sbi], vconc[sbi]); //w = x.*y
      
      // Then we add boundary terms ans other source terms into RHS.
      // RHS = 1.0 * bcvcorr + 1.0 * v_sources_corr + 1.0 * rhs
      VecAXPBYPCZ(vcumulative_corr[sbi], 1.0, 1.0, 1.0, bcvcorr[sbi], v_sources_corr[sbi]);   //z = ax + by + cz
      
      // Then we set the new previous concentration.
      VecCopy(vconc[sbi], vpconc[sbi]); // pconc = conc
      // And finally proceed with transport matrix multiplication.
      MatMultAdd(tm, vpconc[sbi], vcumulative_corr[sbi], vconc[sbi]); // conc=tm*pconc + bc
      END_TIMER("mat mult");
    }
    
    END_TIMER("convection-one step");
}


void ConvectionTransport::set_target_time(double target_time)
{

    //sets target_mark_type (it is fixed) to be met in next_time()
    time_->marks().add(TimeMark(target_time, target_mark_type));

    // This is done every time TOS calls update_solution.
    // If CFL condition is changed, time fixation will change later from TOS.
    
    // Set the same constraint as was set last time.
    time_->set_upper_constraint(time_constraint_cfl,cfl_max_step);
    
    // fixing convection time governor till next target_mark_type (got from TOS or other)
    // may have marks for data changes
    time_->fix_dt_until_mark();
}


//=============================================================================
// CREATE TRANSPORT MATRIX
//=============================================================================
void ConvectionTransport::create_transport_matrix_mpi() {

    START_TIMER("convection_matrix_assembly");

    ElementFullIter el2 = ELEMENT_FULL_ITER_NULL(mesh_);
    ElementFullIter elm = ELEMENT_FULL_ITER_NULL(mesh_);
    struct Edge *edg;
    unsigned int n;
    int s, j, np, rank, new_j, new_i;
    double aij, aii;
        
    MatZeroEntries(tm);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);

    double flux, flux2, edg_flux;


    vector<double> edge_flow(mesh_->n_edges(),0.0);
    for(unsigned int i=0; i < mesh_->n_edges() ; i++) { // calculate edge Qv
        Edge &edg = mesh_->edges[i];
        for( int s=0; s < edg.n_sides; s++) {
            flux = mh_dh->side_flux( *(edg.side(s)) );
            if ( flux > 0)  edge_flow[i]+= flux;
        }
    }

    if (balance_ != nullptr)
    	balance_->start_mass_assembly(subst_idx);

    aii = 0.0;

    for (unsigned int loc_el = 0; loc_el < el_ds->lsize(); loc_el++) {
        elm = mesh_->element(el_4_loc[loc_el]);
        new_i = row_4_el[elm.index()];

        double csection = data_.cross_section.value(elm->centre(), elm->element_accessor());
        double por_m = data_.porosity.value(elm->centre(), elm->element_accessor());

        if (balance_ != nullptr)
        {
        	for (unsigned int sbi=0; sbi<n_subst_; ++sbi)
        		balance_->add_mass_matrix_values(subst_idx[sbi], elm->region().bulk_idx(), {row_4_el[el_4_loc[loc_el]]}, {csection*por_m*elm->measure()} );
        }
        
        FOR_ELEMENT_SIDES(elm,si) {
            // same dim
            flux = mh_dh->side_flux( *(elm->side(si)) );
            if (elm->side(si)->cond() == NULL) {
                 edg = elm->side(si)->edge();
                 edg_flux = edge_flow[ elm->side(si)->edge_idx() ];
                 FOR_EDGE_SIDES(edg,s)
                    // this test should also eliminate sides facing to lower dim. elements in comp. neighboring
                    // These edges on these sides should have just one side
                    if (edg->side(s) != elm->side(si)) {
                        j = ELEMENT_FULL_ITER(mesh_, edg->side(s)->element()).index();
                        new_j = row_4_el[j];

                        flux2 = mh_dh->side_flux( *(edg->side(s)));
                        if ( flux2 > 0.0 && flux <0.0)
                            aij = -(flux * flux2 / ( edg_flux * elm->measure() * csection * por_m) );
                        else aij =0;
                        MatSetValue(tm, new_i, new_j, aij, INSERT_VALUES);
                    }
                if (flux > 0.0)
                    aii -= (flux / (elm->measure() * csection * por_m) );
            } else {
                if (flux > 0.0)
                    aii -= (flux / (elm->measure() * csection * por_m) );
            }
        }  // end same dim     //ELEMENT_SIDES

        FOR_ELM_NEIGHS_VB(elm,n) // comp model
            {
                el2 = ELEMENT_FULL_ITER(mesh_, elm->neigh_vb[n]->side()->element() ); // higher dim. el.
                ASSERT( el2 != elm, "Elm. same\n");
                new_j = row_4_el[el2.index()];
                flux = mh_dh->side_flux( *(elm->neigh_vb[n]->side()) );

                // volume source - out-flow from higher dimension
                if (flux > 0.0)  aij = flux / (elm->measure() * csection * por_m);
                else aij=0;
                MatSetValue(tm, new_i, new_j, aij, INSERT_VALUES);
                // out flow from higher dim. already accounted

                // volume drain - in-flow to higher dimension
                if (flux < 0.0) {
                        aii -= (-flux) / (elm->measure() * csection * por_m);                           // diagonal drain
                        aij = (-flux) / (el2->measure() *
                                        data_.cross_section.value(el2->centre(), el2->element_accessor()) *
                                        data_.porosity.value(el2->centre(), el2->element_accessor()));
                } else aij=0;
                MatSetValue(tm, new_j, new_i, aij, INSERT_VALUES);
            }

        MatSetValue(tm, new_i, new_i, aii, INSERT_VALUES);
        
        cfl_flow_[loc_el] = fabs(aii);
        aii = 0.0;
    } // END ELEMENTS

    if (balance_ != nullptr)
    	balance_->finish_mass_assembly(subst_idx);
    
    MatAssemblyBegin(tm, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(tm, MAT_FINAL_ASSEMBLY);

    is_convection_matrix_scaled = false;
    END_TIMER("convection_matrix_assembly");

    transport_matrix_time = time_->t();
}





//=============================================================================
//      OUTPUT VECTOR GATHER
//=============================================================================
void ConvectionTransport::output_vector_gather() {

    unsigned int sbi;
    IS is;

    ISCreateGeneral(PETSC_COMM_SELF, mesh_->n_elements(), row_4_el, PETSC_COPY_VALUES, &is); //WithArray
    VecScatterCreate(vconc[0], is, out_conc[0].get_data_petsc(), PETSC_NULL, &vconc_out_scatter);
    for (sbi = 0; sbi < n_subst_; sbi++) {
        VecScatterBegin(vconc_out_scatter, vconc[sbi], out_conc[sbi].get_data_petsc(), INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(vconc_out_scatter, vconc[sbi], out_conc[sbi].get_data_petsc(), INSERT_VALUES, SCATTER_FORWARD);
    }
    VecScatterDestroy(&(vconc_out_scatter));
    ISDestroy(&(is));
}


double **ConvectionTransport::get_concentration_matrix() {
	return conc;
}

void ConvectionTransport::get_par_info(int * &el_4_loc_out, Distribution * &el_distribution_out){
	el_4_loc_out = this->el_4_loc;
	el_distribution_out = this->el_ds;
	return;
}

int *ConvectionTransport::get_el_4_loc(){
	return el_4_loc;
}

int *ConvectionTransport::get_row_4_el(){
	return row_4_el;
}



void ConvectionTransport::calculate_cumulative_balance()
{
	for (unsigned int sbi=0; sbi<n_subst_; ++sbi)
	{
		balance_->calculate_cumulative_sources(sbi, vpconc[sbi], time_->dt());
		balance_->calculate_cumulative_fluxes(sbi, vpconc[sbi], time_->dt());
	}
}


void ConvectionTransport::calculate_instant_balance()
{
	for (unsigned int sbi=0; sbi<n_subst_; ++sbi)
	{
		balance_->calculate_mass(sbi, vconc[sbi]);
		balance_->calculate_source(sbi, vconc[sbi]);
		balance_->calculate_flux(sbi, vconc[sbi]);
	}
}


void ConvectionTransport::output_data() {

    if (time_->is_current( time_->marks().type_output() )) {
        output_vector_gather();

		data_.output_fields.set_time(time_->step());
		data_.output_fields.output(output_stream_);

    }
}

void ConvectionTransport::set_balance_object(boost::shared_ptr<Balance> balance)
{
	balance_ = balance;
	subst_idx = balance_->add_quantities(substances_.names());
}
