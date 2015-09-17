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


ConvectionTransport::ConvectionTransport(Mesh &init_mesh, const Input::Record &in_rec)
: TransportBase(init_mesh, in_rec)
{
	START_TIMER("ConvectionTransport");
	this->eq_data_ = &data_;

    //mark type of the equation of convection transport (created in EquationBase constructor) and it is fixed
    time_ = new TimeGovernor(in_rec.val<Input::Record>("time"));
    target_mark_type = time_->equation_fixed_mark_type();

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
    is_convection_matrix_scaled = false;
    need_time_rescaling=true;

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
    VecDestroy(&v_sources_corr);
    MatDestroy(&tm);
    
    VecDestroy(vconc);
    VecDestroy(bcvcorr);
    VecDestroy(vpconc);
    VecDestroy(vcumulative_corr);
    
    for (sbi = 0; sbi < n_subst_; sbi++) {
      //no mpi vectors
      xfree(sources_density[sbi]);
      xfree(sources_conc[sbi]);
      xfree(sources_sigma[sbi]);
      xfree(cumulative_corr[sbi]);
    }
    
    
    xfree(sources_corr);
    
    xfree(sources_density);
    xfree(sources_conc);
    xfree(sources_sigma);
    xfree(cumulative_corr);
    
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

    unsigned int i;
    int sbi, n_subst;
    n_subst = n_subst_;


    sources_corr = new double[el_ds->lsize()];
    sources_density = (double**) xmalloc(n_subst * sizeof(double*));
    sources_conc = (double**) xmalloc(n_subst * sizeof(double*));
    sources_sigma = (double**) xmalloc(n_subst * sizeof(double*));
    
    cumulative_corr = (double**) xmalloc(n_subst * sizeof(double*));
    for (sbi = 0; sbi < n_subst; sbi++) {
      sources_density[sbi] = (double*) xmalloc(el_ds->lsize() * sizeof(double));
      sources_conc[sbi] = (double*) xmalloc(el_ds->lsize() * sizeof(double));
      sources_sigma[sbi] = (double*) xmalloc(el_ds->lsize() * sizeof(double));
      cumulative_corr[sbi] = (double*) xmalloc(el_ds->lsize() * sizeof(double));
    }

    conc = (double**) xmalloc(n_subst * sizeof(double*));
    out_conc.clear();
    out_conc.resize(n_subst);
    for (sbi = 0; sbi < n_subst; sbi++) {
        conc[sbi] = (double*) xmalloc(el_ds->lsize() * sizeof(double));
        out_conc[sbi].resize( el_ds->size() );
        for (i = 0; i < el_ds->lsize(); i++) {
            conc[sbi][i] = 0.0;
        }
    }
}

//=============================================================================
//	ALLOCATION OF TRANSPORT VECTORS (MPI)
//=============================================================================
void ConvectionTransport::alloc_transport_structs_mpi() {

    int sbi, n_subst, rank, np;
    n_subst = n_subst_;

    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);

    bcvcorr = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vconc = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vpconc = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vcumulative_corr = (Vec*) xmalloc(n_subst * (sizeof(Vec)));


    VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), PETSC_DECIDE,
            sources_corr, &v_sources_corr);

    for (sbi = 0; sbi < n_subst; sbi++) {
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

        VecZeroEntries(vcumulative_corr[sbi]);
        VecZeroEntries(out_conc[sbi].get_data_petsc());
    }


    MatCreateAIJ(PETSC_COMM_WORLD, el_ds->lsize(), el_ds->lsize(), mesh_->n_elements(),
            mesh_->n_elements(), 16, PETSC_NULL, 4, PETSC_NULL, &tm);

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

}


//=============================================================================
// COMPUTE SOURCES
//=============================================================================
void ConvectionTransport::compute_concentration_sources(unsigned int sbi) {

  //temporary variables
  unsigned int loc_el;
  double conc_diff, csection, por_m;
  ElementAccessor<3> ele_acc;
  arma::vec3 p;
    
  //TODO: would it be possible to check the change in data for chosen substance? (may be in multifields?)
  
  //checking if the data were changed
    if( (data_.sources_density.changed() )
          || (data_.sources_conc.changed() )
          || (data_.sources_sigma.changed() )
          || (data_.cross_section.changed())
		  || (data_.porosity.changed() ))
      {
        START_TIMER("sources_reinit");
        for (loc_el = 0; loc_el < el_ds->lsize(); loc_el++) 
        {
          ele_acc = mesh_->element_accessor(el_4_loc[loc_el]);
          p = ele_acc.centre();
          
          por_m = data_.porosity.value(p, ele_acc);

          //if(data_.sources_density.changed_during_set_time)
          sources_density[sbi][loc_el] = data_.sources_density.value(p, ele_acc)(sbi)/por_m;
      
          //if(data_.sources_conc.changed_during_set_time)
          sources_conc[sbi][loc_el] = data_.sources_conc.value(p, ele_acc)(sbi);
        
          //if(data_.sources_sigma.changed_during_set_time)
          sources_sigma[sbi][loc_el] = data_.sources_sigma.value(p, ele_acc)(sbi)/por_m;
        }
        END_TIMER("sources_reinit");

        Element *ele;
        if (balance_ != nullptr)
        {
        	START_TIMER("Balance source assembly");
        	balance_->start_source_assembly(sbi);

        	//now computing source concentrations: density - sigma (source_conc - actual_conc)
        	for (loc_el = 0; loc_el < el_ds->lsize(); loc_el++)
            {
        		ele = mesh_->element(el_4_loc[loc_el]);
        		p = ele->centre();

        		csection = data_.cross_section.value(p, ele->element_accessor());
        		por_m = data_.porosity.value(p, ele->element_accessor());

        		balance_->add_source_matrix_values(sbi, ele->region().bulk_idx(), {row_4_el[el_4_loc[loc_el]]}, {sources_sigma[sbi][loc_el]*ele->measure()*por_m*csection});
        		balance_->add_source_rhs_values(sbi, ele->region().bulk_idx(), {row_4_el[el_4_loc[loc_el]]}, {sources_density[sbi][loc_el]*ele->measure()*por_m*csection});
            }

        	balance_->finish_source_assembly(sbi);
        	END_TIMER("Balance source assembly");
        }
      }


    //now computing source concentrations: density - sigma (source_conc - actual_conc)
    START_TIMER("calculate sources_corr");
    for (loc_el = 0; loc_el < el_ds->lsize(); loc_el++) 
        {
          conc_diff = sources_conc[sbi][loc_el] - conc[sbi][loc_el];
          if ( conc_diff > 0.0)
            sources_corr[loc_el] = ( sources_density[sbi][loc_el]
                                     + conc_diff * sources_sigma[sbi][loc_el] );
          else
            sources_corr[loc_el] = sources_density[sbi][loc_el];
        }
    END_TIMER("calculate sources_corr");
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
    	set_boundary_conditions();
    	for (unsigned int sbi=0; sbi<n_subst_; ++sbi)
    		compute_concentration_sources(sbi);

    	calculate_instant_balance();
    }

    // write initial condition
	output_data();
}



void ConvectionTransport::update_solution() {

    START_TIMER("convection-one step");


    START_TIMER("data reinit");
    data_.set_time(time_->step()); // set to the last computed time

    ASSERT(mh_dh, "Null MH object.\n" );
    // update matrix and sources if neccessary


    if (mh_dh->time_changed() > transport_matrix_time  || data_.porosity.changed()) {
        create_transport_matrix_mpi();
        is_convection_matrix_scaled=false;

        // need new fixation of the time step

        time_->set_upper_constraint(cfl_max_step);
        time_->fix_dt_until_mark();

        set_boundary_conditions();
        // scale boundary sources
        for (unsigned int sbi=0; sbi<n_subst_; sbi++) VecScale(bcvcorr[sbi], time_->estimate_dt());

        need_time_rescaling = true;
    } else {
        // possibly read boundary conditions
        if (data_.bc_conc.changed() ) {
            set_boundary_conditions();
            // scale boundary sources
            for (unsigned int  sbi=0; sbi<n_subst_; sbi++) VecScale(bcvcorr[sbi], time_->dt());
        }
    }

    if (need_time_rescaling) {
        if ( is_convection_matrix_scaled ) {
            // rescale matrix
            MatShift(tm, -1.0);
            MatScale(tm, time_->estimate_dt()/time_->dt() );
            MatShift(tm, 1.0);

            for (unsigned int sbi=0; sbi<n_subst_; sbi++) VecScale(bcvcorr[sbi], time_->estimate_dt()/time_->dt());

        } else {
            // scale fresh convection term matrix
            MatScale(tm, time_->estimate_dt());
            MatShift(tm, 1.0);
            is_convection_matrix_scaled = true;

        }
        need_time_rescaling = false;
    }

    END_TIMER("data reinit");


    // proceed to actually computed time
    // explicit scheme use values from previous time and then set then new time
    time_->next_time();


    for (unsigned int sbi = 0; sbi < n_subst_; sbi++) {
      // one step in MOBILE phase
      
      START_TIMER("compute_concentration_sources");
      //sources update  
      compute_concentration_sources(sbi);
     
      VecScale(v_sources_corr, time_->dt());
      //vcumulative_corr[sbi] = 1.0 * bcvcorr[sbi] + v_sources_corr;
      VecWAXPY(vcumulative_corr[sbi],1.0,bcvcorr[sbi],v_sources_corr);
      END_TIMER("compute_concentration_sources");

      START_TIMER("mat mult");
      VecCopy(vconc[sbi], vpconc[sbi]); // pconc = conc
      MatMultAdd(tm, vpconc[sbi], vcumulative_corr[sbi], vconc[sbi]); // conc=tm*pconc + bc
      END_TIMER("mat mult");
    }

    END_TIMER("convection-one step");
}


void ConvectionTransport::set_target_time(double target_time)
{

    //sets target_mark_type (it is fixed) to be met in next_time()
    time_->marks().add(TimeMark(target_time, target_mark_type));

    // make new time step fixation, invalidate scaling
    // same is done when matrix has changed in compute_one_step
    time_->set_upper_constraint(cfl_max_step);
    
    // fixing convection time governor till next target_mark_type (got from TOS or other)
    // may have marks for data changes
    time_->fix_dt_until_mark();
    need_time_rescaling = true;

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
    double max_sum, aij, aii;
        
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

    max_sum = 0.0;
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

        if (fabs(aii) > max_sum)
            max_sum = fabs(aii);
        aii = 0.0;
    } // END ELEMENTS

    if (balance_ != nullptr)
    	balance_->finish_mass_assembly(subst_idx);

    double glob_max_sum;

    MPI_Allreduce(&max_sum,&glob_max_sum,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
    cfl_max_step = 1 / glob_max_sum;
    
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
	Vec vpconc_diff;
	const double *pconc;
	double *pconc_diff;

	for (unsigned int sbi=0; sbi<n_subst_; ++sbi)
	{
		VecDuplicate(vpconc[sbi], &vpconc_diff);
		VecGetArrayRead(vpconc[sbi], &pconc);
		VecGetArray(vpconc_diff, &pconc_diff);
		for (unsigned int loc_el=0; loc_el<el_ds->lsize(); ++loc_el)
		{
			if (pconc[loc_el] < sources_conc[sbi][loc_el])
				pconc_diff[loc_el] = sources_conc[sbi][loc_el] - pconc[loc_el];
			else
				pconc_diff[loc_el] = 0;
		}
		balance_->calculate_cumulative_sources(sbi, vpconc_diff, time_->dt());
		balance_->calculate_cumulative_fluxes(sbi, vpconc[sbi], time_->dt());

		VecRestoreArray(vpconc_diff, &pconc_diff);
		VecRestoreArrayRead(vpconc[sbi], &pconc);
		VecDestroy(&vpconc_diff);
	}
}


void ConvectionTransport::calculate_instant_balance()
{
	Vec vconc_diff;
	const double *conc;
	double *conc_diff;

	for (unsigned int sbi=0; sbi<n_subst_; ++sbi)
	{
		VecDuplicate(vconc[sbi], &vconc_diff);
		VecGetArrayRead(vconc[sbi], &conc);
		VecGetArray(vconc_diff, &conc_diff);
		for (unsigned int loc_el=0; loc_el<el_ds->lsize(); ++loc_el)
		{
			if (conc[loc_el] < sources_conc[sbi][loc_el])
				conc_diff[loc_el] = sources_conc[sbi][loc_el] - conc[loc_el];
			else
				conc_diff[loc_el] = 0;
		}

		balance_->calculate_mass(sbi, vconc[sbi]);
		balance_->calculate_source(sbi, vconc_diff);
		balance_->calculate_flux(sbi, vconc[sbi]);

		VecRestoreArray(vconc_diff, &conc_diff);
		VecRestoreArrayRead(vconc[sbi], &conc);
		VecDestroy(&vconc_diff);
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
