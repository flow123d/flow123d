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

#include "io/output.h"

#include "la/distribution.hh"

#include "la/sparse_graph.hh"
#include <iostream>
#include <iomanip>
#include <string>

// TODO: move partitioning into mesh_ and remove this include
#include "flow/darcy_flow_mh.hh"
#include "flow/old_bcd.hh"
#include "input/accessors.hh"
#include "input/input_type.hh"

#include "coupling/time_governor.hh"

#include "fields/field_algo_base.hh"
#include "fields/field_values.hh"
#include "fields/field_elementwise.hh" 
#include "reaction/isotherm.hh" // SorptionType enum

namespace IT = Input::Type;

IT::Selection ConvectionTransport::EqData::sorption_type_selection = IT::Selection("TransportSorptionType")
    .add_value(none,"none","No sorption considered")
    .add_value(Isotherm::linear,"linear","Linear isotherm described sorption considered.")
    .add_value(Isotherm::freundlich,"freundlich","Freundlich isotherm described sorption considered")
    .add_value(Isotherm::langmuir,"langmuir","Langmuir isotherm described sorption considered")
    .close();


IT::Selection ConvectionTransport::EqData::output_selection =
		EqData().output_fields.make_output_field_selection("ConvectionTransport_Output")
		.close();


ConvectionTransport::EqData::EqData() : TransportBase::TransportEqData()
{
	ADD_FIELD(bc_conc, "Boundary conditions for concentrations.", "0.0");
    	bc_conc.read_field_descriptor_hook = OldBcdInput::trans_conc_hook;
    	bc_conc.units( UnitSI().kg().m(-3) );
	ADD_FIELD(init_conc, "Initial concentrations.", "0.0");
    	init_conc.units( UnitSI().kg().m(-3) );

    output_fields += *this;
    output_fields += conc_mobile.name("conc").units( UnitSI().kg().m(-3) );
}


ConvectionTransport::ConvectionTransport(Mesh &init_mesh, const Input::Record &in_rec)
: TransportBase(init_mesh, in_rec)
{
	this->eq_data_ = &data_;

    //mark type of the equation of convection transport (created in EquationBase constructor) and it is fixed
    time_ = new TimeGovernor(in_rec.val<Input::Record>("time"));
    target_mark_type = time_->equation_fixed_mark_type();

    cfl_max_step = time_->end_time();

    in_rec.val<Input::Array>("substances").copy_to(subst_names_);
    n_subst_ = subst_names_.size();
    INPUT_CHECK(n_subst_ >= 1 ,"Number of substances must be positive.\n");

    Input::Iterator<Input::Record> it = in_rec.find<Input::Record>("mass_balance");
    if (it) mass_balance_ = new MassBalance(this, *it);

    data_.set_n_components(n_subst_);
    data_.set_mesh(init_mesh);
    data_.set_input_list( in_rec.val<Input::Array>("input_fields") );
    data_.set_limit_side(LimitSide::right);


    sub_problem = 0;

    make_transport_partitioning();
    alloc_transport_vectors();
    alloc_transport_structs_mpi();
    transport_matrix_time = -1.0; // or -infty
    is_convection_matrix_scaled = false;
    need_time_rescaling=true;

    // register output vectors
    output_rec = in_rec.val<Input::Record>("output_stream");
	data_.conc_mobile.init(subst_names_);
	data_.conc_mobile.set_mesh(*mesh_);
	data_.output_fields.output_type(OutputTime::ELEM_DATA);

	for (unsigned int sbi=0; sbi<n_subst_; sbi++)
	{
		// create shared pointer to a FieldElementwise and push this Field to output_field on all regions
		std::shared_ptr<FieldElementwise<3, FieldValue<3>::Scalar> > output_field_ptr(new FieldElementwise<3, FieldValue<3>::Scalar>(out_conc[MOBILE][sbi], n_subst_, mesh_->n_elements()));
		data_.conc_mobile[sbi].set_field(mesh_->region_db().get_region_set("ALL"), output_field_ptr, 0);
	}
	data_.output_fields.set_limit_side(LimitSide::right);
	output_stream_ = OutputTime::create_output_stream(output_rec);
	output_stream_->add_admissible_field_names(in_rec.val<Input::Array>("output_fields"), data_.output_selection);
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

    if (mass_balance_ != NULL)
    	delete mass_balance_;

    //Destroy mpi vectors at first
    VecDestroy(&v_sources_corr);
    MatDestroy(&tm);
    
    VecDestroy(vconc);
    VecDestroy(bcvcorr);
    VecDestroy(vpconc);
    VecDestroy(vcumulative_corr);
    VecDestroy(vconc_out);
    
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
    
    delete output_stream_;

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
			conc[MOBILE][sbi][index] = value(sbi);
    }

}

//=============================================================================
//	ALLOCATE OF TRANSPORT VARIABLES (ELEMENT & NODES)
//=============================================================================
void ConvectionTransport::alloc_transport_vectors() {

    unsigned int i;
    int sbi, n_subst, ph;
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

    conc = (double***) xmalloc(MAX_PHASES * sizeof(double**));
    out_conc = (double***) xmalloc(MAX_PHASES * sizeof(double**));
    for (ph = 0; ph < MAX_PHASES; ph++) {
        if ((sub_problem & ph) == ph) {
            conc[ph] = (double**) xmalloc(n_subst * sizeof(double*));
            out_conc[ph] = (double**) xmalloc(n_subst * sizeof(double*));
            for (sbi = 0; sbi < n_subst; sbi++) {
                conc[ph][sbi] = (double*) xmalloc(el_ds->lsize() * sizeof(double));
                out_conc[ph][sbi] = (double*) xmalloc(el_ds->size() * sizeof(double));
                for (i = 0; i < el_ds->lsize(); i++) {
                    conc[ph][sbi][i] = 0.0;
                }
                for (i = 0; i < el_ds->size(); i++) {
                	out_conc[ph][sbi][i] = 0.0;
                }
            }
        } else {
            conc[ph] = NULL;
            out_conc[ph] = NULL;
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


    vconc_out = (Vec*) xmalloc(n_subst * (sizeof(Vec))); // extend to all
    

    VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), PETSC_DECIDE,
            sources_corr, &v_sources_corr);

    for (sbi = 0; sbi < n_subst; sbi++) {
        VecCreateMPI(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(), &bcvcorr[sbi]);
        VecZeroEntries(bcvcorr[sbi]);
        VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), mesh_->n_elements(), conc[MOBILE][sbi],
                &vconc[sbi]);

        VecCreateMPI(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(), &vpconc[sbi]);
        VecZeroEntries(vconc[sbi]);
        VecZeroEntries(vpconc[sbi]);

        // SOURCES
        VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), mesh_->n_elements(),
        		cumulative_corr[sbi],&vcumulative_corr[sbi]);

        VecCreateSeqWithArray(PETSC_COMM_SELF,1, mesh_->n_elements(), out_conc[MOBILE][sbi], &vconc_out[sbi]);

        VecZeroEntries(vcumulative_corr[sbi]);
        VecZeroEntries(vconc_out[sbi]);
    }


    MatCreateAIJ(PETSC_COMM_WORLD, el_ds->lsize(), el_ds->lsize(), mesh_->n_elements(),
            mesh_->n_elements(), 16, PETSC_NULL, 4, PETSC_NULL, &tm);

}


void ConvectionTransport::set_boundary_conditions()
{
    ElementFullIter elm = ELEMENT_FULL_ITER_NULL(mesh_);

    unsigned int sbi, loc_el;
    
    // Assembly bcvcorr vector
    for(sbi=0; sbi < n_subst_; sbi++) VecZeroEntries(bcvcorr[sbi]);


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
                    }
                }
            }

        }
    }

    for (sbi=0; sbi<n_subst_; sbi++)  	VecAssemblyBegin(bcvcorr[sbi]);
    for (sbi=0; sbi<n_subst_; sbi++)   	VecAssemblyEnd(bcvcorr[sbi]);

}


//=============================================================================
// COMPUTE SOURCES
//=============================================================================
void ConvectionTransport::compute_concentration_sources(unsigned int sbi) {

  //temporary variables
  unsigned int loc_el;
  double conc_diff, csection;
  ElementAccessor<3> ele_acc;
  arma::vec3 p;
    
  //TODO: would it be possible to check the change in data for chosen substance? (may be in multifields?)
  
  //checking if the data were changed
    if( (data_.sources_density.changed() )
          || (data_.sources_conc.changed() )
          || (data_.sources_sigma.changed() )
          || (data_.cross_section.changed()))
      {
        START_TIMER("sources_reinit");
        for (loc_el = 0; loc_el < el_ds->lsize(); loc_el++) 
        {
          ele_acc = mesh_->element_accessor(el_4_loc[loc_el]);
          p = ele_acc.centre();
          
          csection = data_.cross_section.value(p, ele_acc);

          //if(data_.sources_density.changed_during_set_time) 
          sources_density[sbi][loc_el] = data_.sources_density.value(p, ele_acc)(sbi)*csection;
      
          //if(data_.sources_conc.changed_during_set_time)
          sources_conc[sbi][loc_el] = data_.sources_conc.value(p, ele_acc)(sbi);
        
          //if(data_.sources_sigma.changed_during_set_time)
          sources_sigma[sbi][loc_el] = data_.sources_sigma.value(p, ele_acc)(sbi)*csection;
        }
      }
    
    //now computing source concentrations: density - sigma (source_conc - actual_conc)
    for (loc_el = 0; loc_el < el_ds->lsize(); loc_el++) 
        {
          conc_diff = sources_conc[sbi][loc_el] - conc[MOBILE][sbi][loc_el];
          if ( conc_diff > 0.0)
            sources_corr[loc_el] = ( sources_density[sbi][loc_el]
                                     + conc_diff * sources_sigma[sbi][loc_el] )
                                   * time_->dt();
          else
            sources_corr[loc_el] = sources_density[sbi][loc_el] * time_->dt();
        }
}

void ConvectionTransport::compute_concentration_sources_for_mass_balance(unsigned int sbi) {

	//temporary variables
	unsigned int loc_el;
	double conc_diff, csection;
	ElementAccessor<3> ele_acc;
	arma::vec3 p;

	double *pconc;
	VecGetArray(vpconc[sbi], &pconc);

	//TODO: would it be possible to check the change in data for chosen substance? (may be in multifields?)

	//checking if the data were changed
	if( (data_.sources_density.changed() )
		  || (data_.sources_conc.changed() )
		  || (data_.sources_sigma.changed() )
		  || (data_.cross_section.changed()))
	{
		START_TIMER("sources_reinit");
		for (loc_el = 0; loc_el < el_ds->lsize(); loc_el++)
		{
			ele_acc = mesh_->element_accessor(el_4_loc[loc_el]);
			p = ele_acc.centre();

			csection = data_.cross_section.value(p, ele_acc);

			//if(data_.sources_density.changed_during_set_time)
			sources_density[sbi][loc_el] = data_.sources_density.value(p, ele_acc)(sbi)*csection;

			//if(data_.sources_conc.changed_during_set_time)
			sources_conc[sbi][loc_el] = data_.sources_conc.value(p, ele_acc)(sbi);

			//if(data_.sources_sigma.changed_during_set_time)
			sources_sigma[sbi][loc_el] = data_.sources_sigma.value(p, ele_acc)(sbi)*csection;
		}
	}

    //now computing source concentrations: density - sigma (source_conc - actual_conc)
    for (loc_el = 0; loc_el < el_ds->lsize(); loc_el++)
    {
    	conc_diff = sources_conc[sbi][loc_el] - pconc[loc_el];
    	if ( conc_diff > 0.0)
    		sources_corr[loc_el] = sources_density[sbi][loc_el] + conc_diff * sources_sigma[sbi][loc_el];
    	else
            sources_corr[loc_el] = sources_density[sbi][loc_el];
    }
}


void ConvectionTransport::zero_time_step()
{
	ASSERT_EQUAL(time_->tlevel(), 0);

	data_.mark_input_times(target_mark_type);
	data_.set_time(*time_);

    set_initial_condition();


    // write initial condition
	output_data();
}



void ConvectionTransport::update_solution() {

    START_TIMER("convection-one step");
    
    unsigned int sbi;
    
    START_TIMER("data reinit");
    data_.set_time(*time_); // set to the last computed time

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
        for (sbi=0; sbi<n_subst_; sbi++) VecScale(bcvcorr[sbi], time_->estimate_dt());

        need_time_rescaling = true;
    } else {
        // possibly read boundary conditions
        if (data_.bc_conc.changed() ) {
            set_boundary_conditions();
            // scale boundary sources
            for (sbi=0; sbi<n_subst_; sbi++) VecScale(bcvcorr[sbi], time_->dt());
        }
    }

    if (need_time_rescaling) {
        if ( is_convection_matrix_scaled ) {
            // rescale matrix
            MatShift(tm, -1.0);
            MatScale(tm, time_->estimate_dt()/time_->dt() );
            MatShift(tm, 1.0);

            for (sbi=0; sbi<n_subst_; sbi++) VecScale(bcvcorr[sbi], time_->estimate_dt()/time_->dt());

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
    time_->next_time(); // explicit scheme use values from previous time and then set then new time


    for (sbi = 0; sbi < n_subst_; sbi++) {
      // one step in MOBILE phase
      
      START_TIMER("compute_concentration_sources");
      //sources update  
      compute_concentration_sources(sbi);  
     
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

    max_sum = 0.0;
    aii = 0.0;

    for (unsigned int loc_el = 0; loc_el < el_ds->lsize(); loc_el++) {
        elm = mesh_->element(el_4_loc[loc_el]);
        new_i = row_4_el[elm.index()];

        double csection = data_.cross_section.value(elm->centre(), elm->element_accessor());
        double por_m = data_.porosity.value(elm->centre(), elm->element_accessor());

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
    VecScatterCreate(vconc[0], is, vconc_out[0], PETSC_NULL, &vconc_out_scatter);
    for (sbi = 0; sbi < n_subst_; sbi++) {
        VecScatterBegin(vconc_out_scatter, vconc[sbi], vconc_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(vconc_out_scatter, vconc[sbi], vconc_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
    }
    VecScatterDestroy(&(vconc_out_scatter));
    ISDestroy(&(is));
}


double ***ConvectionTransport::get_concentration_matrix() {
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



void ConvectionTransport::calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance)
{
    double mass_flux[n_substances()];
    double *pconc[n_substances()];

    for (unsigned int sbi=0; sbi<n_substances(); sbi++)
    	VecGetArray(vpconc[sbi], &pconc[sbi]);

    FOR_BOUNDARIES(mesh_, bcd) {

        // !! there can be more sides per one boundary
        int index = row_4_el[bcd->side()->element().index()];
        if (!el_ds->is_local(index)) continue;
        int loc_index = index-el_ds->begin();

        double por_m = data_.porosity.value(bcd->side()->element()->centre(), bcd->side()->element()->element_accessor() );
        double water_flux = mh_dh->side_flux(*(bcd->side()));
        if (water_flux < 0) {
        	arma::vec bc_conc = data_.bc_conc.value( bcd->element()->centre(), bcd->element_accessor() );
        	for (unsigned int sbi=0; sbi<n_substances(); sbi++)
        		mass_flux[sbi] = water_flux*bc_conc[sbi]*por_m;
        } else {
        	for (unsigned int sbi=0; sbi<n_substances(); sbi++)
        		mass_flux[sbi] = water_flux*pconc[sbi][loc_index]*por_m;
        }

        Region r = bcd->region();
        if (! r.is_valid()) xprintf(Msg, "Invalid region, ele % d, edg: % d\n", bcd->bc_ele_idx_, bcd->edge_idx_);
        unsigned int bc_region_idx = r.boundary_idx();

        for (unsigned int sbi=0; sbi<n_substances(); sbi++)
        {
            bcd_balance[sbi][bc_region_idx] += mass_flux[sbi];

            if (water_flux > 0) bcd_plus_balance[sbi][bc_region_idx] += mass_flux[sbi];
            else bcd_minus_balance[sbi][bc_region_idx] += mass_flux[sbi];
        }
    }

}

void ConvectionTransport::calc_elem_sources(vector<vector<double> > &mass, vector<vector<double> > &src_balance)
{
    for (unsigned int sbi=0; sbi<n_substances(); sbi++)
    {
        compute_concentration_sources_for_mass_balance(sbi);
        double *sources = sources_corr;

        FOR_ELEMENTS(mesh_,elem)
        {
        	int index = row_4_el[elem.index()];
        	if (!el_ds->is_local(index)) continue;
        	ElementAccessor<3> ele_acc = elem->element_accessor();
        	double por_m = data_.porosity.value(elem->centre(), ele_acc);
        	double csection = data_.cross_section.value(elem->centre(), ele_acc);
        	int loc_index = index - el_ds->begin();

        	//temporary fix - mass balance works only when no reaction term is present
			double sum_sol_phases = conc[0][sbi][loc_index];

			mass[sbi][ele_acc.region().bulk_idx()] += por_m*csection*sum_sol_phases*elem->measure();
			src_balance[sbi][ele_acc.region().bulk_idx()] += sources[loc_index]*elem->measure();
        }
    }
}



void ConvectionTransport::output_data() {

    if (time_->is_current( time_->marks().type_output() )) {
        output_vector_gather();

		data_.output_fields.set_time(*time_);
		data_.output_fields.output(output_stream_);


        if (mass_balance_) 	mass_balance_->output(time_->t());

    }
}
