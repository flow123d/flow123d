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

//#include "ppfcs.h"
//#include "btc.h" XX
//#include "reaction.h" XX

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

#include "fields/field_base.hh"
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
		IT::Selection("ConvectionTransport_Output")
		.copy_values(EqData().output_fields.make_output_field_selection())
		.close();


ConvectionTransport::EqData::EqData() : TransportBase::TransportEqData()
{
	ADD_FIELD(bc_conc, "Boundary conditions for concentrations.", "0.0");
	ADD_FIELD(init_conc, "Initial concentrations.", "0.0");
    ADD_FIELD(por_imm, "Porosity material parameter of the immobile zone. Vector, one value for every substance.", "0.0");
    ADD_FIELD(alpha, "Diffusion coefficient of non-equilibrium linear exchange between mobile and immobile zone (dual porosity)."
            " Vector, one value for every substance.", "0.0");
    ADD_FIELD(sorp_type, "Type of sorption isotherm.", "\"none\"");
    sorp_type.input_selection(&sorption_type_selection);
    ADD_FIELD(sorp_coef0, "First parameter of sorption: Scaling of the isothem for all types. Vector, one value for every substance. ", "0.0");
    ADD_FIELD(sorp_coef1, "Second parameter of sorption: exponent( Freundlich isotherm), limit concentration (Langmuir isotherm). "
            "Vector, one value for every substance.", "1.0");
    ADD_FIELD(phi, "Fraction of the total sorption surface exposed to the mobile zone, in interval (0,1). "
            "Used only in combination with dual porosity model. Vector, one value for every substance.", "1.0");

    bc_conc.read_field_descriptor_hook = OldBcdInput::trans_conc_hook;

    output_fields += *this;
    output_fields += conc_mobile.name("mobile_p0").units("M/L^3");
}

/*
RegionSet ConvectionTransport::EqData::read_descriptor_hook(Input::Record rec) {
    // Base method EqDataBase::read_boundary_list_item must be called first!
    RegionSet domain = EqDataBase::read_boundary_list_item(rec);
    FilePath bcd_file;

    // read transport boundary conditions using old file format .tbc
    if (rec.opt_val("old_boundary_file", bcd_file) )
        OldBcdInput::instance()->read_transport(bcd_file, bc_conc);

    return domain;
}*/



ConvectionTransport::ConvectionTransport(Mesh &init_mesh, const Input::Record &in_rec)
: TransportBase(init_mesh, in_rec)
{
    //mark type of the equation of convection transport (created in EquationBase constructor) and it is fixed
    target_mark_type = this->mark_type() | TimeGovernor::marks().type_fixed_time();
    output_mark_type = this->mark_type() | TimeGovernor::marks().type_fixed_time() | time_->marks().type_output();
    time_ = new TimeGovernor(in_rec.val<Input::Record>("time"), target_mark_type);
    time_->marks().add_time_marks(0.0,
        in_rec.val<Input::Record>("output_stream").val<double>("time_step"),
        time_->end_time(), output_mark_type );

    cfl_max_step = time_->end_time();

    in_rec.val<Input::Array>("substances").copy_to(subst_names_);
    n_subst_ = subst_names_.size();
    INPUT_CHECK(n_subst_ >= 1 ,"Number of substances must be positive.\n");

    Input::Iterator<Input::Record> it = in_rec.find<Input::Record>("mass_balance");
    if (it)
    	mass_balance_ = new MassBalance(this, *it);

    data_.init_conc.n_comp(n_subst_);
    data_.bc_conc.n_comp(n_subst_);
    data_.alpha.n_comp(n_subst_);
    data_.sorp_type.n_comp(n_subst_);
    data_.sorp_coef0.n_comp(n_subst_);
    data_.sorp_coef1.n_comp(n_subst_);
    data_.sources_density.n_comp(n_subst_);
    data_.sources_sigma.n_comp(n_subst_);
    data_.sources_conc.n_comp(n_subst_);
    data_.set_mesh(init_mesh);
    data_.set_input_list( in_rec.val<Input::Array>("data") );
    data_.mark_input_times(target_mark_type);

    data_.set_limit_side(LimitSide::right);
    data_.set_time(*time_);


    sorption = in_rec.val<bool>("sorption_enable");
    dual_porosity = in_rec.val<bool>("dual_porosity");
    // reaction_on = in_rec.val<bool>("transport_reactions");


    sub_problem = 0;
    if (dual_porosity == true)
        sub_problem += 1;
    if (sorption == true)
        sub_problem += 2;

    make_transport_partitioning();
    alloc_transport_vectors();
    alloc_transport_structs_mpi();
    transport_matrix_time = -1.0; // or -infty
    set_initial_condition();

    is_convection_matrix_scaled = false;
    need_time_rescaling=true;

    // register output vectors
    output_rec = in_rec.val<Input::Record>("output_stream");
    int ierr, rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    ASSERT(ierr == 0, "Error in MPI_Comm_rank.");
    if (rank == 0)
    {
    	vector<string> output_names(subst_names_);
    	for (vector<string>::iterator it=output_names.begin(); it!=output_names.end(); it++)
    		*it += "_mobile";
    	data_.conc_mobile.init(output_names);
    	data_.conc_mobile.set_mesh(*mesh_);
    	data_.output_fields.output_type(OutputTime::ELEM_DATA);

    	for (int sbi=0; sbi<n_subst_; sbi++)
    	{
    		// create shared pointer to a FieldElementwise and push this Field to output_field on all regions
    		std::shared_ptr<FieldElementwise<3, FieldValue<3>::Scalar> > output_field_ptr(new FieldElementwise<3, FieldValue<3>::Scalar>(out_conc[MOBILE][sbi], n_subst_, mesh_->n_elements()));
    		data_.conc_mobile[sbi].set_field(mesh_->region_db().get_region_set("ALL"), output_field_ptr, 0);
    	}
        data_.output_fields.set_limit_side(LimitSide::right);
        output_stream_ = OutputTime::output_stream(output_rec);
        output_stream_->add_admissible_field_names(in_rec.val<Input::Array>("output_fields"), data_.output_selection);
    }

    // write initial condition
    output_vector_gather();
    if (rank == 0)
    {
    	data_.output_fields.set_time(*time_);
    	data_.output_fields.output(output_stream_);
    }
}


//=============================================================================
// MAKE TRANSPORT
//=============================================================================
void ConvectionTransport::make_transport_partitioning() {

    F_ENTRY;

    //int rank, np, i, j, k, row_MH, a;
    //struct DarcyFlowMH *water=transport->problem->water;
/*
    SparseGraph *ele_graph = new SparseGraphMETIS(mesh_->n_elements()); // graph for partitioning
    Distribution init_ele_ds = ele_graph->get_distr(); // initial distr.
    int *loc_part = new int[init_ele_ds.lsize()]; // partitionig in initial distribution

    make_element_connection_graph(mesh_, ele_graph, true);
    WARN_ASSERT(ele_graph->is_symmetric(),"Attention graph for partitioning is not symmetric!\n");

    ele_graph->partition(loc_part);

    delete ele_graph;
*/
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
    unsigned int sbi, ph;

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
    
    /*
    for (ph = 0; ph < MAX_PHASES; ph++) {
      if ((sub_problem & ph) == ph) {
        for (sbi = 0; sbi < n_subst_; sbi++) {
          xfree(conc[ph][sbi]);
          xfree(out_conc[ph][sbi]);   
        }
        xfree(conc[ph]);
        xfree(out_conc[ph]);
      }
    }
    
    DBGMSG("inner conc vecs freed\n");
    
    xfree(conc);
    xfree(out_conc);
    //*/
}

/*
//=============================================================================
// GET REACTION
//=============================================================================
void ConvectionTransport::get_reaction(int i,oReaction *reaction) {
	react[i] = *reaction;
}
*/
//=============================================================================
// RECOMPUTE MATRICES
//=============================================================================

/*
void ConvectionTransport::set_flow_field_vector(const MH_DofHandler &dh){
    // DBGMSG("set_flow_fieldvec\n");
    mh_dh = &dh;
	create_transport_matrix_mpi();
};
*/

void ConvectionTransport::set_cross_section_field(Field< 3, FieldValue<3>::Scalar >* cross_section) {
    data_.cross_section = cross_section;
}




void ConvectionTransport::set_initial_condition()
{
    FOR_ELEMENTS(mesh_, elem)
    {
    	if (!el_ds->is_local(row_4_el[elem.index()])) continue;

    	unsigned int index = row_4_el[elem.index()] - el_ds->begin();
    	ElementAccessor<3> ele_acc = mesh_->element_accessor(elem.index());
		arma::vec value = data_.init_conc.value(elem->centre(), ele_acc);

		for (int sbi=0; sbi<n_subst_; sbi++)
		{
			conc[MOBILE][sbi][index] = value(sbi);
			//pconc[MOBILE][sbi][index] = value(sbi);
		}
    }

}

//=============================================================================
//	ALLOCATE OF TRANSPORT VARIABLES (ELEMENT & NODES)
//=============================================================================
void ConvectionTransport::alloc_transport_vectors() {

    int i, sbi, n_subst, ph; //, j;
    //ElementIter elm;
    n_subst = n_subst_;


   // printf("%d\t\n",n_substances);
   // getchar();

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
    //pconc = (double***) xmalloc(MAX_PHASES * sizeof(double**));
    out_conc = (double***) xmalloc(MAX_PHASES * sizeof(double**));
    //transport->node_conc = (double****) xmalloc(MAX_PHASES * sizeof(double***));
    for (ph = 0; ph < MAX_PHASES; ph++) {
        if ((sub_problem & ph) == ph) {
            conc[ph] = (double**) xmalloc(n_subst * sizeof(double*)); //(MAX_PHASES * sizeof(double*));
            //pconc[ph] = (double**) xmalloc(n_subst * sizeof(double*));
            out_conc[ph] = (double**) xmalloc(n_subst * sizeof(double*));
            //  transport->node_conc[sbi] = (double***) xmalloc(MAX_PHASES * sizeof(double**));
            //}
            //}
            for (sbi = 0; sbi < n_subst; sbi++) {
                conc[ph][sbi] = (double*) xmalloc(el_ds->lsize() * sizeof(double));
                //pconc[ph][sbi] = (double*) xmalloc(el_ds->lsize() * sizeof(double));
                out_conc[ph][sbi] = (double*) xmalloc(el_ds->size() * sizeof(double));
                // transport->node_conc[sbi][ph] = (double**)xmalloc((mesh__->n_elements() ) * sizeof(double*));
                for (i = 0; i < el_ds->lsize(); i++) {
                    conc[ph][sbi][i] = 0.0;
                    //pconc[ph][sbi][i] = 0.0;

                }
                for (i = 0; i < el_ds->size(); i++) {
                	out_conc[ph][sbi][i] = 0.0;
                }
                /*
                 i = 0;
                 FOR_ELEMENTS(elm){
                 transport->node_conc[sbi][ph][i++]=(double*)xmalloc((elm->n_nodes) * sizeof(double));
                 for(j = 0 ;j < elm->n_nodes ; j++)
                 transport->node_conc[sbi][ph][i-1][j] = 0.0;
                 }*/
            }
        } else {
            conc[ph] = NULL;
            //pconc[ph] = NULL;
            out_conc[ph] = NULL;
            //transport->node_conc[sbi][ph] = NULL;
        }
    }
}

//=============================================================================
//	ALLOCATION OF TRANSPORT VECTORS (MPI)
//=============================================================================
void ConvectionTransport::alloc_transport_structs_mpi() {

    int sbi, n_subst, ierr, rank, np; //, i, j, ph;
    //ElementIter elm;
    n_subst = n_subst_;

    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);

    bcvcorr = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vconc = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vpconc = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vcumulative_corr = (Vec*) xmalloc(n_subst * (sizeof(Vec)));


    // if( rank == 0)
    vconc_out = (Vec*) xmalloc(n_subst * (sizeof(Vec))); // extend to all
    

    ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), PETSC_DECIDE,
            sources_corr, &v_sources_corr);

    for (sbi = 0; sbi < n_subst; sbi++) {
        ierr = VecCreateMPI(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(), &bcvcorr[sbi]);
        VecZeroEntries(bcvcorr[sbi]);
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), mesh_->n_elements(), conc[MOBILE][sbi],
                &vconc[sbi]);

//        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(),
//                pconc[MOBILE][sbi], &vpconc[sbi]);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(), &vpconc[sbi]);
        VecZeroEntries(vconc[sbi]);
        VecZeroEntries(vpconc[sbi]);

        // SOURCES
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1, el_ds->lsize(), mesh_->n_elements(),
        		cumulative_corr[sbi],&vcumulative_corr[sbi]);

        //  if(rank == 0)
        ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1, mesh_->n_elements(), out_conc[MOBILE][sbi], &vconc_out[sbi]);

        VecZeroEntries(vcumulative_corr[sbi]);
        VecZeroEntries(vconc_out[sbi]);
    }


    ierr = MatCreateAIJ(PETSC_COMM_WORLD, el_ds->lsize(), el_ds->lsize(), mesh_->n_elements(),
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
            double csection = data_.cross_section->value(elm->centre(), elm->element_accessor());
            double por_m = data_.por_m.value(elm->centre(), elm->element_accessor());

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


    //VecView(bcvcorr[0],PETSC_VIEWER_STDOUT_SELF);
    //exit(0);
}


//=============================================================================
// COMPUTE SOURCES
//=============================================================================
void ConvectionTransport::compute_concentration_sources(unsigned int sbi) {

  //temporary variables
  unsigned int loc_el;
  double conc_diff;
  ElementAccessor<3> ele_acc;
  arma::vec3 p;
    
  //TODO: would it be possible to check the change in data for chosen substance? (may be in multifields?)
  
  //checking if the data were changed
    if( (data_.sources_density.changed() )
          || (data_.sources_conc.changed() )
          || (data_.sources_sigma.changed() ) )
      {
        START_TIMER("sources_reinit");
        for (loc_el = 0; loc_el < el_ds->lsize(); loc_el++) 
        {
          ele_acc = mesh_->element_accessor(el_4_loc[loc_el]);
          p = ele_acc.centre();
          
          //if(data_.sources_density.changed_during_set_time) 
          sources_density[sbi][loc_el] = data_.sources_density.value(p, ele_acc)(sbi);
      
          //if(data_.sources_conc.changed_during_set_time)
          sources_conc[sbi][loc_el] = data_.sources_conc.value(p, ele_acc)(sbi);
        
          //if(data_.sources_sigma.changed_during_set_time)
          sources_sigma[sbi][loc_el] = data_.sources_sigma.value(p, ele_acc)(sbi);
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
	double conc_diff;
	ElementAccessor<3> ele_acc;
	arma::vec3 p;

	double *pconc;
	VecGetArray(vpconc[sbi], &pconc);

	//TODO: would it be possible to check the change in data for chosen substance? (may be in multifields?)

	//checking if the data were changed
	if( (data_.sources_density.changed() )
		  || (data_.sources_conc.changed() )
		  || (data_.sources_sigma.changed() ) )
	{
		START_TIMER("sources_reinit");
		for (loc_el = 0; loc_el < el_ds->lsize(); loc_el++)
		{
			ele_acc = mesh_->element_accessor(el_4_loc[loc_el]);
			p = ele_acc.centre();

			//if(data_.sources_density.changed_during_set_time)
			sources_density[sbi][loc_el] = data_.sources_density.value(p, ele_acc)(sbi);

			//if(data_.sources_conc.changed_during_set_time)
			sources_conc[sbi][loc_el] = data_.sources_conc.value(p, ele_acc)(sbi);

			//if(data_.sources_sigma.changed_during_set_time)
			sources_sigma[sbi][loc_el] = data_.sources_sigma.value(p, ele_acc)(sbi);
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


void ConvectionTransport::compute_one_step() {

    START_TIMER("convection-one step");
    
    unsigned int loc_el,sbi;
    
    START_TIMER("data reinit");
    data_.set_time(*time_); // set to the last computed time

    ASSERT(mh_dh, "Null MH object.\n" );
    // update matrix and sources if neccessary


    if (mh_dh->time_changed() > transport_matrix_time  || data_.por_m.changed()) {
        DBGMSG("mh time: %f tm: %f por: %d\n", mh_dh->time_changed(), transport_matrix_time, data_.por_m.changed());
        create_transport_matrix_mpi();

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
            for (sbi=0; sbi<n_subst_; sbi++) VecScale(bcvcorr[sbi], time_->estimate_dt());
        }
    }

    if (need_time_rescaling) {
        if ( is_convection_matrix_scaled ) {
            // rescale matrix
            //for (unsigned int sbi=0; sbi<n_substances; sbi++) VecScale(bcvcorr[sbi], time_->dt()/time_->estimate_dt());
            MatShift(tm, -1.0);
            MatScale(tm, time_->estimate_dt()/time_->dt() );
            MatShift(tm, 1.0);
            DBGMSG("rescaling matrix\n");

            for (sbi=0; sbi<n_subst_; sbi++) VecScale(bcvcorr[sbi], time_->estimate_dt()/time_->dt());

        } else {
            // scale fresh convection term matrix
            //for (unsigned int sbi=0; sbi<n_substances; sbi++) VecScale(bcvcorr[sbi], time_->estimate_dt());
            MatScale(tm, time_->estimate_dt());
            MatShift(tm, 1.0);
            is_convection_matrix_scaled = true;

        }
        need_time_rescaling = false;
    }

    // update source vectors
//    for (unsigned int sbi = 0; sbi < n_substances; sbi++) {
//            MatMult(bcm, bcv[sbi], bcvcorr[sbi]);
//            VecView(bcv[sbi],PETSC_VIEWER_STDOUT_SELF);
//            getchar();
//            VecView(bcvcorr[sbi],PETSC_VIEWER_STDOUT_SELF);
//           getchar();
//    }


    END_TIMER("data reinit");


    // proceed to actually computed time
    //time_->view("CONVECTION");
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
      //VecView(vconc[sbi],PETSC_VIEWER_STDOUT_SELF);
      END_TIMER("mat mult");

     //}

     //START_TIMER("dual porosity/old-sorption");
     
    /*if(sorption == true) for(int loc_el = 0; loc_el < el_ds->lsize(); loc_el++)
    {
      for(int i_subst = 0; i_subst < n_subst_; i_subst++)
      {
        //following conditional print is here for comparison of old and new type of sorption input concentrations
        if(i_subst < (n_subst_ - 1)) cout << conc[MOBILE][i_subst][loc_el] << ", ";
          else cout << conc[MOBILE][i_subst][loc_el] << endl;
      }
	}
    for (sbi = 0; sbi < n_subst_; sbi++) {*/
    
    START_TIMER("old_sorp_step");
        if ((dual_porosity == true) || (sorption == true) )
            // cycle over local elements only in any order
            for (loc_el = 0; loc_el < el_ds->lsize(); loc_el++) {

                if (dual_porosity == true)
                    transport_dual_porosity(loc_el, mesh_->element(el_4_loc[loc_el]), sbi);
                if (sorption == true)
                    transport_sorption(loc_el, mesh_->element(el_4_loc[loc_el]), sbi);

            }
        // transport_node_conc(mesh_,sbi,problem->transport_sub_problem);  // vyresit prepocet
      //END_TIMER("dual porosity/old-sorption");
    }
    END_TIMER("old_sorp_step");
    
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

/*
void ConvectionTransport::preallocate_transport_matrix() {

    for (unsigned int loc_el = 0; loc_el < el_ds->lsize(); loc_el++) {
        elm = mesh_->element(el_4_loc[loc_el]);
        new_i = row_4_el[elm.index()];

        n_on_proc=1; n_off_proc=0;
        FOR_ELEMENT_SIDES(elm,si) {
            if (elm->side(si)->cond() == NULL) {
                    edg = elm->side(si)->edge();
                    FOR_EDGE_SIDES(edg,s) {
                           j = ELEMENT_FULL_ITER(mesh_, edg->side(s)->element()).index();
                           new_j = row_4_el[j];
                           if (el_ds->is_local(new_j)) n_on_proc++;
                           else n_off_proc++;
                    }
                    n_on_proc--; // do not count diagonal entry more then once
            }
        }

        FOR_ELM_NEIGHS_VB(elm,n) // comp model
            {
                el2 = ELEMENT_FULL_ITER(mesh_, elm->neigh_vb[n]->side()->element() ); // higher dim. el.
                ASSERT( el2 != elm, "Elm. same\n");
                    if (flux > 0.0) {
                        // volume source - out-flow from higher dimension
                        j = el2.index();
                        new_j = row_4_el[j];
                        MatSetValue(tm, new_i, new_j, aij, INSERT_VALUES);
                        // out flow from higher dim. already accounted
                    }
                    if (flux < 0.0) {
                        // volume drain - in-flow to higher dimension
                        aij = (-flux) / (el2->measure() *
                                        data_.cross_section->value(el2->centre(), el2->element_accessor()) *
                                        data_.por_m.value(el2->centre(), el2->element_accessor()));
                        new_j = row_4_el[el2.index()];
                        MatSetValue(tm, new_j, new_i, aij, INSERT_VALUES);

                        // diagonal drain
                        aii -= (-flux) / (elm->measure() * csection * por_m);
                    }

                //} // end comp model
            }
    } // END ELEMENTS

}*/

//=============================================================================
// CREATE TRANSPORT MATRIX
//=============================================================================
void ConvectionTransport::create_transport_matrix_mpi() {

    DBGMSG("TM assembly\n");
    START_TIMER("convection_matrix_assembly");

    ElementFullIter el2 = ELEMENT_FULL_ITER_NULL(mesh_);
    ElementFullIter elm = ELEMENT_FULL_ITER_NULL(mesh_);
    struct Edge *edg;
    //struct Neighbour *ngh;
    //struct Transport *transport;
    int n, s, j, np, rank, new_j, new_i; //, i;
    double max_sum, aij, aii; //, *solution;
    /*
    DarcyFlow *water;

    water = problem->water;
    solution = water.solution();
    id = water
*/

    /*
     FOR_ELEMENTS(elm)
     FOR_ELEMENT_SIDES(elm,i){
     printf("id: %d side: %d flux:  %f\n",elm->id, i,elm->side(i)->flux);
     getchar();
     }

     getchar();
     */

        
    MatZeroEntries(tm);
//    MatZeroEntries(bcm);


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

        double csection = data_.cross_section->value(elm->centre(), elm->element_accessor());
        double por_m = data_.por_m.value(elm->centre(), elm->element_accessor());

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
//                if (flux < 0.0) {
//                    aij = -(flux / (elm->measure() * csection * por_m) );
//                    j = elm->side(si)->cond_idx() ;
                    // DBGMSG("BCM, i: %d j:%d\n", new_i , j);
//                    MatSetValue(bcm, new_i, j, aij, INSERT_VALUES);
                    // vyresit BC matrix !!!!
                    //   printf("side in elm:%d value:%f\n ",elm->id,svector->val[j-1]);
                    //   printf("%d\t%d\n",elm->id,id2pos(problem,elm->side(si)->id,problem->spos_id,BC));

//                }
                if (flux > 0.0)
                    aii -= (flux / (elm->measure() * csection * por_m) );
            }
        }  // end same dim     //ELEMENT_SIDES

        FOR_ELM_NEIGHS_VB(elm,n) // comp model
            //FOR_NEIGH_ELEMENTS(elm->neigh_vb[n],s)
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
                                        data_.cross_section->value(el2->centre(), el2->element_accessor()) *
                                        data_.por_m.value(el2->centre(), el2->element_accessor()));
                } else aij=0;
                MatSetValue(tm, new_j, new_i, aij, INSERT_VALUES);
            }

        MatSetValue(tm, new_i, new_i, aii, INSERT_VALUES);

        if (fabs(aii) > max_sum)
            max_sum = fabs(aii);
        aii = 0.0;
        //   i++;
    } // END ELEMENTS

    double glob_max_sum;

    MPI_Allreduce(&max_sum,&glob_max_sum,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
    //xprintf(Msg,"CFL: glob_max_sum=%f\n",glob_max_sum);
    cfl_max_step = 1 / glob_max_sum;
    //time_step = 0.9 / glob_max_sum;
    
    DBGMSG("start assembly\n");
    MatAssemblyBegin(tm, MAT_FINAL_ASSEMBLY);
//    MatAssemblyBegin(bcm, MAT_FINAL_ASSEMBLY);
    

    MatAssemblyEnd(tm, MAT_FINAL_ASSEMBLY);
//    MatAssemblyEnd(bcm, MAT_FINAL_ASSEMBLY);
    DBGMSG("end assembly\n");


    // MPI_Barrier(PETSC_COMM_WORLD);

     //MatView(tm,PETSC_VIEWER_STDOUT_SELF);


    is_convection_matrix_scaled = false;
    END_TIMER("convection_matrix_assembly");

    transport_matrix_time = time_->t();
}


//=============================================================================
// COMPUTE CONCENTRATIONS IN THE NODES FROM THE ELEMENTS
//=============================================================================
/*
 void transport_node_conc(struct Transport *transport)
 {
 TNode* nod;
 mesh_* mesh_ = transport->mesh_;
 int i;
 double P,N,scale, *pi;
 int min_elm_dim;
 int sbi,sub;
 pi = transport_aloc_pi(mesh_);

 FOR_NODES(nod)
 {
 P = N = 0;
 min_elm_dim=3;
 FOR_NODE_ELEMENTS(nod,i)
 {
 if (nod->element[i]->dim < min_elm_dim)
 min_elm_dim =  nod->element[i]->dim;
 }
 //printf("%d %d %d\n",nod->id,i,min_elm_dim);
 FOR_NODE_ELEMENTS(nod,i)
 {
 if  (nod->element[i]->dim == min_elm_dim)
 {
 pi[ i ] = sqrt( SQUARE(nod->x - nod->element[i]->centre[0]) +
 SQUARE(nod->y - nod->element[i]->centre[1]) +
 SQUARE(nod->z - nod->element[i]->centre[2]) );
 if (pi[i] > ZERO)
 P += 1 / pi[ i ]; //       {	P += pi[ i ]; N++;}
 }
 //else  printf("vynechano %d ",i);

 }

 //node_init(nod,sbi,sub);

 FOR_NODE_ELEMENTS(nod,i) {
 if ((nod->element[i]->dim == min_elm_dim) && (P != 0)) {
 if (pi[i] > ZERO) { //       if ((pi[i]> ZERO) && (N!=1)){
 scale = 1 / (pi[i] * P); //   scale = (1 - (pi[ i ] / P)) / (N - 1);
 nod->conc[sbi] += nod->element[i]->conc[sbi] * scale;
 if ((sub & 1) == 1)
 nod->conc_immobile[sbi]
 += nod->element[i]->conc_immobile[sbi] * scale;
 if ((sub & 2) == 2)
 nod->conc_sorb[sbi] += nod->element[i]->conc_sorb[sbi]
 * scale;
 if ((sub & 3) == 3)
 nod->conc_immobile_sorb[sbi]
 += nod->element[i]->conc_immobile_sorb[sbi]
 * scale;
 } else {
 nod->conc[sbi] = nod->element[i]->conc[sbi];
 if ((sub & 1) == 1)
 nod->conc_immobile[sbi]
 = nod->element[i]->conc_immobile[sbi];
 if ((sub & 2) == 2)
 nod->conc_sorb[sbi] = nod->element[i]->conc_sorb[sbi];
 if ((sub & 3) == 3)
 nod->conc_immobile_sorb[sbi]
 = nod->element[i]->conc_immobile_sorb[sbi];
 break;
 }
 }
 } //for node elm
 }
 xfree( pi );
 }
 */
//=============================================================================
// COMPUTE DUAL POROSITY  ( ANALYTICAL EQUATION )
//
// assume: elm_pos - position of values of pconc and conc vectors corresponding to a mesh_ element
//         sbi - matter index
//         material - material on corresponding mesh_ element
//=============================================================================
void ConvectionTransport::transport_dual_porosity( int elm_pos, ElementFullIter elem, int sbi) {

    double conc_avg = 0.0;
    unsigned int loc_el;
    //int id;
    //double ***conc = transport->conc;
    //double ***pconc = transport->pconc;
    double cm, pcm, ci, pci, por_m, por_imm, alpha;

    		por_m = data_.por_m.value(elem->centre(), elem->element_accessor());
    		por_imm = data_.por_imm.value(elem->centre(), elem->element_accessor());
    		alpha = data_.alpha.value(elem->centre(), elem->element_accessor())(sbi);
    		pcm = conc[MOBILE][sbi][elm_pos];
    		pci = conc[IMMOBILE][sbi][elm_pos];
    		// ---compute average concentration------------------------------------------
    		conc_avg = ((por_m * pcm) + (por_imm * pci)) / (por_m + por_imm);

    		if ((conc_avg != 0.0) && (por_imm != 0.0)) {
        		// ---compute concentration in mobile area-----------------------------------
        		cm = (pcm - conc_avg) * exp(-alpha * ((por_m + por_imm) / (por_m * por_imm)) * time_->dt()) + conc_avg;

        		// ---compute concentration in immobile area---------------------------------
        		ci = (pci - conc_avg) * exp(-alpha * ((por_m + por_imm) / (por_m * por_imm)) * time_->dt()) + conc_avg;
        		// --------------------------------------------------------------------------
        		//printf("\n%f\t%f\t%f",conc_avg,cm,ci);
//                         DBGMSG("cm: %f  ci: %f  pcm: %f  pci: %f  conc_avg: %f  alpha: %f  por_m: %f  por_imm: %f  time_dt: %f\n",
//                                 cm, ci, pcm, pci, conc_avg, alpha, por_m, por_imm, time_->dt());
        		//getchar();

        		conc[MOBILE][sbi][elm_pos] = cm;
        		conc[IMMOBILE][sbi][elm_pos] = ci;
    		}

}
//=============================================================================
//      TRANSPORT SORPTION
//=============================================================================
void ConvectionTransport::transport_sorption( int elm_pos, ElementFullIter elem, int sbi) {

    double conc_avg = 0.0;
    double conc_avg_imm = 0.0;
    double n, Nm, Nimm;
    //int id;
    double phi = data_.phi.value(elem->centre(), elem->element_accessor());
    double por_m = data_.por_m.value(elem->centre(), elem->element_accessor());
    double por_imm = data_.por_imm.value(elem->centre(), elem->element_accessor());
    arma::Col<unsigned int> sorp_type = data_.sorp_type.value(elem->centre(), elem->element_accessor());
    arma::vec sorp_coef0 = data_.sorp_coef0.value(elem->centre(), elem->element_accessor());
    arma::vec sorp_coef1 = data_.sorp_coef1.value(elem->centre(), elem->element_accessor());

    n = 1 - (por_m + por_imm);
    Nm = por_m;
    Nimm = por_imm;

    conc_avg = conc[MOBILE][sbi][elm_pos] + conc[MOBILE_SORB][sbi][elm_pos] * n / Nm; // cela hmota do poru

    //cout << "input concentration for old sorption is " << conc[MOBILE][sbi][elm_pos] << endl;
    if (conc_avg != 0) {
        compute_sorption(conc_avg, sorp_coef0[sbi], sorp_coef1[sbi], sorp_type[sbi], &conc[MOBILE][sbi][elm_pos],
                &conc[MOBILE_SORB][sbi][elm_pos], Nm / n, n * phi / Nm);

    }
    //printf("\n%f\t%f\t",n * phi / Nm,n * phi / Nm);
    //printf("\n%f\t%f\t",n * phi / Nimm,n * (1 - phi) / Nimm);
    // getchar();

    if ((dual_porosity == true) && (por_imm != 0)) {
        conc_avg_imm = conc[IMMOBILE][sbi][elm_pos] + conc[IMMOBILE_SORB][sbi][elm_pos] * n / Nimm; // cela hmota do poru

        if (conc_avg_imm != 0) {
            compute_sorption(conc_avg_imm, sorp_coef0[sbi], sorp_coef1[sbi], sorp_type[sbi], &conc[IMMOBILE][sbi][elm_pos],
                    &conc[IMMOBILE_SORB][sbi][elm_pos], Nimm / n, n * (1 - phi) / Nimm);
        }
    }

}
//=============================================================================
//      COMPUTE SORPTION
//=============================================================================
void ConvectionTransport::compute_sorption(double conc_avg, double sorp_coef0, double sorp_coef1, unsigned int sorp_type, double *concx, double *concx_sorb, double Nv,
        double N) {
    double Kx = sorp_coef0 * N;
    double parameter;
    double NR, pNR, cz, tcz;

    double ad = 1e4;
    double tolerence = 1e-8;
    int i;

    pNR = 0;

    //if(conc_avg > 1e-20)
    switch (sorp_type) {
    case Isotherm::linear: //linear
        *concx = conc_avg / (1 + Kx);
        //    *concx_sorb = (conc_avg - *concx) * Nv;   // s = Kd *c  [kg/m^3]
        break;
    case Isotherm::freundlich: //freundlich
        parameter = sorp_coef1;
        cz = pow(ad / (Kx * parameter), 1 / (parameter - 1));
        tcz = ad / parameter;
        NR = cz;
        for (i = 0; i < 20; i++) //Newton Raphson iteration cycle
        {
            NR -= (NR + ((NR > cz) ? Kx * pow(NR, parameter) : tcz * NR) - conc_avg) / (1 + ((NR > cz) ? parameter * Kx * pow(NR,
                    parameter - 1) : tcz));
            if ((NR <= cz) || (fabs(NR - pNR) < tolerence * NR ) )
                break;
            pNR = NR;
            //                if(NR < 0) printf("\n%f\t",NR);
        }
        *concx = NR;
        break;
    case Isotherm::langmuir: // langmuir
        parameter = sorp_coef1;
        NR = 0;
        //Kx = sorp_coef0/N;
        for (i = 0; i < 5; i++) //Newton Raphson iteration cycle
        {
            //NR -= (NR + (NR * Kx * parameter) / (1 + NR * Kx) - conc_avg) / (1 + Kx * parameter / pow(1 + NR * Kx, 2));
            NR -= (NR + (N * NR * parameter * sorp_coef0)/( 1 + NR * sorp_coef0 ) - conc_avg)/(1 + N * sorp_coef0 * parameter/pow((1 + NR * sorp_coef0), 2));
            if (fabs(NR - pNR) < tolerence *NR)
                break;
            pNR = NR;
        }
        *concx = NR;
        //   *concx_sorb = (conc_avg - *concx) * Nv;
        break;
    }

    *concx_sorb = (conc_avg - *concx) * Nv;
}
//=============================================================================
//      TIME STEP (RECOMPUTE)
//=============================================================================
//void ConvectionTransport::compute_time_step() {
//    double problem_save_step = OptGetDbl("Global", "Save_step", "1.0");
//    double problem_stop_time = OptGetDbl("Global", "Stop_time", "1.0");
//    save_step = (int) ceil(problem_save_step / time_step); // transport  rev
//    time_step = problem_save_step / save_step;
//    steps = save_step * (int) floor(problem_stop_time / problem_save_step);
//
//}
//=============================================================================
//      TRANSPORT ONE STEP
//=============================================================================

#if 0
//=============================================================================
//      TRANSPORT UNTIL TIME
//=============================================================================
void ConvectionTransport::transport_until_time(double time_interval) {
    	int step = 0;
    	register int t;
    	/*Distribution *distribution;
    	int *el_4_loc;

    	// Chemistry initialization
    	Linear_reaction *decayRad = new Linear_reaction(time_step,mesh_,n_substances,dual_porosity);
    	this->get_par_info(el_4_loc, distribution);
    	decayRad->set_concentration_matrix(pconc, distribution, el_4_loc);
    	Semchem_interface *Semchem_reactions = new Semchem_interface(time_step,mesh_,n_substances,dual_porosity);
    	Semchem_reactions->set_el_4_loc(el_4_loc);
    	Semchem_reactions->set_concentration_matrix(pconc, distribution, el_4_loc);*/

	    for (t = 1; t <= steps; t++) {
	    	time += time_step;
	     //   SET_TIMER_SUBFRAMES("TRANSPORT",t);  // should be in destructor as soon as we have class iteration counter
		START_TIMER("transport_step");
	    	compute_one_step();
		END_TIMER("transport_step");

		     /*/ Calling linear reactions and Semchem together
		    	  for (int loc_el = 0; loc_el < el_ds->lsize(); loc_el++) {
		    		 START_TIMER("decay_step");
		    	   	 (*decayRad).compute_reaction(pconc[MOBILE], loc_el);
		    	     if (dual_porosity == true) {
		    	    	(*decayRad).compute_reaction(pconc[IMMOBILE], loc_el);
		    	     }
		    	     END_TIMER("decay_step");
		    	     START_TIMER("semchem_step");
		    	     if(Semchem_reactions->semchem_on == true){
		    	    	 Semchem_reactions->set_timestep(time_step);
		    	    	 //Semchem_reactions->compute_reaction(dual_porosity, mesh_->element(el_4_loc[loc_el]), loc_el, pconc[MOBILE], pconc[IMMOBILE]);
		    	    	 Semchem_reactions->compute_reaction(dual_porosity, mesh_->element(el_4_loc[loc_el]), loc_el, pconc);
		    	     }
		    	     END_TIMER("semchem_step");
		    	  }*/

	        step++;
	        //&& ((ConstantDB::getInstance()->getInt("Problem_type") != PROBLEM_DENSITY)
	        if ((save_step == step) || (write_iterations)) {
	            xprintf( Msg, "Output\n");
                    output_vector_gather();

                    // Register concentrations data on elements
                    for(int subst_id=0; subst_id<n_subst_; subst_id++) {
                        output_time->register_elem_data(subst_names_[subst_id], "", out_conc[MOBILE][subst_id], mesh_->n_elements());
                    }
                    output_time->write_data(time);
                //  if (ConstantDB::getInstance()->getInt("Problem_type") != STEADY_SATURATED)
                // output_time(problem, t * time_step); // time variable flow field
                step = 0;
            }
	    }
}

#endif

//=============================================================================
//      OUTPUT VECTOR GATHER
//=============================================================================
void ConvectionTransport::output_vector_gather() {

    unsigned int sbi/*, rank, np*/;
    IS is;
    //PetscViewer inviewer;

    //	MPI_Barrier(PETSC_COMM_WORLD);
/*    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);*/


    //ISCreateStride(PETSC_COMM_SELF,mesh_->n_elements(),0,1,&is);
    ISCreateGeneral(PETSC_COMM_SELF, mesh_->n_elements(), row_4_el, PETSC_COPY_VALUES, &is); //WithArray
    VecScatterCreate(vconc[0], is, vconc_out[0], PETSC_NULL, &vconc_out_scatter);
    for (sbi = 0; sbi < n_subst_; sbi++) {
        VecScatterBegin(vconc_out_scatter, vconc[sbi], vconc_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(vconc_out_scatter, vconc[sbi], vconc_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
    }
    //VecView(transport->vconc[0],PETSC_VIEWER_STDOUT_WORLD);
    //VecView(transport->vconc_out[0],PETSC_VIEWER_STDOUT_WORLD);
    VecScatterDestroy(&(vconc_out_scatter));
    ISDestroy(&(is));
}


void ConvectionTransport::get_parallel_solution_vector(Vec &vc){
	return;
};

void ConvectionTransport::get_solution_vector(double* &vector, unsigned int &size){
	return;
};


double ***ConvectionTransport::get_concentration_matrix() {
	return conc;
}

void ConvectionTransport::get_par_info(int * &el_4_loc_out, Distribution * &el_distribution_out){
	el_4_loc_out = this->el_4_loc;
	el_distribution_out = this->el_ds;
	return;
}

bool ConvectionTransport::get_dual_porosity(){
	return dual_porosity;
}

int *ConvectionTransport::get_el_4_loc(){
	return el_4_loc;
}

int *ConvectionTransport::get_row_4_el(){
	return row_4_el;
}

/*
int ConvectionTransport::get_n_substances() {
	return n_subst_;
}
*/


void ConvectionTransport::calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance)
{
    double mass_flux[n_substances()];
    double *pconc[n_substances()];

    for (int sbi=0; sbi<n_substances(); sbi++)
    	VecGetArray(vpconc[sbi], &pconc[sbi]);

    FOR_BOUNDARIES(mesh_, bcd) {

        // !! there can be more sides per one boundary
        int index = row_4_el[bcd->side()->element().index()];
        if (!el_ds->is_local(index)) continue;
        int loc_index = index-el_ds->begin();

        double por_m = data_.por_m.value(bcd->side()->element()->centre(), bcd->side()->element()->element_accessor() );
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
        	double por_m = data_.por_m.value(elem->centre(), ele_acc);
        	double csection = data_.cross_section->value(elem->centre(), ele_acc);
        	int loc_index = index - el_ds->begin();
			double sum_sol_phases = 0;

			for (int ph=0; ph<MAX_PHASES; ph++)
			{
				if ((sub_problem & ph) == ph)
					sum_sol_phases += conc[ph][sbi][loc_index];
			}

			mass[sbi][ele_acc.region().bulk_idx()] += por_m*csection*sum_sol_phases*elem->measure();
			src_balance[sbi][ele_acc.region().bulk_idx()] += sources[loc_index]*elem->measure();
        }
    }
}



void ConvectionTransport::output_data() {

    if (time_->is_current(output_mark_type)) {

        DBGMSG("\nTOS: output time: %f\n", time_->t());
        output_vector_gather();

        int rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        if (rank == 0)
        {
        	// Register fresh output data
			//Input::Record output_rec = this->in_rec_->val<Input::Record>("output");
			//OutputTime::register_data<3, FieldValue<3>::Scalar>(output_rec, OutputTime::ELEM_DATA, &data_.conc_mobile);
			data_.output_fields.set_time(*time_);
			data_.output_fields.output(output_stream_);
        }

        if (mass_balance_)
        	mass_balance_->output(time_->t());

        //for synchronization when measuring time by Profiler
        MPI_Barrier(MPI_COMM_WORLD);
    }
}
