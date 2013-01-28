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
 * @todo
 * - podrobne komentare k funkcim
 * - extremne podrobny komentar ke strukture/objektu
 *   zejmena komentar k prerekvizitam
 *   popis jak transport funguje
 * - udelat z toho objekt
 * - oddelit output -> transport_output
 *   objekt pro transport by jen pocital koncentrace
 *
 */
//#include <cmath>
//#include "system/system.hh"
//#include "xio.h"


#include "system/system.hh"

#include "system/math_fce.h"
#include "mesh/mesh.h"
#include "transport/transport.h"
//#include "transport/sources.hh"
#include "io/output.h"
#include "materials.hh"
#include "io/read_ini.h"
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

#include "coupling/time_governor.hh"


ConvectionTransport::ConvectionTransport(Mesh &init_mesh, MaterialDatabase &material_database, const Input::Record &in_rec)
: EquationBase(init_mesh, material_database, in_rec)
{
    F_ENTRY;

    // [Density]
    /*
    max_dens_it = OptGetInt("Density", "Density_max_iter", "20");
    dens_implicit = OptGetBool("Density", "Density_implicit", "no");
    dens_eps = OptGetDbl("Density", "Eps_iter", "1.0e-5");
    write_iterations = OptGetBool("Density", "Write_iterations", "no");
    dens_step = OptGetInt("Density", "Density_steps", "1");
    */

    //double problem_stop_time = OptGetDbl("Global", "Stop_time", "1.0");

    //mark type of the equation of convection transport (created in EquationBase constructor) and it is fixed
    target_mark_type = this->mark_type() | TimeGovernor::marks().type_fixed_time();
    time_ = new TimeGovernor(in_rec.val<Input::Record>("time"), target_mark_type);
    

    sorption = in_rec.val<bool>("sorption_enable");
    dual_porosity = in_rec.val<bool>("dual_porosity");
    // reaction_on = in_rec.val<bool>("transport_reactions");


    in_rec.val<Input::Array>("substances").copy_to(substance_name);
    n_substances = substance_name.size();
    INPUT_CHECK(n_substances >= 1 ,"Number of substances must be positive.\n");

    pepa=false; reaction_on = false;
    // pepa = OptGetBool("Transport", "Decay", "no"); //PEPA
    // type = OptGetInt("Transport", "Decay_type", "-1"); //PEPA

    sub_problem = 0;
    if (dual_porosity == true)
        sub_problem += 1;
    if (sorption == true)
        sub_problem += 2;

    mat_base->read_transport_materials(dual_porosity, sorption,n_substances);
    make_transport_partitioning();
    alloc_transport_vectors();
    alloc_transport_structs_mpi();
    read_initial_condition(in_rec.val<FilePath>("initial_file"));



    // read times of time dependent boundary condition and check the input files
    
    // bc marks do not influent choice of time step (not fixed_time type)
    bc_mark_type_ = time_->marks().type_bc_change() | equation_mark_type_;
    std::vector<double> bc_times;


    FILE * f;
    string fname;
    bc_fname = in_rec.val<FilePath>("boundary_file");
    DBGMSG("bc file: %s\n", bc_fname.c_str());
    Input::Iterator<Input::Array> it = in_rec.find<Input::Array>("bc_times");
    if ( ! it ) {
        // only one boundary condition, check filename and read bc condition
        bc_time_level = -1;
        fname = make_bc_file_name(-1);
        DBGMSG("bc file: %s\n ", fname.c_str());
        if ( !(f = xfopen(fname.c_str(), "rt")) ) {
            xprintf(UsrErr,"Missing file: %s", fname.c_str());
        }
        xfclose(f);
        read_bc_vector();
    } else {
        it->copy_to(bc_times);

        bc_time_level = 0;
        // set initial bc to zero
        for(unsigned int sbi=0; sbi < n_substances; ++sbi) VecZeroEntries(bcv[sbi]);
        // check files for bc time levels and set TimeMarks
        for(unsigned int i=0;i<bc_times.size();++i) {
            fname = make_bc_file_name(i);
            DBGMSG("bc file_: %s\n ", fname.c_str());
            if ( !(f = xfopen(fname.c_str(), "rt")) ) {
                xprintf(UsrErr,"Missing file: %s", fname.c_str());
            }
            xfclose(f);

            time_->marks().add(TimeMark(bc_times[i],bc_mark_type_));
        }
    }


    is_convection_matrix_scaled = false;

    /*
    if (ConstantDB::getInstance()->getInt("Problem_type") == PROBLEM_DENSITY) {
        sscales = OptGetStr("Transport", "Substances_density_scales", "1.0");
        subst_scales(sscales);
    }*/

   // n_elements = mesh_->n_elements();

/*
    reaction = NULL;
    n_reaction = 0;
*/

    Input::Iterator<FilePath> sources_it = in_rec.find<FilePath>("sources_file");

    if (sources_it) {
        transportsources = new TransportSources( n_substances, *el_ds );
        transportsources->read_concentration_sources( string(* sources_it), row_4_el, mesh_);
    } else {
        transportsources = NULL;
    }



    output_vector_gather();


}


//=============================================================================
// MAKE TRANSPORT
//=============================================================================
void ConvectionTransport::make_transport_partitioning() {

    F_ENTRY;

    int rank, np, i, j, k, row_MH, a;
    //struct DarcyFlowMH *water=transport->problem->water;

    SparseGraph *ele_graph = new SparseGraphMETIS(mesh_->n_elements()); // graph for partitioning
    Distribution init_ele_ds = ele_graph->get_distr(); // initial distr.
    int *loc_part = new int[init_ele_ds.lsize()]; // partitionig in initial distribution

    make_element_connection_graph(mesh_, ele_graph, true);
    WARN_ASSERT(ele_graph->is_symmetric(),"Attention graph for partitioning is not symmetric!\n");

    ele_graph->partition(loc_part);

    delete ele_graph;

    int *id_4_old = (int *) xmalloc(mesh_->n_elements() * sizeof(int));
    i = 0;
    FOR_ELEMENTS(mesh_, ele)
        id_4_old[i] = i, i++;
    id_maps(mesh_->n_elements(), id_4_old, init_ele_ds, (int *) loc_part, el_ds, el_4_loc, row_4_el);

    delete[] loc_part;
    xfree(id_4_old);

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
// TODO
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

void ConvectionTransport::set_flow_field_vector(const MH_DofHandler &dh){
    // DBGMSG("set_flow_fieldvec\n");
    mh_dh = &dh;
	create_transport_matrix_mpi();
};

double ***ConvectionTransport::get_out_conc(){
	return out_conc;
}

vector<string> &ConvectionTransport::get_substance_names(){
	return substance_name;
}



//=============================================================================
//
//=============================================================================
void ConvectionTransport::subst_scales(char *line) {
    int sbi;

    ASSERT(!( (n_substances < 1) || (line == NULL) ),"Bad parameter of the function subst_scales()\n");
    substance_density_scale= (double*) xmalloc( n_substances* sizeof(double));
    for (sbi = 0; sbi < n_substances; sbi++)
    	substance_density_scale[sbi] = atof(strtok(sbi == 0 ? line : NULL, " \t,;"));
}
//=============================================================================
// READ INITIAL CONDITION
//=============================================================================
void ConvectionTransport::read_initial_condition(string concentration_fname) {
        F_ENTRY;

        FILE	*in;		  // input file
		char     line[ LINE_SIZE ]; // line of data file
		unsigned int sbi,index, id, eid, i,n_concentrations, global_idx;

		in = xfopen( concentration_fname, "rt" );

		skip_to( in, "$Concentrations" );
		xfgets( line, LINE_SIZE - 2, in );
		n_concentrations = atoi( xstrtok( line) );
		INPUT_CHECK(!( n_concentrations < 1 ),"Number of concentrations < 1 in function read_concentration_list()\n");
	    INPUT_CHECK(!( mesh_->n_elements() != n_concentrations),"Different number of elements and concentrations\n");


	    for (i = 0; i < n_concentrations; i++) {
	    	xfgets( line, LINE_SIZE - 2, in );
	    	// DBGMSG("%s\n",line);
	    	//printf("%s\t%d\n",line,n_concentrations);
	    	ASSERT(!(line == NULL),"NULL as argument of function parse_concentration_line()\n");
	    	eid    = atoi( xstrtok( line) );
	    	global_idx =row_4_el[mesh_->element.find_id(eid).index()];
	    	if ( el_ds->is_local(global_idx) ) {
	    		index = global_idx - el_ds->begin();
	    		for( sbi = 0; sbi < n_substances; sbi++ ){
	    			conc[MOBILE][ sbi ][index] = atof( xstrtok( NULL) );
	    			pconc[MOBILE][ sbi ][index] = conc[MOBILE][ sbi ][index];
	    		}
	    	}
	    }

		xfclose( in );

}

//=============================================================================
//	ALLOCATE OF TRANSPORT VARIABLES (ELEMENT & NODES)
//=============================================================================
void ConvectionTransport::alloc_transport_vectors() {

    int i, j, sbi, n_subst, ph;
    ElementIter elm;
    n_subst = n_substances;


   // printf("%d\t\n",n_substances);
   // getchar();

     cumulative_corr = (double**) xmalloc(n_subst * sizeof(double*));
     for (sbi = 0; sbi < n_subst; sbi++)
         cumulative_corr[sbi] = (double*) xmalloc(el_ds->lsize() * sizeof(double));
 /*
     for (sbi = 0; sbi < n_subst; sbi++)
         for (i = 0; i < el_ds->size(); i++)
             cumulative_corr[sbi][i] = 0.0;
   */
    conc = (double***) xmalloc(MAX_PHASES * sizeof(double**));
    pconc = (double***) xmalloc(MAX_PHASES * sizeof(double**));
    out_conc = (double***) xmalloc(MAX_PHASES * sizeof(double**));
    //transport->node_conc = (double****) xmalloc(MAX_PHASES * sizeof(double***));
    for (ph = 0; ph < MAX_PHASES; ph++) {
        if ((sub_problem & ph) == ph) {
            conc[ph] = (double**) xmalloc(n_subst * sizeof(double*)); //(MAX_PHASES * sizeof(double*));
            pconc[ph] = (double**) xmalloc(n_subst * sizeof(double*));
            out_conc[ph] = (double**) xmalloc(n_subst * sizeof(double*));
            //  transport->node_conc[sbi] = (double***) xmalloc(MAX_PHASES * sizeof(double**));
            //}
            //}
            for (sbi = 0; sbi < n_subst; sbi++) {
                conc[ph][sbi] = (double*) xmalloc(el_ds->lsize() * sizeof(double));
                pconc[ph][sbi] = (double*) xmalloc(el_ds->lsize() * sizeof(double));
                out_conc[ph][sbi] = (double*) xmalloc(el_ds->size() * sizeof(double));
                // transport->node_conc[sbi][ph] = (double**)xmalloc((mesh__->n_elements() ) * sizeof(double*));
                for (i = 0; i < el_ds->lsize(); i++) {
                    conc[ph][sbi][i] = 0.0;
                    pconc[ph][sbi][i] = 0.0;

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
            pconc[ph] = NULL;
            out_conc[ph] = NULL;
            //transport->node_conc[sbi][ph] = NULL;
        }
    }
}
//=============================================================================
//	ALLOCATE OF TRANSPORT (DENSITY VECTORS)
//=============================================================================
void ConvectionTransport::alloc_density_vectors() {

    int ph, sbi, i, sub;
    int n_subst = n_substances;
    int n_elements = mesh_->n_elements();

    sub = sub_problem;

    scalar_it = (double*) xmalloc(n_elements * sizeof(double)); // Zatim nevyuzito
    prev_conc = (double***) xmalloc(MAX_PHASES * sizeof(double**)); //transport->prev_conc = (double***) xmalloc(n_subst * sizeof(double**));

    for (ph = 0; ph < MAX_PHASES; ph++) {
     if ((sub & ph) == ph) {        
      prev_conc[ph] = (double**) xmalloc(n_subst * sizeof(double*)); //transport->prev_conc[sbi] = (double**) xmalloc(MAX_PHASES * sizeof(double*));

        for (sbi = 0; sbi < n_elements; sbi++)
                prev_conc[ph][sbi] = (double*) xmalloc(n_elements * sizeof(double));

                for (i = 0; i < n_elements; i++)
                    prev_conc[ph][sbi][i] = 0.0;
    } else  prev_conc[ph] = NULL;
    }
}
//=============================================================================
//	ALLOCATION OF TRANSPORT VECTORS (MPI)
//=============================================================================
void ConvectionTransport::alloc_transport_structs_mpi() {

    int i, j, sbi, n_subst, ph, ierr, rank, np;
    ElementIter elm;
    n_subst = n_substances;

    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);

    bcv = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    bcvcorr = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vconc = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vpconc = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
/*
    vsources_density = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vsources_sigma = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vsources_conc = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
  //  vsources_corr = (Vec*) xmalloc((n_subst+1) * (sizeof(Vec)));
    vsources_corr = (Vec*) xmalloc(n_subst * (sizeof(Vec))); */
    vcumulative_corr = (Vec*) xmalloc(n_subst * (sizeof(Vec)));


    // if( rank == 0)
    vconc_out = (Vec*) xmalloc(n_subst * (sizeof(Vec))); // extend to all

    // TODO: should be replaced by Distribution(Block) or better remove whole boundary matrix with these vectors

    lb_col = (int*) xmalloc(np * sizeof(int)); // local rows in BM
    lb_col[0] = 0;
    for (i = 0; i < (np -1); i++) {
        lb_col[i] = (int)floor(mesh_->n_boundaries() / np);
    }
    lb_col[np - 1] = mesh_->n_boundaries() - (np - 1) * lb_col[0];

    /*
     if((transport->sub_problem & 1) == 1)
     transport->vconc_im=(Vec*)xmalloc(n_subst * (sizeof (Vec)));
     else
     transport->vconc_im = NULL;

     if((transport->sub_problem & 2) == 2)
     transport->vconc_so=(Vec*)xmalloc(n_subst * (sizeof (Vec)));
     else
     transport->vconc_so = NULL;

     if((transport->sub_problem & 3) == 3)
     transport->vconc_im_so=(Vec*)xmalloc(n_subst * (sizeof (Vec)));
     else
     transport->vconc_im_so = NULL;

     */


    for (sbi = 0; sbi < n_subst; sbi++) {
        //ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->lb_col[rank],mesh__->n_boundaries(),&transport->bcv[sbi]);
        //ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh__->n_elements(),&transport->bcvcorr[sbi]);
        //ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh__->n_elements(),&transport->vconc[sbi]);
        //ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh__->n_elements(),&transport->vpconc[sbi]);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, lb_col[rank], mesh_->n_boundaries(), &bcv[sbi]);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(), &bcvcorr[sbi]);
        VecZeroEntries(bcv[sbi]);
        VecZeroEntries(bcvcorr[sbi]);
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(), conc[MOBILE][sbi],
                &vconc[sbi]);

        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(),
                pconc[MOBILE][sbi], &vpconc[sbi]);
        VecZeroEntries(vconc[sbi]);
        VecZeroEntries(vpconc[sbi]);

        // SOURCES
        /*
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(),
        		sources_density[sbi],&vsources_density[sbi]);
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(),
        		sources_sigma[sbi],&vsources_sigma[sbi]);
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(),
        		sources_conc[sbi],&vsources_conc[sbi]);
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(),
        		sources_corr[sbi],&vsources_corr[sbi]); */
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(),
        		cumulative_corr[sbi],&vcumulative_corr[sbi]);

        //  if(rank == 0)
        ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, mesh_->n_elements(), out_conc[MOBILE][sbi], &vconc_out[sbi]);

        VecZeroEntries(vcumulative_corr[sbi]);
        VecZeroEntries(vconc_out[sbi]);
        //ierr = VecCreateMPI(PETSC_COMM_SELF ,transport->mesh_->n_elements(),transport->mesh_->n_elements(),&transport->vconc_out[sbi]); /*xx*/
        /*
         if(transport->vconc_im != NULL)
         ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh_->n_elements(),&transport->vconc_im[sbi]);
         if(transport->vconc_so != NULL)
         ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh_->n_elements(),&transport->vconc_so[sbi]);
         if(transport->vconc_im_so != NULL)
         ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh_->n_elements(),&transport->vconc_im_so[sbi]);
         */
    }


  //  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, el_ds->lsize(), mesh_->n_elements(),	/**/
  //          		sources_corr[n_substances],&vsources_corr[n_substances]);				/**/


    //
    //ierr=MatCreateMPIAIJ(PETSC_COMM_WORLD,transport->l_row[rank],transport->l_row[rank],mesh_->n_elements(),mesh_->n_elements(),
    //		PETSC_NULL,transport->d_row[rank],PETSC_NULL,transport->od_row[rank],&transport->tm);

    //ierr=MatCreateMPIAIJ(PETSC_COMM_WORLD,transport->l_row[rank],transport->lb_col[rank],mesh_->n_elements(),mesh_->n_boundaries(),
    //				PETSC_NULL,transport->db_row[rank],PETSC_NULL,transport->odb_row[rank],&transport->bcm);

    ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD, el_ds->lsize(), el_ds->lsize(), mesh_->n_elements(),
            mesh_->n_elements(), 8, PETSC_NULL, 1, PETSC_NULL, &tm);

    ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD, el_ds->lsize(), lb_col[rank], mesh_->n_elements(),
            mesh_->n_boundaries(), 2, PETSC_NULL, 0, PETSC_NULL, &bcm);

    /*
     MatCreateMPIAIJ(comm,m,n,M,N,0,PETSC_NULL,0,PETSC_NULL,&transport->tm);
     MatCreateMPIAIJ(comm,m,n,M,N,0,PETSC_NULL,0,PETSC_NULL,&transport->bcm);

     ierr=MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,transport->l_row[rank],mesh_->n_elements(),mesh_->n_elements(),
     PETSC_NULL,transport->d_row[rank],PETSC_NULL,transport->od_row[rank],&transport->tm);

     ierr=MatCreateMPIAIJ(PETSC_COMM_WORLD,transport->l_row[rank],transport->lb_col[rank],mesh_->n_elements(),mesh_->n_boundaries(),
     PETSC_NULL,transport->db_row[rank],PETSC_NULL,transport->odb_row[rank],&transport->bcm);
     */


}

std::string ConvectionTransport::make_bc_file_name(int level) {

    if (level >= 0 ) {
        std::stringstream name_str;
        name_str << bc_fname << "_" << setfill('0') << setw(3) << level;
        return name_str.str();
    } else {
        return bc_fname;
    }
}
void ConvectionTransport::read_bc_vector() {

    int rank, sbi;


    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0) {
        xprintf(Msg, "CONVECTION: Reading BC (level = %d)\n",bc_time_level);

        int bcd_id, boundary_id, boundary_index;
        double bcd_conc;
        char line[LINE_SIZE]; // line of data file

        // make bc filename

        FILE *in = xfopen(make_bc_file_name(bc_time_level).c_str(), "rt");
        skip_to(in, "$Transport_BCD");
        xfgets(line, LINE_SIZE - 2, in);
        int n_bcd = atoi(xstrtok(line));
        for (int i_bcd = 0; i_bcd < n_bcd; i_bcd++) {
            xfgets(line, LINE_SIZE - 2, in);
            bcd_id = atoi(xstrtok(line)); // scratch transport bcd id
            boundary_id = atoi(xstrtok(NULL));
            //        DBGMSG("transp b. id: %d\n",boundary_id);
            //if (OldBcdInput::instance()->bcd_ids_[i_bcd] != boundary_id) xprintf(UsrErr, "Incorrect flow ID in Transport BCD file: %d\n", boundary_id);
            boundary_index = *(Boundary::id_to_bcd.find_id(bcd_id));
            INPUT_CHECK(boundary_index >= 0,"Wrong boundary index %d for bcd id %d in transport bcd file!", boundary_id, bcd_id);
            for (sbi = 0; sbi < n_substances; sbi++) {
                bcd_conc = atof(xstrtok(NULL));
                DBGMSG(" b idx: %d sunst: %d conc: %f\n", boundary_index, sbi, bcd_conc);
                VecSetValue(bcv[sbi], boundary_index, bcd_conc, INSERT_VALUES);
            }
        }
        xfclose(in);
    }
    for (sbi = 0; sbi < n_substances; sbi++)
        VecAssemblyBegin(bcv[sbi]);
    //for (sbi = 0; sbi < n_substances; sbi++)
    //    VecZeroEntries(bcvcorr[sbi]);
    for (sbi = 0; sbi < n_substances; sbi++)
        VecAssemblyEnd(bcv[sbi]);

    // update source vectors
    // TODO: rather use Lazy dependency
    if (bc_time_level != -1)
        for (unsigned int sbi = 0; sbi < n_substances; sbi++) {
            MatMult(bcm, bcv[sbi], bcvcorr[sbi]);
        }

    // VecView(bcvcorr[0],PETSC_VIEWER_STDOUT_SELF);
    // getchar();
    bc_time_level++;

}

void ConvectionTransport::compute_one_step() {

    MaterialDatabase::Iter material;
    int sbi;

    // possibly read boundary conditions

    if (bc_time_level != -1 && time_->is_current(bc_mark_type_))
        read_bc_vector();

    // proceed to actually computed time
    //time_->view("CONVECTION");
    time_->next_time();

    for (sbi = 0; sbi < n_substances; sbi++) {
        // one step in MOBILE phase
        if (transportsources != NULL) {
            //DBGMSG("component: %d\n", sbi);

            //if (vcumulative_corr[sbi][10] >0) { int i =1;}
            //if (bcvcorr[sbi][10] >0) { int i =1;}
            //if (conc[sbi][10] >0) { int i =1;}
            VecAXPBYPCZ(vcumulative_corr[sbi], 1.0, time_->dt(), 0.0, bcvcorr[sbi],
                    transportsources->compute_concentration_sources(sbi, conc[MOBILE][sbi] )
                    );
        } else {
            VecCopy(bcvcorr[sbi], vcumulative_corr[sbi]);
        }

        //VecView(vpconc[sbi],PETSC_VIEWER_STDOUT_SELF);

        MatMultAdd(tm, vpconc[sbi], vcumulative_corr[sbi], vconc[sbi]); // conc=tm*pconc + bc
        //VecView(vconc[sbi],PETSC_VIEWER_STDOUT_SELF);

        VecCopy(vconc[sbi], vpconc[sbi]); // pconc = conc

        if ((dual_porosity == true) || (sorption == true) || (pepa == true) || (reaction_on == true))
            // cycle over local elements only in any order
            for (int loc_el = 0; loc_el < el_ds->lsize(); loc_el++) {
                material = (mesh_->element(el_4_loc[loc_el]))->material;

                if (dual_porosity == true)
                    transport_dual_porosity(loc_el, material, sbi);
                if (sorption == true)
                    transport_sorption(loc_el, material, sbi);

                // if (reaction_on == true)
                //    transport_reaction(trans, loc_el, material, sbi);


            }
        // transport_node_conc(mesh_,sbi,problem->transport_sub_problem);  // vyresit prepocet
    }
}


void ConvectionTransport::set_target_time(double target_time)
{
    //sets target_mark_type (it is fixed) to be met in next_time()
    time_->marks().add(TimeMark(target_time, target_mark_type));
    
    //returns integer, one can check here whether the constraint has been set or not
    time_->set_upper_constraint(cfl_max_step);
    
    //fixing convection time governor till next target_mark_type (got from TOS or other)
    time_->fix_dt_until_mark();
    
    //time_->view("CONVECTION");    //show convection time governor

    if ( is_convection_matrix_scaled ) {
        // rescale matrix
        MatScale(bcm, time_->dt()/time_->estimate_dt());
        MatShift(tm, -1.0);
        MatScale(tm, time_->dt()/time_->estimate_dt() );
        MatShift(tm, 1.0);
    } else {
        // scale fresh convection term matrix
        MatScale(bcm, time_->estimate_dt());
        MatScale(tm, time_->estimate_dt());
        MatShift(tm, 1.0);
    }

    // update source vectors
    for (unsigned int sbi = 0; sbi < n_substances; sbi++) {
            MatMult(bcm, bcv[sbi], bcvcorr[sbi]);
//            VecView(bcv[sbi],PETSC_VIEWER_STDOUT_SELF);
//            getchar();
//            VecView(bcvcorr[sbi],PETSC_VIEWER_STDOUT_SELF);
//           getchar();
    }

    is_convection_matrix_scaled = true;
}


//=============================================================================
// CREATE TRANSPORT MATRIX
//=============================================================================
void ConvectionTransport::create_transport_matrix_mpi() {

    START_TIMER("transport_matrix_assembly");

    ElementFullIter el2 = ELEMENT_FULL_ITER_NULL(mesh_);
    ElementFullIter elm = ELEMENT_FULL_ITER_NULL(mesh_);
    struct Edge *edg;
    struct Neighbour *ngh;
    //struct Transport *transport;
    int n, s, i, j, np, rank, new_j, new_i;
    double max_sum, aij, aii, *solution;
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
    MatZeroEntries(bcm);


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
    START_TIMER("matrix_assembly_mpi");

    for (int loc_el = 0; loc_el < el_ds->lsize(); loc_el++) {
        elm = mesh_->element(el_4_loc[loc_el]);
        new_i = row_4_el[elm.index()];

        FOR_ELEMENT_SIDES(elm,si) {
            // same dim
            flux = mh_dh->side_flux( *(elm->side(si)) );
            if (elm->side(si)->cond() == NULL) {
                if (flux < 0.0) {
                    edg = elm->side(si)->edge();
                    edg_flux = edge_flow[ elm->side(si)->edge_idx() ];
                    if ( edg_flux > ZERO)
                        FOR_EDGE_SIDES(edg,s)
                            // this test should also eliminate sides facing to lower dim. elements in comp. neighboring
                            // These edges on these sides should have just one side
                            if (edg->side(s) != elm->side(si)) {
                                flux2 = mh_dh->side_flux( *(edg->side(s)));
                                if ( flux2 > 0.0 ) {
                                    aij = -(flux * flux2 / ( edg_flux * elm->volume()
                                           * elm->material->por_m));
                                    j = ELEMENT_FULL_ITER(mesh_, edg->side(s)->element()).index();
                                    new_j = row_4_el[j];
                                    MatSetValue(tm, new_i, new_j, aij, INSERT_VALUES);
                                }
                            }
                }
                if (flux > 0.0)
                    aii -= (flux / (elm->volume() * elm->material->por_m));
            } else {
                if (flux < 0.0) {
                    aij = -(flux / (elm->volume() * elm->material->por_m));
                    j = elm->side(si)->cond_idx() ;
                    DBGMSG("BCM, i: %d j:%d\n", new_i , j);
                    MatSetValue(bcm, new_i, j, aij, INSERT_VALUES);
                    // vyresit BC matrix !!!!
                    //   printf("side in elm:%d value:%f\n ",elm->id,svector->val[j-1]);
                    //   printf("%d\t%d\n",elm->id,id2pos(problem,elm->side(si)->id,problem->spos_id,BC));

                }
                if (flux > 0.0)
                    aii -= (flux / (elm->volume() * elm->material->por_m));
            }
        }  // end same dim     //ELEMENT_SIDES

        FOR_ELM_NEIGHS_VB(elm,n) // comp model
            //FOR_NEIGH_ELEMENTS(elm->neigh_vb[n],s)
            {
                el2 = ELEMENT_FULL_ITER(mesh_, elm->neigh_vb[n]->side()->element() ); // higher dim. el.
                ASSERT( el2 != elm, "Elm. same\n");
                //if (elm.id() != el2.id()) {
                    flux = mh_dh->side_flux( *(elm->neigh_vb[n]->side()) );
                    if (flux > 0.0) {
                        // volume source - out flow from higher dimension
                        aij = flux / (elm->volume() * elm->material->por_m);
                        j = el2.index();
                        new_j = row_4_el[j];
                        MatSetValue(tm, new_i, new_j, aij, INSERT_VALUES);
                        // out flow from higher dim. already accounted
                    }
                    if (flux < 0.0) {
                        // volume drain - in flow to higher dimension
                        // in flow to higher dim.
                        aij = -(flux / (el2->volume() * el2->material->por_m));
                        new_j = row_4_el[el2.index()];
                        MatSetValue(tm, new_j, new_i, aij, INSERT_VALUES);

                        // diagonal drain
                        aii += flux / (elm->volume() * elm->material->por_m);

                    }

                //} // end comp model
            }
	/*
        FOR_ELM_NEIGHS_VV(elm,n) { //non-comp model
            ngh = elm->neigh_vv[n];
            FOR_NEIGH_ELEMENTS(ngh,s) {

                el2 = ELEMENT_FULL_ITER(mesh_, ngh->element[s]);
                if (elm.id() != el2.id()) {
                    flux = ngh->sigma * ngh->geom_factor * (el2->scalar - elm->scalar);
                    if (flux > 0.0) {
                        aij = flux / (elm->volume() * elm->material->por_m); // -=
                        j = el2.index();
                        new_j = row_4_el[j];
                        MatSetValue(tm, new_i, new_j, aij, INSERT_VALUES);
                    }
                    if (flux < 0.0)
                        aii += flux / (elm->volume() * elm->material->por_m);
                }
            }
        } // end non-comp model
      */
        MatSetValue(tm, new_i, new_i, aii, INSERT_VALUES);

        if (fabs(aii) > max_sum)
            max_sum = fabs(aii);
        aii = 0.0;
        //   i++;
    } // END ELEMENTS

    double glob_max_sum;

    MPI_Allreduce(&max_sum,&glob_max_sum,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);
    cfl_max_step = 1 / glob_max_sum;
    //time_step = 0.9 / glob_max_sum;
    
    DBGMSG("start assembly\n");
    MatAssemblyBegin(tm, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(bcm, MAT_FINAL_ASSEMBLY);
    

    MatAssemblyEnd(tm, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(bcm, MAT_FINAL_ASSEMBLY);
    DBGMSG("end assembly\n");


    // MPI_Barrier(PETSC_COMM_WORLD);
    /*
     MatView(transport->tm,PETSC_VIEWER_STDOUT_SELF);
     getchar();
     */
    is_convection_matrix_scaled = false;
    END_TIMER("transport_matrix_assembly");
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
void ConvectionTransport::transport_dual_porosity( int elm_pos, MaterialDatabase::Iter material, int sbi) {

    double conc_avg = 0.0;
    int id;
    //double ***conc = transport->conc;
    //double ***pconc = transport->pconc;
    double cm, pcm, ci, pci, por_m, por_imm, alpha;

    por_m = material->por_m;
    por_imm = material->por_imm;
    alpha = material->alpha[sbi];
    pcm = pconc[MOBILE][sbi][elm_pos];
    pci = pconc[IMMOBILE][sbi][elm_pos];
    // ---compute average concentration------------------------------------------
    conc_avg = ((por_m * pcm) + (por_imm * pci)) / (por_m + por_imm);

    if ((conc_avg != 0.0) && (por_imm != 0.0)) {
        // ---compute concentration in mobile area-----------------------------------
        cm = (pcm - conc_avg) * exp(-alpha * ((por_m + por_imm) / (por_m * por_imm)) * time_->dt()) + conc_avg;

        // ---compute concentration in immobile area---------------------------------
        ci = (pci - conc_avg) * exp(-alpha * ((por_m + por_imm) / (por_m * por_imm)) * time_->dt()) + conc_avg;
        // --------------------------------------------------------------------------
        //printf("\n%f\t%f\t%f",conc_avg,cm,ci);
        //getchar();

        conc[MOBILE][sbi][elm_pos] = cm;
        pconc[MOBILE][sbi][elm_pos] = cm;
        conc[IMMOBILE][sbi][elm_pos] = ci;
        pconc[IMMOBILE][sbi][elm_pos] = ci;
    }

    /*
     // ---compute average concentration------------------------------------------
     conc_avg = (( material->por_m * elm->pconc[sbi] )
     + (material->por_imm * elm->pconc_immobile[sbi] ))
     / ( material->por_m + material->por_imm );

     if((conc_avg != 0) && (material->por_imm != 0))
     {
     // ---compute concentration in mobile area-----------------------------------
     elm->conc[sbi] = ( elm->pconc[sbi] - conc_avg )
     * exp( - material->alpha[sbi] * ((material->por_m + material->por_imm)
     / (material->por_m * material->por_imm)) * transport->time_step )
     + conc_avg;

     // ---compute concentration in immobile area---------------------------------
     elm->conc_immobile[sbi] = ( elm->pconc_immobile[sbi] - conc_avg )
     * exp( - material->alpha[sbi] * ((material->por_m + material->por_imm)
     / (material->por_m * material->por_imm)) * transport->time_step )
     + conc_avg;
     // --------------------------------------------------------------------------
     //printf("\n%f\t%f\t%f",conc_avg,elm->conc[sbi],elm->conc_immobile[sbi]);
     //getchar();

     elm->pconc[sbi] = elm->conc[sbi];
     elm->pconc_immobile[sbi] = elm->conc_immobile[sbi];
     }
     */
}
//=============================================================================
//      TRANSPORT SORPTION
//=============================================================================
void ConvectionTransport::transport_sorption( int elm_pos, MaterialDatabase::Iter mtr, int sbi) {

    double conc_avg = 0.0;
    double conc_avg_imm = 0.0;
    double n, Nm, Nimm;
    int id;
    double phi;

    phi = mtr->phi;

    if ((mtr->sorp_coef[sbi].size() == 0) || (mtr->por_m == 1)) return;

    n = 1 - (mtr->por_m + mtr->por_imm);
    Nm = mtr->por_m;
    Nimm = mtr->por_imm;

    conc_avg = pconc[MOBILE][sbi][elm_pos] + pconc[MOBILE_SORB][sbi][elm_pos] * n / Nm; // cela hmota do poru


    if (conc_avg != 0) {
        compute_sorption(conc_avg, mtr->sorp_coef[sbi], mtr->sorp_type[sbi], &conc[MOBILE][sbi][elm_pos],
                &conc[MOBILE_SORB][sbi][elm_pos], Nm / n, n * phi / Nm);

        pconc[MOBILE][sbi][elm_pos] = conc[MOBILE][sbi][elm_pos];
        pconc[MOBILE_SORB][sbi][elm_pos] = conc[MOBILE_SORB][sbi][elm_pos];
    }
    //printf("\n%f\t%f\t",n * phi / Nm,n * phi / Nm);
    //printf("\n%f\t%f\t",n * phi / Nimm,n * (1 - phi) / Nimm);
    // getchar();

    if ((dual_porosity == true) && (mtr->por_imm != 0)) {
        conc_avg_imm = pconc[IMMOBILE][sbi][elm_pos] + pconc[IMMOBILE_SORB][sbi][elm_pos] * n / Nimm; // cela hmota do poru

        if (conc_avg_imm != 0) {
            compute_sorption(conc_avg_imm, mtr->sorp_coef[sbi], mtr->sorp_type[sbi], &conc[IMMOBILE][sbi][elm_pos],
                    &conc[IMMOBILE_SORB][sbi][elm_pos], Nimm / n, n * (1 - phi) / Nimm);

            pconc[IMMOBILE][sbi][elm_pos] = conc[IMMOBILE][sbi][elm_pos];
            pconc[IMMOBILE_SORB][sbi][elm_pos] = conc[IMMOBILE_SORB][sbi][elm_pos];
        }
    }

}
//=============================================================================
//      COMPUTE SORPTION
//=============================================================================
void ConvectionTransport::compute_sorption(double conc_avg, vector<double> &sorp_coef, int sorp_type, double *concx, double *concx_sorb, double Nv,
        double N) {
    double Kx = sorp_coef[0] * N;
    double parameter;// = sorp_coef[1];
    double NR, pNR, cz, tcz;
    //double lZero = 0.0000001;
    double ad = 1e4;
    int i;

    pNR = 0;

    //if(conc_avg > 1e-20)
    switch (sorp_type) {
    case 1: //linear
        *concx = conc_avg / (1 + Kx);
        //    *concx_sorb = (conc_avg - *concx) * Nv;   // s = Kd *c  [kg/m^3]
        break;
    case 2: //freundlich
        parameter = sorp_coef[1];
        cz = pow(ad / (Kx * parameter), 1 / (parameter - 1));
        tcz = ad / parameter;
        NR = cz;
        for (i = 0; i < 20; i++) //Newton Raphson iteration cycle
        {
            NR -= (NR + ((NR > cz) ? Kx * pow(NR, parameter) : tcz * NR) - conc_avg) / (1 + ((NR > cz) ? parameter * Kx * pow(NR,
                    parameter - 1) : tcz));
            if ((NR <= cz) || (fabs(NR - pNR) < ZERO))
                break;
            pNR = NR;
            //                if(NR < 0) printf("\n%f\t",NR);
        }
        *concx = NR;
        break;
    case 3: // langmuir
        parameter = sorp_coef[1];
        NR = 0;
        for (i = 0; i < 5; i++) //Newton Raphson iteration cycle
        {
            NR -= (NR + (NR * Kx * parameter) / (1 + NR * Kx) - conc_avg) / (1 + Kx * parameter / pow(1 + NR * Kx, 2));
            if (fabs(NR - pNR) < ZERO)
                break;
            pNR = NR;
        }
        *concx = NR;
        //   *concx_sorb = (conc_avg - *concx) * Nv;
        break;
    }
    /*   else{
     *concx_sorb = 0.0;
     *concx = 0.0;
     return;
     }

     if((fabs(conc_avg - *concx) > 1e-20))   */
    *concx_sorb = (conc_avg - *concx) * Nv;
    /*   else{
     *concx_sorb = 0.0;
     if(fabs(*concx) < 1e-20 )
     *concx = 0.0;
     else
     *concx = conc_avg;
     }                               */

    /*
     if(DBL_EQ(conc_avg, *concx + *concx_sorb / Nv) != 1){
     printf("\n%f\t%f\t%f\t%f",conc_avg,*concx + *concx_sorb / Nv,*concx,*concx_sorb / Nv);
     getchar();
     } */
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
                    for(int subst_id=0; subst_id<n_substances; subst_id++) {
                        output_time->register_elem_data(substance_name[subst_id], "", out_conc[MOBILE][subst_id], mesh_->n_elements());
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

    int sbi/*, rank, np*/;
    IS is;
    PetscViewer inviewer;

    //	MPI_Barrier(PETSC_COMM_WORLD);
/*    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);*/


    //ISCreateStride(PETSC_COMM_SELF,mesh_->n_elements(),0,1,&is);
    ISCreateGeneral(PETSC_COMM_SELF, mesh_->n_elements(), row_4_el, PETSC_COPY_VALUES, &is); //WithArray
    VecScatterCreate(vconc[0], is, vconc_out[0], PETSC_NULL, &vconc_out_scatter);
    for (sbi = 0; sbi < n_substances; sbi++) {
        VecScatterBegin(vconc_out_scatter, vconc[sbi], vconc_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(vconc_out_scatter, vconc[sbi], vconc_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
    }
    //VecView(transport->vconc[0],PETSC_VIEWER_STDOUT_WORLD);
    //VecView(transport->vconc_out[0],PETSC_VIEWER_STDOUT_WORLD);
    VecScatterDestroy(&(vconc_out_scatter));
    ISDestroy(&(is));
}

//=============================================================================
//      COMPARE DENSITY ITERATION
//=============================================================================
int ConvectionTransport::compare_dens_iter() {
/*
    ElementIter elm;
    double max_err;
    max_err = 0;
    //	FOR_ELEMENTS( elm )
    //		xprintf(Msg,"%f %f %f\n",elm->scalar , elm->scalar_it,elm->scalar - elm->scalar_it);
    FOR_ELEMENTS(mesh_,  elm ) {
        if (fabs(elm->scalar - elm->scalar_it) > max_err) {
            max_err = fabs(elm->scalar - elm->scalar_it);
            //xprintf(Msg,"%f %f %f\n",elm->scalar , elm->scalar_it, elm->scalar - elm->scalar_it);
        }
    }
    xprintf(Msg,"Maximum pressure difference in iteration: %f10.8\n",max_err);
    if (max_err > dens_eps)
        return 0;
    else
        return 1;*/
}
//=============================================================================
//      RESTART ITERATION CONCENTRATION
//=============================================================================
void ConvectionTransport::restart_iteration_C() {
    int sbi, n_subst, sub, ph;


   // struct Transport *transport = problem->transport;
   // double ***conc, ***prev_conc;

    n_subst = n_substances;
    sub = sub_problem;
    //conc = transport->conc;
    //prev_conc = transport->prev_conc;

    for (sbi = 0; sbi < n_subst; sbi++)
        for (ph = 0; ph < 4; ph++)
            if (conc[sbi][ph] != NULL)
                memcpy(conc[sbi][ph], prev_conc[sbi][ph], mesh_->n_elements() * sizeof(double));

}
//=============================================================================
//      SAVE & RESTART ITERATION OF PRESSURE
//=============================================================================
void ConvectionTransport::save_restart_iteration_H() {
/*
    ElementIter elm;
    FOR_ELEMENTS(mesh_,  elm ) {
        elm->scalar_it = elm->scalar;
    }
    */
}
//=============================================================================
//      SAVE TIME STEP CONCENTRATION
//=============================================================================
void ConvectionTransport::save_time_step_C() {

    int sbi, n_subst, sub, ph;
   // struct Transport *transport = problem->transport;
    double ***conc, ***prev_conc;

    n_subst = n_substances;
    sub = sub_problem;
    //conc = transport->conc;
    //prev_conc = transport->prev_conc;

    for (ph = 0; ph < 4; ph++)
        for (sbi = 0; sbi < n_subst; sbi++)
            if (conc[ph][sbi] != NULL)
                memcpy(prev_conc[ph][sbi], conc[ph][sbi], mesh_->n_elements() * sizeof(double));

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

double ***ConvectionTransport::get_prev_concentration_matrix(){
	return pconc;
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

int ConvectionTransport::get_n_substances() {
	return n_substances;
}
