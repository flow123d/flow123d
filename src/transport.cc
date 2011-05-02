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

#include "constantdb.h"
#include "system/system.hh"
#include "system/math_fce.h"
#include "problem.h"
#include "mesh/mesh.h"
#include "transport.h"
#include "output.h"
#include "materials.hh"
#include "read_ini.h"
#include "ppfcs.h"
//#include "btc.h" XX
//#include "reaction.h" XX
#include "flow/darcy_flow_mh.hh"
#include "system/par_distribution.hh"
#include "mesh/ini_constants_mesh.hh"
#include "sparse_graph.hh"
#include "semchem/semchem_interface.hh"
#include "reaction/linear_reaction.hh"
#include <string.h>

static double *transport_aloc_pi(Mesh*);


//=============================================================================
// MAKE TRANSPORT
//=============================================================================
void ConvectionTransport::make_transport_partitioning() {

    F_ENTRY;

    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);
    int rank, np, i, j, k, row_MH, a;
    //struct DarcyFlowMH *water=transport->problem->water;

    SparseGraph *ele_graph = new SparseGraphMETIS(mesh->n_elements()); // graph for partitioning
    Distribution init_ele_ds = ele_graph->get_distr(); // initial distr.
    int *loc_part = new int[init_ele_ds.lsize()]; // partitionig in initial distribution

    make_element_connection_graph(mesh, ele_graph, true);
    WARN_ASSERT(ele_graph->is_symmetric(),"Attention graph for partitioning is not symmetric!\n");

    ele_graph->partition(loc_part);

    delete ele_graph;

    int *id_4_old = (int *) xmalloc(mesh->n_elements() * sizeof(int));
    i = 0;
    FOR_ELEMENTS(ele)
        id_4_old[i] = i, i++;
    id_maps(mesh->n_elements(), id_4_old, init_ele_ds, (int *) loc_part, el_ds, el_4_loc, row_4_el);

    delete[] loc_part;
    xfree(id_4_old);

}

//ConvectionTransport::ConvectionTransport(struct Problem *problemMaterialDatabase, Mesh *init_mesh)
//: problem(problem), mesh(init_mesh)
ConvectionTransport::ConvectionTransport(MaterialDatabase *material_database, Mesh *init_mesh)
: mat_base(material_database), mesh(init_mesh)
{
    transport_init();
}

ConvectionTransport::~ConvectionTransport()
{

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
// MAKE TRANSPORT
//=============================================================================
void ConvectionTransport::transport_init() {
    //struct Transport *transport = problem->transport;
//	mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);
    char *snames, *sscales;
    F_ENTRY;

    // [Density]
    max_dens_it = OptGetInt("Density", "Density_max_iter", "20");
    dens_implicit = OptGetBool("Density", "Density_implicit", "no");
    dens_eps = OptGetDbl("Density", "Eps_iter", "1.0e-5");
    write_iterations = OptGetBool("Density", "Write_iterations", "no");
    dens_step = OptGetInt("Density", "Density_steps", "1");
    // [Transport]
    transport_on = OptGetBool("Transport", "Transport_on", "no");
    sorption = OptGetBool("Transport", "Sorption", "no");
    dual_porosity = OptGetBool("Transport", "Dual_porosity", "no");
    reaction_on = OptGetBool("Transport", "Reactions", "no");
    concentration_fname = IONameHandler::get_instance()->get_input_file_name(OptGetFileName("Transport", "Concentration", "\\"));
    transport_bcd_fname = IONameHandler::get_instance()->get_input_file_name(OptGetFileName("Transport", "Transport_BCD", "\\"));
    transport_out_fname = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Transport", "Transport_out", "\\"));
    transport_out_im_fname = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Transport", "Transport_out_im", "\\"));
    transport_out_sorp_fname = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Transport", "Transport_out_sorp", "\\"));
    transport_out_im_sorp_fname = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Transport", "Transport_out_im_sorp", "\\"));

    pepa = OptGetBool("Transport", "Decay", "no"); //PEPA
    type = OptGetInt("Transport", "Decay_type", "-1"); //PEPA


    DBGMSG("Transport substances.\n");
    n_substances = OptGetInt("Transport", "N_substances", NULL );
    snames = OptGetStr("Transport", "Substances", "none");
    subst_names(snames);
    if (ConstantDB::getInstance()->getInt("Problem_type") == PROBLEM_DENSITY) {
        sscales = OptGetStr("Transport", "Substances_density_scales", "1.0");
        subst_scales(sscales);
    }

   // n_elements = mesh->n_elements();

/*
    reaction = NULL;
    n_reaction = 0;
*/
    sub_problem = 0;
    if (dual_porosity == true)
        sub_problem += 1;
    if (sorption == true)
        sub_problem += 2;

    frame = 0;
    time = 0;

  //  problem->material_database->read_transport_materials(dual_porosity, sorption,
  //          n_substances);



        make_transport_partitioning();
        alloc_transport_vectors();
        read_initial_condition();
        alloc_transport_structs_mpi();
        fill_transport_vectors_mpi();

        transport_output_init();
        transport_output();

    INPUT_CHECK(!(n_substances < 1 ),"Number of substances must be positive\n");
}
//=============================================================================
//
//=============================================================================
//char ConvectionTransport::**subst_names(int n_subst, char *line) {
void ConvectionTransport::subst_names(char *line) {
    int sbi;

    ASSERT(!( (n_substances < 1) || (line == NULL) ),"Bad parameter of function subst_names()\n");
    substance_name = (char**) xmalloc(n_substances * sizeof(char*));
    for (sbi = 0; sbi < n_substances; sbi++)
        substance_name[sbi] = xstrcpy(strtok((sbi == 0 ? line : NULL), " \t,;"));
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
void ConvectionTransport::read_initial_condition() {
		Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);
		FILE	*in;		  // input file
		char     line[ LINE_SIZE ]; // line of data file
		int sbi,index, id, eid, i,n_concentrations, global_idx;

		xprintf( Msg, "Reading concentrations...")/*orig verb 2*/;
	//	getchar();
	//	printf("%s\n",concentration_fname);
	//	getchar();
		in = xfopen( concentration_fname, "rt" );
		//in = xfopen( "test1.tic", "rt" );

		skip_to( in, "$Concentrations" );
		xfgets( line, LINE_SIZE - 2, in );
		n_concentrations = atoi( xstrtok( line) );
		INPUT_CHECK(!( n_concentrations < 1 ),"Number of concentrations < 1 in function read_concentration_list()\n");
	    INPUT_CHECK(!( mesh->n_elements() != n_concentrations),"Different number of elements and concentrations\n");


	    for (i = 0; i < n_concentrations; i++) {
	    	//printf("%s\n",line);
	    	xfgets( line, LINE_SIZE - 2, in );
	    	ASSERT(!(line == NULL),"NULL as argument of function parse_concentration_line()\n");
	    	id    = atoi( xstrtok( line) );	// TODO: id musi byt >0 nebo >=0 ???
	    	INPUT_CHECK(!( id < 0 ),"Id number of concentration must be > 0\n");
	    	eid    = atoi( xstrtok( NULL) );
	    	int global_idx =row_4_el[mesh->element.find_id(eid).index()];
	    	if ( el_ds->is_local(global_idx) ) {
	    		index = global_idx - el_ds->begin();
	    		for( sbi = 0; sbi < n_substances; sbi++ ){
	    			conc[MOBILE][ sbi ][index] = atof( xstrtok( NULL) );
	    			pconc[MOBILE][ sbi ][index] = conc[MOBILE][ sbi ][index];
	    		}
	    	}
	    }

		xfclose( in );
		xprintf( MsgVerb, " %d concentrations readed. ", n_concentrations )/*orig verb 4*/;
		xprintf( Msg, "O.K.\n")/*orig verb 2*/;


        // TODO:
        // Why explicit application of boundary condition ?
        // Even if it is necessary , there should be i = 0; before cycle
        /*
         FOR_ELEMENTS(elm)
         FOR_ELEMENT_SIDES(elm,s)
         if (elm->side[s]->cond != NULL) {
         transport->conc[MOBILE][sbi][i]
         = elm->side[s]->cond->transport_bcd->conc[sbi];
         transport->pconc[MOBILE][sbi][i++]
         = elm->side[s]->cond->transport_bcd->conc[sbi];
         }*/


}
//=============================================================================
//	ALLOCATE OF TRANSPORT VARIABLES (ELEMENT & NODES)
//=============================================================================
void ConvectionTransport::alloc_transport_vectors() {
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    int i, j, sbi, n_subst, ph;
    ElementIter elm;
    n_subst = n_substances;
    

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
	for(sbi = 0; sbi < n_subst; sbi++){
           conc[ph][sbi] = (double*) xmalloc(el_ds->lsize() * sizeof(double));
           pconc[ph][sbi] = (double*) xmalloc(el_ds->lsize() * sizeof(double));
           out_conc[ph][sbi] = (double*) xmalloc(el_ds->size() * sizeof(double));
           // transport->node_conc[sbi][ph] = (double**)xmalloc((mesh->n_elements() ) * sizeof(double*));
           for (i = 0; i < el_ds->lsize(); i++) {
             conc[ph][sbi][i] = 0.0;
             pconc[ph][sbi][i] = 0.0;
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
      }else {
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
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    int ph, sbi, i, sub;
    int n_subst = n_substances;
    int n_elements = mesh->n_elements();

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
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

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

    // if( rank == 0)
    vconc_out = (Vec*) xmalloc(n_subst * (sizeof(Vec))); // extend to all

    // TODO: should be replaced by Distribution(Block) or better remove whole boundary matrix with these vectors
    lb_col = (int*) xmalloc(np * sizeof(int)); // local rows in BM
    for (i = 0; i < np; i++) {
        lb_col[i] = mesh->n_boundaries() / np;
    }
    lb_col[np - 1] = mesh->n_boundaries() - (np - 1) * lb_col[0];

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
        //ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->lb_col[rank],mesh->n_boundaries(),&transport->bcv[sbi]);
        //ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh->n_elements(),&transport->bcvcorr[sbi]);
        //ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh->n_elements(),&transport->vconc[sbi]);
        //ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh->n_elements(),&transport->vpconc[sbi]);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, lb_col[rank], mesh->n_boundaries(), &bcv[sbi]);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, el_ds->lsize(), mesh->n_elements(), &bcvcorr[sbi]);
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, el_ds->lsize(), mesh->n_elements(), conc[MOBILE][sbi],
                &vconc[sbi]);
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, el_ds->lsize(), mesh->n_elements(),
                pconc[MOBILE][sbi], &vpconc[sbi]);

        //  if(rank == 0)
        ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, mesh->n_elements(), out_conc[MOBILE][sbi], &vconc_out[sbi]);

        //ierr = VecCreateMPI(PETSC_COMM_SELF ,transport->mesh->n_elements(),transport->mesh->n_elements(),&transport->vconc_out[sbi]); /*xx*/
        /*
         if(transport->vconc_im != NULL)
         ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh->n_elements(),&transport->vconc_im[sbi]);
         if(transport->vconc_so != NULL)
         ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh->n_elements(),&transport->vconc_so[sbi]);
         if(transport->vconc_im_so != NULL)
         ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh->n_elements(),&transport->vconc_im_so[sbi]);
         */
    }
    //
    //ierr=MatCreateMPIAIJ(PETSC_COMM_WORLD,transport->l_row[rank],transport->l_row[rank],mesh->n_elements(),mesh->n_elements(),
    //		PETSC_NULL,transport->d_row[rank],PETSC_NULL,transport->od_row[rank],&transport->tm);

    //ierr=MatCreateMPIAIJ(PETSC_COMM_WORLD,transport->l_row[rank],transport->lb_col[rank],mesh->n_elements(),mesh->n_boundaries(),
    //				PETSC_NULL,transport->db_row[rank],PETSC_NULL,transport->odb_row[rank],&transport->bcm);

    ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD, el_ds->lsize(), el_ds->lsize(), mesh->n_elements(),
            mesh->n_elements(), 8, PETSC_NULL, 1, PETSC_NULL, &tm);

    ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD, el_ds->lsize(), lb_col[rank], mesh->n_elements(),
            mesh->n_boundaries(), 2, PETSC_NULL, 0, PETSC_NULL, &bcm);

    /*
     MatCreateMPIAIJ(comm,m,n,M,N,0,PETSC_NULL,0,PETSC_NULL,&transport->tm);
     MatCreateMPIAIJ(comm,m,n,M,N,0,PETSC_NULL,0,PETSC_NULL,&transport->bcm);

     ierr=MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,transport->l_row[rank],mesh->n_elements(),mesh->n_elements(),
     PETSC_NULL,transport->d_row[rank],PETSC_NULL,transport->od_row[rank],&transport->tm);

     ierr=MatCreateMPIAIJ(PETSC_COMM_WORLD,transport->l_row[rank],transport->lb_col[rank],mesh->n_elements(),mesh->n_boundaries(),
     PETSC_NULL,transport->db_row[rank],PETSC_NULL,transport->odb_row[rank],&transport->bcm);
     */
}
//=============================================================================
//	FILL TRANSPORT VECTORS (MPI)
//=============================================================================
void ConvectionTransport::fill_transport_vectors_mpi() {

    int rank, sbi;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0) {

    xprintf( Msg, "Reading transport boundary conditions...");

    int bcd_id, boundary_id, boundary_index;
    double bcd_conc;
    char     line[ LINE_SIZE ]; // line of data file
    const std::string& transport_file_name = IONameHandler::get_instance()->get_input_file_name(OptGetFileName("Transport", "Transport_BCD", "\\"));
    FILE *in = xfopen( transport_file_name, "rt" );

    skip_to( in, "$Transport_BCD" );
    xfgets( line, LINE_SIZE - 2, in );
    int n_bcd = atoi( xstrtok( line) );
    for(int i_bcd=0; i_bcd<n_bcd; i_bcd++) {
        xfgets( line, LINE_SIZE - 2, in );
        bcd_id    = atoi( xstrtok( line) ); // scratch transport bcd id
        boundary_id    = atoi( xstrtok( NULL) );
//        DBGMSG("transp b. id: %d\n",boundary_id);
        boundary_index = mesh->boundary.find_id(boundary_id).index();
        INPUT_CHECK(boundary_index >= 0,"Wrong boundary index %d for bcd id %d in transport bcd file!", boundary_id, bcd_id);
        for( sbi = 0; sbi < n_substances; sbi++ ) {
            bcd_conc = atof( xstrtok( NULL) );
            VecSetValue(bcv[sbi], boundary_index, bcd_conc, INSERT_VALUES);
        }
    }
    xfclose( in );
    xprintf( MsgVerb, " %d transport conditions read. ", n_bcd );
    xprintf( Msg, "O.K.\n");
    }

    for(sbi=0;sbi < n_substances;sbi++) VecAssemblyBegin(bcv[sbi]);
    for(sbi=0;sbi < n_substances;sbi++) VecZeroEntries(bcvcorr[sbi]);
    for(sbi=0;sbi < n_substances;sbi++) VecAssemblyEnd(bcv[sbi]);


    /*
     VecView(transport->bcv[0],PETSC_VIEWER_STDOUT_SELF);
     getchar();
     VecView(transport->bcv[1],PETSC_VIEWER_STDOUT_SELF);
     getchar();
     */

}
//=============================================================================
// CREATE TRANSPORT MATRIX
//=============================================================================
void ConvectionTransport::create_transport_matrix_mpi() {
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    ElementFullIter el2 = ELEMENT_FULL_ITER_NULL;
    ElementFullIter elm = ELEMENT_FULL_ITER_NULL;
    struct Edge *edg;
    struct Neighbour *ngh;
    //struct Transport *transport;
    int n, s, i, j, np, rank, new_j, new_i;
    double max_sum, flux, aij, aii, *solution;
    /*
    DarcyFlow *water;

    water = problem->water;
    solution = water.solution();
    id = water
*/

    /*
     FOR_ELEMENTS(elm)
     FOR_ELEMENT_SIDES(elm,i){
     printf("id: %d side: %d flux:  %f\n",elm->id, i,elm->side[i]->flux);
     getchar();
     }

     getchar();
     */




    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);

    FOR_EDGES(edg) { // calculate edge Qv
        edg->faux = 0;
        FOR_EDGE_SIDES(edg,s)
            if (edg->side[s]->flux > 0)
                edg->faux += edg->side[s]->flux;
    }

    max_sum = 0.0;
    aii = 0.0;
    for (int loc_el = 0; loc_el < el_ds->lsize(); loc_el++) {
        elm = mesh->element(el_4_loc[loc_el]);
        new_i = row_4_el[elm.index()];

        FOR_ELEMENT_SIDES(elm,si) // same dim
            if (elm->side[si]->cond == NULL) {
                if (elm->side[si]->flux < 0.0) {
                    if (elm->side[si]->neigh_bv != NULL) { //comp model
                        aij = -(elm->side[si]->flux / (elm->volume * elm->material->por_m));
                        j = ELEMENT_FULL_ITER(elm->side[si]->neigh_bv->element[0]).index();
                        new_j = row_4_el[j];
                        MatSetValue(tm, new_i, new_j, aij, INSERT_VALUES);
                    }
                    // end comp model
                    else {
                        edg = elm->side[si]->edge;
                        if (edg->faux > ZERO)
                            FOR_EDGE_SIDES(edg,s)
                                if ((edg->side[s]->id != elm->side[si]->id) && (edg->side[s]->flux > 0.0)) {
                                    aij = -(elm->side[si]->flux * edg->side[s]->flux / (edg->faux * elm->volume
                                            * elm->material->por_m));
                                    j = ELEMENT_FULL_ITER(edg->side[s]->element).index();
                                    new_j = row_4_el[j];
                                    MatSetValue(tm, new_i, new_j, aij, INSERT_VALUES);
                                }
                    }
                }
                if (elm->side[si]->flux > 0.0)
                    aii -= (elm->side[si]->flux / (elm->volume * elm->material->por_m));
            } else {
                if (elm->side[si]->flux < 0.0) {
                    aij = -(elm->side[si]->flux / (elm->volume * elm->material->por_m));
                    j = BOUNDARY_FULL_ITER(elm->side[si]->cond).index();
                    MatSetValue(bcm, new_i, j, aij, INSERT_VALUES);
                    // vyresit BC matrix !!!!
                    //   printf("side in elm:%d value:%f\n ",elm->id,svector->val[j-1]);
                    //   printf("%d\t%d\n",elm->id,id2pos(problem,elm->side[si]->id,problem->spos_id,BC));

                }
                if (elm->side[si]->flux > 0.0)
                    aii -= (elm->side[si]->flux / (elm->volume * elm->material->por_m));
            } // end same dim     //ELEMENT_SIDES


        FOR_ELM_NEIGHS_VB(elm,n) // comp model
            FOR_NEIGH_ELEMENTS(elm->neigh_vb[n],s)
                if (elm.id() != ELEMENT_FULL_ITER(elm->neigh_vb[n]->element[s]).id()) {
                    if (elm->neigh_vb[n]->side[s]->flux > 0.0) {
                        aij = elm->neigh_vb[n]->side[s]->flux / (elm->volume * elm->material->por_m);
                        j = ELEMENT_FULL_ITER(elm->neigh_vb[n]->element[s]).index();
                        new_j = row_4_el[j];
                        MatSetValue(tm, new_i, new_j, aij, INSERT_VALUES);
                    }
                    if (elm->neigh_vb[n]->side[s]->flux < 0.0)
                        aii += elm->neigh_vb[n]->side[s]->flux / (elm->volume * elm->material->por_m);
                } // end comp model
        FOR_ELM_NEIGHS_VV(elm,n) { //non-comp model
            ngh = elm->neigh_vv[n];
            FOR_NEIGH_ELEMENTS(ngh,s) {

                el2 = ELEMENT_FULL_ITER(ngh->element[s]);
                if (elm.id() != el2.id()) {
                    flux = ngh->sigma * ngh->geom_factor * (el2->scalar - elm->scalar);
                    if (flux > 0.0) {
                        aij = flux / (elm->volume * elm->material->por_m); // -=
                        j = el2.index();
                        new_j = row_4_el[j];
                        MatSetValue(tm, new_i, new_j, aij, INSERT_VALUES);
                    }
                    if (flux < 0.0)
                        aii += flux / (elm->volume * elm->material->por_m);
                }
            }
        } // end non-comp model

        MatSetValue(tm, new_i, new_i, aii, INSERT_VALUES);

        if (fabs(aii) > max_sum)
            max_sum = fabs(aii);
        aii = 0.0;
        //   i++;
    } // END ELEMENTS

    double glob_max_sum;

    MPI_Allreduce(&max_sum,&glob_max_sum,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);

    MatAssemblyBegin(tm, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(bcm, MAT_FINAL_ASSEMBLY);

    max_step = 1 / glob_max_sum;
    time_step = 0.9 / glob_max_sum;

    MatAssemblyEnd(tm, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(bcm, MAT_FINAL_ASSEMBLY);

    // MPI_Barrier(PETSC_COMM_WORLD);
    /*
     MatView(transport->tm,PETSC_VIEWER_STDOUT_SELF);
     getchar();
     */

}
//=============================================================================
// CREATE TRANSPORT MATRIX STEP MPI
//=============================================================================
void ConvectionTransport::transport_matrix_step_mpi(double time_step) {
    MatScale(bcm, time_step);
    MatScale(tm, time_step);
    MatShift(tm, 1.0);

    /*
     MatView(transport->bcm,PETSC_VIEWER_STDOUT_WORLD);
     getchar();
     */

}
//=============================================================================
// MATRIX VECTOR PRODUCT (MPI)
//=============================================================================
void ConvectionTransport::calculate_bc_mpi() {
    int sbi, n_subst;

    n_subst = n_substances;

    for (sbi = 0; sbi < n_subst; sbi++)
        MatMult(bcm, bcv[sbi], bcvcorr[sbi]);
    /*
     VecView(transport->bcvcorr[0],PETSC_VIEWER_STDOUT_SELF);
     getchar();
     */
}
//=============================================================================
// MATRIX VECTOR PRODUCT (MPI)
//=============================================================================
void ConvectionTransport::transport_step_mpi(Mat *tm, Vec *conc, Vec *pconc, Vec *bc) {

    MatMultAdd(*tm, *pconc, *bc, *conc); // conc=tm*pconc + bc
    VecSwap(*conc, *pconc); // pconc = conc


}
//=============================================================================
// ALLOCATING MEMORY FOR TRANSPORT VARIABLE PI
//=============================================================================
double *transport_aloc_pi(Mesh* mesh) {
    int max_elm;
    double *pi;

    NodeIter nod;

    max_elm = 0;
    FOR_NODES(nod)
        if (max_elm < nod->n_elements)
            max_elm = nod->n_elements;
    pi = (double *) xmalloc(max_elm * sizeof(double));

    return pi;
}
//=============================================================================
// COMPUTE CONCENTRATIONS IN THE NODES FROM THE ELEMENTS
//=============================================================================
/*
 void transport_node_conc(struct Transport *transport)
 {
 TNode* nod;
 Mesh* mesh = transport->mesh;
 int i;
 double P,N,scale, *pi;
 int min_elm_dim;
 int sbi,sub;
 pi = transport_aloc_pi(mesh);

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
// assume: elm_pos - position of values of pconc and conc vectors corresponding to a mesh element
//         sbi - matter index
//         material - material on corresponding mesh element
//=============================================================================
void ConvectionTransport::transport_dual_porosity( int elm_pos, MaterialDatabase::Iter material, int sbi) {
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

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
        cm = (pcm - conc_avg) * exp(-alpha * ((por_m + por_imm) / (por_m * por_imm)) * time_step) + conc_avg;

        // ---compute concentration in immobile area---------------------------------
        ci = (pci - conc_avg) * exp(-alpha * ((por_m + por_imm) / (por_m * por_imm)) * time_step) + conc_avg;
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
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    double conc_avg = 0.0;
    double conc_avg_imm = 0.0;
    double n, Nm, Nimm;
    int id;
    double phi;

    phi = mtr->phi;

    if ((mtr->sorp_coef[sbi][0] == 0) || (mtr->por_m == 1))
        return;

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
void ConvectionTransport::compute_time_step() {

    double problem_save_step = OptGetDbl("Global", "Save_step", "1.0");
    double problem_stop_time = OptGetDbl("Global", "Stop_time", "1.0");
    save_step = (int) ceil(problem_save_step / time_step); // transport  rev
    time_step = problem_save_step / save_step;
    steps = save_step * (int) floor(problem_stop_time / problem_save_step);

}
//=============================================================================
//      TRANSPORT ONE STEP
//=============================================================================
void ConvectionTransport::transport_one_step() {


}
//=============================================================================
//      TRANSPORT UNTIL TIME
//=============================================================================
void ConvectionTransport::transport_one_step(double time_interval) {


}
//=============================================================================
//      CONVECTION
//=============================================================================
void ConvectionTransport::convection() {
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);
    MaterialDatabase::Iter material;

    int step;
    register int t;
    int n_subst,sbi,elm_pos,rank,i,size;
//  	double **reaction_matrix;
  	Linear_reaction *decayRad;


    START_TIMER("TRANSPORT");

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    xprintf( Msg, "Calculating transport...")/*orig verb 2*/;
    n_subst = n_substances;

    //  flow_cs(trans); //DECOVALEX
    // int tst = 1; // DECOVALEX

    create_transport_matrix_mpi();
    compute_time_step();


    /*  }
     else{
     steps = (int)floor(problem->transport->update_dens_time / problem->transport->time_step);
     save_step = steps + 1;
     }*/

    transport_matrix_step_mpi(time_step); //TIME STEP
    calculate_bc_mpi();

     if (rank == 0) {
        printf("  %d computing cycles, %d writing steps\n", steps, ((int) (steps / save_step) + 1))/*orig verb 6*/;
        xprintf( MsgVerb, "  %d computing cycles, %d writing steps\n",steps ,((int)(steps / save_step) + 1) )/*orig verb 6*/;
    }
    //output_FCS(trans); //DECOVALEX


    step = 0;
    //fw_chem = fopen("vystup.txt","w"); fclose(fw_chem); //makes chemistry output file clean, before transport is computed
    for (t = 1; t <= steps; t++) {
    	time += time_step;
     //   SET_TIMER_SUBFRAMES("TRANSPORT",t);  // should be in destructor as soon as we have class iteration counter
        step++;
        for (sbi = 0; sbi < n_subst; sbi++) {
            /*
             if(tst && (particle_test(trans) > 100.0) && tst){ // DECOVALEX
             clear_tbc(trans);
             tst = 0;
             }
             output_AGE(trans,(t-1) * trans->time_step); // DECOVALEX
             */

            transport_step_mpi(&tm, &vconc[sbi], &vpconc[sbi], &bcvcorr[sbi]);

            if ((dual_porosity == true) || (sorption == true) || (pepa == true) || (reaction_on == true))
                // cycle over local elements only in any order
                for (int loc_el = 0; loc_el < el_ds->lsize(); loc_el++) {
                    material = (mesh->element(el_4_loc[loc_el])) -> material;

                    if (dual_porosity == true)
                        transport_dual_porosity(loc_el, material, sbi);
                    if (sorption == true)
                        transport_sorption(loc_el, material, sbi);
                    /*
                     if (reaction_on == true)
                     transport_reaction(trans, loc_el, material, sbi);

                     */
                }
            // transport_node_conc(mesh,sbi,problem->transport_sub_problem);  // vyresit prepocet
        }
        xprintf( Msg, "Time : %f\n",time);
        //======================================
        //              CHEMISTRY
        //======================================
    if(OptGetBool("Semchem_module", "Compute_reactions", "no") == true)
    {
            if (t == 1) { //initial value of t == 1 & it is incremented at the beginning of the cycle
                priprav();
            }
            for (int loc_el = 0; loc_el < el_ds->lsize(); loc_el++) {
                //xprintf(Msg,"\nKrok %f\n",trans->time_step);
                che_vypocetchemie(dual_porosity, time_step, mesh->element(el_4_loc[loc_el]), loc_el, conc[MOBILE], conc[IMMOBILE]);
            }// for cycle running over elements
        }
    //===================================================
    //     RADIOACTIVE DECAY + FIRST ORDER REACTIONS
    //===================================================
    if(OptGetBool("Decay_module", "Compute_decay", "no") == true){
            int rows, cols, dec_nr, nr_of_decay, dec_name_nr = 1;
            //char dec_name[30];

            if (t == 1) {
                decayRad = new Linear_reaction(n_subst, time_step);
    	}
            for (int loc_el = 0; loc_el < el_ds->lsize(); loc_el++) {
    		(*decayRad).compute_reaction(pconc[MOBILE], n_subst, loc_el);
                if (dual_porosity == true) {
    			(*decayRad).compute_reaction(pconc[IMMOBILE], n_subst, loc_el);
                }
            }
    }/*
    else{
            xprintf(Msg,"\nDecay is not computed.\n");
    }*/
        //======================================

        //   save_step == step;
        //&& ((ConstantDB::getInstance()->getInt("Problem_type") != PROBLEM_DENSITY)
        if ((save_step == step) || (write_iterations)) {
            xprintf( Msg, "Output\n");
            //if (size != 1)
            	//time = t * time_step;
            	transport_output();
            //transport_output(trans, t * time_step, ++frame);
          //  if (ConstantDB::getInstance()->getInt("Problem_type") != STEADY_SATURATED)
               // output_time(problem, t * time_step); // time variable flow field
            //	output_transport_time_BTC(trans, t * trans->time_step); // BTC test - spatne vypisuje casy
            //output_transport_time_CS(problem, t * problem->time_step);
            step = 0;
            //sorb_mob_arr = NULL;
        }
    }
    xprintf( Msg, "O.K.\n");
    transport_output_finish();
}

//=============================================================================
//      OUTPUT VECTOR GATHER
//=============================================================================
void ConvectionTransport::output_vector_gather() {
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    int sbi, rank, np;
    IS is;
    PetscViewer inviewer;

    //	MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);


    //ISCreateStride(PETSC_COMM_SELF,mesh->n_elements(),0,1,&is);
    ISCreateGeneral(PETSC_COMM_SELF, mesh->n_elements(), row_4_el, &is); //WithArray
    VecScatterCreate(vconc[0], is, vconc_out[0], PETSC_NULL, &vconc_out_scatter);
    for (sbi = 0; sbi < n_substances; sbi++) {
        VecScatterBegin(vconc_out_scatter, vconc[sbi], vconc_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(vconc_out_scatter, vconc[sbi], vconc_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
    }
    //VecView(transport->vconc[0],PETSC_VIEWER_STDOUT_WORLD);
    //VecView(transport->vconc_out[0],PETSC_VIEWER_STDOUT_WORLD);
    VecScatterDestroy(vconc_out_scatter);
    ISDestroy(is);
}
//=============================================================================
//      TRANSPORT OUTPUT
//=============================================================================
//void transport_output(double ***out_conc,char **subst_names ,int n_subst,double time, int frame, char *transport_out_fname) {
void ConvectionTransport::transport_output() {
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	output_vector_gather();
	if (rank == 0){
		switch (ConstantDB::getInstance()->getInt("Pos_format_id")) {
		case POS_BIN:
			output_transport_time_bin( out_conc, substance_name ,n_substances, time, frame, (char*)transport_out_fname.c_str());
			break;
		case POS_ASCII:
			output_transport_time_ascii(out_conc, substance_name ,n_substances, time, frame, (char*)transport_out_fname.c_str());
			break;
		case VTK_SERIAL_ASCII:
			output_transport_time_vtk_serial_ascii(out_conc, substance_name ,n_substances, time, frame, (char*)transport_out_fname.c_str());
			break;
		case VTK_PARALLEL_ASCII:
			xprintf(UsrErr, "VTK_PARALLEL_ASCII: not implemented yet\n");
			break;
		}
	}
	frame++;
}
//=============================================================================
//      TRANSPORT OUTPUT INIT
//=============================================================================
//void transport_output_init(char *transport_out_fname) {
void ConvectionTransport::transport_output_init() {
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if (rank == 0){
		switch (ConstantDB::getInstance()->getInt("Pos_format_id")) {
		case POS_BIN:
			output_msh_init_bin(mesh, (char*)transport_out_fname.c_str());
			break;
		case POS_ASCII:
			output_msh_init_ascii(mesh, (char*)transport_out_fname.c_str());
			break;
    	case VTK_SERIAL_ASCII:
    		output_msh_init_vtk_serial_ascii( (char*)transport_out_fname.c_str());
    		break;
    	case VTK_PARALLEL_ASCII:
    		xprintf(UsrErr, "VTK_PARALLEL_ASCII: not implemented yet\n");
    		break;
		}
	}
}
//=============================================================================
//      TRANSPORT OUTPUT FINISH
//=============================================================================
void ConvectionTransport::transport_output_finish() {
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if (rank == 0){
		switch (ConstantDB::getInstance()->getInt("Pos_format_id")) {
		case POS_BIN:
			/* There is no need to do anything for this file format */
			break;
		case POS_ASCII:
			/* There is no need to do anything for this file format */
			break;
		case VTK_SERIAL_ASCII:
			output_msh_finish_vtk_serial_ascii((char*)transport_out_fname.c_str());
			break;
		case VTK_PARALLEL_ASCII:
			xprintf(UsrErr, "VTK_PARALLEL_ASCII: not implemented yet\n");
			break;
		}
	}
}
//=============================================================================
//      COMPARE DENSITY ITERATION
//=============================================================================
int ConvectionTransport::compare_dens_iter() {
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    ElementIter elm;
    double max_err;
    max_err = 0;
    //	FOR_ELEMENTS( elm )
    //		xprintf(Msg,"%f %f %f\n",elm->scalar , elm->scalar_it,elm->scalar - elm->scalar_it);
    FOR_ELEMENTS( elm ) {
        if (fabs(elm->scalar - elm->scalar_it) > max_err) {
            max_err = fabs(elm->scalar - elm->scalar_it);
            //xprintf(Msg,"%f %f %f\n",elm->scalar , elm->scalar_it, elm->scalar - elm->scalar_it);
        }
    }
    xprintf(Msg,"Maximum pressure difference in iteration: %f10.8\n",max_err);
    if (max_err > dens_eps)
        return 0;
    else
        return 1;
}
//=============================================================================
//      RESTART ITERATION CONCENTRATION
//=============================================================================
void ConvectionTransport::restart_iteration_C() {
    int sbi, n_subst, sub, ph;

    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

   // struct Transport *transport = problem->transport;
   // double ***conc, ***prev_conc;

    n_subst = n_substances;
    sub = sub_problem;
    //conc = transport->conc;
    //prev_conc = transport->prev_conc;

    for (sbi = 0; sbi < n_subst; sbi++)
        for (ph = 0; ph < 4; ph++)
            if (conc[sbi][ph] != NULL)
                memcpy(conc[sbi][ph], prev_conc[sbi][ph], mesh->n_elements() * sizeof(double));

}
//=============================================================================
//      SAVE & RESTART ITERATION OF PRESSURE
//=============================================================================
void ConvectionTransport::save_restart_iteration_H() {
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    ElementIter elm;
    FOR_ELEMENTS( elm ) {
        elm->scalar_it = elm->scalar;
    }
}
//=============================================================================
//      SAVE TIME STEP CONCENTRATION
//=============================================================================
void ConvectionTransport::save_time_step_C() {
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

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
                memcpy(prev_conc[ph][sbi], conc[ph][sbi], mesh->n_elements() * sizeof(double));

}
