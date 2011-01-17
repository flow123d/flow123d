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

#include "constantdb.h"
#include "system.hh"
#include "math_fce.h"
#include "problem.h"
#include "mesh.h"
#include "transport.h"
#include "transport_bcd.h"
#include "boundaries.h"
#include "concentrations.h"
#include "neighbours.h"
#include "elements.h"
#include "output.h"
#include "materials.hh"
#include "read_ini.h"
#include "ppfcs.h"
#include "btc.h"
#include "reaction.h"
#include "darcy_flow_mh.hh"
#include "par_distribution.hh"
#include "mesh/ini_constants_mesh.hh"
#include "sparse_graph.hh"
#include "profiler.hh"

//void init_transport_vectors_mpi(struct Transport *transport);

static void create_transport_matrix_mpi(struct Transport *transport);
static void alloc_transport_structs_mpi(struct Transport *transport);
//static void create_transport_structures_mpi(struct Transport *transport);
static void fill_transport_vectors_mpi(struct Transport *transport);
static void calculate_bc_mpi(struct Transport *transport);
static void transport_partioning(struct Transport *transport,int np);
static void output_vector_gather(struct Transport *transport);
static void transport_matrix_step_mpi(struct Transport *transport,double time_step);
static void make_transport_partitioning(struct Transport *transport);

static void alloc_transport_matrix(struct Transport *transport);
static void alloc_transport_vectors(struct Transport *transport);
static void alloc_density_vectors(struct Transport *transport);
static void transport_vectors_init(struct Transport *transport);
static void create_transport_matrix(struct Problem *problem);
static void transport_matrix_step(struct TMatrix *tmatrix,double time_step);
static void matvecs(struct TMatrix *tmatrix,double *pconc,double *conc);
static double *transport_aloc_pi(Mesh*);
static void transport_dual_porosity(struct Transport *transport, int elm_pos, int sbi);
static void transport_sorption(struct Transport *transport, int elm_pos, int sbi);
static void compute_sorption(double conc_avg,vector<double> &sorp_coef,int sorp_type,double *concx,
        double *concx_sorb,double Nv, double N);
static char **subst_names( int n_subst, char *line );
static double *subst_scales( int n_subst, char *line );

//PEPA CHUDOBA
static void decay( struct Transport *transport, int elm_pos, int type);


//=============================================================================
// MAKE TRANSPORT
//=============================================================================
void make_transport(struct Transport *transport){
    F_ENTRY;

//	transport_init(transport);

    make_id2pos_list();
    transport->problem->material_database->read_transport_materials(transport->dual_porosity, transport->sorption, transport->n_substances);


	if(transport->mpi != 1){
		alloc_transport_matrix(transport);
		alloc_transport_vectors(transport);
		transport_vectors_init(transport);
		if(ConstantDB::getInstance()->getInt("Problem_type") == PROBLEM_DENSITY)
			alloc_density_vectors(transport);
	}
	else{
        alloc_transport_vectors(transport); 
        transport_vectors_init(transport); 
	make_transport_partitioning(transport);
        alloc_transport_structs_mpi(transport);
        fill_transport_vectors_mpi(transport);
	}

	/*
	FOR_ELEMENTS_IT(i){
		printf("%d\t%d\n",i->id,i - transport->mesh->element->begin());
		getchar();
	}
*/

}
//=============================================================================
// MAKE TRANSPORT
//=============================================================================
void make_transport_partitioning(struct Transport *transport){

  Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);
  int rank,np,i,j,k,row_MH,a;
  //struct DarcyFlowMH *water=transport->problem->water;
  
  SparseGraph *ele_graph = new SparseGraphMETIS(mesh->n_elements());                     // graph for partitioning
  Distribution init_ele_ds = ele_graph->get_distr();  // initial distr.
  int *loc_part = new int[init_ele_ds.lsize()];                                     // partitionig in initial distribution

  make_element_connection_graph(mesh, ele_graph, true);
  WARN_ASSERT(ele_graph->is_symmetric(),"Attention graph for partitioning is not symmetric!\n");

  ele_graph->partition(loc_part);

  delete ele_graph;

  int *id_4_old = (int *) xmalloc(mesh->n_edges * sizeof(int));
  i = 0;
  FOR_ELEMENTS(ele)  id_4_old[i] = i, i++;
  id_maps(mesh->n_elements(), id_4_old, init_ele_ds, (int *) loc_part,
          transport->el_ds, transport->el_4_loc, transport->row_4_el);

  delete loc_part;
  xfree(id_4_old);

}
//=============================================================================
// ALLOCATE TRANSPORT
//=============================================================================
void alloc_transport(struct Problem *problem){
	struct Transport *transport;

	problem->transport = (struct Transport*)xmalloc(sizeof(struct Transport));
	transport = problem->transport;
	transport->problem = problem;
	transport->mpi = 1;
}
//=============================================================================
// MAKE TRANSPORT
//=============================================================================
void transport_init(struct Problem *problem)
{
    struct Transport *transport = problem->transport;
	char *snames, *sscales;

	F_ENTRY;

	// [Density]
    transport->max_dens_it      = OptGetInt( "Density", "Density_max_iter", "20");
    transport->dens_implicit    = OptGetBool( "Density", "Density_implicit", "no");
    transport->dens_eps         = OptGetDbl( "Density", "Eps_iter", "1.0e-5" );
    transport->write_iterations = OptGetBool( "Density", "Write_iterations", "no");
    transport->dens_step        = OptGetInt( "Density", "Density_steps", "1" );
    // [Transport]
    transport -> transport_on   = OptGetBool( "Transport", "Transport_on", "no" );
    transport -> sorption       = OptGetBool( "Transport", "Sorption", "no" );
    transport -> dual_porosity  = OptGetBool( "Transport", "Dual_porosity", "no" );
    transport -> reaction_on    = OptGetBool( "Transport", "Reactions", "no");
    transport -> concentration_fname  = OptGetFileName( "Transport", "Concentration", "\\" );
    transport -> transport_bcd_fname  = OptGetFileName( "Transport", "Transport_BCD", "\\" );
    transport -> transport_out_fname  = OptGetFileName( "Transport", "Transport_out", "\\" );
    transport -> transport_out_im_fname       = OptGetFileName( "Transport", "Transport_out_im", "\\" );
    transport -> transport_out_sorp_fname     = OptGetFileName( "Transport", "Transport_out_sorp", "\\" );
    transport -> transport_out_im_sorp_fname  = OptGetFileName( "Transport", "Transport_out_im_sorp", "\\" );

    transport -> pepa   = OptGetBool( "Transport", "Decay", "no" ); //PEPA
    transport -> type   = OptGetInt( "Transport", "Decay_type", "-1" ); //PEPA

    DBGMSG("Transport substances.\n");
	transport->n_substances     = OptGetInt( "Transport", "N_substances", NULL );
 	snames				        = OptGetStr( "Transport", "Substances", "none" );
 	transport -> substance_name	= subst_names( transport->n_substances, snames );
	if (ConstantDB::getInstance()->getInt("Problem_type") == PROBLEM_DENSITY) {
		sscales  = OptGetStr( "Transport", "Substances_density_scales", "1.0" );
    	transport -> substance_density_scale = subst_scales( transport->n_substances, sscales );
	}

	transport->reaction = NULL;
	transport->n_reaction = 0;

	transport->sub_problem = 0;
    if(transport->dual_porosity == true)
    	transport->sub_problem += 1;
    if(transport->sorption == true)
    	transport->sub_problem += 2;

    INPUT_CHECK(!(transport->n_substances < 1 ),"Number of substances must be positive\n");
}
//=============================================================================
//
//=============================================================================
char **subst_names( int n_subst, char *line )
{
	char **sn;
	int sbi;

	ASSERT(!( (n_subst < 1) || (line == NULL) ),"Bad parameter of function subst_names()\n");
	sn = (char**) xmalloc( n_subst * sizeof( char* ) );
	for( sbi = 0; sbi < n_subst; sbi++ )
		sn[ sbi ] = xstrcpy( strtok( (sbi == 0 ? line : NULL) , " \t,;" ) );
	return sn;
}
//=============================================================================
//
//=============================================================================
double *subst_scales( int n_subst, char *line )
{
    double *ss;
    int sbi;

    ASSERT(!( (n_subst < 1) || (line == NULL) ),"Bad parameter of the function subst_scales()\n");
    ss = (double*) xmalloc( n_subst * sizeof( double ) );
    for( sbi = 0; sbi < n_subst; sbi++ )
        ss[ sbi ] = atof( strtok( sbi == 0 ? line : NULL , " \t,;" ) );
    return ss;
}
//=============================================================================
// ALLOCATE TRANSPORT STRUCTURES
//=============================================================================
void transport_vectors_init(struct Transport *transport)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

//    ElementIter elm;
//    int s;
	int i,sbi,n_subst = transport->n_substances;


	for(sbi=0;sbi<n_subst;sbi++){
	  i = 0;
	  FOR_ELEMENTS(elm){
	                transport->conc[sbi][MOBILE][i] = elm->start_conc->conc[sbi]; // = 0;
	                transport->pconc[sbi][MOBILE][i++] = elm->start_conc->conc[sbi];
	        }

	        FOR_ELEMENTS(elm)
	                FOR_ELEMENT_SIDES(elm,s)
	                        if(elm->side[s]->cond != NULL){
	                                transport->conc[sbi][MOBILE][i] = elm->side[s]->cond->transport_bcd->conc[sbi] ;
	                                transport->pconc[sbi][MOBILE][i++] = elm->side[s]->cond->transport_bcd->conc[sbi];
	                        }
	}
}
//=============================================================================
// ALLOCATE TRANSPORT STRUCTURES
//=============================================================================
void alloc_transport_matrix(struct Transport *transport)
{
  int *n;
  n = max_entry();

      Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

  transport->tmatrix = (struct TMatrix*)xmalloc(sizeof(struct TMatrix));
  transport->tmatrix->lenrow = (int*)xmalloc(mesh->n_elements() * sizeof(int));
  transport->tmatrix->irowst = (int*)xmalloc(mesh->n_elements() * sizeof(int));
  transport->tmatrix->jcn = (int*)xmalloc(n[0] * sizeof(int));
  transport->tmatrix->val = (double*)xmalloc(n[0] * sizeof(double));
  xfree(n);
}
//=============================================================================
//	ALLOCATE OF TRANSPORT VARIABLES (ELEMENT & NODES)
//=============================================================================
void alloc_transport_vectors(struct Transport *transport)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    int i,j,sbi,n_subst,ph;
  ElementIter elm;
  n_subst = transport->n_substances;

  transport->conc = (double***)xmalloc(n_subst * sizeof(double**));
  transport->pconc = (double***)xmalloc(n_subst * sizeof(double**));
  transport->node_conc = (double****)xmalloc(n_subst * sizeof(double***));

  for(sbi=0;sbi<n_subst;sbi++){

	  transport->conc[sbi] = (double**)xmalloc(MAX_PHASES * sizeof(double*));
	  transport->pconc[sbi] = (double**)xmalloc(MAX_PHASES * sizeof(double*));
	  transport->node_conc[sbi] = (double***)xmalloc(MAX_PHASES * sizeof(double**));



      for(ph = 0;ph < 4;ph++)
		  if((transport->sub_problem & ph) == ph){

			  if(ph == 0){
			      transport->conc[sbi][ph] = (double*)xmalloc((mesh->n_elements() + mesh->n_boundaries()) * sizeof(double));
			      transport->pconc[sbi][ph] = (double*)xmalloc((mesh->n_elements() + mesh->n_boundaries()) * sizeof(double));
			  }
			  else{
				  transport->conc[sbi][ph] = (double*)xmalloc((mesh->n_elements() ) * sizeof(double));
				  transport->pconc[sbi][ph] = (double*)xmalloc((mesh->n_elements() ) * sizeof(double));

			  }

			  transport->node_conc[sbi][ph] = (double**)xmalloc((mesh->n_elements() ) * sizeof(double*));



			  for(i=0;i<mesh->n_elements();i++){
				  transport->conc[sbi][ph][i] = 0.0;
				  transport->pconc[sbi][ph][i] = 0.0;
			  }

			  i = 0;
			  FOR_ELEMENTS(elm){
					transport->node_conc[sbi][ph][i++]=(double*)xmalloc((elm->n_nodes) * sizeof(double));
					for(j = 0 ;j < elm->n_nodes ; j++)
						transport->node_conc[sbi][ph][i-1][j] = 0.0;
			  }

		  }
		  else{
			  transport->conc[sbi][ph] = NULL;
			  transport->pconc[sbi][ph] = NULL;
			  transport->node_conc[sbi][ph] = NULL;
		  }
  }
}
//=============================================================================
//	ALLOCATE OF TRANSPORT VARIABLES (ELEMENT & NODES)
//=============================================================================
void alloc_density_vectors(struct Transport *transport)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    int ph,sbi,i,sub;
	int n_subst = transport->n_substances;
	int n_elements = mesh->n_elements();

	sub = transport->sub_problem;

	transport->scalar_it = (double*)xmalloc(n_elements * sizeof(double)); // Zatim nevyuzito
	transport->prev_conc = (double***)xmalloc(n_subst * sizeof(double**));

	 for(sbi=0;sbi<n_subst;sbi++){

		  transport->prev_conc[sbi] = (double**)xmalloc(MAX_PHASES * sizeof(double*));

	      for(ph = 0;ph < 4;ph++)
			  if((sub & ph) == ph){
				  transport->prev_conc[sbi][ph] = (double*)xmalloc(n_elements * sizeof(double));

				  for(i=0;i<n_elements;i++)
					  transport->prev_conc[sbi][ph][i] = 0.0;

			  }
			  else
				  transport->prev_conc[sbi][ph] = NULL;

	 }
}
//=============================================================================
// CREATE TRANSPORT MATRIX
//=============================================================================
void create_transport_matrix(struct Problem *problem)
{
  Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);
  
  ElementIter elm;
  struct Edge *edg;
  struct Transport *transport;
  struct SVector *svector;
  struct TMatrix *tmatrix;
  int i,j,k,l,n,s,itemp,ec,*max_ent;
  double ftemp,max_sum,sum,flux;

  max_sum = 0.0;

  transport = problem->transport;
  tmatrix = problem->transport->tmatrix;


  FOR_EDGES(edg){                                         // calculate edge Qv
        edg->faux = 0;
        FOR_EDGE_SIDES(edg,s)
                if(edg->side[s]->flux > 0)
                        edg->faux += edg->side[s]->flux;
  }

  max_ent = max_entry();

  svector = (struct SVector*)xmalloc(sizeof(struct SVector)); //uvolnit
  svector->pos = (int*)xmalloc(max_ent[1] * sizeof(int));
  svector->val = (double*)xmalloc(max_ent[1] * sizeof(double));

  transport->tmatrix->rows = mesh->n_elements();

  i = -1;
  ec = 0;
  FOR_ELEMENTS(elm){
    i++; // pos
    for(k=0;k<max_ent[1];k++){
        svector->pos[k] = 0;
        svector->val[k] = 0.0;
    }
    svector->pos[0] = i;
    j = 1;
      // prislusne porozity!!!!!!!

    FOR_ELEMENT_SIDES(elm,si)            // same dim
      if(elm->side[si]->cond == NULL){
        if(elm->side[si]->flux < 0.0 ){

            if(elm->side[si]->neigh_bv != NULL){         //comp model
                svector->val[j] = -(elm->side[si]->flux / (elm->volume * elm->material->por_m));
                //svector->pos[j++] = id2pos(mesh,ELEMENT_FULL_ITER(elm->side[si]->neigh_bv->element[0]).id(),mesh->epos_id,ELM);
                svector->pos[j++] = ELEMENT_FULL_ITER(elm->side[si]->neigh_bv->element[0]) - mesh->element.begin();
        //        printf("side in elm:%d value:%f pos:%d\n ",elm->id,svector->val[j-1],svector->pos[j-1]);
        //        getchar();
            }
                                                                   // end comp model
            else{
                edg = elm->side[si]->edge;
                if(edg->faux > ZERO)
                        FOR_EDGE_SIDES(edg,s)
                                if((edg->side[s]->id != elm->side[si]->id) && (edg->side[s]->flux > 0.0)){
                                        svector->val[j] = -(elm->side[si]->flux *  edg->side[s]->flux / (edg->faux * elm->volume * elm->material->por_m));
                                        //svector->pos[j++] = id2pos(mesh,ELEMENT_FULL_ITER().id(),mesh->epos_id,ELM);
                                        svector->pos[j++] = ELEMENT_FULL_ITER(edg->side[s]->element) - mesh->element.begin();
                                    //    printf("side in elm:%d value:%f\n ",elm->id,svector->val[j-1]);
                                }
            }
        }
        if(elm->side[si]->flux > 0.0){
                svector->val[0] -= (elm->side[si]->flux / (elm->volume * elm->material->por_m));
            //    printf("side in elm:%d value:%f\n ",elm->id,svector->val[0]);
        }
      }
      else{
        if(elm->side[si]->flux < 0.0){
                svector->val[j] = -(elm->side[si]->flux / (elm->volume * elm->material->por_m));
               // svector->pos[j++] = id2pos(mesh,elm->side[si]->id,mesh->spos_id,BC);
                svector->pos[j++] = BOUNDARY_FULL_ITER(elm->side[si]->cond) - mesh->boundary.begin() + mesh->n_elements();
             //   printf("side in elm:%d value:%f\n ",elm->id,svector->val[j-1]);
             //   printf("%d\t%d\n",elm->id,id2pos(problem,elm->side[si]->id,problem->spos_id,BC));

        }
        if(elm->side[si]->flux > 0.0){
                svector->val[0] -= (elm->side[si]->flux / (elm->volume * elm->material->por_m));
            //    printf("side in elm:%d value:%f\n ",elm->id,svector->val[0]);
        }
      }         // end same dim     //ELEMENT_SIDES


    FOR_ELM_NEIGHS_VB(elm,n)                    // comp model
      FOR_NEIGH_ELEMENTS(elm->neigh_vb[n],s)
        if (elm.id() != ELEMENT_FULL_ITER(elm->neigh_vb[n]->element[s]).id()){
          if (elm->neigh_vb[n]->side[s]->flux > 0.0){
                svector->val[j] = elm->neigh_vb[n]->side[s]->flux / (elm->volume * elm->material->por_m);
              //  svector->pos[j++] = id2pos(mesh,ELEMENT_FULL_ITER(elm->neigh_vb[n]->element[s]).id(),mesh->epos_id,ELM);
                svector->pos[j++] = ELEMENT_FULL_ITER(elm->neigh_vb[n]->element[s]) - mesh->element.begin();
          }
          if (elm->neigh_vb[n]->side[s]->flux < 0.0){
                svector->val[0] += elm->neigh_vb[n]->side[s]->flux / (elm->volume * elm->material->por_m);
          }
       }                             // end comp model


  Neighbour *ngh;
  ElementFullIter el2 = ELEMENT_FULL_ITER_NULL;
  FOR_ELM_NEIGHS_VV(elm,n) {                     //non-comp model
    ngh=elm->neigh_vv[n];
    FOR_NEIGH_ELEMENTS(ngh,s) {

      el2=ELEMENT_FULL_ITER(ngh->element[s]);
      if (elm.id()!=el2.id())
      {
        flux = ngh->sigma*ngh->geom_factor * (el2->scalar - elm->scalar);
        if (flux > 0.0){
           svector->val[j] = flux / (elm->volume * elm->material->por_m);         // -=
           svector->pos[j++] = ELEMENT_FULL_ITER(el2) - mesh->element.begin(); //id2pos(mesh,el2.id(),mesh->epos_id,ELM);
        }
        if (flux < 0.0){
          svector->val[0] += flux / (elm->volume * elm->material->por_m);
        }
      }
    }
  } // end non-comp model



    if(fabs(svector->val[0]) > max_sum)
        max_sum = fabs(svector->val[0]);


    transport->tmatrix->lenrow[i] = j;

    for(k=0;k<tmatrix->lenrow[i]-1;k++)
        for(l=0;l<transport->tmatrix->lenrow[i]-1;l++)
        if(svector->pos[l] > svector->pos[l+1]){
          itemp = svector->pos[l+1];
          ftemp = svector->val[l+1];
          svector->pos[l+1] = svector->pos[l];
          svector->val[l+1] = svector->val[l];
          svector->pos[l] = itemp;
          svector->val[l] = ftemp;
        }

    tmatrix->irowst[i]=ec;

    sum = 0.0;
    for(k = ec; k < ec + transport->tmatrix->lenrow[i]; k++){
        transport->tmatrix->jcn[k] = svector->pos[k-ec];
        transport->tmatrix->val[k] = svector->val[k-ec];
    }

    ec += j;          // entry count


   }     // ELEMENTS

   xfree(max_ent);

   transport->max_step = 1 / max_sum;
   transport->time_step = 0.9 / max_sum;


   xfree(svector->pos);
   xfree(svector->val);
   xfree(svector);
   /*
   printf("\n");
   k = 0;
   for(i=0;i<problem->mesh->n_elements();i++){
        for(j=0;j<tmatrix->lenrow[i] ;j++){
                printf("%d\t%d\t%f\n",i,tmatrix->jcn[k],tmatrix->val[k]);
                k++;
        }
   }


   printf("\nEntry count: %d\n",ec);
   */
}
//=============================================================================
//	ALLOCATION OF TRANSPORT VECTORS (MPI)
//=============================================================================
void alloc_transport_structs_mpi(struct Transport *transport)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    int i,j,sbi,n_subst,ph,ierr,rank,np;
  ElementIter elm;
  n_subst = transport->n_substances;

	MPI_Barrier(PETSC_COMM_WORLD);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &np);

	transport_partioning(transport,np);

  transport->bcv=(Vec*)xmalloc(n_subst * (sizeof (Vec)));
  transport->bcvcorr=(Vec*)xmalloc(n_subst * (sizeof (Vec)));
  transport->vconc=(Vec*)xmalloc(n_subst * (sizeof (Vec)));
  transport->vpconc=(Vec*)xmalloc(n_subst * (sizeof (Vec)));

 // if( rank == 0)
	  transport->vconc_out=(Vec*)xmalloc(n_subst * (sizeof (Vec))); // extend to all


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

  for(sbi=0;sbi<n_subst;sbi++){
	  //ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->lb_col[rank],mesh->n_boundaries(),&transport->bcv[sbi]);
	  //ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh->n_elements(),&transport->bcvcorr[sbi]);
	  //ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh->n_elements(),&transport->vconc[sbi]);
	  //ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh->n_elements(),&transport->vpconc[sbi]);

	  ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->lb_col[rank],mesh->n_boundaries(),&transport->bcv[sbi]);
	  ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->el_ds->lsize(),mesh->n_elements(),&transport->bcvcorr[sbi]);
	  VecDuplicate(transport->bcvcorr[sbi],&transport->vconc[sbi]);
	  VecDuplicate(transport->bcvcorr[sbi],&transport->vpconc[sbi]);


	//  if(rank == 0)
		  ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,mesh->n_elements(),transport->conc[sbi][MOBILE],&transport->vconc_out[sbi]);

		  //ierr = VecCreateMPI(PETSC_COMM_SELF ,transport->mesh->n_elements(),transport->mesh->n_elements(),&transport->vconc_out[sbi]); /*xx*/

	  if(transport->vconc_im != NULL)
		  ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh->n_elements(),&transport->vconc_im[sbi]);
	  if(transport->vconc_so != NULL)
		  ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh->n_elements(),&transport->vconc_so[sbi]);
	  if(transport->vconc_im_so != NULL)
		  ierr = VecCreateMPI(PETSC_COMM_WORLD,transport->l_row[rank],mesh->n_elements(),&transport->vconc_im_so[sbi]);
  }
  //
	//ierr=MatCreateMPIAIJ(PETSC_COMM_WORLD,transport->l_row[rank],transport->l_row[rank],mesh->n_elements(),mesh->n_elements(),
	//		PETSC_NULL,transport->d_row[rank],PETSC_NULL,transport->od_row[rank],&transport->tm);

	//ierr=MatCreateMPIAIJ(PETSC_COMM_WORLD,transport->l_row[rank],transport->lb_col[rank],mesh->n_elements(),mesh->n_boundaries(),
	//				PETSC_NULL,transport->db_row[rank],PETSC_NULL,transport->odb_row[rank],&transport->bcm);

	ierr=MatCreateMPIAIJ(PETSC_COMM_WORLD,transport->el_ds->lsize(),transport->el_ds->lsize(),mesh->n_elements(),mesh->n_elements(),
			8,PETSC_NULL,1,PETSC_NULL,&transport->tm);

	ierr=MatCreateMPIAIJ(PETSC_COMM_WORLD,transport->el_ds->lsize(),transport->lb_col[rank],mesh->n_elements(),mesh->n_boundaries(),
					2,PETSC_NULL,0,PETSC_NULL,&transport->bcm);





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
void fill_transport_vectors_mpi(struct Transport *transport)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    int i,new_i,sbi,ph,id,n_subst,s,k,rank;
  ElementFullIter elm = ELEMENT_FULL_ITER(NULL);
  BoundaryIter bc;
  double start_conc;
  n_subst = transport->n_substances;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


  for(sbi=0;sbi<n_subst;sbi++){
		for(int loc_el=0;loc_el<transport->el_ds->lsize();loc_el++) {
		    elm=mesh->element(transport->el_4_loc[loc_el]);


		start_conc = elm->start_conc->conc[sbi];
		i = elm.index();
		new_i =  transport->row_4_el[i];
		VecSetValue(transport->vconc[sbi],new_i,start_conc,INSERT_VALUES);
		VecSetValue(transport->vpconc[sbi],new_i,start_conc,INSERT_VALUES);
	}
		k = 0;
		FOR_BOUNDARIES(bc)
			VecSetValue(transport->bcv[sbi],k++,bc->transport_bcd->conc[sbi],INSERT_VALUES);

				 // VecDuplicate(transport->vconc[sbi][ph],&transport->vpconc[sbi][ph]);

	VecZeroEntries(transport->bcvcorr[sbi]);
//	if( rank == 0)
//		VecZeroEntries(transport->vconc_out[sbi]); /*xx*/
	VecAssemblyBegin(transport->bcvcorr[sbi]);
//	if( rank == 0)
//		VecAssemblyBegin(transport->vconc_out[sbi]); /*xx*/
	VecAssemblyBegin(transport->vconc[sbi]);
	VecAssemblyBegin(transport->vpconc[sbi]);
	VecAssemblyBegin(transport->bcv[sbi]);

	if(transport->vconc_im != NULL)
		VecZeroEntries(transport->vconc_im[sbi]);
	if(transport->vconc_so != NULL)
		VecZeroEntries(transport->vconc_so[sbi]);
	if(transport->vconc_im_so != NULL)
		VecZeroEntries(transport->vconc_im_so[sbi]);
  }

  for(sbi=0;sbi<n_subst;sbi++){
//	  if( rank == 0)
//		  VecAssemblyEnd(transport->vconc_out[sbi]); /*xx*/
  	  VecAssemblyEnd(transport->vconc[sbi]);
  	  VecAssemblyEnd(transport->vpconc[sbi]);
  	  VecAssemblyEnd(transport->bcv[sbi]);
  	  VecAssemblyEnd(transport->bcvcorr[sbi]);

	  if(transport->vconc_im != NULL)
		  VecAssemblyEnd(transport->vconc_im[sbi]);
	  if(transport->vconc_so != NULL)
		  VecAssemblyEnd(transport->vconc_so[sbi]);
	  if(transport->vconc_im_so != NULL)
		  VecAssemblyEnd(transport->vconc_im_so[sbi]);
  }
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
void create_transport_matrix_mpi(struct Transport *transport)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

  ElementFullIter el2 = ELEMENT_FULL_ITER_NULL;
  ElementFullIter elm = ELEMENT_FULL_ITER_NULL;
  struct Edge *edg;
  struct Neighbour *ngh;
  //struct Transport *transport;
  int n,s,i,j,np,rank,new_j,new_i;
  double max_sum,flux,aij,aii;

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

  FOR_EDGES(edg){                                         // calculate edge Qv
        edg->faux = 0;
        FOR_EDGE_SIDES(edg,s)
                if(edg->side[s]->flux > 0)
                        edg->faux += edg->side[s]->flux;
  }

  max_sum = 0.0;
  aii = 0.0;
  for(int loc_el=0;loc_el<transport->el_ds->lsize();loc_el++) {
      elm=mesh->element(transport->el_4_loc[loc_el]);
      new_i =  transport->row_4_el[elm.index()];



    FOR_ELEMENT_SIDES(elm,si)            // same dim
      if(elm->side[si]->cond == NULL){
        if(elm->side[si]->flux < 0.0 ){
            if(elm->side[si]->neigh_bv != NULL){         //comp model
                aij = -(elm->side[si]->flux / (elm->volume * elm->material->por_m));
                j = ELEMENT_FULL_ITER(elm->side[si]->neigh_bv->element[0]).index();
                new_j = transport->row_4_el[j];
                MatSetValue (transport->tm,new_i,new_j,aij,INSERT_VALUES);
            }
                                                                   // end comp model
            else{
                edg = elm->side[si]->edge;
                if(edg->faux > ZERO)
                        FOR_EDGE_SIDES(edg,s)
                                if((edg->side[s]->id != elm->side[si]->id) && (edg->side[s]->flux > 0.0)){
                                        aij = -(elm->side[si]->flux *  edg->side[s]->flux / (edg->faux * elm->volume * elm->material->por_m));
                                        j = ELEMENT_FULL_ITER(edg->side[s]->element).index();
                                        new_j = transport->row_4_el[j];
                                        MatSetValue (transport->tm,new_i,new_j,aij,INSERT_VALUES);
                                }
            }
        }
        if(elm->side[si]->flux > 0.0)
                aii -= (elm->side[si]->flux / (elm->volume * elm->material->por_m));
      }
      else{
        if(elm->side[si]->flux < 0.0){
                aij = -(elm->side[si]->flux / (elm->volume * elm->material->por_m));
                j = BOUNDARY_FULL_ITER(elm->side[si]->cond).index();
                MatSetValue (transport->bcm,new_i,j,aij,INSERT_VALUES);
                // vyresit BC matrix !!!!
             //   printf("side in elm:%d value:%f\n ",elm->id,svector->val[j-1]);
             //   printf("%d\t%d\n",elm->id,id2pos(problem,elm->side[si]->id,problem->spos_id,BC));

        }
        if(elm->side[si]->flux > 0.0)
                aii -= (elm->side[si]->flux / (elm->volume * elm->material->por_m));
      }         // end same dim     //ELEMENT_SIDES


    FOR_ELM_NEIGHS_VB(elm,n)                    // comp model
      FOR_NEIGH_ELEMENTS(elm->neigh_vb[n],s)
        if (elm.id() != ELEMENT_FULL_ITER(elm->neigh_vb[n]->element[s]).id()){
          if (elm->neigh_vb[n]->side[s]->flux > 0.0){
                aij = elm->neigh_vb[n]->side[s]->flux / (elm->volume * elm->material->por_m);
                j = ELEMENT_FULL_ITER(elm->neigh_vb[n]->element[s]).index();
                new_j = transport->row_4_el[j];
                MatSetValue (transport->tm,new_i,new_j,aij,INSERT_VALUES);
          }
          if (elm->neigh_vb[n]->side[s]->flux < 0.0)
                aii += elm->neigh_vb[n]->side[s]->flux / (elm->volume * elm->material->por_m);
       }                             // end comp model
  FOR_ELM_NEIGHS_VV(elm,n) {                     //non-comp model
    ngh=elm->neigh_vv[n];
    FOR_NEIGH_ELEMENTS(ngh,s) {

      el2=ELEMENT_FULL_ITER(ngh->element[s]);
      if (elm.id()!=el2.id())
      {
        flux = ngh->sigma*ngh->geom_factor * (el2->scalar - elm->scalar);
        if (flux > 0.0){
           aij = flux / (elm->volume * elm->material->por_m);         // -=
           j = el2.index();
           new_j = transport->row_4_el[j];
           MatSetValue (transport->tm,new_i,new_j,aij,INSERT_VALUES);
        }
        if (flux < 0.0)
          aii += flux / (elm->volume * elm->material->por_m);
      }
    }
  } // end non-comp model

	MatSetValue (transport->tm,new_i,new_i,aii,INSERT_VALUES);

    if(fabs(aii) > max_sum)
        max_sum = fabs(aii);
    aii = 0.0;
 //   i++;
   }     // END ELEMENTS

  double glob_max_sum;

  MPI_Allreduce(&max_sum,&glob_max_sum,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);

  MatAssemblyBegin(transport->tm,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(transport->bcm,MAT_FINAL_ASSEMBLY);

  transport->max_step = 1 / glob_max_sum;
  transport->time_step = 0.9 / glob_max_sum;

  MatAssemblyEnd(transport->tm,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(transport->bcm,MAT_FINAL_ASSEMBLY);






 // MPI_Barrier(PETSC_COMM_WORLD);
  /*
  MatView(transport->tm,PETSC_VIEWER_STDOUT_SELF);
  getchar();
  */

 }
//=============================================================================
// CREATE TRANSPORT MATRIX STEP MPI
//=============================================================================
void transport_matrix_step_mpi(struct Transport *transport,double time_step)
{
	  MatScale(transport->bcm,time_step);
	  MatScale(transport->tm,time_step);
	  MatShift(transport->tm,1.0);

/*
		MatView(transport->bcm,PETSC_VIEWER_STDOUT_WORLD);
		getchar();
*/

}
//=============================================================================
// MATRIX VECTOR PRODUCT (MPI)
//=============================================================================
void calculate_bc_mpi(struct Transport *transport)
{
	int sbi,n_subst;

	n_subst = transport->n_substances;

	for(sbi = 0;sbi < n_subst;sbi++)
		MatMult(transport->bcm,transport->bcv[sbi],transport->bcvcorr[sbi]);
/*
	  VecView(transport->bcvcorr[0],PETSC_VIEWER_STDOUT_SELF);
	  getchar();
*/
}
//=============================================================================
// MATRIX VECTOR PRODUCT (MPI)
//=============================================================================
void transport_step_mpi(Mat *tm,Vec *conc,Vec *pconc,Vec *bc)
{



   MatMultAdd(*tm,*pconc,*bc,*conc);  // conc=tm*pconc + bc
   VecSwap(*conc,*pconc); // pconc = conc


}
//=============================================================================
// CREATE TRANSPORT MATRIX STEP
//=============================================================================
void transport_matrix_step(struct TMatrix *tmatrix,double time_step)
{
  int i,j,k;

  k=0;
  for(i = 0; i < tmatrix->rows; i++)
     for(j = 0;j<tmatrix->lenrow[i];j++){
        tmatrix->val[k++] *= time_step;
        if(tmatrix->jcn[k-1] == i)
                tmatrix->val[k-1] += 1;
  }
   /*
   printf("\n");
   k = 0;
   for(i=0;i<problem->mesh->n_elements();i++){
        for(j=0;j<tmatrix->lenrow[i] ;j++){
                printf("%d\t%d\t%f\n",i,tmatrix->jcn[k],tmatrix->val[k]);
                k++;
        }
   }
   getchar();
    */

}
//=============================================================================
// MATRIX VECTOR PRODUCT
//=============================================================================
void matvecs(struct TMatrix *tmatrix,double *pconc,double *conc)
{
   double temp;
   int i,j,k;

  // for(k=0;k<tmatrix->rows;k++)
  //    printf("%f\n",conc[k]);

   //for(m = 0; m < n; m++){
   k=0;
   //printf("%d\n",m);
        for(i=0;i<tmatrix->rows;i++){ // row
                temp = 0;
                for(j=k;j<(tmatrix->lenrow[i]+k);j++){
                //        getchar();
               //         printf("%d\t%d\t%f\t%f\t%d\t%f\n",i,j,transport->conc[i],transport->pconc[i],tmatrix->jcn[j], transport->pconc[tmatrix->jcn[j]]);
                        temp += tmatrix->val[j] * pconc[tmatrix->jcn[j]];
              //          printf("\n%d\t%f",i,transport->conc[i]);
                }
                conc[i] = temp;
                k += tmatrix->lenrow[i];
                //getchar();
        }
        memcpy(pconc,conc,tmatrix->rows * sizeof(double));     //!!!
   //}
  // getchar();
}
//=============================================================================
// ALLOCATING MEMORY FOR TRANSPORT VARIABLE PI
//=============================================================================
double *transport_aloc_pi(Mesh* mesh)
{
  int max_elm;
  double *pi;

  NodeIter nod;

  max_elm = 0;
  FOR_NODES(nod)
    if (max_elm < nod->n_elements) max_elm = nod->n_elements;
  pi = (double *) xmalloc ( max_elm * sizeof( double ));

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
//=============================================================================
void transport_dual_porosity(struct Transport *transport, int elm_pos, int sbi)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    double conc_avg = 0.0;
   MaterialDatabase::Iter material;
   int id;
   double ***conc=transport->conc,***pconc=transport->pconc,cm,pcm,ci,pci,por_m,por_imm,alpha;

  // TODO: remove epos_id and find_id
//  DBGMSG("el pos: %d\n",elm_pos);
  id = mesh->epos_id[elm_pos];
  material =  mesh->element.find_id(id)->material;

  por_m = material->por_m;
  por_imm = material->por_imm;
  alpha = material->alpha[sbi];
  pcm = pconc[sbi][MOBILE][elm_pos];
  pci = pconc[sbi][IMMOBILE][elm_pos];

  // ---compute average concentration------------------------------------------
  conc_avg = (( por_m * pcm )
           + (por_imm * pci ))
           / ( por_m + por_imm );

  if((conc_avg != 0.0) && (por_imm != 0.0))
  {
  // ---compute concentration in mobile area-----------------------------------
  cm = ( pcm - conc_avg )
                 * exp( - alpha * ((por_m + por_imm)
                 / (por_m * por_imm)) * transport->time_step )
                 + conc_avg;

  // ---compute concentration in immobile area---------------------------------
  ci = ( pci - conc_avg )
                 * exp( - alpha * ((por_m + por_imm)
                 / (por_m * por_imm)) * transport->time_step )
                 + conc_avg;
  // --------------------------------------------------------------------------
  //printf("\n%f\t%f\t%f",conc_avg,cm,ci);
  //getchar();

  conc[sbi][MOBILE][elm_pos] = cm;
  pconc[sbi][MOBILE][elm_pos] = cm;
  conc[sbi][IMMOBILE][elm_pos] = ci;
  pconc[sbi][IMMOBILE][elm_pos] = ci;
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
void transport_sorption(struct Transport *transport, int elm_pos, int sbi)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    double conc_avg = 0.0;
   double conc_avg_imm = 0.0;
   double n, Nm, Nimm;
   MaterialDatabase::Iter mtr;
   int id;
   double phi;
   double ***conc=transport->conc,***pconc=transport->pconc;

   // TODO: remove epos_id and find_id
   id = mesh->epos_id[elm_pos];
   mtr =  mesh->element.find_id(id)->material;



   phi = mtr->phi;

   if((mtr->sorp_coef[sbi][0] == 0) || (mtr->por_m == 1)) return;

   n = 1-(mtr->por_m + mtr->por_imm);
   Nm = mtr->por_m;
   Nimm = mtr->por_imm;

   conc_avg = pconc[sbi][MOBILE][elm_pos] + pconc[sbi][MOBILE_SORB][elm_pos]  * n/Nm; // cela hmota do poru



   if(conc_avg !=0)
   {
   compute_sorption(conc_avg,mtr->sorp_coef[sbi],mtr->sorp_type[sbi],
           &conc[sbi][MOBILE][elm_pos],&conc[sbi][MOBILE_SORB][elm_pos], Nm/n , n * phi / Nm);

   pconc[sbi][MOBILE][elm_pos] = conc[sbi][MOBILE][elm_pos];
   pconc[sbi][MOBILE_SORB][elm_pos] = conc[sbi][MOBILE_SORB][elm_pos];
   }
   //printf("\n%f\t%f\t",n * phi / Nm,n * phi / Nm);
   //printf("\n%f\t%f\t",n * phi / Nimm,n * (1 - phi) / Nimm);
  // getchar();

   if( (transport->dual_porosity == true) && (mtr->por_imm != 0) )
   {
   conc_avg_imm = pconc[sbi][IMMOBILE][elm_pos] + pconc[sbi][IMMOBILE_SORB][elm_pos]  * n/Nimm; // cela hmota do poru

   if(conc_avg_imm !=0)
        {
        compute_sorption(conc_avg_imm,mtr->sorp_coef[sbi],mtr->sorp_type[sbi],
        &conc[sbi][IMMOBILE][elm_pos],&conc[sbi][IMMOBILE_SORB][elm_pos], Nimm/n , n * (1 - phi) / Nimm);

        pconc[sbi][IMMOBILE][elm_pos] = conc[sbi][IMMOBILE][elm_pos];
        pconc[sbi][IMMOBILE_SORB][elm_pos] = conc[sbi][IMMOBILE_SORB][elm_pos];
        }
   }

}
//=============================================================================
//      COMPUTE SORPTION
//=============================================================================
void compute_sorption(double conc_avg,vector<double> &sorp_coef,int sorp_type,double *concx,double *concx_sorb,
double Nv, double N)
{
double Kx = sorp_coef[0] * N;
double parameter;// = sorp_coef[1];
double NR,pNR,cz,tcz;
//double lZero = 0.0000001;
double ad = 1e4;
int i;

        pNR = 0;

//if(conc_avg > 1e-20)
   switch(sorp_type){
        case 1:         //linear
                *concx = conc_avg / (1 + Kx);
            //    *concx_sorb = (conc_avg - *concx) * Nv;   // s = Kd *c  [kg/m^3]
                break;
        case 2:         //freundlich
                parameter = sorp_coef[1];
                cz = pow(ad /(Kx * parameter),1 / (parameter - 1));
                tcz = ad / parameter;
                NR = cz;
                for(i=0; i<20; i++)  //Newton Raphson iteration cycle
                {
                NR -= ( NR +  ((NR > cz ) ? Kx * pow(NR, parameter) : tcz * NR)  - conc_avg) /
               (1 + ((NR > cz ) ?  parameter * Kx * pow(NR, parameter - 1) : tcz ));
                if ((NR <= cz) || (fabs(NR - pNR) < ZERO ))break;
                pNR = NR;
//                if(NR < 0) printf("\n%f\t",NR);
                }
                *concx = NR;
                break;
        case 3:         // langmuir
                parameter = sorp_coef[1];
                NR = 0;
                for(i=0; i<5;i++)  //Newton Raphson iteration cycle
                {
                NR -= (NR + (NR*Kx*parameter)/(1+NR*Kx) - conc_avg) /
                        (1 + Kx*parameter/pow(1+NR*Kx,2) );
                if (fabs(NR - pNR) < ZERO ) break;
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
//      CONVECTION
//=============================================================================
void convection(struct Transport *trans)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    int steps,step,save_step,frame = 0;
  register int t;
  int n_subst,sbi,elm_pos,rank;
  struct Problem *problem = trans->problem;
  struct TMatrix *tmatrix=problem->transport->tmatrix;
  double ***pconc;
  double ***conc;

  START_TIMER("TRANSPORT");
  // int tst = 1; // DECOVALEX

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);



  xprintf( Msg, "Calculating transport...")/*orig verb 2*/;
  n_subst = mesh->n_substances;

//  flow_cs(trans); //DECOVALEX


  if(trans->mpi != 1){
	  create_transport_matrix(problem);
  }
  else
  create_transport_matrix_mpi(trans);

 // MatView(trans->tm,PETSC_VIEWER_STDOUT_WORLD);


//  if(problem->type != PROBLEM_DENSITY){
	  save_step = (int)ceil(problem->save_step / trans->time_step);      // transport  rev
	  trans->time_step = problem->save_step / save_step;
	  steps = save_step * (int)floor(problem->stop_time / problem->save_step);
/*  }
 else{
	  steps = (int)floor(problem->transport->update_dens_time / problem->transport->time_step);
	  save_step = steps + 1;
  }*/

  if(trans->mpi != 1)
	  transport_matrix_step(tmatrix,trans->time_step);
  else{
	  transport_matrix_step_mpi(trans,trans->time_step);
	  calculate_bc_mpi(trans);
	}

  if(rank == 0){
	  printf("  %d computing cycles, %d writing steps\n",steps ,((int)(steps / save_step) + 1) )/*orig verb 6*/;
	  xprintf( MsgVerb, "  %d computing cycles, %d writing steps\n",steps ,((int)(steps / save_step) + 1) )/*orig verb 6*/;
	}
  pconc = trans->pconc;
  conc = trans->conc;

   //output_FCS(trans); //DECOVALEX


  for(t = 1; t <= steps; t++)
  {
    step ++;
    for(sbi=0;sbi<n_subst;sbi++)
    {
/*
      if(tst && (particle_test(trans) > 100.0) && tst){ // DECOVALEX
        	  		clear_tbc(trans);
        	  		tst = 0;
      }
      output_AGE(trans,(t-1) * trans->time_step); // DECOVALEX
 */



  if(trans->mpi != 1)
    matvecs(trans->tmatrix,trans->pconc[sbi][MOBILE],trans->conc[sbi][MOBILE]);
  else
    transport_step_mpi(&trans->tm,&trans->vconc[sbi],&trans->vpconc[sbi],&trans->bcvcorr[sbi]);

      if((trans->dual_porosity == true) || (trans->sorption == true) || (trans->pepa == true) || (trans->reaction_on == true))
    	  for(elm_pos = 0;elm_pos < mesh->n_elements();elm_pos++){
    		  if(trans->dual_porosity == true)
    			  transport_dual_porosity(trans,elm_pos,sbi);
    		  if(trans->sorption == true)
    			  transport_sorption(trans,elm_pos,sbi);
    		  if(trans->pepa == true)
    			  decay(trans,elm_pos,trans->type); // pepa chudoba
    	      if(trans->reaction_on == true)
    	    	  transport_reaction(trans,elm_pos,sbi);

           }
     // transport_node_conc(mesh,sbi,problem->transport_sub_problem);  // vyresit prepocet
    }
    xprintf( Msg, "Time : %f\n",trans->time_step*t);

    if ((save_step == step) && ((ConstantDB::getInstance()->getInt("Problem_type") != PROBLEM_DENSITY) || (trans-> write_iterations))){
    	if(trans->mpi == 1)
    		output_vector_gather(trans);

    	if(rank == 0)
    		transport_output(trans, t * trans->time_step,++frame);
        if (ConstantDB::getInstance()->getInt("Problem_type") != STEADY_SATURATED)
        	output_time(problem, t * trans->time_step);    // time variable flow field
	//	output_transport_time_BTC(trans, t * trans->time_step); // BTC test - spatne vypisuje casy
	  //output_transport_time_CS(problem, t * problem->time_step);
      step = 0;
    }
  }
  xprintf( Msg, "O.K.\n");

}

//=============================================================================
//      OUTPUT VECTOR GATHER
//=============================================================================
void output_vector_gather(struct Transport *transport)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    int sbi,rank,np;
	IS is;
	PetscViewer inviewer;

//	MPI_Barrier(PETSC_COMM_WORLD);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &np);

	//if(np >1){
		//ISCreateStride(PETSC_COMM_SELF,mesh->n_elements(),0,1,&is);
		ISCreateGeneral(PETSC_COMM_SELF,mesh->n_elements(),transport->row_4_el,&is); //WithArray
		VecScatterCreate(transport->vconc[0],is,transport->vconc_out[0],PETSC_NULL,&transport->vconc_out_scatter);
		for(sbi = 0; sbi < transport->n_substances; sbi++){
			VecScatterBegin(transport->vconc_out_scatter,transport->vconc[sbi],transport->vconc_out[sbi],INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd(transport->vconc_out_scatter,transport->vconc[sbi],transport->vconc_out[sbi],INSERT_VALUES,SCATTER_FORWARD);
			}
		//VecView(transport->vconc[0],PETSC_VIEWER_STDOUT_WORLD);
		//VecView(transport->vconc_out[0],PETSC_VIEWER_STDOUT_WORLD);
		VecScatterDestroy(transport->vconc_out_scatter);
		ISDestroy(is);
		/*
	}
	else{
		for(sbi = 0; sbi < transport->n_substances; sbi++)
			VecCopy(transport->vconc[sbi],transport->vconc_out[sbi]);
	}*/
}
//=============================================================================
//      TRANSPORT OUTPUT
//=============================================================================
void transport_output(struct Transport *transport,double time,int frame)
{
    switch(ConstantDB::getInstance()->getInt("Pos_format_id")){
    case POS_BIN:
        output_transport_time_bin(transport, time, frame, transport->transport_out_fname);
        break;
    case POS_ASCII:
        output_transport_time_ascii(transport, time, frame, transport->transport_out_fname);
        break;
    case VTK_SERIAL_ASCII:
        output_transport_time_vtk_serial_ascii(transport, time, frame, transport->transport_out_fname);
        break;
    case VTK_PARALLEL_ASCII:
        xprintf(UsrErr, "VTK_PARALLEL_ASCII: not implemented yet\n");
        break;
    }
}
//=============================================================================
//      TRANSPORT OUTPUT INIT
//=============================================================================
void transport_output_init(struct Transport *transport)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    switch(ConstantDB::getInstance()->getInt("Pos_format_id")){
    case POS_BIN:
        output_msh_init_bin(mesh, transport->transport_out_fname);
        break;
    case POS_ASCII:
        output_msh_init_ascii(mesh, transport->transport_out_fname);
        break;
    case VTK_SERIAL_ASCII:
        output_msh_init_vtk_serial_ascii(transport->problem, transport->transport_out_fname);
        break;
    case VTK_PARALLEL_ASCII:
        xprintf(UsrErr, "VTK_PARALLEL_ASCII: not implemented yet\n");
        break;
    }
}

//=============================================================================
//      TRANSPORT OUTPUT FINISH
//=============================================================================
void transport_output_finish(struct Transport *transport)
{
    switch(ConstantDB::getInstance()->getInt("Pos_format_id")){
    case POS_BIN:
        /* There is no need to do anything for this file format */
        break;
    case POS_ASCII:
        /* There is no need to do anything for this file format */
        break;
    case VTK_SERIAL_ASCII:
        output_msh_finish_vtk_serial_ascii(transport->problem, transport->transport_out_fname);
        break;
    case VTK_PARALLEL_ASCII:
        xprintf(UsrErr, "VTK_PARALLEL_ASCII: not implemented yet\n");
        break;
    }
}

//=============================================================================
//      COMPARE DENSITY ITERATION
//=============================================================================
int compare_dens_iter(struct Problem *problem)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    ElementIter elm;
	double max_err;
	max_err=0;
//	FOR_ELEMENTS( elm )
//		xprintf(Msg,"%f %f %f\n",elm->scalar , elm->scalar_it,elm->scalar - elm->scalar_it);
	FOR_ELEMENTS( elm ){
		if (fabs(elm->scalar - elm->scalar_it) > max_err) {
			 max_err = fabs(elm->scalar - elm->scalar_it);
			 //xprintf(Msg,"%f %f %f\n",elm->scalar , elm->scalar_it, elm->scalar - elm->scalar_it);
			 }
		}
	xprintf(Msg,"Maximum pressure difference in iteration: %f10.8\n",max_err);
	if (max_err > problem->transport->dens_eps) return 0; else return 1;
}
//=============================================================================
//      RESTART ITERATION CONCENTRATION
//=============================================================================
void restart_iteration_C(struct Problem *problem)
{
    int sbi,n_subst,sub,ph;

  Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    struct Transport *transport = problem->transport;
    double ***conc, ***prev_conc;

    n_subst = problem->transport->n_substances;
    sub = problem->transport->sub_problem;
    conc = transport->conc;
    prev_conc = transport->prev_conc;


    for(sbi =0; sbi < n_subst; sbi++)
    	for(ph = 0; ph < 4; ph++)
			if(conc[sbi][ph] != NULL)
				memcpy(conc[sbi][ph],prev_conc[sbi][ph],mesh->n_elements() * sizeof(double));

}
//=============================================================================
//      SAVE & RESTART ITERATION OF PRESSURE
//=============================================================================
void save_restart_iteration_H(struct Problem *problem)
{
        Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

        ElementIter elm;
	FOR_ELEMENTS( elm )
	{
		elm->scalar_it = elm->scalar;
	}
}
//=============================================================================
//      SAVE TIME STEP CONCENTRATION
//=============================================================================
void save_time_step_C(struct Problem *problem)
{
  Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    int sbi,n_subst,sub,ph;
    struct Transport *transport = problem->transport;
    double ***conc, ***prev_conc;

    n_subst = problem->transport->n_substances;
    sub = problem->transport->sub_problem;
    conc = transport->conc;
    prev_conc = transport->prev_conc;


    for(sbi =0; sbi < n_subst; sbi++)
    	for(ph = 0; ph < 4; ph++)
			if(conc[sbi][ph] != NULL)
				memcpy(prev_conc[sbi][ph],conc[sbi][ph],mesh->n_elements() * sizeof(double));

}
//
//=============================================================================
//      PEPA CHUDOBA
//=============================================================================
void decay( struct Transport *transport, int elm_pos, int type)
{
	int ph = 0;
	double delta_t = transport->time_step;
	double ***conc,***pconc;
	int i;
	conc = transport->conc;
	pconc = transport->pconc;


	switch(type){
	case 1:
		conc[1][ph][elm_pos]= pconc[1][ph][elm_pos]*(exp(-0.000000459*delta_t));     //Be 10
		conc[2][ph][elm_pos]= pconc[2][ph][elm_pos]*(exp(-0.000122*delta_t));        //C 14
		conc[3][ph][elm_pos]= pconc[3][ph][elm_pos]*(exp(-0.0000023*delta_t));        //Cl 36
		conc[4][ph][elm_pos]= pconc[4][ph][elm_pos]*(exp(-0.0000068*delta_t));        //Ca 41
		conc[5][ph][elm_pos]= pconc[5][ph][elm_pos]*(exp(-0.00000912*delta_t));        //Ni 59
		conc[6][ph][elm_pos]= pconc[6][ph][elm_pos]*(exp(-0.00702*delta_t));        //Ni 63
		conc[7][ph][elm_pos]= pconc[7][ph][elm_pos]*(exp(-0.00000195*delta_t));        //Se 79
		conc[8][ph][elm_pos]= pconc[8][ph][elm_pos]*(exp(-0.0241*delta_t));        //Sr 90
		conc[9][ph][elm_pos]= pconc[9][ph][elm_pos]*(exp(-0.000173*delta_t));        //Mo 93
		conc[10][ph][elm_pos]= pconc[10][ph][elm_pos]*(exp(-0.000000453*delta_t));    //Zr 93
		conc[11][ph][elm_pos]= pconc[11][ph][elm_pos]*(exp(-0.0000347*delta_t));    //Nb 94
		conc[12][ph][elm_pos]= pconc[12][ph][elm_pos]*(exp(-0.00000324*delta_t));    //Tc 99
		conc[13][ph][elm_pos]= pconc[13][ph][elm_pos]*(exp(-0.000000107*delta_t));    //Pd 107
		conc[14][ph][elm_pos]= pconc[14][ph][elm_pos]*(exp(-0.00158*delta_t));        //Ag 108m
		conc[15][ph][elm_pos]= pconc[15][ph][elm_pos]*(exp(-0.00000299*delta_t));    //Sn 126
		conc[16][ph][elm_pos]= pconc[16][ph][elm_pos]*(exp(-0.0000000431*delta_t));    //I 129
		conc[17][ph][elm_pos]= pconc[17][ph][elm_pos]*(exp(-0.000000301*delta_t));    //Cs 135
		conc[18][ph][elm_pos]= pconc[18][ph][elm_pos]*(exp(-0.0231*delta_t));        //Cs 137
		conc[19][ph][elm_pos]= pconc[19][ph][elm_pos]*(exp(-0.0077*delta_t));        //Sm 151
		conc[20][ph][elm_pos]= pconc[20][ph][elm_pos]*(exp(-0.000578*delta_t));        //Ho 166m
		conc[21][ph][elm_pos]= pconc[21][ph][elm_pos]*(exp(-0.000433*delta_t));        //Ra 226
		conc[22][ph][elm_pos]= pconc[22][ph][elm_pos]*(exp(-0.0000944*delta_t));    //Th 229
		conc[23][ph][elm_pos]= pconc[23][ph][elm_pos]*(exp(-0.0000092*delta_t));    //Th 230
		conc[24][ph][elm_pos]= pconc[24][ph][elm_pos]*(exp(-0.0000212*delta_t));    //Pa 231
		conc[25][ph][elm_pos]= pconc[25][ph][elm_pos]*(exp(-0.0000000000493*delta_t));    //Th 232
		conc[26][ph][elm_pos]= pconc[26][ph][elm_pos]*(exp(-0.00000435*delta_t));    //U 233
		conc[27][ph][elm_pos]= pconc[27][ph][elm_pos]*(exp(-0.00000282*delta_t));    //U 234
		conc[28][ph][elm_pos]= pconc[28][ph][elm_pos]*(exp(-0.000000000985*delta_t));    //U 235
		conc[29][ph][elm_pos]= pconc[29][ph][elm_pos]*(exp(-0.0000000296*delta_t));    //U 236
		conc[30][ph][elm_pos]= pconc[30][ph][elm_pos]*(exp(-0.000000324*delta_t));    //Np 237
		conc[31][ph][elm_pos]= pconc[31][ph][elm_pos]*(exp(-0.0079*delta_t));        //Pu 238
		conc[32][ph][elm_pos]= pconc[32][ph][elm_pos]*(exp(-0.000000000155*delta_t));    //U 238
		conc[33][ph][elm_pos]= pconc[33][ph][elm_pos]*(exp(-0.0000288*delta_t));    //Pu 239
		conc[34][ph][elm_pos]= pconc[34][ph][elm_pos]*(exp(-0.000106*delta_t));        //Pu 240
		conc[35][ph][elm_pos]= pconc[35][ph][elm_pos]*(exp(-0.0016*delta_t));        //Am 241
		conc[36][ph][elm_pos]= pconc[36][ph][elm_pos]*(exp(-0.00492*delta_t));        //Am 242m
		conc[37][ph][elm_pos]= pconc[37][ph][elm_pos]*(exp(-0.00000186*delta_t));    //Pu 242
		conc[38][ph][elm_pos]= pconc[38][ph][elm_pos]*(exp(-0.000094*delta_t));        //Am 243
		conc[39][ph][elm_pos]= pconc[39][ph][elm_pos]*(exp(-0.0383*delta_t));        //Cm 244
		conc[40][ph][elm_pos]= pconc[40][ph][elm_pos]*(exp(-0.0000817*delta_t));    //Cm 245
		conc[41][ph][elm_pos]= pconc[41][ph][elm_pos]*(exp(-0.000146*delta_t));        //Cm 246
		for(i=0;i<42;i++){
			pconc[i][ph][elm_pos] = conc[i][ph][elm_pos];
		}
		break;
	case 2:
		conc[1][ph][elm_pos]= pconc[1][ph][elm_pos]*(exp(-0.000433*delta_t));        //Ra 226
		conc[2][ph][elm_pos]= pconc[2][ph][elm_pos]*(exp(-0.0000944*delta_t));    //Th 229
		conc[3][ph][elm_pos]= pconc[3][ph][elm_pos]*(exp(-0.0000092*delta_t));    //Th 230
		conc[4][ph][elm_pos]= pconc[4][ph][elm_pos]*(exp(-0.0000212*delta_t));    //Pa 231
		conc[5][ph][elm_pos]= pconc[5][ph][elm_pos]*(exp(-0.0000000000493*delta_t));    //Th 232
		conc[6][ph][elm_pos]= pconc[6][ph][elm_pos]*(exp(-0.00000435*delta_t));    //U 233
		conc[7][ph][elm_pos]= pconc[7][ph][elm_pos]*(exp(-0.00000282*delta_t));    //U 234
		conc[8][ph][elm_pos]= pconc[8][ph][elm_pos]*(exp(-0.000000000985*delta_t));    //U 235
		conc[9][ph][elm_pos]= pconc[9][ph][elm_pos]*(exp(-0.0000000296*delta_t));    //U 236
		conc[10][ph][elm_pos]= pconc[10][ph][elm_pos]*(exp(-0.000000324*delta_t));    //Np 237
		conc[11][ph][elm_pos]= pconc[11][ph][elm_pos]*(exp(-0.0079*delta_t));        //Pu 238
		conc[12][ph][elm_pos]= pconc[12][ph][elm_pos]*(exp(-0.000000000155*delta_t));    //U 238
		conc[13][ph][elm_pos]= pconc[13][ph][elm_pos]*(exp(-0.0000288*delta_t));    //Pu 239
		conc[14][ph][elm_pos]= pconc[14][ph][elm_pos]*(exp(-0.000106*delta_t));        //Pu 240
		conc[15][ph][elm_pos]= pconc[15][ph][elm_pos]*(exp(-0.0016*delta_t));        //Am 241
		conc[16][ph][elm_pos]= pconc[16][ph][elm_pos]*(exp(-0.00492*delta_t));        //Am 242m
		conc[17][ph][elm_pos]= pconc[17][ph][elm_pos]*(exp(-0.00000186*delta_t));    //Pu 242
		conc[18][ph][elm_pos]= pconc[18][ph][elm_pos]*(exp(-0.000094*delta_t));        //Am 243
		conc[19][ph][elm_pos]= pconc[19][ph][elm_pos]*(exp(-0.0383*delta_t));        //Cm 244
		conc[20][ph][elm_pos]= pconc[20][ph][elm_pos]*(exp(-0.0000817*delta_t));    //Cm 245
		conc[21][ph][elm_pos]= pconc[21][ph][elm_pos]*(exp(-0.000146*delta_t));        //Cm 246
		for(i=0;i<22;i++){
			pconc[i][ph][elm_pos] = conc[i][ph][elm_pos];
		}
		break;
	case 3:
		//conc[0][ph][elm_pos]= pconc[0][ph][elm_pos]*(exp(-0*delta_t));//
		conc[1][ph][elm_pos]= pconc[1][ph][elm_pos]*(exp(-0.000000459*delta_t));     //Be 10
		conc[2][ph][elm_pos]= pconc[2][ph][elm_pos]*(exp(-0.000122*delta_t));        //C 14
		conc[3][ph][elm_pos]= pconc[3][ph][elm_pos]*(exp(-0.0000023*delta_t));        //Cl 36
		conc[4][ph][elm_pos]= pconc[4][ph][elm_pos]*(exp(-0.0000068*delta_t));        //Ca 41
		conc[5][ph][elm_pos]= pconc[5][ph][elm_pos]*(exp(-0.00000912*delta_t));        //Ni 59
		conc[6][ph][elm_pos]= pconc[6][ph][elm_pos]*(exp(-0.00702*delta_t));        //Ni 63
		conc[7][ph][elm_pos]= pconc[7][ph][elm_pos]*(exp(-0.000173*delta_t));        //Mo 93
		conc[8][ph][elm_pos]= pconc[8][ph][elm_pos]*(exp(-0.000000453*delta_t));    //Zr 93
		conc[9][ph][elm_pos]= pconc[9][ph][elm_pos]*(exp(-0.0000347*delta_t));    //Nb 94
		conc[10][ph][elm_pos]= pconc[10][ph][elm_pos]*(exp(-0.00158*delta_t));        //Ag 108m
		conc[11][ph][elm_pos]= pconc[11][ph][elm_pos]*(exp(-0.000578*delta_t));        //Ho 166m
		for(i=0;i<12;i++){
			pconc[i][ph][elm_pos] = conc[i][ph][elm_pos];
		}
		break;
	case 4:
		//conc[0][ph][elm_pos]= pconc[0][ph][elm_pos]*(exp(-0*delta_t));//
		conc[1][ph][elm_pos]= pconc[1][ph][elm_pos]*(exp(-0.00000195*delta_t));        //Se 79
		conc[2][ph][elm_pos]= pconc[2][ph][elm_pos]*(exp(-0.0241*delta_t));        //Sr 90
		conc[3][ph][elm_pos]= pconc[3][ph][elm_pos]*(exp(-0.000000453*delta_t));    //Zr 93
		conc[4][ph][elm_pos]= pconc[4][ph][elm_pos]*(exp(-0.00000324*delta_t));    //Tc 99
		conc[5][ph][elm_pos]= pconc[5][ph][elm_pos]*(exp(-0.000000107*delta_t));    //Pd 107
		conc[6][ph][elm_pos]= pconc[6][ph][elm_pos]*(exp(-0.00000299*delta_t));    //Sn 126
		conc[7][ph][elm_pos]= pconc[7][ph][elm_pos]*(exp(-0.0000000431*delta_t));    //I 129
		conc[8][ph][elm_pos]= pconc[8][ph][elm_pos]*(exp(-0.000000301*delta_t));    //Cs 135
		conc[9][ph][elm_pos]= pconc[9][ph][elm_pos]*(exp(-0.0231*delta_t));        //Cs 137
		conc[10][ph][elm_pos]= pconc[10][ph][elm_pos]*(exp(-0.0077*delta_t));        //Sm 151
	/*	if(pconc[2][ph][elm_pos] != 0){
		printf("\nA: %f\t A+: %f",pconc[0][ph][elm_pos],conc[0][ph][elm_pos]);
		printf("\nC: %f\t C+: %f",pconc[2][ph][elm_pos],conc[2][ph][elm_pos]);
		getchar();
		} */
		for(i=0;i<11;i++){
			pconc[i][ph][elm_pos] = conc[i][ph][elm_pos];
		}
		break;
	default:
		break;
		}
}
//=============================================================================
// RETURN MAX NUMBER OF ENTRIES IN THE ROW
//=============================================================================

void transport_partioning(struct Transport *transport,int np)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

  struct Edge *edge;
  int i,j,k,pos,proc,id,llim,ulim,n_pos,l_pos,llimb,ulimb;



  transport->l_row = (int*)xmalloc(np*sizeof(int));  // local rows in TM
  transport->lb_col = (int*)xmalloc(np*sizeof(int)); // local rows in BM
  for(i=0;i< np;i++){
	  transport->l_row[i] = mesh->n_elements() / np;
	  transport->lb_col[i] = mesh->n_boundaries() / np;
  }

  //transport->l_row[i] = (int)floor( (double)mesh->n_elements() / (double)np);

  transport->l_row[np-1] = mesh->n_elements() - (np - 1 ) * transport->l_row[0];
  transport->lb_col[np-1] = mesh->n_boundaries() - (np - 1 ) * transport->lb_col[0];


  transport->d_row =(int**)xmalloc(np*sizeof(int*));
  transport->od_row =(int**)xmalloc(np*sizeof(int*));
  transport->db_row =(int**)xmalloc(np*sizeof(int*));
  transport->odb_row =(int**)xmalloc(np*sizeof(int*));
  for(i=0; i< np; i++){
	  transport->d_row[i] =(int*)xmalloc(transport->l_row[i]*sizeof(int));
	  transport->od_row[i] =(int*)xmalloc(transport->l_row[i]*sizeof(int));
	  transport->db_row[i] =(int*)xmalloc(transport->l_row[i]*sizeof(int));
	  transport->odb_row[i] =(int*)xmalloc(transport->l_row[i]*sizeof(int));
  }

  /*   l_row[0]  !!!!  */


  FOR_ELEMENTS( elm )
  {
	  pos = ELEMENT_FULL_ITER(elm) - mesh->element.begin();
	  //pos = id2pos(mesh,elm.id(),mesh->epos_id,ELM);


	  proc = pos/transport->l_row[0];
	  if(proc == np)
		  proc = np -1;

	  //if((pos - proc *transport-> l_row[0])) l_pos = 0;

	  l_pos = pos - proc * transport->l_row[0];

//	  printf("\nl_pos:%d\tl_row[%d]%d",l_pos,proc,transport->l_row[proc]);

	  transport->d_row[proc][l_pos] = 1;
	  transport->od_row[proc][l_pos] = 0;
	  transport->db_row[proc][l_pos] = 0;
	  transport->odb_row[proc][l_pos] = 0;

	  llim = proc * transport->l_row[0];
	  if(proc < (np -1))
		ulim = proc * transport->l_row[0] - 1;
	  else
		ulim = mesh->n_elements() - 1;

	  llimb = proc * transport->lb_col[0];
	  if(proc < (np -1))
		ulimb = proc * transport->lb_col[0] - 1;
	  else
		ulimb = mesh->n_boundaries() - 1;

	//  printf("\nllim:%d\tulim:%d\t\tD:%d\tOD:%d",llim,ulim,transport->d_row[proc][l_pos],transport->od_row[proc][l_pos]);

    FOR_ELEMENT_SIDES(elm,si){                         //same dim
      if (elm->side[si]->cond != NULL){  //BC

    	//n_pos = id2pos(mesh,elm->side[si]->id,mesh->spos_id,BC) - mesh->n_elements();
    	 n_pos = BOUNDARY_FULL_ITER(elm->side[si]->cond) - mesh->boundary.begin();
		  if((n_pos < llimb) || (n_pos > ulimb))
			  transport->odb_row[proc][l_pos]++;
		  else
			  transport->db_row[proc][l_pos]++;
      }
      else{
		  edge = elm->side[ si ]->edge;
		  FOR_EDGE_SIDES(edge,j){
			 // id = ELEMENT_FULL_ITER(edge->side[j]->element).id();
			 // n_pos = id2pos(mesh,id,mesh->epos_id,ELM);
			  n_pos = ELEMENT_FULL_ITER(edge->side[j]->element) - mesh->element.begin();

			  if(n_pos != pos)
			  {
				  if((n_pos < llim) || (n_pos > ulim))
					  transport->od_row[proc][l_pos]++;
				  else
					  transport->d_row[proc][l_pos]++;
			  }
		  }
		}
		  if(elm->side[si]->neigh_bv != NULL){   // comp model
			  //id = ELEMENT_FULL_ITER(elm->side[si]->neigh_bv->element[0]).id();
			  //n_pos = id2pos(mesh,id,mesh->epos_id,ELM);
			  n_pos = ELEMENT_FULL_ITER(elm->side[si]->neigh_bv->element[0]) - mesh->element.begin();
			  if(n_pos != pos)
			  {
				  if((n_pos < llim) || (n_pos > ulim))
					  transport->od_row[proc][l_pos]++;
				  else
					  transport->d_row[proc][l_pos]++;
			  }
		  }

      }                                            // end same dim


   //   printf("SD id:%d,size:%d\n",elm->id,size);
     //

   FOR_ELM_NEIGHS_VB(elm,i){                    // comp model
       id = ELEMENT_FULL_ITER(elm->neigh_vb[i]->element[1]).id();
       n_pos = id2pos(mesh,id,mesh->epos_id,ELM);
       if(n_pos != pos)
       {
			  if((n_pos < llim) || (n_pos > ulim))
				  transport->od_row[proc][l_pos]++;
			  else
				  transport->d_row[proc][l_pos]++;
       }


  //     printf("VB id:%d,size:%d\n",elm->id,size);
      }                                        // end comp model


   FOR_ELM_NEIGHS_VB(elm,i){                    // non-comp model

	   id = ELEMENT_FULL_ITER(elm->neigh_vb[i]->element[0]).id();
	   if(id == elm.id())
		   id = ELEMENT_FULL_ITER(elm->neigh_vb[i]->element[1]).id();

       n_pos = id2pos(mesh,id,mesh->epos_id,ELM);
       if(n_pos != pos)
       {
			  if((n_pos < llim) || (n_pos > ulim))
				  transport->od_row[proc][l_pos]++;
			  else
				  transport->d_row[proc][l_pos]++;
       }

   }
  //  l_pos++;
  }

/*
  getchar();
  for(i=0; i< transport->l_row[0]; i++)
	  printf("\ni:%d\tD:%d\tOD:%d",i,transport->d_row[proc][i],transport->od_row[proc][i]);
  getchar();
*/

}
