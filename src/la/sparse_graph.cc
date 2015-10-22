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
 *
 * @file
 * @ingroup la
 * @brief  Construction and partitioning of a sparse graph
 *
 */

#include "system/global_defs.h"
#include "system/system.hh"
#include "la/distribution.hh"
#include "sparse_graph.hh"

#include "boost/lambda/lambda.hpp"
#include <algorithm>

extern "C" {
#include <metis.h>
}


/*****************************************************************************************
 SPARSE GRAPH
 *****************************************************************************************/


SparseGraph::SparseGraph(const Distribution &distr)
    : vtx_distr(distr),
      rows(NULL),
      adj(NULL),
      adj_weights(NULL),
      part_to_check(NULL),
      adj_of_proc( vtx_distr.np() )
{
    // positions only of local vertexes
    vtx_XYZ= new float[vtx_distr.lsize()+1];
    vtx_weights= new int[vtx_distr.lsize()+1];
    for(unsigned int i=0; i<vtx_distr.lsize(); i++) vtx_weights[i]=1;
}



SparseGraph::SparseGraph(int loc_size, MPI_Comm comm)
    : vtx_distr(loc_size, comm),
      rows(NULL),
      adj(NULL),
      adj_weights(NULL),
      adj_of_proc( vtx_distr.np() )
{
    // positions only of local vertexes
    vtx_XYZ= new float[vtx_distr.lsize()+1];
    vtx_weights= new int[vtx_distr.lsize()+1];
    for(unsigned int i=0; i<vtx_distr.lsize(); i++) vtx_weights[i]=1;
}

void SparseGraph::set_edge(const int a, const int b,int weight)
{
    Edge e={a,b, weight};

    adj_of_proc[vtx_distr.get_proc(a)].push(e);
}




void SparseGraph::set_vtx_position(const int vtx, const float xyz[3],int weight)
{
    ASSERT(vtx_distr.is_local(vtx),"Can not set vertex position for nonlocal vertex %d.\n",vtx);
    int loc_index=vtx-vtx_distr.begin();
    memcpy(vtx_XYZ+3*loc_index, xyz, 3*sizeof(float));
    vtx_weights[loc_index]=weight;
}


/**
 *   Edge comparison. For lexical sorting of local edges.
 */

bool operator <(const SparseGraph::Edge &a,const SparseGraph::Edge  &b)
{
    return (a.from < b.from ||
            ( a.from == b.from && a.to < b.to));
}


// TODO: split into smaller functions
// check trivial case : Localized distribution
void SparseGraph::finalize()
{
   ASSERT( adj==NULL, "Graph is already finalized\n");

   unsigned int proc;
   int total_size;
   vector< stack<Edge> >::iterator s;
   unsigned int edge_size=3;   // 3 = number of integers in Edge to send

   /////////////////////////////////////
   // communicate inserted edges

   // unroll data in adj_of_proc into send buffer
   int * sdispls = (int *) xmalloc( vtx_distr.np() * sizeof(int) );
   int * scounts = (int *) xmalloc( vtx_distr.np() * sizeof(int) );

   total_size=0;
   for(proc=0, s=adj_of_proc.begin(); s!=adj_of_proc.end();++s, ++proc)
   {
       sdispls[proc] = total_size;
       scounts[proc] = edge_size * (s)->size(); 
       total_size += edge_size * ((s)->size)();
   }
   int * sendbuf = (int *) xmalloc( (total_size+1) * sizeof(int ) );

   Edge edge;
   int buf_pos=0;
   for(proc=0, s = adj_of_proc.begin(); s!=adj_of_proc.end(); ++s, ++proc)
   {
       ASSERT(sdispls[proc] == buf_pos,
               "Mismatch between displacement %d and buffer position %d. \n", sdispls[proc], buf_pos );
       while ( ! (s)->empty() ) {
           edge=(s)->top();
           sendbuf[buf_pos++]=edge.from;
           sendbuf[buf_pos++]=edge.to;
           sendbuf[buf_pos++]=edge.weight;
           (s)->pop();
       }
   }

   // communicate send sizes
   int * rcounts = (int *) xmalloc( vtx_distr.np() * sizeof(int) );
   int * rdispls = (int *) xmalloc( vtx_distr.np() * sizeof(int) );
   MPI_Alltoall( scounts, 1, MPI_INT, rcounts, 1, MPI_INT, vtx_distr.get_comm());

   // prepare receiving buffers
   total_size=0;
   for(proc=0; proc<vtx_distr.np(); proc++)
   {
       rdispls[proc] = total_size;
       total_size += rcounts[proc];
   }

   int * recvbuf = (int *) xmalloc( (total_size+1) * sizeof(int ) );

   MPI_Alltoallv (
           sendbuf, scounts, sdispls, MPI_INT,
           recvbuf,  rcounts, rdispls, MPI_INT,
           vtx_distr.get_comm() );

   xfree(sendbuf);
   xfree(scounts);
   xfree(sdispls);
   xfree(rcounts);
   xfree(rdispls);

   /////////////////////////////////////
   // construct local adj and rows arrays
   Edge *edges= (Edge *) recvbuf;
   int size=total_size/edge_size;
   std::sort(edges, edges + size);

   allocate_sparse_graph(vtx_distr.lsize() + 1, size+1);
   rows[0]=0;

   if (size != 0) {

   int i;
   int row=0;
   int i_adj=0;
   int loc_from;
   Edge unknown_edge={-1,-1,0};
   Edge *last_edge=&unknown_edge;

   for(i=0;i<size;i++) {
       if (! (*last_edge < edges[i]) ) continue; // skip equivalent edges
       last_edge=edges+i;

       ASSERT(vtx_distr.is_local(edges[i].from),
               "Received non-local edge: %d %d at position %d\n",edges[i].from, edges[i].to,i);

       loc_from=edges[i].from-vtx_distr.begin();
       ASSERT( row <= loc_from, "Decrease in sorted edges at %d\n",i);

       while ( row < loc_from ) rows[++row]=i;
       adj[i_adj]=edges[i].to;
       adj_weights[i_adj]=edges[i].weight;
       i_adj++;
   }
   rows[++row]=i; // i==size (of adj array)

   }

   xfree(recvbuf);
}


/**
 * Check if the subgraphs of the partitions are connected.
 */
bool SparseGraph::check_subgraph_connectivity(int *part)
{
    ASSERT( vtx_distr.lsize(0)==vtx_distr.size() , "Check of graph continuity not yet implemented for paralel case.\n");
    if (vtx_distr.myp()!=0) return(true);

    part_to_check=part;
    checked_vtx.resize(vtx_distr.size(), 0);
    std::vector<bool> checked_proc(vtx_distr.np(), false);

    unsigned int n_proc=0;
    for(unsigned int vtx=0; n_proc<vtx_distr.np() && vtx<vtx_distr.size(); vtx++) {
        if (checked_vtx[vtx] != 2) {
            proc_to_check=part_to_check[vtx];
            // check if the processor is still unvisited
            if ( checked_proc[proc_to_check] )
                xprintf(Warn, "Disconnected subgraph %d detected at vertex %d.\n",part_to_check[vtx],vtx);
            // DFS unvisited vertex
            checked_vtx[vtx]=1;
            DFS(vtx);
            checked_proc[proc_to_check]=true;
            n_proc++;
        }
    }
    checked_vtx.clear();

    DBGMSG("Connectivity of subgraphs is OK.\n");
    return (true);
}

/**
 * Recursive part of Deep First Search algorithm of the connectivity check.
 */
void SparseGraph::DFS(int vtx)
{
    ASSERT( vtx>=0 && vtx< (int) vtx_distr.size(),"Invalid entry vertex %d in DFS.\n",vtx);
    int neighbour;
    for(int i_neigh=rows[vtx]; i_neigh< rows[vtx+1];i_neigh++) {
        neighbour = adj[i_neigh];
        if (part_to_check[neighbour] == proc_to_check && checked_vtx[neighbour] == 0) {
            checked_vtx[neighbour]=1;
            DFS(neighbour);
        }
    }
    checked_vtx[vtx]=2;
}



void SparseGraph::view()
{
    ASSERT( adj,"Can not view non finalized graph.\n");
    int row,col;
    xprintf(Msg,"SparseGraph\n");
    for(row=0; row < (int) vtx_distr.lsize(); row++) {
        xprintf(Msg,"edges from this vertex: %d\n",rows[row+1]);
        for(col=rows[row]; col<rows[row+1]; col++) {
            xprintf(Msg,"edge (v1, v2): %d %d\n",row+vtx_distr.begin(), adj[col]);
        }
    }
}


bool SparseGraph::is_symmetric()
{
    ASSERT( rows && adj, "Graph is not yet finalized.");

    int loc_row, row, row_pos;
    int col_pos,col,loc_col;

    for(loc_row=0;loc_row< (int) vtx_distr.lsize();loc_row++) {
        row=loc_row+vtx_distr.begin();
        for(row_pos=rows[loc_row];row_pos<rows[loc_row+1];row_pos++) {
            col=adj[row_pos];
            if (vtx_distr.is_local(col)) {
                // find the local column
                loc_col=col-vtx_distr.begin();
                for(col_pos=rows[loc_col]; col_pos<rows[loc_col+1]; col_pos++)
                  if (adj[col_pos]==row) break;
                if (adj[col_pos]!=row) return false;
            }
        }
    }

    return true;
}


/**
 * Free all non-null pointers.
 */
SparseGraph::~SparseGraph()
{
   delete[] vtx_XYZ;
   delete[] vtx_weights;
}

/*!
 * @brief     Output a sparse graph.
 * @return    Changed ostream.
 */
ostream & operator <<(ostream & out, const SparseGraph &fcall)
{
    // TODO
    return out; //dodelat...
}


/*****************************************************************************************
 SPARSE GRAPH PETSC
 *****************************************************************************************/

void SparseGraphPETSC::allocate_sparse_graph(int lsize_vtxs, int lsize_adj)
{
    PetscMalloc(lsize_vtxs*sizeof(int), &rows);
    PetscMalloc(lsize_adj*sizeof(int),&adj);
    PetscMalloc(lsize_adj*sizeof(int),&adj_weights);
}


void SparseGraphPETSC::partition(int *loc_part)
{
    ASSERT( adj && rows,"Can not make partition of non finalized graph.\n");

    MatCreateMPIAdj(vtx_distr.get_comm(), vtx_distr.lsize(),vtx_distr.size(),
            rows, adj,adj_weights, &petsc_adj_mat);
    MatPartitioningCreate(vtx_distr.get_comm(),& petsc_part);
    MatPartitioningSetAdjacency(petsc_part, petsc_adj_mat);
    MatPartitioningSetFromOptions(petsc_part);
    MatPartitioningApply(petsc_part,&part_IS);

    // copy the partition array
    const PetscInt *is_array;
    ISGetIndices(part_IS, &is_array);
    memcpy(loc_part, is_array, vtx_distr.lsize()*sizeof(int));
    ISRestoreIndices(part_IS, &is_array);
}


SparseGraphPETSC::~SparseGraphPETSC()
{
    ISDestroy(&part_IS);
    MatPartitioningDestroy(&petsc_part);
    MatDestroy(&petsc_adj_mat);

    if (adj) PetscFree(adj);
    if (rows) PetscFree(rows);
    if (adj_weights) PetscFree(adj_weights);

    // automaticaly calls base class destructor
}


/*****************************************************************************************
 SPARSE GRAPH METIS
 *****************************************************************************************/


void SparseGraphMETIS::allocate_sparse_graph(int lsize_vtxs, int lsize_adj)
{
    rows = new int[lsize_vtxs];
    adj = new int[lsize_adj];
    adj_weights = new int [lsize_adj];
}

void SparseGraphMETIS::partition(int *part)
{
    ASSERT( vtx_distr.lsize(0)==vtx_distr.size(),
            "METIS could be used only with localized distribution.\n");

    if (vtx_distr.np()==1) {
        for(unsigned int i=0;i<vtx_distr.size();i++) part[i]=0;
        return;
    } else {
        if (vtx_distr.myp()==0) {
                  int n_vtx=vtx_distr.size();
                  int n_proc=vtx_distr.np();
                  int num_flag=0;     // indexing style (0-C, 1-Fortran)
                  int edgecut;

/***********************************************************************************
 *  SETTING OPTIONS
 */        
#if (METIS_VER_MAJOR >= 5)
                  if ( sizeof(idx_t) != sizeof(int) ) {
                    printf("ERROR in GRAPH_DIVIDE_C: Wrong type of integers for METIS.\n");
                    abort();
                  }
                  int ncon = 1;
                  real_t ubvec[1];
                  ubvec[0] = 1.001;
                  int options[METIS_NOPTIONS];

                  for (unsigned int i = 0;i < METIS_NOPTIONS;i++) options[i] = -1;
                  
                  options[METIS_OPTION_OBJTYPE]   = METIS_OBJTYPE_CUT;
                  options[METIS_OPTION_CTYPE]     = METIS_CTYPE_RM;
                  options[METIS_OPTION_IPTYPE]    = METIS_IPTYPE_GROW;
                  options[METIS_OPTION_RTYPE]     = METIS_RTYPE_GREEDY;
                  options[METIS_OPTION_NCUTS]     = 1;
                  options[METIS_OPTION_NSEPS]     = 1;
                  options[METIS_OPTION_NUMBERING] = num_flag;
                  options[METIS_OPTION_NITER]     = 10;
                  options[METIS_OPTION_SEED]      = 12345;
                  options[METIS_OPTION_MINCONN]   = 1;
                  options[METIS_OPTION_CONTIG]    = 0;
                  options[METIS_OPTION_COMPRESS]  = 0;
                  options[METIS_OPTION_CCORDER]   = 0;
                  options[METIS_OPTION_UFACTOR]   = 0;
                  /*options[METIS_OPTION_DBGLVL]    = METIS_DBG_INFO;*/
                  options[METIS_OPTION_DBGLVL]    = 0;
#else
                  /* weights */
                  int wgtflag=3;
                  int options[5];
                  for (unsigned int  i = 0; i < 5; i++ )  options[i] = 0;
              
#endif
                  

                  /* Initialize parts */
                  for (unsigned int  i = 0; i < vtx_distr.lsize(); i++ )   part[i] = num_flag;

/***********************************************************************************
 *  CALL METIS using optimal algorithm depending on the number of partitions
 */        

                  if (n_proc  <= 8) {

#if (METIS_VER_MAJOR >= 5)
                      options[METIS_OPTION_PTYPE]     = METIS_PTYPE_RB;
                      options[METIS_OPTION_UFACTOR]   = 1;
                      METIS_PartGraphRecursive(&n_vtx, &ncon, rows, adj,
                                                vtx_weights, NULL, adj_weights, &n_proc, NULL,
                                                ubvec, options, &edgecut,part);
#else
                      // has to be checked
                      METIS_PartGraphRecursive(&n_vtx, rows, adj, 
                                                vtx_weights, adj_weights, &wgtflag, &num_flag, 
                                                &n_proc, options, &edgecut, part);
#endif
                  } else {
    
#if (METIS_VER_MAJOR >= 5)
                      options[METIS_OPTION_PTYPE]     = METIS_PTYPE_KWAY;
                      options[METIS_OPTION_UFACTOR]   = 30;
                      METIS_PartGraphKway(&n_vtx, &ncon, rows, adj,
                                          vtx_weights, NULL, adj_weights, &n_proc, NULL,
                                          ubvec, options, &edgecut, part);
#else
                      METIS_PartGraphKway(&n_vtx,rows,adj, //vtx distr, local vtx begins, edges of local vtxs
                                  vtx_weights,adj_weights,&wgtflag,&num_flag, // vertex, edge weights, ...
                                  &n_proc,options,&edgecut,part);
#endif
                  }     
        }
    }
}

SparseGraphMETIS::~SparseGraphMETIS()
{
    if (adj) delete[](adj);
    if (rows) delete[](rows);
    if (adj_weights) delete[](adj_weights);

    // automaticaly calls base class destructor
}
