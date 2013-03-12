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

#include "system/system.hh"
#include "system/par_distribution.hh"
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
      adj_of_proc( vtx_distr.np() ),
      adj(NULL),
      rows(NULL),
      adj_weights(NULL)
{
    F_ENTRY;

    // positions only of local vertexes
    vtx_XYZ= new float[vtx_distr.lsize()+1];
    vtx_weights= new int[vtx_distr.lsize()+1];
    for(int i=0; i<vtx_distr.lsize(); i++) vtx_weights[i]=1;    
}



SparseGraph::SparseGraph(int loc_size)
    : vtx_distr(loc_size),
      adj_of_proc( vtx_distr.np() ),
      adj(NULL),
      rows(NULL),
      adj_weights(NULL)
{
    F_ENTRY;

    // positions only of local vertexes
    vtx_XYZ= new float[vtx_distr.lsize()+1];
    vtx_weights= new int[vtx_distr.lsize()+1];
    for(int i=0; i<vtx_distr.lsize(); i++) vtx_weights[i]=1;    
}

void SparseGraph::set_edge(const int a, const int b,int weight)
{
    F_ENTRY;

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
   F_ENTRY;
   ASSERT( adj==NULL, "Graph is already finalized\n");

   int proc;
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
   //DBGMSG("construct adj\n");
   Edge *edges= (Edge *) recvbuf;
   int size=total_size/edge_size;
   std::sort(edges, edges + size);

   allocate_sparse_graph(vtx_distr.lsize() + 1, size+1);
   adj = (int *) xmalloc( (size+1) * sizeof(int) );
   rows = (int *) xmalloc( (vtx_distr.lsize() + 1) * sizeof(int) );
   adj_weights = (int *) xmalloc( (size+1) * sizeof(int) );
   
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
    F_ENTRY;

    ASSERT( vtx_distr.lsize(0)==vtx_distr.size() , "Check of graph continuity not yet implemented for paralel case.\n");
    if (vtx_distr.myp()!=0) return(true);

    part_to_check=part;
    checked_vtx.resize(vtx_distr.size(), 0);
    std::vector<bool> checked_proc(vtx_distr.np(), false);

    int n_proc=0;
    for(int vtx=0; n_proc<vtx_distr.np() && vtx<vtx_distr.size(); vtx++) {
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
    ASSERT( vtx>=0 && vtx<vtx_distr.size(),"Invalid entry vertex %d in DFS.\n",vtx);
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
    ASSERT(NONULL(adj),"Can not view non finalized graph.\n");
    int row,col;
    xprintf(Msg,"SparseGraph\n");
    for(row=0; row < vtx_distr.lsize(); row++) {
        xprintf(Msg,"edges from this vertex: %d\n",rows[row+1]);
        for(col=rows[row]; col<rows[row+1]; col++) {
            xprintf(Msg,"edge (v1, v2): %d %d\n",row+vtx_distr.begin(), adj[col]);
        }
    }
}


bool SparseGraph::is_symmetric()
{
    ASSERT( NONULL(rows) && NONULL(adj), "Graph is not yet finalized.");

    int loc_row, row, row_pos;
    int col_pos,col,loc_col;

    for(loc_row=0;loc_row<vtx_distr.lsize();loc_row++) {
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
    ASSERT(NONULL(adj) && NONULL(rows),"Can not make partition of non finalized graph.\n");

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
    F_ENTRY;

    ASSERT( vtx_distr.lsize(0)==vtx_distr.size(),
            "METIS could be used only with localized distribution.\n");

    if (vtx_distr.np()==1) {
        for(int i=0;i<vtx_distr.size();i++) part[i]=0;

    } else {


        int n_vtx=vtx_distr.size();
        int n_proc=vtx_distr.np();
        int wght_flag=3;    // which wghts are given (0-none,1-edges,2-vtxs,3-both)
        int num_flag=0;     // indexing style (0-C, 1-Fortran)
        int options[8];
        int edgecut;


        if (vtx_distr.myp()==0) {
            options[0]=0;
            options[4]=255;  //dbg_lvl

            DBGMSG("METIS call\n");
            METIS_PartGraphKway(&n_vtx,rows,adj, //vtx distr, local vtx begins, edges of local vtxs
                vtx_weights,adj_weights,&wght_flag,&num_flag, // vertex, edge weights, ...
                &n_proc,options,&edgecut,part);

            DBGMSG("Graph edge cut: %d\n",edgecut);
        }
        MPI_Bcast( part, n_vtx, MPI_INT, 0, vtx_distr.get_comm() );
    }
}

SparseGraphMETIS::~SparseGraphMETIS()
{
    if (adj) delete[](adj);
    if (rows) delete[](rows);
    if (adj_weights) delete[](adj_weights);

    // automaticaly calls base class destructor
}
