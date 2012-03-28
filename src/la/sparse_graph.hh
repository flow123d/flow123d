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
 * @brief Distributed sparse graphs, partitioning.
 *
 *
 */
#ifndef SPARSE_GRAPH_HH_
#define SPARSE_GRAPH_HH_


#include <vector>
#include <stack>
#include <ostream>
#include <petscmat.h>

#include "la/distribution.hh"
using namespace std;





/**
 * @brief  Virtual class for construction and partitioning of a distributed sparse graph.
 *
 * An empty graph is created in constructor (collective). You has to provide either Distribution of the vertices or
 * local number of vertices and communicator. Than you can add edges with set_edge method. Added edges are
 * stored in the buffer on calling  processor. Than you have to call finalize method (collective) to make sparse graph.
 * Graph partitioning is either directly through METIS or through PETSC (ParMETIS, and others). There is derived class for
 * particular partitioning method.
 *
 * Actual implementation is not fully scalable. row sizes (vtx degree)
 * are accumulated in global arrays, then communicated to all processors and finally
 * only local part is used. There should be only local part of the final array with aditional buffer for ex-processor
 * values.
 *
 * Other issue is that we hav eproblem using Parmetis so acctualy only METIS can be used on
 * a parallel graph located on the first processor.
 *
 * TODO: use Boost parallel graphs and write interface form them to METIS, Parmetis or PETSC
 *
 */
class SparseGraph {
public:

    /**
     * Construct an empty graph form given Distribution of vertices.
     */
    SparseGraph(const Distribution &distr);

    /**
     * Construct an empty graph for given number of local vertices of the graph.
     * Make its own distribution object.
     */
    SparseGraph(int loc_size);

    /**
     * Store an edge. We just store all inserted edges and count support arrays to
     * communicate edges to the right processors
     *
     * @param[in] a - starting vertex of the edge
     * @param[in] b - terminal vertex of the edge
     * @param[in] weight - optional weight of the edge (integer)
     *
     */
    void set_edge(const int a, const int b, int weight=1);

    /**
     * Set position and possibly weight of a local vertex. Assume that vtx is an index of local vertex.
     * Positions are used for initial distribution when using ParMETIS.
     *
     * @param[in] vtx - global vertex index (from zero)
     * @param[in] xyz - coordinates of vetrex position
     * @parem[in] weight - optional weight of the vertex
     *
     */
    void set_vtx_position(const int vtx, const float xyz[3], int weight=1);

    /**
     * @brief  Make sparse graph structures: rows, adj
     *
     * 1) send edges to the owner of .from vertex
     * 2) sort local edges
     * 3) fill rows, adj;   remove duplicities
     * @param[in] vtx - global vertex index (from zero)
     * @param[in] xyz - coordinates of vetrex position
     *
     * Assume that vtx is an index of local vertex.
     */
    void finalize();

    /**
     * Fills an array of integers by the partitioning of the vertices.
     * The array size is number of local verices according to the distribution of the graph.
     * see get_distr() method. The array has to be PREALLOCATED by user to the correct size.
     *
     */
    virtual void partition(int *loc_part) = 0;

    /**
     * Check if the subgraphs of the given partitioning are connected.
     * Works only for fully localized graph (on one processor).
     */
    bool check_subgraph_connectivity(int *part);
    void DFS(int vtx);

    /**
     * Returns reference to the distribution of the graph.
     */
    Distribution get_distr() {return vtx_distr;}

    /**
     * Simple graph output. Has to be already finalized.
     */
    void view();

    /**
     * Check symmetry of the local part of the graph. Has to be already finalized.
     */
    bool is_symmetric();
    virtual ~SparseGraph();

    //friend ostream & operator <<(ostream & out, const SparseGraph &sg);
protected:
    // Auxiliary Edge type.
    typedef struct {
       int from, to, weight;          ///< Edge weights for communication (optional).
    } Edge;


    Distribution vtx_distr;       ///< distribution of vertexes
    int *rows;                    ///< starts of vtx rows in adj array,
                                  ///< vtx i have edge to vertexes adj[rows[i]] .. adj[rows[i+1]-1]
    float * vtx_XYZ;              ///< optional vertex coordinates (global array)
    int *vtx_weights;             ///< Vertex weights for computations (optional).
    int *adj;                     ///< sparse adjency
    int *adj_weights;

    int *part_to_check;                    ///< created partitioning used through check of connectivity
    int proc_to_check;                     ///< subgraph to check
    std::vector<int> checked_vtx;          ///< coloring of DFS algorithm

    // temporary arrays
    vector< stack<Edge> > adj_of_proc;      ///< storage for graph edges for individual processors

    /**
     * Use virtual method for allocation since for PETSC interface we have to use PetscMalloc.
     */
    virtual void allocate_sparse_graph(int lsize_vtxs, int lsize_adj)=0;

    friend bool operator <(const Edge &a,const Edge  &b);

};

ostream & operator <<(ostream & out, const SparseGraph &sg);



/**
 * Sparse Graph that use PETSC for partitioning.
 */
class SparseGraphPETSC : public SparseGraph
{


public:
    /**
     *  Construct empty graph only from number of vertices, use default initial distribution.
     */
    SparseGraphPETSC(int n_vtxs)
        : SparseGraph(Distribution(Distribution::Block, n_vtxs)), petsc_adj_mat(0), petsc_part(0), part_IS(0) {}

    /**
     *  Construct empty graph with given distribution of vertices.
     */
    SparseGraphPETSC(const Distribution &distr)
        : SparseGraph(distr), petsc_adj_mat(0), petsc_part(0), part_IS(0) {}
    /**
     * Implementation of partitioning using PETSC interface to various partitioning software.
     */
    virtual void partition(int *loc_part);

private:
    Mat petsc_adj_mat;
    MatPartitioning petsc_part;
    IS part_IS;

    virtual void allocate_sparse_graph(int lsize_vtxs, int lsize_adj);
    virtual ~SparseGraphPETSC();

};


/**
 * Sparse Graph that use METIS for partitioning.
 */
class SparseGraphMETIS : public SparseGraph
{
public:
    /**
     *  Construct empty graph only from number of vertices,
     *  use localized distribution in order to use sequential METIS library.
     */
    SparseGraphMETIS(int n_vtxs)
        : SparseGraph(Distribution(Distribution::Localized, n_vtxs)) {}

    /**
     * Implementation of partitioning using METIS.
     */
    virtual void partition(int *loc_part);

private:
    virtual void allocate_sparse_graph(int lsize_vtxs, int lsize_adj);
    virtual ~SparseGraphMETIS();
};


#endif /* SPARSE_GRAPH_HH_ */
