/*!
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file  Distributed sparse graphs, partitioning.
 *
 *
 */
#ifndef SPARSE_GRAPH_HH_
#define SPARSE_GRAPH_HH_

#include <par_distribution.hh>
#include <vector>
#include <stack>
#include <ostream>
#include <petscmat.h>

using namespace std;


// Auxiliary Edge type.
typedef struct {
   int from, to, weight;          ///< Edge weights for communication (optional).
} SparseGraphEdge;



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
protected:
    typedef SparseGraphEdge Edge;

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
    vector< stack<Edge> > adj_of_proc;      // storage for graph edges for individual processors
    //int *degree;                  // global array to count vtx degrees, temporary
    //int *num_edges_of_proc;       // count local inserted edges of all processors

    virtual void allocate_sparse_graph(int lsize_vtxs, int lsize_adj)=0;


public:
    SparseGraph(const Distribution &distr);            // construct by distribution
    SparseGraph(int loc_size);   // construct by local size and comunicator

    void set_edge(const int a, const int b, int weight=1);
    void set_vtx_position(const int vtx, const float xyz[3], int weight=1);
    void finalize();
//    const Distribution get_vtx_distr() { return vtx_distr; }

    virtual void partition(int *loc_part) = 0;
    //void make_partition_PETSC(IS *part, bool check_subgraph_continuity=false);
    //void make_partition_Metis(int *part,  bool check_subgraph_continuity=false);
    bool check_subgraph_connectivity(int *part);
    void DFS(int vtx);

    const int * get_partition();
    Distribution get_distr() {return vtx_distr;}
    void view();
    bool is_symmetric();
    virtual ~SparseGraph();

    friend ostream & operator <<(ostream & out, const SparseGraph &sg);
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
