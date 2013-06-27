/*
 * partition.hh
 *
 *  Created on: May 3, 2013
 *      Author: jb
 */

#ifndef PARTITIONING_HH_
#define PARTITIONING_HH_

#include "input/input_type.hh"
#include "input/accessors.hh"

class Mesh;
class SparseGraph;
class Distribution;

/**
 * @brief Class for the mesh partitioning.
 * This should provide:
 * - construction of cenectivity graph from the mesh
 * - call a partitioner
 * - provide new element numbering and new element distribution
 *   (global numbers for the local part - at least for ghost part)
 *   (+ old indices of elements for the local part of the mesh - in order to do not change current implementation)
 * TODO:
 * - deal with partitioining of boundary mesh
 */
class Partitioning {
public:

    /// Input specification objects.
    static Input::Type::Selection graph_type_sel;
    static Input::Type::Selection tool_sel;
    static Input::Type::Record input_type;

    /**
     *  Constructor. A pointer to the mesh and accessor to an input record have to be provided.
     */
    Partitioning(Mesh *mesh, Input::Record in);

    /**
     * # set partitioning into the elements of the mesh
     */
//    void set_partition_to_mesh();
    /**
     * Get initial distribution.
     */
    const Distribution *get_init_distr() const;
    /**
     * Get local part of mesh partition.
     */
    const int *get_loc_part() const;

    /// Destructor.
    ~Partitioning();

private:
    /**
     * Types of partitioning algorithms.
     */
    enum PartitionTool {
        PETSc,      ///< Use PETSc interface to various partitioing tools.
        METIS       ///< Use direct interface to Metis.
    };

    /**
     * Types of weights used for element partitioning.
     */
    enum PartitionGraphType {
        any_neighboring,               ///< Add edge for any pair of neighboring elements
        any_wight_lower_dim_cuts,      ///< Same as before and assign higher weight to cuts of lower dimension in order to make them stick to one face
        same_dimension_neghboring,     ///< Add edge for any pair of neighboring elements of same dimension (bad for matrix multiply)
    };

    /// The input mesh
    Mesh        *mesh_;
    /// Input Record accessor.
    Input::Record in_;

    /// Graph used to partitioning the mesh.
    SparseGraph *graph_;
    /// Partition numbers for local elements in original distribution of elements given be @p init_el_ds_.
    int *       loc_part_;
    /// Original distribution of elements. Depends on type of partitioner
    Distribution *init_el_ds_;

    /**
     * Creates sparse parallel graph from the mesh (using algorithm given by the key "graph_type" of the input record accessor @p in_
     */
    void make_element_connection_graph();
    /**
     * # Creates sparse parallel graph from the mesh (using algorithm given by the key "graph_type" of the input record accessor @p in_
     * # calls partitioning tool given by the key "tool" of the input record accessor @p in_
     * # result is local part of partition
     */
    void make_partition();

};


#endif /* PARTITIONING_HH_ */
