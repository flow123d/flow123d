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
    static Input::Type::Selection & get_graph_type_sel();
    static Input::Type::Selection & get_tool_sel();
    static Input::Type::Record & get_input_type();

    /**
     *  Constructor. A pointer to the mesh and accessor to an input record have to be provided.
     */
    Partitioning(Mesh *mesh, Input::Record in);

    /**
     * Get initial distribution.
     */
    const Distribution *get_init_distr() const;
    /**
     * Get local part of mesh partition.
     */
    const int *get_loc_part() const;

    /**
     * Creates and returns vector with element partitioning for output.
     */
    shared_ptr< vector<int> > subdomain_id_field_data();

    /**
     * Obsolete see source file for doc.
     */
    void id_maps(int n_ids, int *id_4_old,
                    Distribution * &new_ds, int * &id_4_loc, int * &new_4_id);


    static void id_maps(int n_ids, int *id_4_old,
            const Distribution &old_ds, int *loc_part,
            Distribution * &new_ds, int * &id_4_loc, int * &new_4_id);

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
        any_weight_lower_dim_cuts,      ///< Same as before and assign higher weight to cuts of lower dimension in order to make them stick to one face
        same_dimension_neighboring,     ///< Add edge for any pair of neighboring elements of same dimension (bad for matrix multiply)
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
    /// Sequential partitioning for output.
    shared_ptr< vector<int> > seq_part_;

    /**
     * Creates sparse parallel graph from the mesh (using algorithm given by the key "graph_type" of the input record accessor @p in_
     */
    void make_element_connection_graph();
    /**
     * Creates sparse parallel graph from the mesh (using algorithm given by the key "graph_type" of the input record accessor @p in_)
     * calls partitioning tool given by the key "tool" of the input record accessor @p in_)
     * result is local part of the partitioning. Can be retrieved by @p get_loc_part().
     */
    void make_partition();


};


#endif /* PARTITIONING_HH_ */
