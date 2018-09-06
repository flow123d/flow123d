/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    dofhandler.hh
 * @brief   Declaration of class which handles the ordering of degrees of freedom (dof) and mappings between local and global dofs.
 * @author  Jan Stebel
 */

#ifndef DOFHANDLER_HH_
#define DOFHANDLER_HH_

#include <vector>              // for vector
#include "mesh/side_impl.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
// #include "mesh/mesh_types.hh"  // for ElementFullIter
#include "mesh/long_idx.hh"    // for LongIdx
#include "fem/discrete_space.hh" // for DiscreteSpace
#include "petscvec.h"          // for Vec

template<unsigned int dim> class FiniteElement;
class Mesh;
class Distribution;
class Dof;


/**
 * Class DOFHandlerBase provides an abstract interface for various dof handlers:
 * - basic handler for a given spatial dimension
 * - multi-dimensional handler for all spatial dimensions (1D-3D)
 * - handler for a specific region / part of mesh
 */
class DOFHandlerBase {
public:

    /**
     * @brief Constructor.
     * @param _mesh The mesh.
     */
    DOFHandlerBase(Mesh &_mesh)
    : n_global_dofs_(0), lsize_(0), loffset_(0), max_elem_dofs_(0), mesh_(&_mesh), dof_ds_(0) {};

    /**
     * @brief Getter for the number of all mesh dofs required by the given
     * finite element.
     */
    const unsigned int n_global_dofs() const { return n_global_dofs_; }
    
    /**
     * @brief Returns the number of dofs on the current process.
     */
    const unsigned int lsize() const { return lsize_; }

    /**
     * @brief Returns the offset of the local part of dofs.
     */
    const unsigned int loffset() const { return loffset_; }
    
    /**
     * @brief Returns max. number of dofs on one element.
     */
    const unsigned int max_elem_dofs() const { return max_elem_dofs_; }

    Distribution *distr() const { return dof_ds_; }

    /**
     * @brief Returns the mesh.
     */
    Mesh *mesh() const { return mesh_; }

    /**
     * @brief Fill vector of the global indices of dofs associated to the @p cell.
     *
     * @param cell The cell.
     * @param indices Vector of dof indices on the cell.
     */
    virtual unsigned int get_dof_indices(const ElementAccessor<3> &cell, std::vector<LongIdx> &indices) const = 0;

    /**
     * @brief Fill vector of the indices of dofs associated to the @p cell on the local process.
     *
     * @param cell The cell.
     * @param indices Vector of dof indices on the cell.
     */
    virtual unsigned int get_loc_dof_indices(const ElementAccessor<3> &cell, std::vector<LongIdx> &indices) const =0;
    
    /**
     * @brief Compute hash value of DOF handler.
     */
    virtual std::size_t hash() const =0;

    /// Destructor.
    virtual ~DOFHandlerBase();

protected:

    /**
     * @brief Number of global dofs assigned by the handler.
     */
    unsigned int n_global_dofs_;
    
    /**
     * @brief Number of dofs associated to local process.
     */
    unsigned int lsize_;

    /**
     * @brief Index of the first dof on the local process.
     */
    unsigned int loffset_;

    /// Max. number of dofs per element.
    unsigned int max_elem_dofs_;

    /**
     * @brief Pointer to the mesh to which the dof handler is associated.
     */
    Mesh *mesh_;

    /**
     * @brief Distribution of dofs associated to local process.
     */
     Distribution *dof_ds_;

};




/**
 * @brief Provides the numbering of the finite element degrees of freedom
 * on the computational mesh.
 *
 * Class DOFHandlerMultiDim distributes the degrees of freedom (dof) for
 * a particular triplet of 1d, 2d and 3d finite elements on the computational mesh
 * and provides mappings between local and global dofs.
 * The template parameter @p dim denotes the spatial dimension of
 * the reference finite element.
 *
 * Currently the functionality is restricted to finite elements with internal and nodal dofs,
 * i.e. the neighboring elements can share only dofs on nodes.
 */
class DOFHandlerMultiDim : public DOFHandlerBase {
public:

    /**
     * @brief Constructor.
     * @param _mesh The mesh.
     */
    DOFHandlerMultiDim(Mesh &_mesh);


    /**
     * @brief Distributes degrees of freedom on the mesh needed
     * for the given discrete space.
     *
     * By default, the dof handler is parallel, meaning that each
     * processor has access to dofs on the local elements and on one
     * layer of ghost elements (owned by neighbouring elements).
     * This can be changed by setting @p sequential to true.
     *
     * @param ds         The discrete space consisting of finite elements for each mesh element.
     * @param sequential If true then each processor will have information about all dofs.
     */
    void distribute_dofs(std::shared_ptr<DiscreteSpace> ds,
                         bool sequential = false);

    /**
     * @brief Returns the global indices of dofs associated to the @p cell.
     *
     * @param cell The cell.
     * @param indices Array of dof indices on the cell.
     */
    unsigned int get_dof_indices(const ElementAccessor<3> &cell,
                                 std::vector<LongIdx> &indices) const override;
    
    /**
     * @brief Returns the indices of dofs associated to the @p cell on the local process.
     *
     * @param cell The cell.
     * @param indices Array of dof indices on the cell.
     */
    unsigned int get_loc_dof_indices(const ElementAccessor<3> &cell,
                                     std::vector<LongIdx> &indices) const override;

    /**
     * @brief Returns the global index of local element.
     *
     * @param loc_el Local index of element.
     */
    inline int el_index(int loc_el) const { return el_4_loc[loc_el]; }

    /**
     * @brief Returns the global index of local edge.
     *
     * @param loc_edg Local index of edge.
     */
    inline LongIdx edge_index(int loc_edg) const { return edg_4_loc[loc_edg]; }

    /**
	 * @brief Returns the global index of local neighbour.
	 *
	 * @param loc_nb Local index of neighbour.
	 */
	inline LongIdx nb_index(int loc_nb) const { return nb_4_loc[loc_nb]; }
	
	/**
     * @brief Return number of dofs on given cell.
     * 
     * @param cell Cell accessor.
     */
	unsigned int n_dofs(ElementAccessor<3> cell) const;

	/**
	 * @brief Returns number of local edges.
	 */
    inline unsigned int n_loc_edges() const { return edg_4_loc.size(); }

    /**
     * @brief Returns number of local neighbours.
     */
    inline unsigned int n_loc_nb() const { return nb_4_loc.size(); }

    /**
     * Returns true if element is on local process.
     * @param index Global element index.
     */
    bool el_is_local(int index) const;

    /**
     * @brief Returns finite element object for given space dimension.
     * 
     * @param cell Cell accessor.
     */
    template<unsigned int dim>
    FiniteElement<dim> *fe(const ElementAccessor<3> &cell) const { return ds_->fe<dim>(cell); }
    
    /**
     * @brief Return dof on a given cell.
     * @param cell Mesh cell.
     * @param idof Number of dof on the cell.
     */
    const Dof &cell_dof(ElementAccessor<3> cell,
                        unsigned int idof) const;

    /**
     * Implements @p DOFHandlerBase::hash.
     */
    std::size_t hash() const override;

    /// Destructor.
    ~DOFHandlerMultiDim() override;
    
    
    

private:

    /**
     * @brief Prepare parallel distribution of elements, edges and neighbours.
     */
    void make_elem_partitioning();
    
    /**
     * @brief Initialize vector of starting indices for elements.
     */
    void init_cell_starts();
    
    /**
     * @brief Initialize auxiliary vector of starting indices of nodal dofs.
     * 
     * @param node_dof_starts Vector of starting indices (output).
     */
    void init_node_dof_starts(std::vector<LongIdx> &node_dof_starts);
    
    /**
     * @brief Initialize node_status.
     * 
     * Set VALID_NODE for nodes owned by local elements and
     * INVALID_NODE for nodes owned by ghost elements.
     * 
     * @param node_status Vector of nodal status (output).
     */
    void init_node_status(std::vector<short int> &node_status);
    
    /**
     * @brief Obtain dof numbers on ghost elements from other processor.
     * @param proc  Neighbouring processor.
     * @param dofs  Array where dofs are stored (output).
     */
    void receive_ghost_dofs(unsigned int proc,
                            std::vector<LongIdx> &dofs);

    /**
     * @brief Send dof numbers to other processor.
     * @param proc  Neighbouring processor.
     */    
    void send_ghost_dofs(unsigned int proc);
    
    /** 
     * @brief Update dofs on local elements from ghost element dofs.
     * 
     * @param proc            Neighbouring processor.
     * @param update_cells    Vector of global indices of elements which need to be updated
     *                        from ghost elements.
     * @param dofs            Vector of dof indices on ghost elements from processor @p proc.
     * @param node_dof_starts Vector of starting indices of nodal dofs.
     * @param node_dofs       Vector of nodal dof indices (output).
     */
    void update_local_dofs(unsigned int proc,
                           const std::vector<bool> &update_cells,
                           const std::vector<LongIdx> &dofs,
                           const std::vector<LongIdx> &node_dof_starts,
                           std::vector<LongIdx> &node_dofs
                          );
    
    /**
     * @brief Communicate local dof indices to all processors.
     */
    void create_sequential();

    
    /**
     * Flags used during distribution of dofs to mark node and dof status.
     */
    static const int INVALID_NODE  = 1;
    static const int VALID_NODE    = 2;
    static const int ASSIGNED_NODE = 3;
    static const int INVALID_DOF   = -1;
    
    
    /// Pointer to the discrete space for which the handler distributes dofs.
    std::shared_ptr<DiscreteSpace> ds_;

    /**
     * @brief Starting indices for element dofs (parallel version).
     * 
     * E.g. dof_indices[cell_starts[idx]] = dof number for first dof on the
     * cell with index idx within the paralle structure. To use with element
     * accessor use the following:
     * 
     *   ElementAccessor<3> cell;
     *   ...
     *   // i-th dof number on the cell
     *   dof_indices[cell_starts[row_4_el[cell.idx()]]+i] = ...
     * 
     * Only local and ghost elements are stored, but the vector has size mesh_->n_elements()+1.
     */
    std::vector<LongIdx> cell_starts;
    
    /**
     * @brief Dof numbers on local and ghost elements (parallel version).
     * 
     * Dofs are ordered accordingly with cell_starts and local dof order
     * given by the finite element. See cell_starts for more description.
     */
    std::vector<LongIdx> dof_indices;
    
    /**
     * @brief Starting indices for element dofs (sequential version).
     * 
     * This vector stores information about all mesh elements.
     * See cell_starts for paralle version.
     */
    std::vector<LongIdx> cell_starts_seq;
    
    /**
     * @brief Dof numbers on mesh elements (sequential version).
     * 
     * This vector stores information about all mesh elements.
     * See dof_indices for paralle version.
     */
    std::vector<LongIdx> dof_indices_seq;


	/// Global element index -> index according to partitioning
    LongIdx *row_4_el;
    
    /// Local element index -> global element index
    LongIdx *el_4_loc;
    
    /// Distribution of elements
    Distribution *el_ds_;

    /// Local edge index -> global edge index
    vector<LongIdx> edg_4_loc;

    /// Local neighbour index -> global neighbour index
    vector<LongIdx> nb_4_loc;
    
    /// Indices of local nodes in mesh tree.
    vector<LongIdx> node_4_loc;
    
    /// Indices of ghost cells (neighbouring with local elements).
    vector<LongIdx> ghost_4_loc;
    
    /// Processors of ghost elements.
    set<unsigned int> ghost_proc;
    
    /// Arrays of ghost cells for each neighbouring processor.
    map<unsigned int, vector<LongIdx> > ghost_proc_el;

};




#endif /* DOFHANDLER_HH_ */
