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

#include <vector>                   // for vector
#include <unordered_map>            // for unordered_map
#include "system/index_types.hh"    // for LongIdx
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "tools/general_iterator.hh"
#include "fem/discrete_space.hh" // for DiscreteSpace
#include "la/vector_mpi.hh"       // for VectorMPI
#include "petscvec.h"          // for Vec


template<unsigned int dim> class FiniteElement;
class DHCellAccessor;
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
    : n_global_dofs_(0), lsize_(0), loffset_(0), max_elem_dofs_(0), mesh_(&_mesh), dof_ds_(0) {}

    /**
     * @brief Getter for the number of all mesh dofs required by the given
     * finite element.
     */
    unsigned int n_global_dofs() const { return n_global_dofs_; }
    
    /**
     * @brief Returns the number of dofs on the current process.
     */
    unsigned int lsize() const { return lsize_; }

    /**
     * @brief Returns max. number of dofs on one element.
     */
    unsigned int max_elem_dofs() const { return max_elem_dofs_; }

    std::shared_ptr<Distribution> distr() const { return dof_ds_; }

    /**
     * @brief Returns the mesh.
     */
    Mesh *mesh() const { return mesh_; }

    /**
     * @brief Compute hash value of DOF handler.
     */
    virtual std::size_t hash() const =0;

    /// Destructor.
    virtual ~DOFHandlerBase();

protected:
    
    /**
     * @brief Fill vector of the global indices of dofs associated to the @p cell.
     *
     * @param cell The cell.
     * @param indices Vector of dof indices on the cell.
     */
    virtual unsigned int get_dof_indices(const DHCellAccessor &cell, std::vector<LongIdx> &indices) const = 0;

    /**
     * @brief Returns a vector of the indices of dofs associated to the @p cell on the local process.
     *
     * @param loc_ele_idx local element index.
     */
    virtual LocDofVec get_loc_dof_indices(unsigned int loc_ele_idx) const =0;

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
     std::shared_ptr<Distribution> dof_ds_;

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
     * @param make_elem_part Allow switch off make_element_partitioning, necessary for boundary DOF handler.
     */
    DOFHandlerMultiDim(Mesh &_mesh, bool make_elem_part = true);


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
    void distribute_dofs(std::shared_ptr<DiscreteSpace> ds);

    /** @brief Returns sequential version of the current dof handler.
     * 
     * Collective on all processors.
     */
    std::shared_ptr<DOFHandlerMultiDim> sequential();
    
    /**
     * @brief Returns scatter context from parallel to sequential vectors.
     * 
     * For sequential dof handler it returns null pointer.
     * Collective on all processors.
     */
    std::shared_ptr<VecScatter> sequential_scatter();
    
    /**
     * @brief Allocates PETSc vector according to the dof distribution.
     */
    virtual VectorMPI create_vector();
    
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
	 * @brief Returns number of local edges.
	 */
    inline unsigned int n_loc_edges() const { return edg_4_loc.size(); }

    /**
     * @brief Returns number of local neighbours.
     */
    inline unsigned int n_loc_nb() const { return nb_4_loc.size(); }

    /// Output structure of dof handler.
    void print() const;

    /**
     * Implements @p DOFHandlerBase::hash.
     */
    std::size_t hash() const override;

    /// Returns range of DOF handler cells (only range of own without ghost cells)
    Range<DHCellAccessor> own_range() const;

    /// Returns range over own and ghost cells of DOF handler
    Range<DHCellAccessor> local_range() const;

    /// Returns range over ghosts DOF handler cells
    Range<DHCellAccessor> ghost_range() const;

    /// Return size of own range (number of own cells)
    inline unsigned int own_size() const {
        return el_ds_->lsize();
    }

    /// Return size of local range (number of local cells)
    inline unsigned int local_size() const {
        return el_ds_->lsize()+ghost_4_loc.size();
    }

    /// Return size of ghost range (number of ghost cells)
    inline unsigned int ghost_size() const {
        return ghost_4_loc.size();
    }

    /// Return DHCellAccessor appropriate to ElementAccessor of given idx
    const DHCellAccessor cell_accessor_from_element(unsigned int elm_idx) const;

    /// Return pointer to discrete space for which the handler distributes dofs.
    std::shared_ptr<DiscreteSpace> ds() const { return ds_; }

    /// Get the map between local dof indices and the global ones.
    const std::vector<LongIdx> & get_local_to_global_map() const { return local_to_global_dof_idx_; }

    /// Destructor.
    ~DOFHandlerMultiDim() override;
    
    
    
    friend class DHCellAccessor;
    friend class DHCellSide;
    friend class DHNeighbSide;
    friend class SubDOFHandlerMultiDim;

protected:

    /**
     * Returns true if element is on local process.
     * @param index Global element index.
     */
    bool el_is_local(int index) const;

    /**
     * @brief Prepare parallel distribution of elements, edges and neighbours.
     */
    void make_elem_partitioning();
    
    /**
     * @brief Initialize vector of starting indices for elements.
     */
    void init_cell_starts();
    
    /**
     * @brief Initialize auxiliary vector of starting indices of nodal/edge dofs.
     * 
     * @param node_dof_starts Vector of nodal starting indices (output).
     * @param edge_dof_starts Vector of edge starting indices (output).
     */
    void init_dof_starts(std::vector<LongIdx> &node_dof_starts,
                         std::vector<LongIdx> &edge_dof_starts);
    
    /**
     * @brief Initialize node_status and edge_status.
     * 
     * Set VALID_NFACE for nodes/edges owned by local elements and
     * INVALID_NFACE for nodes/edges owned by ghost elements.
     * 
     * @param node_status Vector of nodal status (output).
     * @param edge_status Vector of edge status (output).
     */
    void init_status(std::vector<short int> &node_status,
                     std::vector<short int> &edge_status);
    
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
     * @param edge_dof_starts Vector of starting indices of edge dofs.
     * @param edge_dofs       Vector of edge dof indices (output).
     */
    void update_local_dofs(unsigned int proc,
                           const std::vector<bool> &update_cells,
                           const std::vector<LongIdx> &dofs,
                           const std::vector<LongIdx> &node_dof_starts,
                           std::vector<LongIdx> &node_dofs,
                           const std::vector<LongIdx> &edge_dof_starts,
                           std::vector<LongIdx> &edge_dofs);
    
    /**
     * @brief Communicate local dof indices to all processors and create new sequential dof handler.
     *
     * Collective on all processors.
     */
    void create_sequential();
    
    /**
     * @brief Returns the global indices of dofs associated to the @p cell.
     *
     * @param cell The cell.
     * @param indices Array of dof indices on the cell.
     */
    unsigned int get_dof_indices(const DHCellAccessor &cell,
                                 std::vector<LongIdx> &indices) const override;
    
    /**
     * @brief Returns the indices of dofs associated to the @p cell on the local process.
     *
     * @param loc_ele_idx local element index.
     */
    inline LocDofVec get_loc_dof_indices(unsigned int loc_ele_idx) const override
    {
        unsigned int ndofs = cell_starts[loc_ele_idx+1]-cell_starts[loc_ele_idx];
        // create armadillo vector on top of existing array
        // vec(ptr_aux_mem, number_of_elements, copy_aux_mem = true, strict = false)
        Idx* mem_ptr = const_cast<Idx*>(&(dof_indices[cell_starts[loc_ele_idx]]));
        return LocDofVec(mem_ptr, ndofs, false, false);
    }


    
    /**
     * Flags used during distribution of dofs to mark n-face and dof status.
     * 
     * INVALID_NFACE  means that on this node/edge/side/cell the current processor
     *                will not distribute dofs.
     * VALID_NFACE    means that on this node/edge/side/cell the current processor
     *                will distribute dofs.
     * ASSIGNED_NFACE means that dofs on this n-face are already distributed.
     * 
     * INVALID_DOF    marks dofs whose global number has to be set by neighbouring processor.
     */
    static const int INVALID_NFACE;
    static const int VALID_NFACE;
    static const int ASSIGNED_NFACE;
    static const int INVALID_DOF;
    
    
    /// Pointer to the discrete space for which the handler distributes dofs.
    std::shared_ptr<DiscreteSpace> ds_;
    
    /// Indicator for parallel/sequential dof handler.
    bool is_parallel_;
    
    /// Sequential dof handler associated to the current (parallel) one.
    std::shared_ptr<DOFHandlerMultiDim> dh_seq_;
    
    /// Scatter context for parallel to sequential vectors.
    std::shared_ptr<VecScatter> scatter_to_seq_;

    /**
     * @brief Starting indices for local (owned+ghost) element dofs.
     * 
     * E.g. dof_indices[cell_starts[idx]] = dof number for first dof on the
     * cell with local index idx within the parallel structure. To use with
     * DHCellAccessor use the following:
     * 
     *   DHCellAccessor<3> cell;
     *   ...
     *   // i-th dof number on the cell
     *   dof_indices[cell_starts[cell.local_idx()]+i] = ...
     * 
     * For sequential dof handler, dofs from all elements are stored and
     * elements are ordered as in mesh_->get_row_4_el().
     */
    std::vector<LongIdx> cell_starts;
    
    /**
     * @brief Dof numbers on local and ghost elements.
     * 
     * Dofs are ordered accordingly with cell_starts and local dof order
     * given by the finite element. See cell_starts for more description.
     */
    std::vector<Idx> dof_indices;
    
    /**
     * @brief Maps local and ghost dof indices to global ones.
     * 
     * First lsize_ entries correspond to dofs owned by local processor,
     * the remaining entries are ghost dofs sorted by neighbouring processor id.
     */
    std::vector<LongIdx> local_to_global_dof_idx_;
    
    /**
     * @brief Maps global element index into local/ghost index (obsolete).
     */
    std::unordered_map<LongIdx,LongIdx> global_to_local_el_idx_;
    
    /// Distribution of elements
    Distribution *el_ds_;

    /// Local edge index -> global edge index
    vector<LongIdx> edg_4_loc;

    /// Local neighbour index -> global neighbour index
    vector<LongIdx> nb_4_loc;
    
    /// Indices of ghost cells (neighbouring with local elements).
    vector<LongIdx> ghost_4_loc;
    
    /// Processors of ghost elements.
    set<unsigned int> ghost_proc;
    
    /// Arrays of ghost cells for each neighbouring processor.
    map<unsigned int, vector<LongIdx> > ghost_proc_el;

};


class SubDOFHandlerMultiDim : public DOFHandlerMultiDim
{
public:
    
    /** @brief Creates a new dof handler for a component of FESystem.
     * 
     * The @p component_idx indicates the index of finite-element
     * from @p dh for which the new sub- handler is made.
     * The numbering of dofs in sub-handler is compatible
     * with the original handler.
     */
    SubDOFHandlerMultiDim(std::shared_ptr<DOFHandlerMultiDim> dh, unsigned int component_idx);
    
    /** @brief Update values in subvector from parent vector.
     * 
     * @p vec    Vector aligned with the parent dof handler.
     * @p subvec Vctor aligned with the current sub-handler.
     */
    void update_subvector(const VectorMPI &vec, VectorMPI &subvec);
    
    /** @brief Update values in parent vector from values of subvector.
     * 
     * @p vec    Vector aligned with parent dof handler.
     * @p subvec Vector aligned with the current sub-handler.
     */
    void update_parent_vector(VectorMPI &vec, const VectorMPI &subvec);
    
    /// Local indices in the parent handler.
    const std::vector<LongIdx> &parent_indices() { return parent_dof_idx_; }
    
    
private:

    /// Get global dof indices of ghost dofs for sub-handler.
    void receive_sub_ghost_dofs(unsigned int proc, vector<LongIdx> &dofs);
    
    /// Send global indices of dofs that are ghost on other processors.
    void send_sub_ghost_dofs(unsigned int proc, const map<LongIdx,LongIdx> &global_to_local_dof_idx);

    /// Parent dof handler.
    std::shared_ptr<DOFHandlerMultiDim> parent_;
    
    /// Index of FE in parent FESystem.
    unsigned int fe_idx_;
    
    /// Local indices in the parent handler.
    std::vector<LongIdx> parent_dof_idx_;
};




#endif /* DOFHANDLER_HH_ */
