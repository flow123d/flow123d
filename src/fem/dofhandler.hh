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
#include "mesh/mesh_types.hh"  // for ElementFullIter
#include "mesh/long_idx.hh"    // for LongIdx
#include "petscvec.h"          // for Vec

template<unsigned int dim> class FiniteElement;
class Mesh;
class Distribution;


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
    DOFHandlerBase(Mesh &_mesh) : global_dof_offset(0), n_dofs(0), lsize_(0), mesh_(&_mesh) {};

    /**
     * @brief Alias for iterator over cells.
     *
     * TODO: Notation to be fixed: element or cell
     * TODO: Iterator goes through cells of all dimensions, but
     * should go only through dim-dimensional ones.
     */
    typedef ElementAccessor<3> CellIterator;

    /**
     * @brief Getter for the number of all mesh dofs required by the given
     * finite element.
     */
    const unsigned int n_global_dofs() const { return n_dofs; }

    /**
     * @brief Returns the number of the first global dof handled by this
     * DOFHandler.
     */
    const unsigned int offset() const { return global_dof_offset; }

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

    /**
     * @brief Returns the parallel distribution of dofs.
     */
    Distribution *distr() const { return ds_; }

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
    virtual unsigned int get_dof_indices(const CellIterator &cell, std::vector<LongIdx> &indices) const = 0;

    /**
     * @brief Fill vector of the indices of dofs associated to the @p cell on the local process.
     *
     * @param cell The cell.
     * @param indices Vector of dof indices on the cell.
     */
    virtual unsigned int get_loc_dof_indices(const CellIterator &cell, std::vector<LongIdx> &indices) const =0;
    
    /**
     * @brief Returns the dof values associated to the @p cell.
     *
     * @param cell The cell.
     * @param values The global vector of values.
     * @param local_values Array of values at local dofs.
     */
    virtual void get_dof_values(const CellIterator &cell, const Vec &values,
            double local_values[]) const = 0;

    /**
     * @brief Compute hash value of DOF handler.
     */
    virtual std::size_t hash() const =0;

    /// Destructor.
    virtual ~DOFHandlerBase() {};

protected:

    /**
     * @brief Index of first global dof.
     *
     * Positive value indicates that the first @p global_dof_offset
     * entries in the global dof vector are reserved for a different
     * DOFHandler.
     */
    unsigned int global_dof_offset;

    /**
     * @brief Number of global dofs assigned by the handler.
     */
    unsigned int n_dofs;

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
    Distribution *ds_;

};




///**
// * @brief Provides the numbering of the finite element degrees of freedom
// * on the computational mesh.
// *
// * Class DOFHandler distributes the degrees of freedom (dof) for
// * a particular finite element on the computational mesh
// * and provides mappings between local and global dofs.
// * The template parameter @p dim denotes the spatial dimension of
// * the reference finite element.
// *
// * Currently the functionality is restricted to discontinuous
// * finite elements, i.e. when the neighboring elements do not
// * share any common dof.
// */
//template<unsigned int dim, unsigned int spacedim>
//class DOFHandler : public DOFHandlerBase {
//public:
//
//    /**
//     * @brief Constructor.
//     * @param _mesh The mesh.
//     */
//    DOFHandler(Mesh &_mesh);
//
//    /**
//     * @brief Alias for iterator over cells.
//     *
//     * TODO: Notation to be fixed: element or cell
//     * TODO: Iterator goes through cells of all dimensions, but
//     * should go only through dim-dimensional ones.
//     */
//    typedef ElementAccessor<3> CellIterator;
//
//    /**
//     * @brief Distributes degrees of freedom on the mesh needed
//     * for the given finite element.
//     *
//     * The additional parameter @p offset allows to reserve space
//     * for another finite element dofs in the beginning of the
//     * global dof vector.
//     *
//     * @param fe The finite element.
//     * @param offset The offset.
//     */
//    void distribute_dofs(FiniteElement<dim,spacedim> &fe, const unsigned int offset = 0);
//
//    /**
//     * @brief Getter for the number of dofs at a single cell.
//     *
//     * This value depends on the given finite element.
//     */
//    const unsigned int n_local_dofs();
//
//    /**
//     * @brief Returns the global indices of dofs associated to the @p cell.
//     *
//     * @param cell The cell.
//     * @param indices Array of dof indices on the cell.
//     */
//    void get_dof_indices(const CellIterator &cell, unsigned int indices[]);
//
//    /**
//     * @brief Returns the dof values associated to the @p cell.
//     *
//     * @param cell The cell.
//     * @param values The global vector of values.
//     * @param local_values Array of values at local dofs.
//     */
//    void get_dof_values(const CellIterator &cell, const Vec &values,
//            double local_values[]);
//
//    /// Destructor.
//    ~DOFHandler();
//
//private:
//
//    /**
//     * @brief Pointer to the finite element class for which the handler
//     * distributes dofs.
//     */
//    FiniteElement<dim,spacedim> *finite_element;
//
//    /**
//     * @brief Number of dofs associated to geometrical entities.
//     *
//     * Global numbers of dofs associated to nodes (object_dofs[0]),
//     * 1D edges (object_dofs[1]), 2D faces (object_difs[2]) and
//     * volumes (object_dofs[3]).
//     */
//    int ***object_dofs;
//
//};



class DOFHandlerMultiDim : public DOFHandlerBase {
public:

    /**
     * @brief Constructor.
     * @param _mesh The mesh.
     */
    DOFHandlerMultiDim(Mesh &_mesh, bool make_elem_part = true);

    /**
     * @brief Alias for iterator over cells.
     *
     * TODO: Notation to be fixed: element or cell
     */
    typedef ElementAccessor<3> CellIterator;

    /**
     * @brief Distributes degrees of freedom on the mesh needed
     * for the given finite elements.
     *
     * The additional parameter @p offset allows to reserve space
     * for another finite element dofs in the beginning of the
     * global dof vector.
     *
     * @param fe1d The 1D finite element.
     * @param fe2d The 2D finite element.
     * @param fe3d The 3D finite element.
     * @param offset The offset.
     */
    void distribute_dofs(FiniteElement<1> &fe1d,
    		FiniteElement<2> &fe2d,
    		FiniteElement<3> &fe3d,
    		const unsigned int offset = 0);

    /**
     * @brief Returns the global indices of dofs associated to the @p cell.
     *
     * @param cell The cell.
     * @param indices Array of dof indices on the cell.
     */
    unsigned int get_dof_indices(const CellIterator &cell, std::vector<LongIdx> &indices) const override;
    
    /**
     * @brief Returns the indices of dofs associated to the @p cell on the local process.
     *
     * @param cell The cell.
     * @param indices Array of dof indices on the cell.
     */
    unsigned int get_loc_dof_indices(const CellIterator &cell, std::vector<LongIdx> &indices) const override;

    /**
     * @brief Returns the dof values associated to the @p cell.
     *
     * @param cell The cell.
     * @param values The global vector of values.
     * @param local_values Array of values at local dofs.
     */
    void get_dof_values(const CellIterator &cell, const Vec &values,
            double local_values[]) const override;

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

    /// Returns finite element object for given space dimension.
    template<unsigned int dim>
    FiniteElement<dim> *fe() const;

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
     * @brief Pointer to the finite element class for which the handler
     * distributes dofs.
     */
    FiniteElement<1> *fe1d_;
    FiniteElement<2> *fe2d_;
    FiniteElement<3> *fe3d_;

    /**
     * @brief Number of dofs associated to geometrical entities.
     *
     * Global numbers of dofs associated to nodes (object_dofs[0]),
     * 1D edges (object_dofs[1]), 2D faces (object_difs[2]) and
     * volumes (object_dofs[3]).
     */
    LongIdx ***object_dofs;


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

};




#endif /* DOFHANDLER_HH_ */
