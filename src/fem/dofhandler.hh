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
 * @brief Declaration of class which handles the ordering of degrees of freedom (dof) and mappings between local and global dofs.
 * @author Jan Stebel
 */

#ifndef DOFHANDLER_HH_
#define DOFHANDLER_HH_

#include <map>
#include <petscmat.h>
#include "mesh/mesh_types.hh"
#include "mesh/elements.h"


template<unsigned int dim, unsigned int spacedim> class FiniteElement;
class Mesh;

/**
 * @brief Provides the numbering of the finite element degrees of freedom
 * on the computational mesh.
 *
 * Class DOFHandler distributes the degrees of freedom (dof) for
 * a particular finite element on the computational mesh
 * and provides mappings between local and global dofs.
 * The template parameter @p dim denotes the spatial dimension of
 * the reference finite element.
 *
 * Currently the functionality is restricted to discontinuous
 * finite elements, i.e. when the neighboring elements do not
 * share any common dof.
 */
template<unsigned int dim, unsigned int spacedim>
class DOFHandler {
public:

    /**
     * @brief Constructor.
     * @param _mesh The mesh.
     */
    DOFHandler(Mesh &_mesh);

    /**
     * @brief Alias for iterator over cells.
     *
     * TODO: Notation to be fixed: element or cell
     * TODO: Iterator goes through cells of all dimensions, but
     * should go only through dim-dimensional ones.
     */
    typedef ElementFullIter CellIterator;

    /**
     * @brief Distributes degrees of freedom on the mesh needed
     * for the given finite element.
     *
     * The additional parameter @p offset allows to reserve space
     * for another finite element dofs in the beginning of the
     * global dof vector.
     *
     * @param fe The finite element.
     * @param offset The offset.
     */
    void distribute_dofs(FiniteElement<dim,spacedim> &fe, const unsigned int offset = 0);

    /**
     * @brief Getter for the number of dofs at a single cell.
     *
     * This value depends on the given finite element.
     */
    const unsigned int n_local_dofs();

    /**
     * @brief Getter for the number of all mesh dofs required by the given
     * finite element.
     */
    const unsigned int n_global_dofs();

    /**
     * @brief Returns the number of the first global dof handled by this
     * DOFHandler.
     */
    const unsigned int offset();

    /**
     * @brief Returns the global indices of dofs associated to the @p cell.
     *
     * @param cell The cell.
     * @param indices Array of dof indices on the cell.
     */
    void get_dof_indices(const CellIterator &cell, unsigned int indices[]);

    /**
     * @brief Returns the dof values associated to the @p cell.
     *
     * @param cell The cell.
     * @param values The global vector of values.
     * @param local_values Array of values at local dofs.
     */
    void get_dof_values(const CellIterator &cell, const Vec &values,
            double local_values[]);

    /**
     * @brief Returns the global number of a local dof.
     *
     * Returns the index of the dof (specified by @p local_dof_id and
     * the iterator to the @p cell) in the global vector of dofs.
     *
     * @param cell The cell.
     * @param local_dof_id Id of the local dof.
     */
    const unsigned int global_dof_id(const CellIterator &cell,
            const unsigned int local_dof_id);

    /**
     * @brief Return the iterator to the first cell.
     */
    CellIterator begin_cell() const;

    /**
     * @brief Return the iterator to the last cell.
     */
    CellIterator end_cell() const;

    /// Destructor.
    ~DOFHandler();

private:

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
     * @brief Pointer to the mesh to which the dof handler is associated.
     */
    Mesh *mesh;

    /**
     * @brief Pointer to the finite element class for which the handler
     * distributes dofs.
     */
    FiniteElement<dim,spacedim> *finite_element;

    /**
     * @brief Number of dofs associated to geometrical entities.
     *
     * Global numbers of dofs associated to nodes (object_dofs[0]),
     * 1D edges (object_dofs[1]), 2D faces (object_difs[2]) and
     * volumes (object_dofs[3]).
     */
    map<void*, int*> object_dofs[dim + 1];

};

#endif /* DOFHANDLER_HH_ */
