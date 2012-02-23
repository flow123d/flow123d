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

using namespace std;

template<unsigned int dim, unsigned int spacedim> class FiniteElement;
class Mesh;

/**
 * Class DOFHandler distributes the degrees of freedom (dof) for
 * a particular finite element on the computational mesh
 * and provides mappings between local and global dofs.
 * The template parameter @p dim denotes the spatial dimension of
 * the reference finite element.
 */
template<unsigned int dim, unsigned int spacedim>
class DOFHandler {
public:

    /**
     * Constructor.
     */
    DOFHandler(Mesh &_mesh);

    /**
     * Alias for iterator over cells
     *
     * TODO: Notation to be fixed: element or cell
     * TODO: Iterator goes through cells of all dimensions, but
     * should go only through dim-dimensional ones.
     */
    typedef ElementFullIter CellIterator;

    /**
     * Distributes degrees of freedom on the mesh needed for the
     * given finite element. The additional parameter @p offset
     * allows to reserve space for another finite element dofs in the
     * beginning of the global dof vector.
     *
     */
    void distribute_dofs(FiniteElement<dim,spacedim> &fe, const unsigned int offset = 0);

    /**
     * Getter for the number of dofs at a single cell. This value
     * depends on the given finite element.
     */
    const unsigned int n_local_dofs();

    /**
     * Getter for the number of all mesh dofs required by the given
     * finite element.
     */
    const unsigned int n_global_dofs();

    /**
     * Returns the number of the first global dof handled by this
     * DOFHandler.
     * @return
     */
    const unsigned int offset();

    /**
     * Returns the global indices of dofs associated to the @p cell.
     */
    void get_dof_indices(const CellIterator &cell, unsigned int indices[]);

    /**
     * Returns the dof values associated to the @p cell.
     */
    void get_dof_values(const CellIterator &cell, const Vec &values,
            double local_values[]);

    /**
     * Returns the index of the dof (specified by @p local_dof_id and
     * the iterator to the @p cell) in the global vector of dofs.
     */
    const unsigned int global_dof_id(const CellIterator &cell,
            const unsigned int local_dof_id);

    /**
     * Return the iterator to the first cell.
     */
    CellIterator begin_cell() const;

    /**
     * Return the iterator to the last cell.
     */
    CellIterator end_cell() const;

    ~DOFHandler();

private:

//    template<unsigned int obj_dim> void get_object_dof_indices(const CellIterator &cell, vector<unsigned int> &indices);

    /*
     * Index of first global dof. Positive value indicates that the
     * first @p global_dof_offset entries in the global dof vector
     * are reserved for a different DOFHandler.
     */
    unsigned int global_dof_offset;

    /*
     * Number of global dofs assigned by the handler.
     */
    unsigned int n_dofs;

    /*
     * Pointer to the mesh to which the dof handler is associated.
     */
    Mesh *mesh;

    /*
     * Pointer to the finite element class for which the handler
     * distributes dofs.
     */
    FiniteElement<dim,spacedim> *finite_element;

    /*
     * Global numbers of dofs associated to nodes (object_dofs[0]),
     * 1D edges (object_dofs[1]), 2D faces (object_difs[2]) and
     * volumes (object_dofs[3]).
     */
    map<void*, int*> object_dofs[dim + 1];

};

#endif /* DOFHANDLER_HH_ */
