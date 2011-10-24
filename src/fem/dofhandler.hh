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
 * $Id: quadrature.hh 1352 2011-09-23 14:14:47Z jan.stebel $
 * $Revision: 1352 $
 * $LastChangedBy: jan.stebel $
 * $LastChangedDate: 2011-09-23 16:14:47 +0200 (Fri, 23 Sep 2011) $
 *
 * @file
 * @brief Declaration of class which handles the ordering of degrees of freedom (dof) and mappings between local and global dofs.
 *  @author Jan Stebel
 */

#ifndef DOFHANDLER_HH_
#define DOFHANDLER_HH_

#include <vector>
#include "mesh/mesh_types.hh"

using namespace std;

template <int> class FiniteElement;
class Mesh;

template <int dim>
class DOFHandler
{
public:

    DOFHandler(Mesh &_mesh);

    /*
     * Alias for iterator over cells
     *
     * TODO: Notation to be fixed: element or cell
     * TODO: Iterator goes through cells of all dimensions, but should go only through dim-dimensional ones.
     */
    typedef ElementFullIter CellIterator;

    /*
     * Go through the mesh and distribute degrees of freedom needed for the given finite element.
     * The additional parameter @p offset allows to reserve space for another finite element dofs in the beginning of the global dof vector.
     *
     */
    void distribute_dofs(const FiniteElement<dim> &fe, const unsigned int offset = 0);

    /*
     * Getter for the number of dofs at a single cell. This value depends on the given finite element.
     */
    const unsigned int n_local_dofs();

    /*
     * Getter for the number of all mesh dofs required by the given finite element.
     */
    const unsigned int n_global_dofs();

    /*
     * Returns the index of the dof (specified by @p local_dof_id and the iterator to the @p cell) in the global vector of dofs.
     */
    const unsigned int global_dof_id(const CellIterator &cell, const unsigned int local_dof_id);

//    /*
//     * Returns the index of the node to which a dof is associated.
//     */
//    const unsigned int node_id(const unsigned int global_dof_id);

    /*
     * Return the iterator to the first cell.
     */
    CellIterator begin_cell() const;

    /*
     * Return the iterator to the last cell.
     */
    CellIterator end_cell() const;

    ~DOFHandler();

private:

    /*
     * Index of first global dof. Positive value indicates that the first @p global_dof_offset entries in the global dof vector are reserved for a different DOFHandler.
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
     * Pointer to the finite element class for which the handler distributes dofs.
     */
    FiniteElement<dim> *finite_element;

    /*
     * Global numbers of cell dofs (cell_dof_ids[icell][idof] = id of global dof corresponding to idof'th local dof on cell icell).
     */
    vector<int*> cell_dof_ids;

};





#endif /* DOFHANDLER_HH_ */
