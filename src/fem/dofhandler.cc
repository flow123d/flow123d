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



#include "fem/dofhandler.hh"
#include "fem/finite_element.hh"
#include "mesh/mesh.h"

#include "fem/fe_p.hh"

FE_P<3,0> fe;



template<unsigned int dim> inline DOFHandler<dim>::DOFHandler(Mesh & _mesh)
: mesh(&_mesh),
  n_dofs(0)
{
}



template<unsigned int dim> inline void DOFHandler<dim>::distribute_dofs(const FiniteElement<dim> & fe, const unsigned int offset)
{
    unsigned int next_free_dof = offset;
    // remember ids of global dofs assigned to nodes
    int node_dof_ids[mesh->node_vector.size()];

    for (int i=0; i<mesh->node_vector.size(); i++) node_dof_ids[i] = -1;

    // TODO: Maybe check if dofs are not yet distributed?

    finite_element = &fe;
    global_dof_offset = offset;

    FOR_ELEMENTS(mesh,cell)
    {
        // skip cells of different dimension
        if (cell->dim != dim) continue;

        cell_dof_ids[cell.id()] = new int[fe->n_dofs()];

        // distribute dofs
        for (int d=0; d<fe.n_dofs(); d++)
        {
            // check whether new dof has to be allocated or we can use an existing one
            // (this applies only to "continuous" finite elements)
            if (fe.dof_is_continuous(d))
            {
                switch (fe.dof_type(d))
                {
                case FE_OBJECT_POINT:
                    // dof is sitting on a node (vertex)
                    // TODO: Class Node has to be modified so that mesh nodes are assigned unique id numbers.
                    int v_id = cell->node[fe.dof_object_id(d)]->id();
                    if (node_dof_ids[v_id] == -1)
                    {
                        node_dof_ids[v_id] = next_free_dof++;
                    }
                    cell_dof_ids[cell.id()][d] = node_dof_ids[v_id];
                    break;
                }
            }
            else
            {
                cell_dof_ids[cell.id()][d] = next_free_dof++;
            }
        }
    }

    n_dofs = next_free_dof - offset;
}



template<unsigned int dim> inline const unsigned int DOFHandler<dim>::n_local_dofs()
{
    return finite_element->n_dofs();
}



template<unsigned int dim> inline const unsigned int DOFHandler<dim>::n_global_dofs()
{
    return n_dofs;
}



template<unsigned int dim> inline const unsigned int DOFHandler<dim>::global_dof_id(const CellIterator &cell, const unsigned int local_dof_id)
{
    return cell_dof_ids[cell.id()][local_dof_id];
}



template<unsigned int dim> inline typename DOFHandler<dim>::CellIterator DOFHandler<dim>::begin_cell() const
{
    return mesh->element.begin();
}



template<unsigned int dim> inline typename DOFHandler<dim>::CellIterator DOFHandler<dim>::end_cell() const
{
    return mesh->element.end();
}

template<unsigned int dim> inline DOFHandler<dim>::~DOFHandler()
{
    for (vector<int*>::iterator icell = cell_dof_ids.begin(); icell != cell_dof_ids.end(); icell++)
        delete[] icell;
}

