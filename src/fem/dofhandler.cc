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


#include "fem/dofhandler.hh"
#include "fem/finite_element.hh"
#include "mesh/mesh.h"
//#include "fem/simplex.hh"










template<unsigned int dim, unsigned int spacedim> inline
DOFHandler<dim,spacedim>::DOFHandler(Mesh & _mesh)
: mesh(&_mesh),
  n_dofs(0),
  global_dof_offset(0),
  finite_element(0)
{
}



template<unsigned int dim, unsigned int spacedim> inline
void DOFHandler<dim,spacedim>::distribute_dofs(FiniteElement<dim,spacedim> & fe, const unsigned int offset)
{
    unsigned int next_free_dof = offset;
    unsigned int n_obj_dofs[dim+1];

    // TODO: Maybe check if dofs are not yet distributed?

    finite_element = &fe;
    global_dof_offset = offset;

    for (int dm=0; dm <= dim; dm++)
    {
        n_obj_dofs[dm] = 0;
        for (unsigned int m=0; m<dof_multiplicities.size(); m++)
            n_obj_dofs[dm] += fe.n_object_dofs(dm, dof_multiplicities[m])*dof_multiplicities[m];
    }

    FOR_ELEMENTS(mesh,cell)
    {
        // skip cells of different dimension
        if (cell->dim != dim) continue;

        // distribute dofs
        // TODO: For the moment we distribute only dofs associated to the cell
        //       In the future we want to distribute dofs on vertices, lines,
        //       and triangles as well.
        object_dofs[dim][cell] = new int[n_obj_dofs[dim]];
        for (int i=0; i<n_obj_dofs[dim]; i++)
           object_dofs[dim][cell][i] = next_free_dof++;
    }


//    FOR_ELEMENTS(mesh,cell)
//    {
//        // skip cells of different dimension
//        if (cell->dim != dim) continue;
//
//        // distribute dofs
//        for (int dm=0; dm<=dim; dm++)
//        {
//            for (int i=0; i<cell->n_sides_by_dim(dm); i++)
//            {
//                void *side = cell->side_by_dim(dm, i);
//                // check if node has already assigned dofs, otherwise
//                // distribute
//                if (object_dofs[dm].find(side) == object_dofs[dm].end())
//                {
//                    object_dofs[dm][side] = new int[n_obj_dofs[dm]];
//                    for (int i=0; i<n_obj_dofs[dm]; i++)
//                        object_dofs[dm][side][i] = next_free_dof++;
//                }
//            }
//        }
//    }

    n_dofs = next_free_dof - offset;
}



template<unsigned int dim, unsigned int spacedim> inline
const unsigned int DOFHandler<dim,spacedim>::n_local_dofs()
{
    return finite_element->n_dofs();
}



template<unsigned int dim, unsigned int spacedim> inline
const unsigned int DOFHandler<dim,spacedim>::n_global_dofs()
{
    return n_dofs;
}


template<unsigned int dim, unsigned int spacedim>
void DOFHandler<dim,spacedim>::get_dof_indices(const CellIterator &cell, unsigned int indices[])
{
    void *side;
    unsigned int offset, pid;

    for (int k=0; k<finite_element->n_object_dofs(dim,DOF_SINGLE); k++)
            indices[k] = object_dofs[dim][cell][k];

//    indices.clear();
//
//    get_object_dof_indices<0>(cell, indices);
//    get_object_dof_indices<1>(cell, indices);
//    get_object_dof_indices<2>(cell, indices);
//    get_object_dof_indices<3>(cell, indices);
}

//template<unsigned int dim> template<unsigned int obj_dim> inline void DOFHandler<dim>::get_object_dof_indices(const CellIterator &cell, unsigned int indices[])
//{
    // TODO: implement for lower dimensional objects

//    void *side;
//    unsigned int offset, pid;
//
//    // loop over cell points/lines/triangles/tetrahedra
//    for (int i=0; i<n_simplex_objects<dim>(obj_dim); i++)
//    {
//        side   = cell->side_by_dim(obj_dim,i);
//        pid    = permutation_id<dim,obj_dim>(cell,i);
//        offset = 0;
//        // loop over dof multiplicities (single dofs, pairs, triples, sextuples)
//        for (vector<unsigned int>::iterator m=dof_multiplicities.begin(); m!=dof_multiplicities.end(); m++)
//        {
//            // loop over particular single dofs/dof pairs/triples/sextuples
//            for (int j=0; j<finite_element.n_object_dofs(obj_dim,*m); j++)
//            {
//                // loop over particular dofs (the single dof/2 dofs in the pair etc.)
//                for (int k=0; k<*m; k++)
//                    indices.push_back(object_dofs[obj_dim][side][offset+Simplex<obj_dim>::pair_permutations[pid][k]]);
//
//                offset += *m;
//            }
//        }
//    }
//}

template<unsigned int dim, unsigned int spacedim> inline
void DOFHandler<dim,spacedim>::get_dof_values(const CellIterator &cell, const Vec &values, double local_values[])
{
    unsigned int indices[finite_element->n_dofs()];

    get_dof_indices(cell, indices);
    VecGetValues(values, finite_element->n_dofs(), (PetscInt *)indices, local_values);
}


template<unsigned int dim, unsigned int spacedim> inline
const unsigned int DOFHandler<dim,spacedim>::global_dof_id(const CellIterator &cell, const unsigned int local_dof_id)
{
    ASSERT(local_dof_id<n_dofs, "Number of local dof is out of range.");
    unsigned int count_dofs = 0;
    for (int dm=0; dm<=dim; dm++)
    {
        for (int i=0; i<cell->n_sides_by_dim(dm); i++)
        {
            int side_dof = count_dofs + cell->n_sides_by_dim(dm) - local_dof_id;
            if (side_dof > 0)
            {
                return object_dofs[dm][cell->side_by_dim(dm,i)][side_dof];
            }
            else
            {
                count_dofs += cell->n_sides_by_dim(dm);
            }
        }
    }
}



template<unsigned int dim, unsigned int spacedim> inline
typename DOFHandler<dim,spacedim>::CellIterator DOFHandler<dim,spacedim>::begin_cell() const
{
    return mesh->element.begin();
}




template<unsigned int dim, unsigned int spacedim> inline
typename DOFHandler<dim,spacedim>::CellIterator DOFHandler<dim,spacedim>::end_cell() const
{
    return mesh->element.end();
}

template<unsigned int dim, unsigned int spacedim> inline
DOFHandler<dim,spacedim>::~DOFHandler()
{
    for (int dm=0; dm<=dim; dm++) object_dofs[dm].clear();
}


template class DOFHandler<1,3>;
template class DOFHandler<2,3>;
template class DOFHandler<3,3>;

