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

#include "fem/discrete_space.hh"
#include "fem/finite_element.hh"
#include "mesh/mesh.h"
#include "mesh/duplicate_nodes.h"



unsigned int EqualOrderDiscreteSpace::n_elem_dofs(const ElementAccessor<3> &cell) const
{
    unsigned int n_dofs = 0;
    switch (cell->dim())
    {
        case 1:
            for (unsigned int d=0; d<fe1_->n_dofs(); d++)
                if (fe1_->dof(d).dim == 1)
                    n_dofs++;
            break;
        case 2:
            for (unsigned int d=0; d<fe2_->n_dofs(); d++)
                if (fe2_->dof(d).dim == 2)
                    n_dofs++;
            break;
        case 3:
            for (unsigned int d=0; d<fe3_->n_dofs(); d++)
                if (fe3_->dof(d).dim == 3)
                    n_dofs++;
            break;
    }
    return n_dofs;
}


unsigned int EqualOrderDiscreteSpace::n_node_dofs(unsigned int nid) const
{
    unsigned int n_dofs = 0;
    unsigned int dim = mesh_->tree->node_dim()[nid];
    switch (dim)
    {
        case 1:
            for (unsigned int d=0; d<fe1_->n_dofs(); d++)
                if (fe1_->dof(d).dim == 0 && fe1_->dof(d).n_face_idx == 0)
                    n_dofs++;
            break;
        case 2:
            for (unsigned int d=0; d<fe2_->n_dofs(); d++)
                if (fe2_->dof(d).dim == 0 && fe2_->dof(d).n_face_idx == 0)
                    n_dofs++;
            break;
        case 3:
            for (unsigned int d=0; d<fe3_->n_dofs(); d++)
                if (fe3_->dof(d).dim == 0 && fe3_->dof(d).n_face_idx == 0)
                    n_dofs++;
            break;
    }
    return n_dofs;
}



template<> FiniteElement<1> *DiscreteSpace::fe(const ElementAccessor<3> &cell) const { return fe1d(cell); }
template<> FiniteElement<2> *DiscreteSpace::fe(const ElementAccessor<3> &cell) const { return fe2d(cell); }
template<> FiniteElement<3> *DiscreteSpace::fe(const ElementAccessor<3> &cell) const { return fe3d(cell); }





