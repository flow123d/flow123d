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
 * @brief Abstract class for description of finite elements.
 * @author Jan Stebel
 */

#ifndef FINITE_ELEMENT_HH_
#define FINITE_ELEMENT_HH_

#include <armadillo>

/*
 * Describes the type of the geometric entity to which the degrees of freedom are associated.
 * For example, in the case of the usual P1 Lagrangean finite element all degrees of freedom
 * are the values at the vertices of the line/triangle/tetrahedron and their corresponding FE_dof_object
 * is FE_OBJECT_POINT.
 */
enum FE_dof_object
{
    /*
     * Dof is associated to a vertex.
     */
    FE_OBJECT_POINT = 0,

    /*
     * Dof is associated to an internal point of a 1D edge.
     */
    FE_OBJECT_LINE = 1,

    /*
     * Dof is associated to an internal point of a 2D face.
     */
    FE_OBJECT_TRIANGLE = 2,

    /*
     * Dof is associated to an internal point of a 3D cell.
     */
    FE_OBJECT_TETRAHEDRON = 3
};

///*
// * Describes the type of continuity at interfaces (vertices, edges, faces) of the functions
// * constructed from the given finite element space at each cell.
// */
//typedef enum Conformity
//{
//    /*
//     * Indicates incompatible continuities of the space.
//     */
//    unknown,
//
//    /*
//     * Discontinuous finite element (no continuity at interfaces).
//     */
//    L2,
//
//    /*
//     * Continuous finite element (continuity at interfaces).
//     */
//    H1
//};

template <unsigned int dim>
class FiniteElement
{
public:

    FiniteElement(unsigned int _n_dofs, const bool *_dof_continuity, const FE_dof_object *_dof_objs, const unsigned int *_dof_obj_ids);

    /*
     * Returns the number of degrees of freedom needed by the finite element.
     */
    const unsigned int n_dofs();

    /*
     * Returns information on the conformity of the global function space.
     */
    const bool dof_is_continuous(int dof);

    /*
     * Returns the type of the geometric object to which the @p dof is associated.
     */
    const FE_dof_object dof_object(int dof);

    /*
     * Returns the number of the geometric entity (node, edge, face, cell) of the @p dof.
     */
    const unsigned int dof_object_id(int dof);

    /*
     * Calculates the value of the @p i-th shape function at the point @p p on the reference element
     */
    virtual double shape_value(const unsigned int i, const arma::vec::fixed<dim> &p) const = 0;

    /*
     * Calculates the gradient of the @p i-th shape function at the point @p p on the reference element.
     * The gradient components are relative to the reference cell coordinate system.
     */
    virtual arma::vec::fixed<dim> shape_grad(const unsigned int i, const arma::vec::fixed<dim> &p) const = 0;

protected:

    const unsigned int number_of_dofs;

    const bool *dof_continuity;

    const FE_dof_object *dof_objects;

    const unsigned int *dof_object_ids;

};






template<unsigned int dim>
FiniteElement<dim>::FiniteElement(unsigned int _n_dofs, const bool *_dof_continuity, const FE_dof_object *_dof_objs, const unsigned int *_dof_obj_ids) :
        number_of_dofs(_n_dofs),
        dof_continuity(_dof_continuity),
        dof_objects(_dof_objs),
        dof_object_ids(_dof_obj_ids)
{}


template<unsigned int dim> inline const unsigned int FiniteElement<dim>::n_dofs()
{
    return number_of_dofs;
}


template<unsigned int dim> inline const bool FiniteElement<dim>::dof_is_continuous(int dof)
{
    ASSERT(dof>=0 && dof<dim, "Dof number is out of range.");
    return dof_continuity[dof];
}


template<unsigned int dim> inline const FE_dof_object FiniteElement<dim>::dof_object(int dof)
{
    ASSERT(dof>=0 && dof<dim, "Dof number is out of range.");
    return dof_objects[dof];
}


template<unsigned int dim> inline const unsigned int FiniteElement<dim>::dof_object_id(int dof)
{
    ASSERT(dof>=0 && dof<dim, "Dof number is out of range.");
    return dof_object_ids[dof];
}


#endif /* FINITE_ELEMENT_HH_ */
