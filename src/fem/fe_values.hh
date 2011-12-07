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
 * @brief Class FEValues calculates finite element data on the actual
 *        cells such as shape function values, gradients, Jacobian of
 *        the mapping from the reference cell etc.
 * @author Jan Stebel
 */

#ifndef FE_VALUES_HH_
#define FE_VALUES_HH_

#include <armadillo>
#include <vector>
#include "fem/mapping.hh"
#include "fem/finite_element.hh"

template<unsigned int dim> class DOFHandler;
template<unsigned int dim> class Quadrature;
//class UpdateFlags;

using namespace arma;
using namespace std;

/**
 * FEValues takes care of the calculation of finite element data on
 * the actual cell such as values of shape functions at quadrature
 * points, gradients of shape functions, Jacobians of the mapping
 * from the reference cell etc.
 * @param dim      dimension of the reference cell
 * @param spacedim dimension of the Euclidean space where the actual
 *                 cell lives
 *
 * TODO: Implement UpdateFlags to indicate which data should be
 * recomputed at each cell.
 */
template<unsigned int dim, unsigned int spacedim>
class FEValues
{
public:

    /**
     * Constructor. Initializes structures and calculates
     * cell-independent data.
     */
    FEValues(Mapping<dim,spacedim> &_mapping,
             const Quadrature<dim> &_quadrature,
             FiniteElement<dim> &_fe//,
//             const UpdateFlags _update_flags
             );

    /**
     * Update cell-dependent data (gradients, Jacobians etc.)
     */
    void reinit(const typename DOFHandler<dim>::CellIterator &cell);

    /**
     * Return the value of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     */
    const double shape_value(const unsigned int function_no, const unsigned int point_no);

    /**
     * Return the gradient of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     */
    const vec::fixed<spacedim> &shape_grad(const unsigned int function_no, const unsigned int point_no);

//    /*
//     * Compute the values of a FE function given by the dof values @p global_dofs at the quadrature points.
//     */
//    void get_function_values(const Vec &global_dofs, vector<double> &values);

private:

    /**
     * The mapping from the reference cell to the actual cell.
     */
    Mapping<dim,spacedim> *mapping;

    /**
     * The quadrature rule used to calculate integrals.
     */
    const Quadrature<dim> *quadrature;

    /**
     * The used finite element.
     */
    FiniteElement<dim> *fe;

//    const UpdateFlags update_flags;

    /**
     * Values at quadrature points of the product of the Jacobian
     * determinant of the mapping and the weight at the particular
     * quadrature point.
     */
    vector<double> JxW_values;

    /**
     * Jacobians of the mapping at the quadrature points.
     */
    vector<mat::fixed<dim,spacedim> > jacobians;

    /**
     * Inverse Jacobians at the quadrature points.
     */
    vector<mat > inverse_jacobians;

    /**
     * Normal vectors at the quadrature points.
     */
    vector<vec::fixed<spacedim> > normal_vectors;

    /**
     * Shape functions evaluated at the quadrature points.
     */
    vector<vec> shape_values;

    /**
     * Gradients of shape functions evaluated at the quadrature points.
     * Each row of the matrix contains the gradient of one shape function.
     */
    vector<mat> shape_gradients;
};



template<unsigned int dim,unsigned int spacedim> inline
FEValues<dim,spacedim>::FEValues(Mapping<dim,spacedim> & _mapping,
                                 const Quadrature<dim> & _quadrature,
                                 FiniteElement<dim> & _fe//,
//                                 const UpdateFlags _update_flags
                                 )
        : mapping(&_mapping),
          quadrature(&_quadrature),
          fe(&_fe)//,
//          update_flags(_update_flags)
{
    jacobians.resize(quadrature->size());
    JxW_values.resize(quadrature->size());
    inverse_jacobians.resize(quadrature->size());
    shape_values.resize(quadrature->size());
    shape_gradients.resize(quadrature->size());
}



template<unsigned int dim,unsigned int spacedim> inline
void FEValues<dim,spacedim>::reinit(const typename DOFHandler<dim>::CellIterator & cell)
{
    // calculate Jacobian of mapping, JxW, inverse Jacobian, normal vector(s)
    mapping->fill_fe_values(cell, *quadrature, jacobians, JxW_values, inverse_jacobians, normal_vectors);

    fe->fill_fe_values(*quadrature, inverse_jacobians, shape_values, shape_gradients);

}


template<unsigned int dim, unsigned int spacedim> inline
const double FEValues<dim,spacedim>::shape_value(const unsigned int function_no, const unsigned int point_no)
{
    return shape_values[point_no][function_no];
}


template<unsigned int dim, unsigned int spacedim> inline
const vec::fixed<spacedim> &FEValues<dim,spacedim>::shape_grad(const unsigned int function_no, const unsigned int point_no)
{
    return shape_gradients[point_no].row(function_no);
}





#endif /* FE_VALUES_HH_ */
