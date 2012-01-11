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
#include "fem/update_flags.hh"

template<unsigned int dim> class DOFHandler;
template<unsigned int dim> class Quadrature;
template<unsigned int dim, unsigned int spacedim> class FEValuesBase;

using namespace arma;
using namespace std;


/**
 * Class FEValuesData holds the arrays of data computed by
 * Mapping and FiniteElement.
 */
template<unsigned int dim, unsigned int spacedim>
class FEValuesData
{
public:

    /**
     * Resize the data arrays.
     */
    void allocate(unsigned int size, UpdateFlags flags);



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
     * Determinants of Jacobians at quadrature points.
     */
    vector<double> determinants;

    /**
     * Inverse Jacobians at the quadrature points.
     */
    vector<mat> inverse_jacobians;

    /**
     * Coordinates of quadrature points in the actual cell coordinate system.
     */
    vector<vec::fixed<spacedim> > points;

    /**
     * Shape functions evaluated at the quadrature points.
     */
    vector<vec> shape_values;

    /**
     * Gradients of shape functions evaluated at the quadrature points.
     * Each row of the matrix contains the gradient of one shape function.
     */
    vector<mat> shape_gradients;

    /**
     * Normal vectors to the element at the quadrature points lying
     * on a side.
     */
    vector<vec::fixed<spacedim> > normal_vectors;

    /**
     * Flags that indicate which finite element quantities are to be computed.
     */
    UpdateFlags update_flags;

    /**
    * Iterator to the last reinit-ed cell.
    */
    typename DOFHandler<dim>::CellIterator *present_cell;

};

/**
* Base class for FEValues and FESideValues
*/
template<unsigned int dim, unsigned int spacedim>
class FEValuesBase
{
public:

    /**
     * Default constructor
     */
    FEValuesBase();

    /**
     * Allocates space for computed data.
     */
    void allocate(Mapping<dim,spacedim> &_mapping,
            Quadrature<dim> &_quadrature,
            FiniteElement<dim> &_fe,
            UpdateFlags flags);
    
    /**
    * Determine quantities to be recomputed on each cell.
    */
    UpdateFlags update_each(UpdateFlags flags);

    /**
     * Return the value of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     */
    const double shape_value(const unsigned int function_no, const unsigned int point_no);

    /**
     * Return the gradient of the @p function_no-th shape function at
     * the @p point_no-th quadrature point.
     */
    const vec::fixed<spacedim> shape_grad(const unsigned int function_no, const unsigned int point_no);

    /**
     * Return the relative volume change of the cell (Jacobian determinant).
     * If dim==spacedim then the sign may be negative, otherwise the
     * result is a positive number.
     */
    const double determinant(const unsigned int point_no);

    /**
     * Return the product of Jacobian determinant and the quadrature
     * weight at given quadrature point.
     */
    const double JxW(const unsigned int point_no);

    /**
     * Return coordinates of the quadrature point in the actual cell system.
     */
    const vec::fixed<spacedim> point(const unsigned int point_no);

    /**
     * Returns the normal vector to a side at given quadrature point.
     */
    const vec::fixed<spacedim> normal_vector(unsigned int point_no);


protected:

    /**
     * The mapping from the reference cell to the actual cell.
     */
    Mapping<dim,spacedim> *mapping;

    /**
     * The quadrature rule used to calculate integrals.
     */
    Quadrature<dim> *quadrature;

    /**
     * The used finite element.
     */
    FiniteElement<dim> *fe;
    
    /**
     * Precomputed maping data.
     */
    MappingInternalData *mapping_data;

    /**
     * Precomputed finite element data.
     */
    FEInternalData *fe_data;

    /**
     * Data computed by the mapping and finite element.
     */
    FEValuesData<dim,spacedim> data;
};




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
class FEValues : public FEValuesBase<dim,spacedim>
{
public:

    FEValues(Mapping<dim,spacedim> &_mapping,
             Quadrature<dim> &_quadrature,
             FiniteElement<dim> &_fe,
             UpdateFlags _flags);

    /**
     * Update cell-dependent data (gradients, Jacobians etc.)
     */
    void reinit(typename DOFHandler<dim>::CellIterator &cell);


};




/**
 * FESideValues takes care of the calculation of finite element data on
 * a side of the actual cell such as values of shape functions at quadrature
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
class FESideValues : public FEValuesBase<dim,spacedim>
{

public:

    /**
     * Constructor. Initializes structures and calculates
     * cell-independent data.
     */
    FESideValues(Mapping<dim,spacedim> &_mapping,
             Quadrature<dim-1> &_sub_quadrature,
             FiniteElement<dim> &_fe,
             UpdateFlags flags);

    ~FESideValues();

    /**
     * Update cell-dependent data (gradients, Jacobians etc.)
     */
    void reinit(typename DOFHandler<dim>::CellIterator &cell,
                Side *side);


private:

    /**
     * Quadrature for the integration on the element sides.
     */
    const Quadrature<dim-1> *sub_quadrature;

};







template<unsigned int dim, unsigned int spacedim> inline
void FEValuesData<dim,spacedim>::allocate(unsigned int size, UpdateFlags flags)
{
    update_flags = flags;

    // resize the arrays of computed quantities
    if (update_flags & update_jacobians)
        jacobians.resize(size);

    if (update_flags & update_volume_elements)
        determinants.resize(size);

    if ((update_flags & update_JxW_values) |
        (update_flags & update_side_JxW_values))
        JxW_values.resize(size);

    if (update_flags & update_inverse_jacobians)
        inverse_jacobians.resize(size);

    if (update_flags & update_values)
        shape_values.resize(size);

    if (update_flags & update_gradients)
        shape_gradients.resize(size);

    if (update_flags & update_quadrature_points)
        points.resize(size);

    if (update_flags & update_normal_vectors)
        normal_vectors.resize(size);
}





template<unsigned int dim,unsigned int spacedim> inline
FEValuesBase<dim,spacedim>::FEValuesBase()
{
}

template<unsigned int dim, unsigned int spacedim>
void FEValuesBase<dim,spacedim>::allocate(Mapping<dim,spacedim> & _mapping,
        Quadrature<dim> & _quadrature,
        FiniteElement<dim> & _fe,
        UpdateFlags _flags)

{
    mapping = &_mapping;
    quadrature = &_quadrature;
    fe = &_fe;

    // add flags required by the finite element or mapping
    data.allocate(quadrature->size(), update_each(_flags));
}

template<unsigned int dim, unsigned int spacedim> inline
UpdateFlags FEValuesBase<dim,spacedim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags | fe->update_each(flags);
    f |= mapping->update_each(f);
    return f;
}

template<unsigned int dim, unsigned int spacedim> inline
const double FEValuesBase<dim,spacedim>::shape_value(const unsigned int function_no, const unsigned int point_no)
{
    return data.shape_values[point_no][function_no];
}


template<unsigned int dim, unsigned int spacedim> inline
const vec::fixed<spacedim> FEValuesBase<dim,spacedim>::shape_grad(const unsigned int function_no, const unsigned int point_no)
{
    return trans(data.shape_gradients[point_no].row(function_no));
}


template<unsigned int dim, unsigned int spacedim> inline
const double FEValuesBase<dim,spacedim>::determinant(const unsigned int point_no)
{
    return data.determinants[point_no];
}

template<unsigned int dim, unsigned int spacedim> inline
const double FEValuesBase<dim,spacedim>::JxW(const unsigned int point_no)
{
    return data.JxW_values[point_no];
}

template<unsigned int dim, unsigned int spacedim> inline
const vec::fixed<spacedim> FEValuesBase<dim,spacedim>::point(const unsigned int point_no)
{
    return data.points[point_no];
}

template<unsigned int dim,unsigned int spacedim> inline
const vec::fixed<spacedim> FEValuesBase<dim,spacedim>::normal_vector(unsigned int point_no)
{
    return data.normal_vectors[point_no];
}













template<unsigned int dim, unsigned int spacedim>
FEValues<dim,spacedim>::FEValues(Mapping<dim,spacedim> &_mapping,
         Quadrature<dim> &_quadrature,
         FiniteElement<dim> &_fe,
         UpdateFlags _flags)
{
    this->allocate(_mapping, _quadrature, _fe, _flags);

    // precompute the maping data and finite element data
    this->mapping_data = this->mapping->initialize(*this->quadrature, this->data.update_flags);
    this->fe_data = this->fe->initialize(*this->quadrature, this->data.update_flags);
}

template<unsigned int dim,unsigned int spacedim> inline
void FEValues<dim,spacedim>::reinit(typename DOFHandler<dim>::CellIterator & cell)
{
    this->data.present_cell = &cell;

    // calculate Jacobian of mapping, JxW, inverse Jacobian, normal vector(s)
    this->mapping->fill_fe_values(cell,
                            *this->quadrature,
                            *this->mapping_data,
                            this->data);

    this->fe->fill_fe_values(*this->quadrature,
                             *this->fe_data,
                             this->data.inverse_jacobians,
                             this->data.shape_values,
                             this->data.shape_gradients,
                             this->data.update_flags);
}









template<unsigned int dim,unsigned int spacedim> inline
FESideValues<dim,spacedim>::FESideValues(Mapping<dim,spacedim> & _mapping,
                                 Quadrature<dim-1> & _sub_quadrature,
                                 FiniteElement<dim> & _fe,
                                 const UpdateFlags _flags)
{
    sub_quadrature = &_sub_quadrature;
    Quadrature<dim> *q = new Quadrature<dim>(_sub_quadrature.size());
    this->allocate(_mapping, *q, _fe, _flags);
}

template<unsigned int dim,unsigned int spacedim>
FESideValues<dim,spacedim>::~FESideValues()
{
    // Since quadrature is an auxiliary internal variable allocated
    // by the constructor, it must be destroyed here.
    delete this->quadrature;
}



template<unsigned int dim,unsigned int spacedim> inline
void FESideValues<dim,spacedim>::reinit(typename DOFHandler<dim>::CellIterator & cell,
                                        Side *side)
{
    this->data.present_cell = &cell;

    // transform the side quadrature points to the cell quadrature points
    this->mapping->transform_subquadrature(cell, *this->quadrature, *side, *sub_quadrature);

    // compute the mapping and finite element data
    this->mapping_data = this->mapping->initialize(*this->quadrature, this->data.update_flags);
    this->fe_data = this->fe->initialize(*this->quadrature, this->data.update_flags);

    // calculate Jacobian of mapping, JxW, inverse Jacobian, normal vector(s)
    this->mapping->fill_fe_side_values(cell,
                                 *side,
                                 *this->quadrature,
                                 *this->mapping_data,
                                 this->data);

    // calculation of finite element data
    this->fe->fill_fe_values(*this->quadrature,
                             *this->fe_data,
                             this->data.inverse_jacobians,
                             this->data.shape_values,
                             this->data.shape_gradients,
                             this->data.update_flags);
}





#endif /* FE_VALUES_HH_ */
