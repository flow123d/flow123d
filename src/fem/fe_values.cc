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


#include "fem/mapping.hh"
#include "quadrature/quadrature.hh"
#include "fem/finite_element.hh"
#include "fem/fe_values.hh"



using namespace arma;
using namespace std;










template<unsigned int dim, unsigned int spacedim> inline
void FEValuesData<dim,spacedim>::allocate(unsigned int size, UpdateFlags flags, bool is_scalar)
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
    {
        if (is_scalar)
        {
            shape_values.resize(size);
        }
        else
        {
            shape_vectors.resize(size);
        }
    }

    if (update_flags & update_gradients)
    {
        if (is_scalar)
        {
            shape_gradients.resize(size);
        }
        else
        {
            shape_grad_vectors.resize(size);
        }
    }

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
        FiniteElement<dim,spacedim> & _fe,
        UpdateFlags _flags)

{
    mapping = &_mapping;
    quadrature = &_quadrature;
    fe = &_fe;

    // add flags required by the finite element or mapping
    data.allocate(quadrature->size(), update_each(_flags), fe->is_scalar());
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
const vec::fixed<spacedim> FEValuesBase<dim,spacedim>::shape_vector(const unsigned int function_no, const unsigned int point_no)
{
    return data.shape_vectors[point_no][function_no];
}

template<unsigned int dim, unsigned int spacedim> inline
const mat::fixed<spacedim,spacedim> FEValuesBase<dim,spacedim>::shape_grad_vector(const unsigned int function_no, const unsigned int point_no)
{
    return data.shape_grad_vectors[point_no][function_no];
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

template<unsigned int dim, unsigned int spacedim> inline
const unsigned int FEValuesBase<dim,spacedim>::n_points()
{
    return quadrature->size();
}

template<unsigned int dim, unsigned int spacedim> inline
const unsigned int FEValuesBase<dim,spacedim>::n_dofs()
{
    return fe->n_dofs();
}












template<unsigned int dim, unsigned int spacedim>
FEValues<dim,spacedim>::FEValues(Mapping<dim,spacedim> &_mapping,
         Quadrature<dim> &_quadrature,
         FiniteElement<dim,spacedim> &_fe,
         UpdateFlags _flags)
{
    this->allocate(_mapping, _quadrature, _fe, _flags);

    // precompute the maping data and finite element data
    this->mapping_data = this->mapping->initialize(*this->quadrature, this->data.update_flags);
    this->fe_data = this->fe->initialize(*this->quadrature, this->data.update_flags);
}

template<unsigned int dim,unsigned int spacedim> inline
void FEValues<dim,spacedim>::reinit(typename DOFHandler<dim,spacedim>::CellIterator & cell)
{
    this->data.present_cell = &cell;

    // calculate Jacobian of mapping, JxW, inverse Jacobian, normal vector(s)
    this->mapping->fill_fe_values(cell,
                            *this->quadrature,
                            *this->mapping_data,
                            this->data);

    this->fe->fill_fe_values(*this->quadrature,
                             *this->fe_data,
                             this->data);
}









template<unsigned int dim,unsigned int spacedim> inline
FESideValues<dim,spacedim>::FESideValues(Mapping<dim,spacedim> & _mapping,
                                 Quadrature<dim-1> & _sub_quadrature,
                                 FiniteElement<dim,spacedim> & _fe,
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
void FESideValues<dim,spacedim>::reinit(typename DOFHandler<dim,spacedim>::CellIterator & cell,
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
                             this->data);
}





template class FEValuesBase<0,3>;
template class FEValuesBase<1,3>;
template class FEValuesBase<2,3>;
template class FEValuesBase<3,3>;

template class FEValues<0,3>;
template class FEValues<1,3>;
template class FEValues<2,3>;
template class FEValues<3,3>;

template class FESideValues<1,3>;
template class FESideValues<2,3>;
template class FESideValues<3,3>;

