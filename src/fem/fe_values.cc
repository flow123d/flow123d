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
 * @file    fe_values.cc
 * @brief   Class FEValues calculates finite element data on the actual
 *          cells such as shape function values, gradients, Jacobian of
 *          the mapping from the reference cell etc.
 * @author  Jan Stebel
 */

#include "fem/mapping.hh"
#include "quadrature/quadrature.hh"
#include "fem/finite_element.hh"
#include "fem/fe_values.hh"



using namespace arma;
using namespace std;










template<unsigned int dim, unsigned int spacedim>
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





template<unsigned int dim,unsigned int spacedim>
FEValuesBase<dim,spacedim>::FEValuesBase()
: mapping(NULL), quadrature(NULL), fe(NULL), mapping_data(NULL), fe_data(NULL)
{
}



template<unsigned int dim,unsigned int spacedim>
FEValuesBase<dim,spacedim>::~FEValuesBase() {
    if (mapping_data) delete mapping_data;
    if (fe_data) delete fe_data;
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



template<unsigned int dim, unsigned int spacedim>
UpdateFlags FEValuesBase<dim,spacedim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags | fe->update_each(flags);
    f |= mapping->update_each(f);
    return f;
}






template<unsigned int dim, unsigned int spacedim>
FEValues<dim,spacedim>::FEValues(Mapping<dim,spacedim> &_mapping,
         Quadrature<dim> &_quadrature,
         FiniteElement<dim,spacedim> &_fe,
         UpdateFlags _flags)
:FEValuesBase<dim, spacedim>()
{
    this->allocate(_mapping, _quadrature, _fe, _flags);

    // precompute the maping data and finite element data
    this->mapping_data = this->mapping->initialize(*this->quadrature, this->data.update_flags);
    this->fe_data = this->fe->initialize(*this->quadrature, this->data.update_flags);
}



template<unsigned int dim,unsigned int spacedim>
void FEValues<dim,spacedim>::reinit(ElementFullIter & cell)
{
	OLD_ASSERT_EQUAL( dim, cell->dim() );
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









template<unsigned int dim,unsigned int spacedim>
FESideValues<dim,spacedim>::FESideValues(Mapping<dim,spacedim> & _mapping,
                                 Quadrature<dim-1> & _sub_quadrature,
                                 FiniteElement<dim,spacedim> & _fe,
                                 const UpdateFlags _flags)
:FEValuesBase<dim,spacedim>()
{
    sub_quadrature = &_sub_quadrature;
    Quadrature<dim> *q = new Quadrature<dim>(_sub_quadrature.size());
    this->allocate(_mapping, *q, _fe, _flags);

    for (unsigned int sid = 0; sid < RefElement<dim>::n_sides; sid++)
    {
    	for (unsigned int pid = 0; pid < RefElement<dim>::n_side_permutations; pid++)
    	{
    		// transform the side quadrature points to the cell quadrature points
    		this->mapping->transform_subquadrature(sid, pid, *sub_quadrature, side_quadrature[sid][pid]);
    		side_mapping_data[sid][pid] = this->mapping->initialize(side_quadrature[sid][pid], this->data.update_flags);
    		side_fe_data[sid][pid] = this->fe->initialize(side_quadrature[sid][pid], this->data.update_flags);
    	}
    }
}



template<unsigned int dim,unsigned int spacedim>
FESideValues<dim,spacedim>::~FESideValues()
{
	for (unsigned int sid=0; sid<RefElement<dim>::n_sides; sid++)
	{
		for (unsigned int pid=0; pid<RefElement<dim>::n_side_permutations; pid++)
		{
			delete side_mapping_data[sid][pid];
			delete side_fe_data[sid][pid];
		}
	}

    // Since quadrature is an auxiliary internal variable allocated
    // by the constructor, it must be destroyed here.
    delete this->quadrature;
}


template<unsigned int dim,unsigned int spacedim>
void FESideValues<dim,spacedim>::reinit(ElementFullIter & cell,
		unsigned int sid)
{
    this->data.present_cell = &cell;

    unsigned int pid = cell->permutation_idx_[sid];

    // calculate Jacobian of mapping, JxW, inverse Jacobian, normal vector(s)
    this->mapping->fill_fe_side_values(cell,
                                 sid,
                                 side_quadrature[sid][pid],
                                 *side_mapping_data[sid][pid],
                                 this->data);

    // calculation of finite element data
    this->fe->fill_fe_values(side_quadrature[sid][pid],
                             *side_fe_data[sid][pid],
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

