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
 * @file    mapping_p1.cc
 * @brief   Class MappingP1 implements the affine transformation of
 *          the unit cell onto the actual cell.
 * @author  Jan Stebel
 */

#include "fem/mapping_p1.hh"
#include "quadrature/quadrature.hh"



using namespace std;



template<unsigned int dim, unsigned int spacedim>
UpdateFlags MappingP1<dim,spacedim>::update_each(UpdateFlags flags)
{
    UpdateFlags f = flags;
    
    if (flags & update_normal_vectors)
        f |= update_inverse_jacobians;

    if ((flags & update_volume_elements) |
        (flags & update_JxW_values) |
        (flags & update_inverse_jacobians))
        f |= update_jacobians;

    return f;
}



template<unsigned int dim, unsigned int spacedim>
auto MappingP1<dim,spacedim>::element_map(ElementAccessor<3> elm) -> ElementMap
{
    ElementMap coords;
    for (unsigned int i=0; i<dim+1; i++)
        coords.col(i) = *elm.node(i);
    return coords;
}


template<unsigned int dim, unsigned int spacedim>
auto MappingP1<dim,spacedim>::jacobian(const ElementMap &coords) -> arma::mat::fixed<spacedim,dim>
{
    arma::mat::fixed<spacedim,dim> jac;
    for (unsigned int i=0; i<spacedim; i++)
        for (unsigned int j=0; j<dim; j++)
            jac(i,j) = coords(i,j+1) - coords(i,0);
    return jac;
}


template<unsigned int dim, unsigned int spacedim>
auto MappingP1<dim,spacedim>::project_real_to_unit(const RealPoint &point, const ElementMap &map) -> BaryPoint
{
    arma::mat::fixed<3, dim> A = map.cols(1,dim);
    for(unsigned int i=0; i < dim; i++ ) {
        A.col(i) -= map.col(0);
    }
    
    arma::mat::fixed<dim, dim> AtA = A.t()*A;
    arma::vec::fixed<dim> Atb = A.t()*(point - map.col(0));
    arma::vec::fixed<dim+1> bary_coord;
    bary_coord.rows(1, dim) = arma::solve(AtA, Atb);
    bary_coord( 0 ) = 1.0 - arma::sum( bary_coord.rows(1,dim) );
    return bary_coord;
}

template<unsigned int dim, unsigned int spacedim>
auto MappingP1<dim,spacedim>::project_unit_to_real(const BaryPoint &point, const ElementMap &map) -> RealPoint
{
    return map * point;
}

template<unsigned int dim, unsigned int spacedim>
auto MappingP1<dim,spacedim>::clip_to_element(BaryPoint &barycentric) -> BaryPoint{
    return RefElement<dim>::clip(barycentric);
}

template <unsigned int dim, unsigned int spacedim>
bool MappingP1<dim,spacedim>::contains_point(arma::vec point, ElementAccessor<3> elm)
{
	arma::vec projection = project_real_to_unit(point, element_map(elm));
	return (projection.min() >= -BoundingBox::epsilon);
}



template class MappingP1<0,3>;  // Only for compilation of DG transport assemble methods, do not use this instance!
template class MappingP1<1,3>;
template class MappingP1<2,3>;
template class MappingP1<3,3>;




