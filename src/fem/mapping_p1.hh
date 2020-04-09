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
 * @file    mapping_p1.hh
 * @brief   Class MappingP1 implements the affine transformation of
 *          the unit cell onto the actual cell.
 * @author  Jan Stebel
 */

#ifndef MAPPING_P1_HH_
#define MAPPING_P1_HH_

#include <string.h>                            // for memcpy
#include <algorithm>                           // for min, swap
#include <cmath>                               // for abs, pow, fabs, log10
#include <ostream>                             // for operator<<
#include <string>                              // for operator<<
#include <armadillo>
#include "fem/dofhandler.hh"                   // for DOFHandlerBase, DOFHan...
#include "fem/mapping.hh"                      // for MappingInternalData (p...
#include "fem/update_flags.hh"                 // for operator&, operator|
#include "mesh/elements.h"                     // for Element::side
#include "mesh/accessors.hh"                   // for Side::node
template <unsigned int dim, unsigned int spacedim> class FEValuesData;
class Quadrature;


/**
 * Auxiliary templates to resolve nonzero matrix sizes for small dimensions. 
 * ( some compilers do not accept ternary operator in template parameters.)
 */ 
template<unsigned int dim>
class MatrixSizes {
public:  
  static const unsigned int dim_minus_one = dim-1;
};  

template<>
class MatrixSizes<0> {
public:  
  static const unsigned int dim_minus_one = 0;
};  


/**
 * @brief Affine mapping between reference and actual cell.
 *
 * Class MappingP1 implements the affine transformation of
 * the reference cell onto the actual cell.
 *
 * @param dim Dimension of the cells.
 * @param spacedim Dimension of the Euclidean space.
 */
template<unsigned int dim, unsigned int spacedim>
class MappingP1 : public Mapping<dim,spacedim>
{
public:

    typedef arma::vec::fixed<dim+1> BaryPoint;
    typedef arma::vec::fixed<spacedim> RealPoint;
    typedef arma::mat::fixed<spacedim, dim+1> ElementMap;
    
    /**
     * @brief Constructor.
     */
    MappingP1();

    /**
     * @brief Initializes the structures and computes static data.
     *
     * @param q Quadrature rule.
     * @param flags Update flags.
     * @return The computed mapping data.
     */
    MappingInternalData *initialize(const Quadrature &q, UpdateFlags flags);

    /**
     * @brief Determines which additional quantities have to be computed.
     *
     * @param flags Update flags for required quantities.
     * @return All necessary flags.
     */
    UpdateFlags update_each(UpdateFlags flags);

    /**
     * @brief Calculates the mapping data on the actual cell.
     *
     * @param cell The actual cell.
     * @param q Quadrature rule.
     * @param data Precomputed mapping data.
     * @param fv_data Data to be computed.
     */
    void fill_fe_values(const ElementAccessor<3> &cell,
                            const Quadrature &q,
                            MappingInternalData &data,
                            FEValuesData<dim,spacedim> &fv_data);

    /**
     * @brief Calculates the mapping data on a side of a cell.
     *
     * @param cell The actual cell.
     * @param sid  Number of the side.
     * @param q The quadrature rule with points on the side.
     * @param data Precomputed mapping data.
     * @param fv_data Data to be computed.
     */
    void fill_fe_side_values(const ElementAccessor<3> &cell,
                            unsigned int sid,
                            const Quadrature &q,
                            MappingInternalData &data,
                            FEValuesData<dim,spacedim> &fv_data);

   
    /**
     * Map from reference element (barycentric coords) to global coord system.
     * Matrix(3, dim+1) M: x_real = M * x_bary;
     * M columns are real coordinates of nodes.
     */
    ElementMap element_map(ElementAccessor<3> elm) const;

    /**
     * Project given point in real coordinates to reference element (barycentic coordinates).
     * Result vector have dimension dim()+1.
     * Use RefElement<dim>::bary_to_local() to get local coordinates.
     */
    BaryPoint project_real_to_unit(const RealPoint &point, const ElementMap &map) const;
    
    /**
     * Project given point from reference element (barycentic coordinates) to real coordinates.
     * Use RefElement<dim>::local_to_bary() to get barycentric coordinates in input.
     */
    RealPoint project_unit_to_real(const BaryPoint &point, const ElementMap &map) const;

    /**
     * Clip a point given by barycentric cocordinates to the element.
     * If the point is out of the element the closest point
     * projection to the element surface is used.
     */
    BaryPoint clip_to_element(BaryPoint &barycentric);

    /// Test if element contains given point.
    bool contains_point(arma::vec point, ElementAccessor<3> elm);

private:

    /**
     * @brief Auxiliary matrix of gradients of shape functions (used for
     * computation of the Jacobian).
     */
    arma::mat::fixed<dim+1,dim> grad;

};








#endif /* MAPPING_P1_HH_ */
