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
 * @file    patch_fe_values.hh
 * @brief   Class FEValues calculates finite element data on the actual
 *          cells such as shape function values, gradients, Jacobian of
 *          the mapping from the reference cell etc.
 * @author  Jan Stebel, David Flanderka
 */

#ifndef PATCH_FE_VALUES_HH_
#define PATCH_FE_VALUES_HH_


#include <string.h>                           // for memcpy
#include <algorithm>                          // for swap
#include <new>                                // for operator new[]
#include <string>                             // for operator<<
#include <vector>                             // for vector
#include "fem/element_values.hh"              // for ElementValues
#include "fem/fe_values.hh"                   // for FEValuesBase
#include "fem/fe_values_views.hh"             // for FEValuesViews
#include "mesh/ref_element.hh"                // for RefElement
#include "mesh/accessors.hh"
#include "fem/update_flags.hh"                // for UpdateFlags
#include "quadrature/quadrature_lib.hh"
#include "fields/eval_subset.hh"

class PatchFEValues;



using Scalar = double;
using Vector = arma::vec3;
using Tensor = arma::mat33;

template <class ValueType>
class ElQ {
public:
    /// Constructor
    ElQ(PatchFeValues *fe_values, uint begin)
    : fe_values_(fe_values), begin_(begin) {}

    ValueType operator()(BulkPoint point);

    ValueType operator()(SidePoint point);

private:
    // attributes:
    PatchFeValues *fe_values_;
    uint begin_;    /// Index of the first component of the Quantity. Size is given by ValueType
}

template <class ValueType>
class FeQ<ValueType> {
public:
    // Class similar to current FeView

    ValueType operator()(uint shape_idx, BulkPoint point);
    ValueType operator()(uint shape_idx, SidePoint point);
    // Implementation for EdgePoint, SidePoint, and JoinPoint shoud have a common implementation
    // resolving to side values

private:
    // attributes:
    PatchFeValues *fe_values_;
    uint begin_;    /// Index of the first component of the Quantity. Size is given by ValueType
}



template<unsigned int spacedim = 3>
class PatchFEValues {
public:

    PatchFEValues()

    /**
     * @brief Return the product of Jacobian determinant and the quadrature
     * weight at given quadrature point.
     *
     * @param quads List of quadratures.
     */
    inline ElQ<Scalar> JxW(std::vector<Quadrature *> quads)
    {
        uint begin = this->n_columns_;
        n_columns_++; // scalar needs one column
        // TODO store to map?? JxW will be pre-computed to column 'begin'.
        return ElQ<Scalar>(&this, begin);
    }

private:
    uint n_columns_;  ///< Number of columns

};


#endif /* PATCH_FE_VALUES_HH_ */
