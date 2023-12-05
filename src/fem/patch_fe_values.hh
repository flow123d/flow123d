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

template<unsigned int spacedim> class PatchFEValues;



using Scalar = double;
using Vector = arma::vec3;
using Tensor = arma::mat33;

template <class ValueType>
class ElQ {
public:
    /// Constructor
    ElQ()
    : fe_values_(nullptr), begin_(0) {}

    /// Constructor
    ElQ(PatchFEValues<3> *fe_values, unsigned int begin)
    : fe_values_(fe_values), begin_(begin) {}

    ValueType operator()(const BulkPoint &point)
    {
    	return 0.0;
    }


    ValueType operator()(const SidePoint &point)
    {
    	return 0.0;
    }


private:
    // attributes:
    PatchFEValues<3> *fe_values_;
    unsigned int begin_;    /// Index of the first component of the Quantity. Size is given by ValueType
};

template <class ValueType>
class FeQ {
public:
    /// Constructor
    FeQ()
    : fe_values_(nullptr), begin_(0) {}

    // Class similar to current FeView
    FeQ(PatchFEValues<3> *fe_values, unsigned int begin)
    : fe_values_(fe_values), begin_(begin) {}


    ValueType operator()(unsigned int shape_idx, const BulkPoint &point)
    {
    	return 0.0;
    }

    ValueType operator()(unsigned int shape_idx, const SidePoint &point)
    {
    	return 0.0;
    }

    // Implementation for EdgePoint, SidePoint, and JoinPoint shoud have a common implementation
    // resolving to side values

private:
    // attributes:
    PatchFEValues<3> *fe_values_;
    unsigned int begin_;    /// Index of the first component of the Quantity. Size is given by ValueType
};



template<unsigned int spacedim = 3>
class PatchFEValues {
public:

    PatchFEValues(unsigned int n_quad_points=0)
    : n_columns_(0) {}

    /**
	 * @brief Initialize structures and calculates cell-independent data.
	 *
	 * @param _quadrature The quadrature rule for the cell associated
     *                    to given finite element or for the cell side.
	 * @param _fe The finite element.
	 * @param _flags The update flags.
	 */
    template<unsigned int DIM>
    void initialize(Quadrature &_quadrature,
                    FiniteElement<DIM> &_fe,
                    UpdateFlags _flags)
    {}

    /// Reinit data.
    void reinit(PatchElementsList patch_elements)
    {}

    /**
     * @brief Return the product of Jacobian determinant and the quadrature
     * weight at given quadrature point.
     *
     * @param quad_list List of quadratures.
     */
    inline ElQ<Scalar> JxW(std::vector<Quadrature *> quad_list)
    {
        uint begin = this->n_columns_;
        n_columns_++; // scalar needs one column
        // TODO store to map?? JxW will be pre-computed to column 'begin'.
        return ElQ<Scalar>(&this, begin);
    }

    /**
     * @brief Return the value of the @p function_no-th shape function at
     * the @p p quadrature point.
     *
     * @param quad_list List of quadratures.
     * @param function_no Number of the shape function.
     */
    inline FeQ<Scalar> scalar_shape(std::vector<Quadrature *> quad_list, unsigned int n_comp)
    {
        uint begin = this->n_columns_;
        n_columns_ += n_comp; // scalar needs one column x n_comp
        // TODO store to map?? shape_value will be pre-computed to column 'begin'.
        return FeQ<Scalar>(&this, begin);
    }

private:
    uint n_columns_;  ///< Number of columns

};


#endif /* PATCH_FE_VALUES_HH_ */
