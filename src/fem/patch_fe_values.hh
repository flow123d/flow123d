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



template<unsigned int spacedim = 3>
class PatchFEValues : public FEValuesBase<PatchFEValues<spacedim>, spacedim> {
public:
    /// Default constructor
	PatchFEValues()
    : FEValuesBase<PatchFEValues<spacedim>, spacedim>(), used_size_(0) {}

    /// Set new value of used_size_
    void resize(unsigned int new_size);

	inline unsigned int used_size() const {
	    return used_size_;
	}

	inline unsigned int max_size() const {
	    return element_data_.size();
	}

protected:
    class ElementFEData
    {
    public:
        ElementFEData() {}

        /// Shape functions evaluated at the quadrature points.
        std::vector<std::vector<double> > shape_values_;

        /// Gradients of shape functions evaluated at the quadrature points.
        /// Each row of the matrix contains the gradient of one shape function.
        std::vector<std::vector<arma::vec::fixed<spacedim> > > shape_gradients_;

        /// Auxiliary object for calculation of element-dependent data.
        std::shared_ptr<ElementValues<spacedim> > elm_values_;

    };

    /// Implement @p FEValuesBase::allocate_in
    void allocate_in() override;

    /// Implement @p FEValuesBase::initialize_in
    void initialize_in(Quadrature &q, unsigned int dim) override;

//    /// Auxiliary storage of FEValuesViews accessors.
//    ViewsCache views_cache_;

    /// Data of elements on patch
    std::vector<ElementFEData> element_data_;

    /// Number of elements on patch. Must be less or equal to size of element_data vector
    unsigned int used_size_;

};

#endif /* PATCH_FE_VALUES_HH_ */
