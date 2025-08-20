/*!
 *
? * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
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
 * @file    integral_data.hh
 * @brief
 * @author  David Flanderka
 */

#ifndef INTEGRAL_DATA_HH_
#define INTEGRAL_DATA_HH_

#include <vector>
#include <memory>
#include <armadillo>
#include <unordered_set>
#include <boost/functional/hash.hpp>      // for boost::hash_value
#include "fem/dh_cell_accessor.hh"
#include "tools/revertable_list.hh"
#include "mesh/range_wrapper.hh"

template <int spacedim> class ElementAccessor;



/// Define Integral hash function - helper struct of IntegralPtrSet
template<typename Integral>
struct IntegralPtrHash {
    std::size_t operator()(const std::shared_ptr<Integral>& key) const {
        return boost::hash_value( std::make_tuple(key->dim(), key->quad()->size()) );
    }
};

/// Content-based equality for shared_ptr<Integral> - helper struct of IntegralPtrSet
template<typename Integral>
struct IntegralPtrEqual {
    bool operator()(const std::shared_ptr<Integral>& lhs, const std::shared_ptr<Integral>& rhs) const {
        return *lhs == *rhs;
    }
};

/// Alias for unordered_set of shared_ptr<Integral> with custom hash
template<typename Integral>
using IntegralPtrSet = std::unordered_set<std::shared_ptr<Integral>, IntegralPtrHash<Integral>, IntegralPtrEqual<Integral>>;


#endif /* INTEGRAL_DATA_HH_ */
