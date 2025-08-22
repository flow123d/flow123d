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
#include <unordered_map>
#include <boost/functional/hash.hpp>      // for boost::hash_value
#include "fem/dh_cell_accessor.hh"
#include "tools/revertable_list.hh"
#include "mesh/range_wrapper.hh"

template <int spacedim> class ElementAccessor;



/// Define Integral Tuple hash function - helper struct of IntegralPtrMap
struct IntegralTplHash {
    std::size_t operator()(std::tuple<uint, uint> tpl) const {
        return boost::hash_value( tpl );
    }

    /// Create tuple from dimennsion and size of Quadrature
    static std::tuple<uint, uint> integral_tuple(uint dim, uint quad_size) {
        return std::make_tuple(dim, quad_size);
    }

};

/// Alias for unordered_map of shared_ptr<Integral> with custom hash
template<typename Integral>
using IntegralPtrMap = std::unordered_map<std::tuple<uint, uint>, std::shared_ptr<Integral>, IntegralTplHash>;



/// Define Operation hash function
template<typename Operation>
struct OperationPtrHash {
    std::size_t operator()(const Operation* key) const {
        std::size_t h1 = std::hash<std::string>()(typeid(*key).name());
        std::size_t h2 = std::hash<std::size_t>()(key->quad()->size());
        return h1 ^ (h2 << 1);
    }
};

/// Content-based equality for Operation*
template<typename Operation>
struct OperationPtrEqual {
    bool operator()(const Operation* lhs, const Operation* rhs) const {
        return *lhs == *rhs;
    }
};

/// Alias for unordered_set of shared_ptr<Integral> with custom hash
template<typename Operation>
using OperationSet = std::unordered_set<Operation *, OperationPtrHash<Operation>, OperationPtrEqual<Operation>>;



#endif /* INTEGRAL_DATA_HH_ */
