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
 * @file    field_model.hh
 * @brief
 */

#ifndef FIELD_MODEL_HH_
#define FIELD_MODEL_HH_

#include <armadillo>
#include <iostream>
#include <tuple>
#include <string>
#include <vector>
#include <utility>
#include <type_traits>

#include "fields/field.hh"
#include "fields/field_common.hh"
#include "fields/field_algo_base.hh"
#include "fields/field_values.hh"
#include "fields/field_value_cache.hh"

template <int spacedim> class ElementAccessor;


/**
 * Wrapper for resolution of the overloaded functions passed as a parameter to the FieldModel.
 *
 * Example:
 @code
   double product(double x, double y) {...}
   Vec product(double x, Vec y) {...}
   Model<3, FieldValue<3>::VectorFixed>::create(wrapper_overload(product), f_scal, f_vec);
 @endcode
 *
 * Should automaticaly resolve the second model_cache_item::eval function.
 */
#define wrap_overload(func) [](auto&&... ps){ return func( std::forward<decltype(ps)>(ps)... ); }

namespace detail
{
    /**
     * base case for building up arguments for the function call
     */
    template< typename CALLABLE, typename TUPLE, int INDEX >
    struct model_cache_item
    {
        template< typename... Vs >
        static auto eval(int i_cache, CALLABLE f, TUPLE t, Vs&&... args) -> decltype(auto)
        {
            return model_cache_item<CALLABLE, TUPLE, INDEX - 1>::eval(
                i_cache,
                f,
                std::forward<decltype(t)>(t),
                std::get<INDEX - 1>(std::forward<decltype(t)>(t))[i_cache],
                std::forward<Vs>(args)...
            );
        }
    };

    /**
     * terminal case - do the actual function call
     */
    template< typename CALLABLE, typename TUPLE >
    struct model_cache_item< CALLABLE, TUPLE, 0 >
    {
        template< typename... Vs >
        static auto eval(int i_cache, CALLABLE f, TUPLE t, Vs&&... args) -> decltype(auto)
        {
            return f(std::forward<Vs>(args)...);
        };
    };


}


/**
 * Class representing field computing form results of other fields.
 *
 * Example of usage:
   @code
    // Functor class with defined operator ()
    class FnProduct {
    public:
        Vector operator() (Scalar a, Vector v) {
             return a(0) * v;
        }
    };

    ...

    // Definition of fields
    Field<3, FieldValue<3>::Scalar > f_scal;
    Field<3, FieldValue<3>::VectorFixed > f_vec;
    Field<3, FieldValue<3>::VectorFixed > result;
    ... // fill data to fields f_scal, f_vec

    // create instance FieldModel class, use helper method Model::create to simply passsing of parameters
  	auto f_product = Model<3, FieldValue<3>::VectorFixed>::create(FnProduct(), f_scal, f_vec);
  	// set field on all regions
  	result.set_field(mesh->region_db().get_region_set("ALL"), f_product);

  	// cache_update
  	FieldValueCache<double> &fvc = result.value_cache();
    f_product.cache_update(fvc, 0, fvc.size(), element_set);
   @endcode
 *
 */
template<int spacedim, class Value, typename Fn, class ... InputFields>
class FieldModel : public FieldAlgorithmBase<spacedim, Value>
{
private:
	Fn* fn;
	typedef std::tuple<InputFields...> FieldsTuple;
    FieldsTuple inputs;

public:
    typedef typename FieldAlgorithmBase<spacedim, Value>::Point Point;

    FieldModel(Fn* func, InputFields... args)
    : fn(func), inputs( std::make_tuple(std::forward<InputFields>(args)...) )
    {}




    void cache_update(FieldValueCache<typename Value::element_type> &data_cache,
				ElementCacheMap &cache_map, unsigned int region_idx)  {
        auto update_cache_data = cache_map.update_cache_data();
        std::unordered_map<unsigned int, typename ElementCacheMap::RegionData>::iterator reg_elm_it =
                update_cache_data.region_cache_indices_map_.find(region_idx);
        unsigned int region_in_cache = reg_elm_it->second.pos_;
        unsigned int i_cache_el_begin = update_cache_data.region_cache_indices_range_[region_in_cache];
        unsigned int i_cache_el_end = update_cache_data.region_cache_indices_range_[region_in_cache+1];
        for(unsigned int i_cache=i_cache_el_begin; i_cache<i_cache_el_end; ++i_cache) {
            data_cache.data().set(i_cache) =
                detail::model_cache_item<
                    Fn,
                    decltype(inputs),
                    std::tuple_size<FieldsTuple>::value
                >::eval(i_cache, fn, inputs);
    	}
    }

    /// Implementation of virtual method
    typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm) override {
        ASSERT(false).error("Forbidden method!\n");
        return this->r_value_;
    }

    /// Implementation of virtual method
    void value_list(const Armor::array &point_list, const ElementAccessor<spacedim> &elm,
                std::vector<typename Value::return_type>  &value_list) override {
        ASSERT(false).error("Forbidden method!\n");
    }

};


/**
 * Auxiliary class to avoid explicit specification of constructor template parameters.
 */
template<int spacedim, class Value>
class Model {
public:

    template<typename Fn, class ... InputFields>
    static auto create(Fn *fn,  InputFields... inputs) -> decltype(auto) {
        return FieldModel<spacedim, Value, Fn, InputFields...>(fn, std::forward<InputFields>(inputs)...);
    }
};










#endif /* FIELD_MODEL_HH_ */
