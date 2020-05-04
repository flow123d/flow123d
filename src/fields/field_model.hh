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
#include "fields/field.hh"
#include "fields/multi_field.hh"

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
     * Extract cached values on index i_cache from the input fields and call given function on these values.
     */

    /// base case for building up arguments for the function call
    template< typename CALLABLE, typename FIELD_TUPLE, int INDEX >
    struct model_cache_item
    {
        template< typename... Vs >
        static auto eval(int i_cache, CALLABLE f, FIELD_TUPLE fields, Vs&&... args) -> decltype(auto) {
            const auto &single_field = std::get < INDEX - 1 > (std::forward<decltype(fields)>(fields));
            return model_cache_item<CALLABLE, FIELD_TUPLE, INDEX - 1>::eval(
                    i_cache, f, std::forward<decltype(fields)>(fields),
                    single_field[i_cache], std::forward<Vs>(args)...);

        }
    };

    /// terminal case - do the actual function call
    template< typename CALLABLE, typename FIELD_TUPLE >
    struct model_cache_item< CALLABLE, FIELD_TUPLE, 0 >
    {
        template< typename... Vs >
        static auto eval(int i_cache, CALLABLE f, FIELD_TUPLE fields, Vs&&... args) -> decltype(auto)
        {
            return f(std::forward<Vs>(args)...);
        };
    };

    /**
     * Check common number of components of the input fields/multifields.
     * Return number of components.
     * Return 0 for no multifields.
     * Throw for different number of components.
     */
    template<typename FIELD_TUPLE, int INDEX >
    struct n_components {
        static uint eval(FIELD_TUPLE fields, uint n_comp) {
            const auto &single_field = std::get < INDEX - 1 > (std::forward<decltype(fields)>(fields));
            uint n_comp_new = single_field.n_comp();
            if (n_comp == 0) {
                n_comp = n_comp_new;
            } else {
                ASSERT_DBG(n_comp == n_comp_new);
            }
            return n_components<FIELD_TUPLE, INDEX - 1>::eval(std::forward<decltype(fields)>(fields), n_comp);
        };
    };

    template<typename FIELD_TUPLE>
    struct n_components<FIELD_TUPLE, 0> {
        static uint eval(FIELD_TUPLE fields, uint n_comp)
        {
            return n_comp;
        };
    };


    /**
     * Return component 'i_comp' of the multifield 'f' or the field 'f'.
     */
    /*template<class FIELD>
    auto field_component(FIELD f, uint i_comp) -> decltype(auto)
    {
        if (f.is_multifield()) {
            return f[i_comp];
        } else {
            return f;
        }
    }*/

    /**
     * Return component 'i_comp' of the multifield 'f'.
     */
    template<int spacedim, class Value>
    auto field_component(const MultiField<spacedim, Value> &f, uint i_comp) -> decltype(auto)
    {
        ASSERT(f.is_multifield());
        return f[i_comp];
    }

    /**
     * Return the field 'f'. Variant to previous method.
     */
    template<int spacedim, class Value>
    auto field_component(const Field<spacedim, Value> &f, uint i_comp) -> decltype(auto)
    {
        ASSERT(!f.is_multifield());
        return f;
    }

    /**
     * For given tuple of fields/multifields and given component index
     * return  the tuple of the selected components.
     */
    template<typename FIELD_TUPLE, int INDEX>
    struct get_components {
        static auto eval(FIELD_TUPLE fields, uint i_comp) -> decltype(auto)
        {
            const auto &single_field = std::get < INDEX - 1 > (std::forward<decltype(fields)>(fields));
            return std::tuple_cat(
                    std::forward_as_tuple(field_component(single_field, i_comp)),
                    get_components<FIELD_TUPLE, INDEX - 1>::eval(
                            std::forward<decltype(fields)>(fields), i_comp)
                    );
        };
    };

    template<typename FIELD_TUPLE>
    struct get_components<FIELD_TUPLE, 0> {
        static auto eval(FIELD_TUPLE fields, uint n_comp) -> decltype(auto)
        {
            return std::forward_as_tuple<>();
        };
    };




    template<typename Function, typename Tuple>
    auto call(Function f, Tuple t)
    {
    }
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

    // Create ElementCacheMap
    std::shared_ptr<EvalPoints> eval_points = std::make_shared<EvalPoints>();
    eval_poinzs->add_bulk<3>( quad ); // initialize EvalPoints
    ElementCacheMap elm_cache_map;
    elm_cache_map.init(eval_points);

    // Definition of fields
    Field<3, FieldValue<3>::Scalar > f_scal;
    Field<3, FieldValue<3>::VectorFixed > f_vec;
    Field<3, FieldValue<3>::VectorFixed > result;
    ... // fill data to fields f_scal, f_vec

    // create instance FieldModel class, use helper method Model::create to simply passsing of parameters
  	auto f_product = Model<3, FieldValue<3>::VectorFixed>::create(FnProduct(), f_scal, f_vec);
  	// set field on all regions
    result.set_mesh( *mesh );
  	result.set_field(mesh->region_db().get_region_set("ALL"), f_product);
    result.cache_allocate(eval_points);
    result.set_time(tg.step(), LimitSide::right);

  	// cache_update
    result.cache_update(elm_cache_map);
   @endcode
 *
 */
template<int spacedim, class Value, typename Fn, class ... InputFields>
class FieldModel : public FieldAlgorithmBase<spacedim, Value>
{
private:
	Fn* fn;
	typedef std::tuple<InputFields...> FieldsTuple;
    FieldsTuple input_fields;

public:
    typedef typename FieldAlgorithmBase<spacedim, Value>::Point Point;

    FieldModel(Fn* func, InputFields... args)
    : fn(func), input_fields( std::forward_as_tuple((args)...) )
    {}




    void cache_update(FieldValueCache<typename Value::element_type> &data_cache,
				ElementCacheMap &cache_map, unsigned int region_idx) override {
        auto update_cache_data = cache_map.update_cache_data();
        unsigned int region_in_cache = update_cache_data.region_cache_indices_range_.find(region_idx)->second;
        unsigned int i_cache_el_begin = update_cache_data.region_value_cache_range_[region_in_cache];
        unsigned int i_cache_el_end = update_cache_data.region_value_cache_range_[region_in_cache+1];
        for(unsigned int i_cache=i_cache_el_begin; i_cache<i_cache_el_end; ++i_cache) {
            data_cache.data().set(i_cache) =
                detail::model_cache_item<
                    Fn,
                    decltype(input_fields),
                    std::tuple_size<FieldsTuple>::value
                >::eval(i_cache, fn, input_fields);
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
    typedef FieldAlgorithmBase<spacedim, Value> FieldBaseType;
    typedef std::shared_ptr< FieldBaseType > FieldBasePtr;

    template<typename Fn, class ... InputFields>
    static auto create(Fn *fn,  InputFields&&... inputs) -> decltype(auto)
    {
        return std::make_shared<FieldModel<spacedim, Value, Fn, InputFields...>>(fn, std::forward<InputFields>(inputs)...);
    }



    template<typename Function, typename Tuple, size_t ... I>
    static auto call_create(Function f, Tuple t, std::index_sequence<I ...>)
    {
        return create(f, std::get<I>(t) ...);
    }

    template<typename Fn, class ... InputFields>
    static auto create_multi(Fn *fn,  InputFields&&... inputs) -> decltype(auto)
    {
        typedef std::tuple<InputFields...> FieldTuple;
        FieldTuple field_tuple = std::forward_as_tuple((inputs)...);
        constexpr uint n_inputs = sizeof...(InputFields);
        uint n_comp = detail::n_components< FieldTuple, n_inputs>::eval(field_tuple, 0);
        ASSERT_DBG(n_comp > 0);
        std::vector<FieldBasePtr> result_components;
        for(uint i=0; i<n_comp; i++) {
            const auto & component_of_inputs = detail::get_components< FieldTuple, n_inputs>::eval(field_tuple, i);
            //const auto & all_args = std::tuple_cat(std::make_tuple(fn), component_of_inputs);
            //FieldBasePtr component_field = detail::call(create<Fn, InputFields&&...>, all_args);

            FieldBasePtr component_field = call_create(fn, component_of_inputs, std::make_index_sequence<n_inputs>{});
            result_components.push_back(component_field);
        }

        return result_components;
    }
};










#endif /* FIELD_MODEL_HH_ */
