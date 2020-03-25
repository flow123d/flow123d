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
 * support for expanding tuples into an overloaded function call's arguments
 */
#define wrap_overload(func) [](auto&&... ps){ return func( std::forward<decltype(ps)>(ps)... ); }

namespace detail
{
//    //
//    // base case for building up arguments for the function call
//    //
//    template< typename CALLABLE, typename TUPLE, int INDEX >
//    struct tuple_into_callable_n
//    {
//        template< typename... Vs >
//        static auto apply(CALLABLE f, TUPLE t, Vs&&... args) -> decltype(auto)
//        {
//            return tuple_into_callable_n<CALLABLE, TUPLE, INDEX - 1>::apply(
//                f,
//                std::forward<decltype(t)>(t),
//                std::get<INDEX - 1>(std::forward<decltype(t)>(t)),
//                std::forward<Vs>(args)...
//            );
//        }
//    };
//
//    //
//    // terminal case - do the actual function call
//    //
//    template< typename CALLABLE, typename TUPLE >
//    struct tuple_into_callable_n< CALLABLE, TUPLE, 0 >
//    {
//        template< typename... Vs >
//        static auto apply(CALLABLE f, TUPLE t, Vs&&... args) -> decltype(auto)
//        {
//            return f(std::forward<Vs>(args)...);
//        };
//    };


    //
    // base case for building up arguments for the function call
    //
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

    //
    // terminal case - do the actual function call
    //
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

//template< typename FUNC, typename TUPLE >
//auto tuple_into_callable(FUNC f, TUPLE&& t) -> decltype(auto)
//{
//    return
//        detail::tuple_into_callable_n<
//            FUNC,
//            decltype(t),
//            std::tuple_size< std::remove_reference_t<TUPLE> >::value
//        >::apply(f, std::forward<decltype(t)>(t) );
//}


template<int spacedim, class Value, class Fn, class ... InputFields>
class FieldModel : public FieldAlgorithmBase<spacedim, Value>
{
private:
	Fn fn;
	//using FnResult = typename std::result_of<Fn()>::type;
	typedef std::tuple<InputFields...> FieldsTuple;
    FieldsTuple inputs;

public:
    typedef typename FieldAlgorithmBase<spacedim, Value>::Point Point;

    FieldModel(Fn functor, InputFields... args)
    : fn(functor), inputs( std::make_tuple(std::forward<InputFields>(args)...) )
    {
//        static_assert( std::is_same<typename Value::return_type, FnResult>::value,
//                "Non-convertible functor type!");
    }




    void cache_update(FieldValueCache<typename Value::element_type> &data_cache,
                unsigned int i_cache_el_begin, unsigned int i_cache_el_end,
                const std::vector< ElementAccessor<spacedim> > &element_set)  {
        for(unsigned int i_cache=i_cache_el_begin; i_cache<i_cache_el_end; ++i_cache) {
            data_cache.data().template mat<Value::NRows_, Value::NCols_>(i_cache) =
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

    template<class Fn, class ... InputFields>
    static auto create(Fn fn,  InputFields... inputs) -> decltype(auto) {
        return FieldModel<spacedim, Value, Fn, InputFields...>(fn, std::forward<InputFields>(inputs)...);
    }
};










#endif /* FIELD_MODEL_HH_ */
