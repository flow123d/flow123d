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
#include <utility>
#include <type_traits>

#include "fields/field.hh"
#include "fields/field_common.hh"
#include "fields/field_values.hh"
#include "fields/field_value_cache.hh"


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


/**
 * Parent class of Field and FieldModel.
 *
 * Will be descendant of FieldCommon in future (needs implementation of pure virtual method).
 */
template<int spacedim, class Value>
class FieldCached /*: public FieldCommon*/ {
public:
	/// Constructor
	FieldCached()
	: fvc(Value::NRows_, Value::NCols_) {}

protected:
    FieldValueCache<typename Value::element_type> fvc;
};


template<int spacedim, class Value, class Fn, class ... InputFields>
class FieldModel : FieldCached<spacedim, Value> {
private:
	Fn fn;
	//using FnResult = typename std::result_of<Fn()>::type;
	typedef std::tuple<InputFields...> FieldsTuple;
    FieldsTuple inputs;

public:
    FieldModel(Fn functor, InputFields... args)
    : fn(functor), inputs( std::make_tuple(std::forward<InputFields>(args)...) )
    {
//        static_assert( std::is_same<typename Value::return_type, FnResult>::value,
//                "Non-convertible functor type!");
    }




    void cache_update()  {
        for(unsigned int i_cache=0; i_cache<this->fvc.size(); ++i_cache) {
            this->fvc.data().template mat<Value::NRows_, Value::NCols_>(i_cache) =
                    //fn( std::get<0>(inputs)[i_cache], std::get<1>(inputs)[i_cache]);
                detail::model_cache_item<
                    Fn,
                    decltype(inputs),
                    std::tuple_size<FieldsTuple>::value
                >::eval(i_cache, fn, inputs);
    	}
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
