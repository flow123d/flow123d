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


/*
 * Expand of std::tuple - functional solution:
 * https://stackoverflow.com/questions/687490/how-do-i-expand-a-tuple-into-variadic-template-functions-arguments
 */
// ------------- UTILITY---------------
template<int...> struct index_tuple{};

template<int I, typename IndexTuple, typename... Types>
struct make_indexes_impl;

template<int I, int... Indexes, typename T, typename ... Types>
struct make_indexes_impl<I, index_tuple<Indexes...>, T, Types...>
{
    typedef typename make_indexes_impl<I + 1, index_tuple<Indexes..., I>, Types...>::type type;
};

template<int I, int... Indexes>
struct make_indexes_impl<I, index_tuple<Indexes...> >
{
    typedef index_tuple<Indexes...> type;
};

template<typename ... Types>
struct make_indexes : make_indexes_impl<0, index_tuple<>, Types...>
{};

// ----------UNPACK TUPLE AND APPLY TO FUNCTION ---------
template<class Ret, class Result, class... Args, int... Indexes >
Ret apply_helper( Ret (*pf)(Args...), Result &res, index_tuple< Indexes... >, std::tuple<Args...>&& tup)
{
    return pf( res, std::forward<Args>( std::get<Indexes>(tup))... );
}

template<class Ret, class Result, class ... Args>
Ret apply(Ret (*pf)(Args...), Result &res, const std::tuple<Args...>&  tup)
{
    return apply_helper(pf, res, typename make_indexes<Args...>::type(), std::tuple<Args...>(tup));
}

template<class Ret, class Result, class ... Args>
Ret apply(Ret (*pf)(Args...), Result &res, std::tuple<Args...>&&  tup)
{
    return apply_helper(pf, res, typename make_indexes<Args...>::type(), std::forward<tuple<Args...>>(tup));
}


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


template<int spacedim, class Value, class Fn, class ... Args>
class FieldModel : FieldCached<spacedim, Value> {
private:
	Fn fn;
    std::tuple<Args...> inputs;

public:
    FieldModel(Fn functor, Args... args)
    : fn(functor), inputs( std::make_tuple(std::forward<Args>(args)...) )
    { static_assert( std::is_same<typename Value::return_type, typename Fn::Result>::value, "Non-convertible functor type!"); }

    void cache_update() {
        auto f0 = std::get<0>(inputs);
        auto f1 = std::get<1>(inputs);

        for(unsigned int i_cache=0; i_cache<this->fvc.size(); ++i_cache) {
            this->fvc.data().template mat<Value::NRows_, Value::NCols_>(i_cache) =
                    fn( f0.data().template mat<Value::NRows_, Value::NCols_>(i_cache),
    	                f1.data().template mat<Value::NRows_, Value::NCols_>(i_cache) );
    	} // */
        // apply(fn.compute, this->fvc, inputs);
        // Solutions for call apply method:
        // - https://www.tutorialspoint.com/function-pointer-to-member-function-in-cplusplus
        // - https://stackoverflow.com/questions/1485983/calling-c-class-methods-via-a-function-pointer
        // - https://stackoverflow.com/questions/54651448/c-passing-overloaded-operator-of-class-as-function-pointer
    }

};


using Scalar = typename arma::Col<double>::template fixed<1>;
using Vector = arma::vec3;
using Tensor = arma::mat33;


class Fn {
public:
    typedef Vector Result;
    typedef Scalar Param0;
    typedef Vector Param1;
    typedef std::tuple< std::shared_ptr<FieldCached<3, Param0>>, std::shared_ptr<FieldCached<3, Param1>> > DepFields;

    static int compute(FieldValueCache<double> &res, FieldValueCache<double> &param0, FieldValueCache<double> &param1) {
        for(unsigned int i_cache=0; i_cache<res.size(); ++i_cache) {
            res.data().template mat<3, 1>(i_cache) =
                    param1.data().template mat<3, 1>(i_cache) * param0.data().template mat<1, 1>(i_cache);
        }
        return 0;
    }

    Result operator() (Param0 a, Param1 v) {
        return a * v;
    }
};




#endif /* FIELD_MODEL_HH_ */
