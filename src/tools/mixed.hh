/*
 * mixed.hh
 *
 *  Created on: Apr 8, 2019
 *      Author: jb
 */

#ifndef SRC_TOOLS_MIXED_HH_
#define SRC_TOOLS_MIXED_HH_

#include <tuple>

const int __spacedim = 3;



template<template<int> class T>
using _MixedBase = std::tuple<T<0>, T<1>, T<2>, T<3>>;

template< template<int dim> class T>
class Mixed : public _MixedBase<T> {
/**
 * Template to simplify storing and passing tuples of instances of dimension parametrized templates.
 * Currently instances for dim = 0,1,2,3 are created. So we assume spacedim=3.
 * Usage:
 *
 * fe_order = 10;
 * gauss = Mixed<QGauss>(q_order);
 * mapping = Mixed<Mapping>(); // specialization for Type<dim, spacedim> also provided.
 * ...
 * fe_values = FEValues<dim>(mapping.get<dim>(), gauss.get<dim>(), ...)
 *
 */
public:
    template<typename... Args>
    Mixed(Args&&... args)
    : _MixedBase<T>(
            T<0>(std::forward<Args>(args)...),
            T<1>(std::forward<Args>(args)...),
            T<2>(std::forward<Args>(args)...),
            T<3>(std::forward<Args>(args)...))
    {}

    template<int i_dim>
    T<i_dim> &get() {
        return std::get<i_dim>(*this);
    }

//    template<typename ...Args, typename R>
//    R forward( R ()) {
//        return
//    }
};




template< template<int, int> class T>
struct FixSpaceDim {
public:
    template<int dim>
    using type =  T<dim, __spacedim>;
};

//template< template<int, int> class T>
//using FixSpaceDim = typename _FixedSpaceDim<T>::template type;


template< template<int, int> class T>
using MixedSpaceDim = Mixed< FixSpaceDim<T>::template type >;



/**
 * Specialization for Type<dim, spacedim> templates.
 * Use fixed spacedim = __spacedim = 3.
 */


#endif /* SRC_TOOLS_MIXED_HH_ */
