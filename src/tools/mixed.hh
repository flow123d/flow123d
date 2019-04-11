/*
 * mixed.hh
 *
 *  Created on: Apr 8, 2019
 *      Author: jb
 */

#ifndef SRC_TOOLS_MIXED_HH_
#define SRC_TOOLS_MIXED_HH_

#include <tuple>
#include <memory>
#include <iostream>

const int __spacedim = 3;



template<template<int> class T>
using _MixedBase = std::tuple<T<0>, T<1>, T<2>, T<3>>;

template< template<int dim> class T>
class Mixed : public _MixedBase<T> {
/**
 * Template to simplify storing and passing tuples of instances of dimension parametrized templates.
 * Currently instances for dim = 0,1,2,3 are created. We assume spacedim=3.
 * Usage:
 *
 * fe_order = 10;
 * auto gauss = Mixed<QGauss>(q_order);      // template <int dim> class QGauss;
 * auto mapping = MixedSpaceDim<Mapping>();  // template <int dim, int spacedim> class Mapping;
 * ...
 * fe_values = FEValues<dim>(mapping.get<dim>(), gauss.get<dim>(), ...)
 *
 * All parameters of the Mixed<T>(...) constructor are forwarder (by value) to the T(...) constructor for every dimension.
 * Only forwarding by value is supported due to problems with MixedPtr.
 * Can not resolve copy constructor correctly.
 *
 */
public:
    template<typename... Args>
    Mixed(Args... args)
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

    template<int i_dim>
    const T<i_dim> &get() const {
        return std::get<i_dim>(*this);
    }

    // Possible collective methods must be implemented in MixedPtr childs.
};


template< template<int, int> class T>
struct FixSpaceDim {
/**
 * Partial template resolution (parameter binding).
 * See implementation of MixedSpaceDim for the usage. Can be used for
 * possible child classes of the Mixed template.
 */
public:
    template<int dim>
    using type =  T<dim, __spacedim>;
};


template< template<int, int> class T>
using MixedSpaceDim = Mixed< FixSpaceDim<T>::template type >;






template<template<int> class T>
using _MixedPtrBase = std::tuple<
        std::shared_ptr<T<0>>, std::shared_ptr<T<1>>,
        std::shared_ptr<T<2>>, std::shared_ptr<T<3>>>;



template< template<int dim> class T>
class MixedPtr : public _MixedPtrBase<T> {
/**
 * Template to simplify storing and passing tuples of shatred_ptr to instances of dimension parametrized templates.
 * Currently instances for dim = 0,1,2,3 are created. We assume spacedim=3.
 * Usage:
 *
 * fe_order = 10;
 * auto gauss = MixedPtr<QGauss>(q_order);      // template <int dim> class QGauss;
 * auto mapping = MixedSpaceDimPtr<Mapping>();  // template <int dim, int spacedim> class Mapping;
 * ...
 * fe_values = FEValues<dim>(*mapping.get<dim>(), *gauss.get<dim>(), ...)
 *
 * All parameters of the Mixed<T>(...) constructor are forwarder (by value) to the T(...) constructor for every dimension.
 * Only forwarding by value is supported due to problems with MixedPtr.
 * We are unable to give precedence to the copy constructor over the prefect forwarding constructor MixedPtr(Args&& ...).
 */
public:

    template < template<int dim> class TT>
    MixedPtr(const MixedPtr<TT> &other)
    : _MixedPtrBase<T>(
            other.get<0>(),
            other.get<1>(),
            other.get<2>(),
            other.get<3>())
    {}

    template<typename... Args>
    MixedPtr(Args... args)
    : _MixedPtrBase<T>(
            std::make_shared<T<0>>(std::forward<Args>(args)...),
            std::make_shared<T<1>>(std::forward<Args>(args)...),
            std::make_shared<T<2>>(std::forward<Args>(args)...),
            std::make_shared<T<3>>(std::forward<Args>(args)...))
    {}


    template<int i_dim>
    std::shared_ptr<T<i_dim>> get() {
        return std::get<i_dim>(*this);
    }

    template<int i_dim>
    const std::shared_ptr<T<i_dim>> get() const {
        return std::get<i_dim>(*this);
    }

    Mixed<T> operator *() {
        return Mixed<T>(
                *(this->get<0>()),
                *(this->get<1>()),
                *(this->get<2>()),
                *(this->get<3>()));
    }

    // Possible collective methods must be implemented in MixedPtr childs.

};

template< template<int, int> class T>
using MixedSpaceDimPtr = MixedPtr< FixSpaceDim<T>::template type >;



#endif /* SRC_TOOLS_MIXED_HH_ */
