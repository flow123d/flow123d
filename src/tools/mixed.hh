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
#include <type_traits>

const int __spacedim = 3;

using Dim = unsigned int;

template<template<Dim ...> class T>
using _MixedBase0 = std::tuple<T<0>, T<1>, T<2>, T<3>>;

template<template<Dim ...> class T>
using _MixedBase1 = std::tuple<T<1>, T<2>, T<3>>;


template< template<Dim ...> class T, int lower_dim=0>
class Mixed : public _MixedBase0<T> {
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
    template < template<Dim...> class TT>
    Mixed(const Mixed<TT> &other)
    : _MixedBase0<T>(
    		T<0>(other.get<0>() ),
			T<1>(other.get<1>() ),
			T<2>(other.get<2>() ),
			T<3>(other.get<3>() ) )
    { static_assert(std::is_convertible<TT<0>, T<0>>::value, "Non-convertible types!"); }

    Mixed(const T<0> &p0,const T<1> &p1,const T<2> &p2,const T<3> &p3)
    : _MixedBase0<T>(p0, p1, p2, p3)
    {}

    Mixed(T<0> &&p0, T<1> &&p1, T<2> &&p2, T<3> &&p3)
    : _MixedBase0<T>(p0, p1, p2, p3)
    {}


    template<typename... Args>
    Mixed(Args... args)
    : _MixedBase0<T>(
            T<0>(std::forward<Args>(args)...),
            T<1>(std::forward<Args>(args)...),
            T<2>(std::forward<Args>(args)...),
            T<3>(std::forward<Args>(args)...))
    {}

    template<Dim i_dim>
    T<i_dim> &get() {
        return std::get<i_dim>(*this);
    }

    template<Dim i_dim>
    const T<i_dim> &get() const {
        return std::get<i_dim>(*this);
    }

    template< template<Dim...> class TParent, typename std::enable_if<std::is_convertible<TParent<0>, T<0>>::value, T<0> >::type >
    operator Mixed<TParent> () const {
    	//ASSERT(std::is_base_of<TParent, T>::value);
        return Mixed<TParent>(
        		TParent<0>(this->get<0>()),
				TParent<1>(this->get<1>()),
				TParent<2>(this->get<2>()),
				TParent<3>(this->get<3>()));
    }

    // Possible collective methods must be implemented in MixedPtr childs.
};




template< template<Dim ...> class T>
class Mixed<T, 1> : public _MixedBase1<T> {
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
    template < template<Dim...> class TT>
    Mixed(const Mixed<TT, 1> &other)
    : _MixedBase1<T>(
            T<1>(other.get<1>() ),
            T<2>(other.get<2>() ),
            T<3>(other.get<3>() ) )
    { static_assert(std::is_convertible<TT<1>, T<1>>::value, "Non-convertible types!"); }

    Mixed(const T<1> &p1,const T<2> &p2,const T<3> &p3)
    : _MixedBase1<T>( p1, p2, p3)
    {}

    Mixed( T<1> &&p1, T<2> &&p2, T<3> &&p3)
    : _MixedBase1<T>(p1, p2, p3)
    {}


    template<typename... Args>
    Mixed(Args... args)
    : _MixedBase1<T>(
            T<1>(std::forward<Args>(args)...),
            T<2>(std::forward<Args>(args)...),
            T<3>(std::forward<Args>(args)...))
    {}

    template<Dim i_dim>
    T<i_dim> &get() {
        return std::get<i_dim-1>(*this);
    }

    template<Dim i_dim>
    const T<i_dim> &get() const {
        return std::get<i_dim-1>(*this);
    }

    template< template<Dim...> class TParent, typename std::enable_if<std::is_convertible<TParent<1>, T<1>>::value, T<1> >::type >
    operator Mixed<TParent> () const {
        //ASSERT(std::is_base_of<TParent, T>::value);
        return Mixed<TParent>(
                TParent<1>(this->get<1>()),
                TParent<2>(this->get<2>()),
                TParent<3>(this->get<3>()));
    }

    // Possible collective methods must be implemented in MixedPtr childs.
};




//template< template<Dim, Dim> class T>
//struct FixSpaceDim {
///**
// * Partial template resolution (parameter binding).
// * See implementation of MixedSpaceDim for the usage. Can be used for
// * possible child classes of the Mixed template.
// */
//public:
//    template<Dim dim>
//    using type =  T<dim, __spacedim>;
//};
//
//
//template< template<Dim, Dim> class T>
//using MixedSpaceDim = Mixed< FixSpaceDim<T>::template type >;






template<template<Dim...> class T>
using _MixedPtrBase0 = std::tuple<
        std::shared_ptr<T<0>>, std::shared_ptr<T<1>>,
        std::shared_ptr<T<2>>, std::shared_ptr<T<3>>>;

template<template<Dim...> class T>
using _MixedPtrBase1 = std::tuple< std::shared_ptr<T<1>>,
        std::shared_ptr<T<2>>, std::shared_ptr<T<3>>>;


template< template<Dim...> class T, int lower_dim=0>
class MixedPtr : public _MixedPtrBase0<T> {
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
    template <Dim dim>
    using TPtr = std::shared_ptr<T<dim>>;

    template < template<Dim...> class TT>
    MixedPtr(const MixedPtr<TT> &other)
    : _MixedPtrBase0<T>(
            other.get<0>(),
            other.get<1>(),
            other.get<2>(),
            other.get<3>())
    {}

    MixedPtr(TPtr<0> p0, TPtr<1> p1, TPtr<2> p2, TPtr<3> p3)
    : _MixedPtrBase0<T>(p0, p1, p2, p3)
    {}

    template<typename... Args>
    MixedPtr(Args... args)
    : _MixedPtrBase0<T>(
            std::make_shared<T<0>>(std::forward<Args>(args)...),
            std::make_shared<T<1>>(std::forward<Args>(args)...),
            std::make_shared<T<2>>(std::forward<Args>(args)...),
            std::make_shared<T<3>>(std::forward<Args>(args)...))
    {}


    template<Dim i_dim>
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


template< template<Dim...> class T>
class MixedPtr<T, 1> : public _MixedPtrBase1<T> {
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
    template <Dim dim>
    using TPtr = std::shared_ptr<T<dim>>;

    template < template<Dim...> class TT>
    MixedPtr(const MixedPtr<TT, 1> &other)
    : _MixedPtrBase1<T>(
            other.get<1>(),
            other.get<2>(),
            other.get<3>())
    {}

    MixedPtr(TPtr<1> p1, TPtr<2> p2, TPtr<3> p3)
    : _MixedPtrBase1<T>(p1, p2, p3)
    {}

    template<typename... Args>
    MixedPtr(Args... args)
    : _MixedPtrBase1<T>(
            std::make_shared<T<1>>(std::forward<Args>(args)...),
            std::make_shared<T<2>>(std::forward<Args>(args)...),
            std::make_shared<T<3>>(std::forward<Args>(args)...))
    {}


    template<Dim i_dim>
    std::shared_ptr<T<i_dim>> get() {
        return std::get<i_dim-1>(*this);
    }

    template<int i_dim>
    const std::shared_ptr<T<i_dim>> get() const {
        return std::get<i_dim-1>(*this);
    }

    Mixed<T, 1> operator *() {
        return Mixed<T, 1>(
                *(this->get<1>()),
                *(this->get<2>()),
                *(this->get<3>()));
    }

    // Possible collective methods must be implemented in MixedPtr childs.

};

//template< template<Dim, Dim> class T>
//using MixedSpaceDimPtr = MixedPtr< FixSpaceDim<T>::template type >;



//template< template <int> class Fn, class ...Args>
//auto dim_switch(Dim dim, Args&&... args) -> decltype( builder.makeObject() )
//{
//    switch (dim)
//    {
//        case 0:
//            return Fn<0>(std::forward<Args>(args)...);
//        case 0:
//            return Fn<1>(std::forward<Args>(args)...);
//        case 0:
//            return Fn<2>(std::forward<Args>(args)...);
//        case 0:
//            return Fn<3>(std::forward<Args>(args)...);
//    }
//
//}


#endif /* SRC_TOOLS_MIXED_HH_ */
