// --------------------------------------------------- GeneralIterator ---------------------------------
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
 * @file    general_iterator.hh
 * @brief   Template Iter serves as general template for internal iterators.
 */

#ifndef GENERAL_ITERATOR_HH_
#define GENERAL_ITERATOR_HH_

/** @brief General iterator template.
 * Provides iterator over objects of type ObjectIn in some container.
 *
 * Operators '*' and '->' returns objects of type ObjectOut
 * 
 * Requires the template object to implement:
 * - comparison operator==()
 * - increment operator++()
 */
template<class ObjectIn, class ObjectOut>
class IterConvert
{
public:
//     IterConvert();

	IterConvert(const ObjectIn& object);

    /// equal operator
    bool operator==(const IterConvert& other);
    /// non-equal operator
    bool operator!=(const IterConvert& other);

    ///  * dereference operator
    const ObjectOut& operator*() const;

    /// -> dereference operator
    const ObjectOut* operator->() const;

    /// prefix increment
    IterConvert& operator++();

private:
    /// Output element of the output mesh.
    ObjectIn object_;
    mutable ObjectOut out_;
};


/**
 * @brief General iterator template.
 *
 * Same as previous but doesn't provide specialization of operators '*' and '->'.
 */
template<class Object>
using Iter = IterConvert<Object, Object>;


/**
 * Create iterator from object
 */
template<class Object>
Iter<Object> make_iter(Object obj) {
	return Iter<Object>(obj);
}

/**
 * Create convertible iterator from object
 */
template<class ObjectIn, class ObjectOut>
IterConvert<ObjectIn, ObjectOut> make_iter(ObjectIn obj) {
	return IterConvert<ObjectIn, ObjectOut>(obj);
}


// --------------------------------------------------- Iter INLINE implementation -----------
// inline IterConvert::IterConvert()
// {}

template<class ObjectIn, class ObjectOut>
inline IterConvert<ObjectIn, ObjectOut>::IterConvert(const ObjectIn& object)
: object_(object)
{}

template<class ObjectIn, class ObjectOut>
inline bool IterConvert<ObjectIn, ObjectOut>::operator==(const IterConvert& other)
{
    return (object_ == other.object_);
}

template<class ObjectIn, class ObjectOut>
inline bool IterConvert<ObjectIn, ObjectOut>::operator!=(const IterConvert& other)
{
    return !( *this == other);
}

template<class ObjectIn, class ObjectOut>
inline const ObjectOut& IterConvert<ObjectIn, ObjectOut>::operator*() const
{
    out_ = (ObjectOut)object_;
    return out_;
}

template<class ObjectIn, class ObjectOut>
inline const ObjectOut* IterConvert<ObjectIn, ObjectOut>::operator->() const
{
    out_ = (ObjectOut)object_;
    return &out_;
}

template<class ObjectIn, class ObjectOut>
inline IterConvert<ObjectIn, ObjectOut>& IterConvert<ObjectIn, ObjectOut>::operator++()
{
    object_.inc();
    return (*this);
}

#endif // GENERAL_ITERATOR_HH_
