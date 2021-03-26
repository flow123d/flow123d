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
 * Provides iterator over objects of type Object in some container.
 *
 * Operators '*' and '->' returns objects of type ObjectOut
 * 
 * Requires the template object to implement:
 * - comparison operator==()
 * - increment operator++()
 */
template<class Object>
class Iter
{
public:
//     IterConvert();

	Iter(const Object& object);

    /// equal operator
    bool operator==(const Iter& other);
    /// non-equal operator
    bool operator!=(const Iter& other);

    ///  * dereference operator
    const Object& operator*() const;

    /// -> dereference operator
    const Object* operator->() const;

    /// prefix increment
    Iter& operator++();

protected:
    Object object_;
};


/** @brief General iterator template.
 * Provides iterator over objects of type ObjectIn in some container.
 * Same as previous but allows conversion of output to type ObjectOut.
 */
template<class ObjectIn, class ObjectOut>
class IterConvert : public Iter<ObjectIn>
{
public:
//     IterConvert();

	IterConvert(const ObjectIn& object);

    ///  * dereference operator
    const ObjectOut& operator*() const;

    /// -> dereference operator
    const ObjectOut* operator->() const;

private:
    mutable ObjectOut out_;
};


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

template<class Object>
inline Iter<Object>::Iter(const Object& object)
: object_(object)
{}

template<class Object>
inline bool Iter<Object>::operator==(const Iter& other)
{
    return (object_ == other.object_);
}

template<class Object>
inline bool Iter<Object>::operator!=(const Iter& other)
{
    return !( *this == other);
}

template<class Object>
inline const Object& Iter<Object>::operator*() const
{
    return object_;
}

template<class Object>
inline const Object* Iter<Object>::operator->() const
{
    return &object_;
}

template<class Object>
inline Iter<Object>& Iter<Object>::operator++()
{
    object_.inc();
    return (*this);
}

template<class ObjectIn, class ObjectOut>
inline IterConvert<ObjectIn, ObjectOut>::IterConvert(const ObjectIn& object)
: Iter<ObjectIn>(object)
{}

template<class ObjectIn, class ObjectOut>
inline const ObjectOut& IterConvert<ObjectIn, ObjectOut>::operator*() const
{
    out_ = (ObjectOut)this->object_;
    return out_;
}

template<class ObjectIn, class ObjectOut>
inline const ObjectOut* IterConvert<ObjectIn, ObjectOut>::operator->() const
{
    out_ = (ObjectOut)this->object_;
    return &out_;
}

#endif // GENERAL_ITERATOR_HH_
