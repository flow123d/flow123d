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
 * @brief   Template GeneralIterator serves as general template for internal iterators.
 */

#ifndef GENERAL_ITERATOR_HH_
#define GENERAL_ITERATOR_HH_

/** @brief General iterator template.
 * Provides iterator over objects in some container.
 * 
 * Requires the template object to implement:
 * - comparison operator==()
 * - increment operator++()
 */
template<class Object>
class GeneralIterator
{
public:
//     GeneralIterator();

    GeneralIterator(const Object& object);

    /// equal operator
    bool operator==(const GeneralIterator& other);
    /// non-equal operator
    bool operator!=(const GeneralIterator& other);

    ///  * dereference operator
    const Object& operator*() const;

    /// -> dereference operator
    const Object* operator->() const;

    /// prefix increment
    GeneralIterator& operator++();

private:
    /// Output element of the output mesh.
    Object object_; 
};


// --------------------------------------------------- GeneralIterator INLINE implementation -----------
// inline GeneralIterator::GeneralIterator()
// {}

template<class Object>
inline GeneralIterator<Object>::GeneralIterator(const Object& object)
: object_(object)
{}

template<class Object>
inline bool GeneralIterator<Object>::operator==(const GeneralIterator& other)
{
    return (object_ == other.object_);
}

template<class Object>
inline bool GeneralIterator<Object>::operator!=(const GeneralIterator& other)
{
    return !( *this == other);
}

template<class Object>
inline const Object& GeneralIterator<Object>::operator*() const
{
    return object_;
}

template<class Object>
inline const Object* GeneralIterator<Object>::operator->() const
{
    return &object_;
}

template<class Object>
inline GeneralIterator<Object>& GeneralIterator<Object>::operator++()
{
    object_.inc();
    return (*this);
}

#endif // GENERAL_ITERATOR_HH_