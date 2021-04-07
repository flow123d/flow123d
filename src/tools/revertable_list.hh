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
 * @file    revertable_list.hh
 * @brief
 */

#ifndef REVARTABLE_LIST_HH_
#define REVARTABLE_LIST_HH_


#include <new>
#include "system/asserts.hh"


/**
 * @brief Struct is a container that encapsulates variable size arrays.
 *
 * Allows to:
 *  1. Add new items to container (method push_back, items are stored as temporary.
 *  2. Mark block of temporary items as final (use method make_permanent)
 *     or cancelled temporary item (method revert_temporary)
 *
 * This algorith allows to add blocks of data, evaluates external condition
 * and possibly reverts unfinished block if condition is not met.
 *
 * Reserved (maximal) size is set in constructor. This size can be enlarged manually
 * through method resize or constructor accepts parameter enlarged_by. If this value
 * is greater than 0 size of container is automatically enlarged during call push_back
 * if container is full.
 */
template<class Type>
struct RevertableList {
public:
    /// Constructor, create new instance with reserved size
	RevertableList(std::size_t reserved_size, std::size_t enlarged_by = 0)
    : temporary_size_(0), permanent_size_(0), enlardeg_by_(enlarged_by)
    {
        data_.resize(reserved_size);
    }

    /// Copy constructor
	RevertableList(const RevertableList& other)
    : data_(other.data_), temporary_size_(other.temporary_size_), permanent_size_(other.permanent_size_), enlardeg_by_(other.enlardeg_by_)
    {}

    /**
     * Resize to new reserved size.
     *
     * New size must be higher than actual size!
     */
    void resize(std::size_t new_size)
    {
    	ASSERT_GT(new_size, reserved_size());
    	data_.resize(new_size);
    }

    /// Return permanent size of list.
    inline std::size_t permanent_size() const
    {
        return permanent_size_;
    }

    /// Return temporary size of list (full size of stored data).
    inline std::size_t temporary_size() const
    {
        return temporary_size_;
    }

    /// Return reserved (maximal) size.
    inline std::size_t reserved_size() const
    {
        return data_.size();
    }

    /**
     * Add new item of list.
     *
     * New item is added to end of list and temporary size value is incremented.
     * Method is equivalent with std::vector::push_back().
     * This method needs to create copy of passed Type.
     */
    inline std::size_t push_back(const Type &t)
    {
        ASSERT_DBG((enlardeg_by_ > 0) || (temporary_size_ < reserved_size())).error("Data array overflowed!\n");
        if (temporary_size_ == reserved_size()) { // enlarge reserved size
        	this->resize( this->reserved_size() + enlardeg_by_ );
        }
        data_[temporary_size_] = t;
        temporary_size_++;
        return temporary_size_;
    }

    /**
     * Create new item in list.
     *
     * New item is created at the end of list and temporary size value is incremented.
     * Method is equivalent with std::vector::emplace_back().
     * Method passes argumets of Type constructor.
     */
    template<class... Args>
    inline std::size_t emplace_back(Args&&... args)
    {
        ASSERT_DBG((enlardeg_by_ > 0) || (temporary_size_ < reserved_size())).error("Data array overflowed!\n");
        if (temporary_size_ == reserved_size()) { // enlarge reserved size
        	this->resize( this->reserved_size() + enlardeg_by_ );
        }
        data_[temporary_size_] = Type( std::forward<Args>(args)... );
        temporary_size_++;
        return temporary_size_;
    }

    /// Finalize temporary part of data.
    inline std::size_t make_permanent()
    {
    	permanent_size_ = temporary_size_;
        return temporary_size_;
    }

    /// Erase temporary part of data.
    inline std::size_t revert_temporary()
    {
    	temporary_size_ = permanent_size_;
        return temporary_size_;
    }

    /// Clear the list.
    inline void reset()
    {
    	temporary_size_ = 0;
    	permanent_size_ = 0;
    }

    inline typename std::vector<Type>::iterator begin()
    {
    	return data_.begin();
    }

    inline typename std::vector<Type>::iterator end()
    {
    	return data_.begin() + permanent_size_;
    }

    /// Return item on given position
    const Type &operator[](std::size_t pos) const {
        ASSERT_LT_DBG(pos, temporary_size_).error("Position is out of data size!\n");
        return data_[pos];
    }

private:
    std::vector<Type> data_;      ///< Vector of items.
    std::size_t temporary_size_;  ///< Temporary size (full size of used data).
    std::size_t permanent_size_;  ///< Final size of data (part of finalize data).
    std::size_t enlardeg_by_;     ///< Allow to enlarge list dynamically during call push_back if reserved size is full
};

#endif /* REVARTABLE_LIST_HH_ */
