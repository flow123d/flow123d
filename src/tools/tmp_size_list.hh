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
 * @file    tmp_size_list.hh
 * @brief
 */

#ifndef TMP_SIZE_LIST_HH_
#define TMP_SIZE_LIST_HH_


#include <new>
#include "system/asserts.hh"


/**
 * @brief Struct is a container that encapsulates variable size arrays.
 *
 * Allows to:
 *  1. Add new items to container, items are stored as temporary.
 *  2. Mark block of temporary items as final (use method finalize_tmp)
 *     or cancelled temporary item (method revert_tmp)
 *
 * This algorith allows to add blocks of data, evaluates external condition
 * and possibly reverts unfinished block if condition is not met.
 */
template<class Type>
struct TmpSizeList {
public:
    /// Constructor, create new instance with reserved size
    TmpSizeList(std::size_t reserved_size)
    : tmp_size_(0), final_size_(0)
    {
        data_.resize(reserved_size);
    }

    /**
     * Resize to new reserved size.
     *
     * New size must be higher than actual size!
     */
    void resize(std::size_t new_size)
    {
    	ASSERT_GT(new_size, max_size());
    	data_.resize(new_size);
    }

    /// Return final size (part of finalize data).
    inline std::size_t final_size() const
    {
        return final_size_;
    }

    /// Return temporary size (full size of stored data).
    inline std::size_t tmp_size() const
    {
        return tmp_size_;
    }

    /// Return maximal (reserved) size.
    inline std::size_t max_size() const
    {
        return data_.size();
    }

    /// Add new item to list.
    inline std::size_t add(const Type &t)
    {
        ASSERT_LT_DBG(tmp_size_, max_size()).error("Data array overflowed!\n");
        data_[tmp_size_] = t;
    	tmp_size_++;
        return tmp_size_;
    }

    /// Finalize temporary part of data.
    inline std::size_t finalize_tmp()
    {
        final_size_ = tmp_size_;
        return tmp_size_;
    }

    /// Erase temporary part of data.
    inline std::size_t revert_tmp()
    {
        tmp_size_ = final_size_;
        return tmp_size_;
    }

    /// Reset stored data.
    inline void reset()
    {
        tmp_size_ = 0;
        final_size_ = 0;
    }

    /// Return item on given position
    const Type &operator[](std::size_t pos) const {
        ASSERT_LT_DBG(pos, tmp_size_).error("Position is out of data size!\n");
        return data_[pos];
    }

private:
    std::vector<Type> data_;  ///< Vector of items.
    std::size_t tmp_size_;    ///< Temporary size (full size of used data).
    std::size_t final_size_;  ///< Final size of data (part of finalize data).
};

#endif /* TMP_SIZE_LIST_HH_ */
