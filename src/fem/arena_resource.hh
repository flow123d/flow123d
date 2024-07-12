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
 * @file    arena_resource.hh
 */

#ifndef ARENA_RESOURCE_HH_
#define ARENA_RESOURCE_HH_

#include <memory_resource>
#include <vector>
#include <iostream>
#include <new>
#include <stdexcept>   // !! Use Flow exception mechanism

#include "system/asserts.hh"


// Final proposal of Arena
template <class Resource>
class ArenaResource : public std::pmr::memory_resource {
protected:
    /// Inner constructor, used only in construction of child arena
    ArenaResource(void *buffer, size_t buffer_size, size_t simd_alignment)
    : upstream_(std::pmr::get_default_resource()),
      buffer_(buffer),
      buffer_size_(buffer_size),
      used_size_(0),
      resource_(buffer_, buffer_size, upstream_),
      simd_alignment_(simd_alignment),
      full_data_(false)
    {
        ASSERT_PERMANENT_EQ( (buffer_size%simd_alignment), 0 );
    }

    /// Returns different upstream resource in debug / release mode
	static inline std::pmr::memory_resource* upstream_resource() {
#ifdef FLOW123D_DEBUG
    	return std::pmr::null_memory_resource();
#else
    	return std::pmr::get_default_resource();
#endif
    }

public:
    /// Constructor. Creates assembly arena
    ArenaResource(size_t buffer_size, size_t simd_alignment, std::pmr::memory_resource* upstream = ArenaResource<Resource>::upstream_resource())
    : upstream_(upstream), // TODO needs use buffer and resource of upstream if default upstream is not used
      buffer_( upstream_->allocate(buffer_size, simd_alignment) ),
      buffer_size_(buffer_size),
      used_size_(0),
      resource_(buffer_, buffer_size, upstream_),
      simd_alignment_(simd_alignment),
      full_data_(false)
    {
        ASSERT_PERMANENT_EQ( (buffer_size%simd_alignment), 0 );
    }


    ~ArenaResource() = default;

    /**
     * Create and return child arena.
     *
     * Child arena is created in free space of actual arena.
     * Actual arena is marked as full (flag full_data_) and cannot allocate new data.
     */
    ArenaResource *get_child_arena() {
        void *p = this->raw_allocate(1, simd_alignment_);
        size_t used_size = (char *)p - (char *)buffer_;
        size_t free_space = buffer_size_ - used_size;
        full_data_ = true;
        return new ArenaResource(p, free_space, simd_alignment_);
    }


    /// Getter for resource
    Resource &resource() {
    	return resource_;
    }

    /// Allocate and return data pointer of n_item array of type T (alignment to length 8 bytes)
    template <class T>
    T* allocate_8(size_t n_items) {
        size_t bytes = sizeof(T) * n_items;
        return (T*)this->raw_allocate(bytes, 8);
    }

    /// Allocate and return data pointer of n_item array of type T (alignment to length given by simd_alignment constructor argument)
    template <class T>
    T* allocate_simd(size_t n_items) {
        size_t bytes = sizeof(T) * n_items;
        return (T*)this->raw_allocate(bytes, simd_alignment_);
    }

    // Reset allocated data
    void reset() {
        resource_.release();
        used_size_ = 0;
#ifdef FLOW123D_DEBUG
    	char *c_buffer = (char *)buffer_;
    	for (size_t i=0; i<buffer_size_; ++i)
    	    c_buffer[i] = 0;
#endif
    }

protected:
    void* raw_allocate(size_t bytes, size_t alignment) {
        ASSERT(!full_data_).error("Allocation of new data is not possible because child arena was created.");
        ASSERT_EQ(buffer_size_%alignment, 0);

        void* p = resource_.allocate(bytes, alignment);
        used_size_ += (bytes + alignment - 1) / alignment * alignment;
        if (p == nullptr) {  // test only in Debug when null_pointer_resource is in use
            throw std::bad_alloc();
        }
        return p;
    }

    /// Override do_allocate to handle allocation logic
    void* do_allocate(size_t bytes, size_t alignment) override {
        return raw_allocate(bytes, alignment);
    }

    /// Override do_deallocate (no-op for monotonic buffer)
    void do_deallocate(FMT_UNUSED void* p, FMT_UNUSED size_t bytes, FMT_UNUSED size_t alignment) override {
        // No-op
    }

    /// Override do_is_equal for memory resource comparison
    bool do_is_equal(const std::pmr::memory_resource& other) const noexcept override {
        return this == &other;
    }

private:
    std::pmr::memory_resource* upstream_;   ///< Pointer to upstream
    void* buffer_;                          ///< Pointer to buffer
    size_t buffer_size_;                    ///< Size of buffer
    size_t used_size_;                      ///< Temporary data member, will be removed
    Resource resource_;                     ///< Resource of arena
    size_t simd_alignment_;                 ///< Size of SIMD alignment
    bool full_data_;                        ///< Flag signs full data (child arena is created)
};


using AssemblyArena = ArenaResource<std::pmr::monotonic_buffer_resource>;
using PatchArena = AssemblyArena;


#endif /* ARENA_RESOURCE_HH_ */
