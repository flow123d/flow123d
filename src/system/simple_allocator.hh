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
 * @file    simple_allocator.cc
 * @ingroup system
 * @brief   Declaring internal allocator class which does not monitor memory
 *          consumption
 */


namespace internal {

    /**
     * Simple allocator which uses default malloc and free functions. Fields allocated
     *   with this allocator will not be included in overall memory consumption.
     */
    template<class T>
    class SimpleAllocator {
    public:
        typedef size_t size_type;
        typedef ptrdiff_t difference_type;
        typedef T *pointer;
        typedef const T *const_pointer;
        typedef T &reference;
        typedef const T &const_reference;
        typedef T value_type;

        SimpleAllocator() { }

        SimpleAllocator(const SimpleAllocator &) { }


        pointer allocate(size_type n, const void * = 0) {
            T *t = (T *) malloc(n * sizeof(T));
            return t;
        }

        void deallocate(void *p, size_type) {
            if (p) {
                free(p);
            }
        }

        pointer address(reference x) const { return &x; }

        const_pointer address(const_reference x) const { return &x; }

        SimpleAllocator<T> &operator=(const SimpleAllocator &) { return *this; }

        void construct(pointer p, const T &val) { new((T *) p) T(val); }

        void destroy(pointer p) { p->~T(); }

        size_type max_size() const { return size_t(-1); }

        template<class U>
        struct rebind {
            typedef SimpleAllocator<U> other;
        };

        template<class U>
        SimpleAllocator(const SimpleAllocator<U> &) { }

        template<class U>
        SimpleAllocator &operator=(const SimpleAllocator<U> &) { return *this; }
    };
}