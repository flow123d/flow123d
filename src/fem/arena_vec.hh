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
 * @file    arena_vec.hh
 */

#ifndef ARENA_VEC_HH_
#define ARENA_VEC_HH_

#include "fem/arena_resource.hh"
#include "system/asserts.hh"

#include <Eigen/Core>
#include <Eigen/Dense>



template<class T> class ArenaOVec;


/**
 * Define vector allocated in ArenaResource and aligned to SIMD size.
 */
template<class T>
class ArenaVec {
public:
    /// Type definition
	typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VecData;
	typedef Eigen::Array<T, Eigen::Dynamic, 1> ArrayData;

	/// Default constructor, set invalid data pointer
    ArenaVec()
    : data_ptr_(nullptr), data_size_(0), arena_(nullptr), scalar_val_( (T)0 ) {}

    /**
     * Constructor. Set scalar value
     */
    ArenaVec(T scalar_val)
    : data_ptr_(nullptr), data_size_(0), arena_(nullptr), scalar_val_(scalar_val) {}

    /**
     * Constructor. Set sizes and allocate data pointer
     */
    ArenaVec(size_t data_size, PatchArena &arena)
    : data_ptr_(nullptr), data_size_(data_size), arena_(&arena), scalar_val_( (T)0 ) {
        data_ptr_ = arena_->allocate_simd<T>( data_size_ );
    }

    /// Copy constructor
    ArenaVec(const ArenaVec<T> &other)
    : data_ptr_(other.data_ptr_), data_size_(other.data_size_),
      arena_(other.arena_), scalar_val_(other.scalar_val_)
    {}

    /**
     * Maps data pointer to Eigen Map of dimensions given data_size_ and returns it.
     */
    inline Eigen::Map<VecData> eigen_map() {
        ASSERT_PTR(data_ptr_);
        return Eigen::Map<VecData>(data_ptr_, data_size_, 1);
    }

    /// Smae as previous but with const modifier
    inline const Eigen::Map<VecData> eigen_map() const {
        ASSERT_PTR(data_ptr_);
        return Eigen::Map<VecData>(data_ptr_, data_size_, 1);
    }

    /**
     * Maps data pointer to Eigen Map of dimensions given data_size_ and returns it.
     */
    inline Eigen::Map<ArrayData> array_map() {
        ASSERT_PTR(data_ptr_);
        return Eigen::Map<ArrayData>(data_ptr_, data_size_, 1);
    }

    /// Smae as previous but with const modifier
    inline const Eigen::Map<ArrayData> array_map() const {
        ASSERT_PTR(data_ptr_);
        return Eigen::Map<ArrayData>(data_ptr_, data_size_, 1);
    }

    /// Return data pointer (development method)
    T* data_ptr() {
        return data_ptr_;
    }

    /// Smae as previous but return const pointer
    const T* data_ptr() const {
        return data_ptr_;
    }

    /// Getter for data_size_
    inline size_t data_size() const {
    	return data_size_;
    }

    /// Getter for arena_
    PatchArena &arena() {
        return *arena_;
    }

    /// Getter for scalar_val_
    T scalar_val() const {
        return scalar_val_;
    }

    /// Set pointer to PatchArena
    inline void set_patch_arena(PatchArena &arena) {
        ASSERT_PTR(arena_);
        this->arena_ = &arena;
    }

    /// Returns copied vector of square root values
    inline ArenaVec<T> sqrt() const {
        ASSERT_PTR(data_ptr_);
        ArenaVec<T> res(data_size_, *arena_);
        Eigen::Map<ArrayData> result_map = res.array_map();
        result_map = this->array_map().sqrt();
        return res;
    }

    /// Returns copied vector of inverse values
    inline ArenaVec<T> inverse() const {
        ASSERT_PTR(data_ptr_);
        ArenaVec<T> res(data_size_, *arena_);
        Eigen::Map<ArrayData> result_map = res.array_map();
        result_map = this->array_map().inverse();
        return res;
    }

    /// Returns copied vector of absolute values
    inline ArenaVec<T> abs() const {
        ASSERT_PTR(data_ptr_);
        ArenaVec<T> res(data_size_, *arena_);
        Eigen::Map<ArrayData> result_map = res.array_map();
        result_map = this->array_map().abs();
        return res;
    }

    /// For development only. TODO remove
    inline T & operator()(std::size_t item) {
        if (data_ptr_ == nullptr) {
            return scalar_val_;
        }
        ASSERT_LT(item, data_size_);
        return data_ptr_[item];
    }

    /// For development only. TODO remove
    inline const T & operator()(std::size_t item) const {
        if (data_ptr_ == nullptr) {
            return scalar_val_;
        }
        ASSERT_LT(item, data_size_);
        return data_ptr_[item];
    }

    /// Assignment operator
    inline ArenaVec<T> &operator=(const ArenaVec<T> &other) {
        data_ptr_ = other.data_ptr_;
        data_size_ = other.data_size_;
        arena_ = other.arena_;
        scalar_val_ = other.scalar_val_;
        return *this;
    }

    /**
     * Addition operator
     *
     * Sums two ArenaVec objects
     * TODO If ASSERT_PTR is thrown needs to implement tests of scalar values (see multiplication operator).
     */
    inline ArenaVec<T> operator+(const ArenaVec<T> &other) const {
        ASSERT_PTR(data_ptr_);
        ASSERT_PTR(other.data_ptr());
        ASSERT_EQ(data_size_, other.data_size());
        ArenaVec<T> res(data_size_, *arena_);
        Eigen::Map<VecData> result_map = res.eigen_map();
        result_map = this->eigen_map() + other.eigen_map();
        return res;
    }

    /**
     * Subtraction operator
     *
     * Subtracts two ArenaVec objects
     * TODO If ASSERT_PTR is thrown needs to implement tests of scalar values (see multiplication operator).
     */
    inline ArenaVec<T> operator-(const ArenaVec<T> &other) const {
        ASSERT_PTR(data_ptr_);
        ASSERT_PTR(other.data_ptr());
        ASSERT_EQ(data_size_, other.data_size());
        ArenaVec<T> res(data_size_, *arena_);
        Eigen::Map<ArrayData> result_map = res.array_map();
        result_map = this->array_map() - other.array_map();
        return res;
    }

    /**
     * Multiplication operator
     *
     * Product Scalar x ArenaVec
     */
    inline ArenaVec<T> operator*(T multi) const {
        ASSERT_PTR(data_ptr_);
        ArenaVec<T> res(data_size_, *arena_);
        Eigen::Map<ArrayData> result_map = res.array_map();
        result_map = this->array_map() * multi;
        return res;
    }

    /**
     * Multiplication operator
     *
     * Product of two ArenaVec objects, checks if one ofe them is define as scalar
     */
    inline ArenaVec<T> operator*(const ArenaVec<T> &other) const {
        if (data_ptr_ == nullptr)
            return other * scalar_val_;
        if (other.data_ptr() == nullptr)
            return this->operator *(other.scalar_val());
        ASSERT_EQ(data_size_, other.data_size());
        ArenaVec<T> res(data_size_, *arena_);
        Eigen::Map<ArrayData> result_map = res.array_map();
        result_map = this->array_map() * other.array_map();
        return res;
    }

    /**
     * Division operator
     *
     * Divides ArenaVec / Scalar
     */
    inline ArenaVec<T> operator/(T div_by) const {
        ASSERT_PTR(data_ptr_);
        ArenaVec<T> res(data_size_, *arena_);
        Eigen::Map<ArrayData> result_map = res.array_map();
        result_map = this->array_map() / div_by;
        return res;
    }

    /**
     * Division operator
     *
     * Divides two ArenaVec objects
     * TODO If ASSERT_PTR is thrown needs to implement tests of scalar values (see multiplication operator).
     */
    inline ArenaVec<T> operator/(const ArenaVec<T> &other) const {
        ASSERT_PTR(data_ptr_);
        ASSERT_PTR(other.data_ptr());
        ASSERT_EQ(data_size_, other.data_size());
        ArenaVec<T> res(data_size_, *arena_);
        Eigen::Map<ArrayData> result_map = res.array_map();
        result_map = this->array_map() / other.array_map();
        return res;
    }

protected:
    /// Constructor. Allows create ArenaVec from ArenaOVec
    ArenaVec(T* data_ptr, size_t data_size, PatchArena &arena)
    : data_ptr_(data_ptr), data_size_(data_size), arena_(&arena) {}

    T* data_ptr_;            ///< Pointer to data array
    size_t data_size_;       ///< Length of data array
    PatchArena *arena_;      ///< Pointer to Arena where intermediate calculations and results are stored, should be changed by set_patch_arena
    T scalar_val_;           ///< Scalar value of T type

    friend class ArenaOVec<T>;
};



/**
 * Define vector allocated in ArenaResource based on ArenaVec with overwrite
 * multiplication operator that executes outer product.
 *
 * Example of usage with conversions between ArenaVec and ArenaOVec:
   @code
    ArenaVec<double> outer_product(ArenaVec<double> a, ArenaVec<double> b) {
        // Convert ArenaVec inputs to ArenaOVec variables
        ArenaOVec<double> a_ovec(a);
        ArenaOVec<double> b_ovec(b);

        // performs outer product
        ArenaOVec<double> result_ovec = a * b;

        // reverse conversion to ArenaVec
        return result_ovec.get_vec();
    }
   @endcode
 *
 * If we consider that size of input vector 'a' is 'M' and size of input vector 'b' is 'N'
 * then size of returned vector is 'M*N'.
 */
template<class T>
class ArenaOVec : public ArenaVec<T> {
public:
    /// Default constructor
    ArenaOVec()
    : ArenaVec<T>() {}

    /// Copy constructor
    ArenaOVec(const ArenaOVec<T> &other)
    : ArenaVec<T>(other) {}

    /// Constructor. Set scalar_val
    ArenaOVec(T scalar_val)
    : ArenaVec<T>(scalar_val) {}

    /**
     * Constructor creates ArenaOVec on data of ArenaVec
     */
    ArenaOVec(const ArenaVec<T> &vec) {
        ASSERT_PTR(vec.data_ptr());

        this->data_ptr_ = vec.data_ptr_;
        this->data_size_ = vec.data_size_;
        this->arena_ = vec.arena_;
    }

    /// Convert ArenaOVec to ArenaVec and its
    ArenaVec<T> get_vec() const {
        return ArenaVec<T>(*this);
    }

    /// Assignment operator
    inline ArenaOVec<T> &operator=(const ArenaOVec<T> &other) {
        ArenaVec<T>::operator=(other);
        return *this;
    }

    /// Addition operator
    inline ArenaOVec<T> operator+(const ArenaOVec<T> &other) const {
        // Test of valid data_ptr is in constructor
        ASSERT_EQ(this->data_size_, other.data_size());
        ArenaVec<T> res_vec(this->data_size_, *this->arena_);
        ArenaOVec<T> res(res_vec);
        Eigen::Map<typename ArenaVec<T>::VecData> result_map = res.eigen_map();
        result_map = this->eigen_map() + other.eigen_map();
        return res;
    }

    /// Multiplication operator
    inline ArenaOVec<T> operator*(const ArenaOVec<T> &other) const {
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatData;

        ArenaVec<T> res_vec(this->data_size_*other.data_size(), *this->arena_);
        ArenaOVec<T> res(res_vec);
        Eigen::Map<MatData> result_map = Eigen::Map<MatData>(res.data_ptr(), this->data_size_, other.data_size());
        result_map = this->eigen_map() * other.eigen_map().transpose();
        return res;
    }
};


#endif /* ARENA_VEC_HH_ */
