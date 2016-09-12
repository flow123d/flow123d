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
 * @file    vec_seq_double.hh
 * @brief   
 */

#ifndef VECTOR_SEQ_DOUBLE_HH_
#define VECTOR_SEQ_DOUBLE_HH_

#include <vector>
#include <memory>
#include "system/global_defs.h"

#include <petscvec.h>

#include "fields/field_elementwise.hh"


/**
 * Auxiliary class for output elementwise concentration vectors
 * in convection transport, sorptions, dual porosity etc.
 *
 * Stores data in two formats:
 *  - shared pointer to std::vector of double
 *  - pointer to PETSC vector that use same data
 *
 * Allows the following functionalities:
 *  - return shared pointer to std::vector of double
 *  - return pointer to PETSC vector
 *  - create shared pointer to FieldElementwise object corresponding with std::vector of double
 */
class VectorSeqDouble {
public:
    typedef typename std::shared_ptr< std::vector<double> > VectorSeq;

	/// Create shared pointer and PETSC vector with given size.
	void resize(unsigned int size)
	{
		data_ptr_ = std::make_shared< std::vector<double> >(size);
		chkerr(VecCreateSeqWithArray(PETSC_COMM_SELF, 1, size, &((*data_ptr_)[0]), &data_petsc_));
		chkerr(VecZeroEntries( data_petsc_ ));
	}

	/// Getter for shared pointer of output data.
	VectorSeq get_data_ptr()
	{
		return data_ptr_;
	}



	/// Getter for PETSC vector of output data (e.g. can be used by scatters).
	Vec &get_data_petsc()
	{
		return data_petsc_;
	}

	/// Create and return shared pointer to FieldElementwise object
	template <int spacedim, class Value>
	std::shared_ptr<FieldElementwise<spacedim, Value> > create_field(unsigned int n_comp)
	{
		std::shared_ptr<FieldElementwise<spacedim, Value> > field_ptr(
		          new FieldElementwise<spacedim, Value>( data_ptr_, n_comp ));
		return field_ptr;
	}

	/// Destructor.
	~VectorSeqDouble()
	{
		if (data_ptr_) chkerr(VecDestroy(&data_petsc_));
	}

    /**
     * Access to the vector element on index @p idx.
     */
    inline double &operator[](unsigned int idx)
    {
    	ASSERT_DBG(idx < data_ptr_->size()) (idx) (data_ptr_->size());
    	return (*data_ptr_)[idx];
    }

private:
	/// shared pointer to vector of data
	VectorSeq data_ptr_;
	/// stored vector of data in PETSC format
	Vec data_petsc_;
};

/**
 * Like VectorSeqDouble but for MPI PETSC vectors. Have acces to local part.
 */
class VectorMPI {
public:
    typedef typename std::vector<double> VectorData;
    typedef typename std::shared_ptr< VectorData > VectorDataPtr;

    VectorMPI() {
    }

    /// Create shared pointer and PETSC vector with given size. COLLECTIVE.
    VectorMPI(unsigned int local_size) {
        resize(local_size);
    }

    /**
     * Resize the vector to given local size. Operation is allowed only if this object is
     * a unique vector object pointing to the actual data.
     */
    void resize(unsigned int local_size) {
        if (data_ptr_.use_count() ==0) {
            data_ptr_ = std::make_shared< std::vector<double> >(local_size);
        } else {
            ASSERT_DBG( data_ptr_.use_count() ==  1 ) ( data_ptr_.use_count() ).error("Object referenced by other pointer. Can not resize.");
            chkerr(VecDestroy(&data_petsc_));
            data_ptr_->resize(local_size);
        }
        chkerr(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, local_size, PETSC_DECIDE, &((*data_ptr_)[0]), &data_petsc_));
        chkerr(VecZeroEntries( data_petsc_ ));
    }

    /// Return new vector with same parallel structure.
    void duplicate(VectorMPI other) {
        this->resize(other.data().size());
    }

    /// Getter for shared pointer of output data.
    VectorDataPtr data_ptr()
    {
        return data_ptr_;
    }

    /// Getter for PETSC vector of output data (e.g. can be used by scatters).
    Vec &petsc_vec()
    {
        return data_petsc_;
    }

    void zero_entries() {
        chkerr(VecZeroEntries( data_petsc_ ));
    }

    VectorData &data()
    {
        ASSERT_DBG(data_ptr_);
        return *data_ptr_;
    }

    /*
    /// Create and return shared pointer to FieldElementwise object
    template <int spacedim, class Value>
    std::shared_ptr<FieldElementwise<spacedim, Value> > create_field(unsigned int n_comp)
    {
        std::shared_ptr<FieldElementwise<spacedim, Value> > field_ptr(
                  new FieldElementwise<spacedim, Value>( data_ptr_, n_comp ));
        return field_ptr;
    }*/

    void swap(VectorMPI &other) {
        OLD_ASSERT_EQUAL(this->data_ptr_->size(), other.data_ptr_->size());
        uint size = this->data_ptr_->size();
        std::swap(this->data_ptr_, other.data_ptr_);
        chkerr(VecDestroy(&data_petsc_));
        chkerr(VecDestroy(&other.data_petsc_));
        chkerr(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, size, PETSC_DECIDE, &((*data_ptr_)[0]), &data_petsc_));
        chkerr(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, size, PETSC_DECIDE, &((*other.data_ptr_)[0]), &other.data_petsc_));
    }


    void copy(VectorMPI &other) {
        OLD_ASSERT_EQUAL(this->data_ptr_->size(), other.data_ptr_->size());
        chkerr(VecCopy(other.data_petsc_, data_petsc_));
    }

    /// Destructor.
    ~VectorMPI()
    {
        if (data_ptr_.use_count() == 1) chkerr(VecDestroy(&data_petsc_));
    }

    /**
     * Access to the vector element on index @p idx.
     */
    inline double &operator[](unsigned int idx)
    {
        ASSERT_DBG(data_ptr_);
        ASSERT_DBG(idx < data_ptr_->size()) (idx) (data_ptr_->size());
        return (*data_ptr_)[idx];
    }

private:

    /// shared pointer to vector of data
    VectorDataPtr data_ptr_;
    /// stored vector of data in PETSC format
    Vec data_petsc_;
};



#endif /* VECTOR_SEQ_DOUBLE_HH_ */
