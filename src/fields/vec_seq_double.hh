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
		VecCreateSeqWithArray(PETSC_COMM_SELF, 1, size, &((*data_ptr_)[0]), &data_petsc_);
		VecZeroEntries( data_petsc_ );
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
		if (data_ptr_) VecDestroy(&data_petsc_);
	}

    /**
     * Access to the vector element on index @p idx.
     */
    inline double &operator[](unsigned int idx)
    {
    	ASSERT(idx < data_ptr_->size(), "Index is out of range.\n");
    	return (*data_ptr_)[idx];
    }

private:
	/// shared pointer to vector of data
	VectorSeq data_ptr_;
	/// stored vector of data in PETSC format
	Vec data_petsc_;
};


#endif /* VECTOR_SEQ_DOUBLE_HH_ */
