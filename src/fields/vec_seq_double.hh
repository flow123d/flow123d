/*
 * vec_seq_double.hh
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
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
 *  - pointer to PETSC vector
 *
 * Allows the following functionalities:
 *  - return shared pointer to std::vector of double
 *  - return pointer to PETSC vector
 *  - create shared pointer to FieldElementwise object corresponding with std::vector of double
 */
class VectorSeqDouble {
public:
	typedef typename std::shared_ptr< std::vector<double> > VectorSeq;

	/// Constructor.
	VectorSeqDouble(unsigned int size)
	{
		data_ptr_ = std::make_shared< std::vector<double> >(size);
		data_petsc_ = (Vec *)&data_ptr_;
	}

	/// Return shared pointer to output data.
	VectorSeq get_data_ptr()
	{
		return data_ptr_;
	}

	/// Return shared pointer to PETSC vector.
	Vec *get_data_petsc()
	{
		return data_petsc_;
	}

	/// Create and return shared pointer to FieldElementwise object
	template <int spacedim, class Value> //3, FieldValue<3>::Scalar
	std::shared_ptr<FieldElementwise<spacedim, Value> > create_field(unsigned int n_comp)
	{
		// TODO: use constructor FieldElementwise(shared_ptr, unsigned int) after the implementation of data cache
		std::shared_ptr<FieldElementwise<spacedim, Value> > field_ptr(
		          new FieldElementwise<spacedim, Value>( *(data_ptr_.get()), n_comp ));
		return field_ptr;
	}

    /**
     * Returns reference to the sub-field (component) of given index @p idx.
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
	Vec *data_petsc_;
};


#endif /* VECTOR_SEQ_DOUBLE_HH_ */
