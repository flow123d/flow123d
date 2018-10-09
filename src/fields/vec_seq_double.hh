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
#include "system/system.hh"
#include "system/global_defs.h"

#include "fields/fe_value_handler.hh"
#include "fields/field_fe.hh"
#include "fem/fe_p.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_system.hh"
#include "fem/dofhandler.hh"
#include "fem/finite_element.hh"

#include <petscvec.h>

class Mesh;
template <int spacedim, class Value> class FieldFE;


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
 *  - create shared pointer to FieldFE object corresponding with std::vector of double
 */
class VectorSeqDouble {
public:
    typedef typename std::shared_ptr< std::vector<double> > VectorSeq;

    /// Constructor.
    VectorSeqDouble() {}

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

	/// Getter for shared pointer of output data.
	unsigned int size()
	{
		ASSERT_PTR(data_ptr_).error("Uninitialized data vector.\n");
		return data_ptr_->size();
	}


	/// Fill all values of data vector with given value.
	void fill(double value)
	{
		ASSERT_PTR(data_ptr_).error("Uninitialized data vector.\n");
		std::fill(data_ptr_->begin(), data_ptr_->end(), value);
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
 * Method creates FieldFE of existing VectorSeqDouble that represents elementwise field.
 *
 * It's necessary to create new VectorSeqDouble of FieldFE, because DOF handler has generally
 * a different ordering than mesh.
 * Then is need to call fill_output_data method.
 *
 * Temporary solution that will be remove after solving issue 995.
 */
template <int spacedim, class Value>
std::shared_ptr<FieldFE<spacedim, Value> > create_field(VectorSeqDouble & vec_seq, Mesh & mesh, unsigned int n_comp)
{
	static MappingP1<1,3> map1;
	static MappingP1<2,3> map2;
	static MappingP1<3,3> map3;

	std::shared_ptr<DOFHandlerMultiDim> dh; // DOF handler object allow create FieldFE
	FiniteElement<0> *fe0; // Finite element objects (allow to create DOF handler)
	FiniteElement<1> *fe1;
	FiniteElement<2> *fe2;
	FiniteElement<3> *fe3;

	switch (n_comp) { // prepare FEM objects for DOF handler by number of components
		case 1: { // scalar
			fe0 = new FE_P_disc<0>(0);
			fe1 = new FE_P_disc<1>(0);
			fe2 = new FE_P_disc<2>(0);
			fe3 = new FE_P_disc<3>(0);
			break;
		}
		case 3: { // vector
			std::shared_ptr< FiniteElement<0> > fe0_ptr = std::make_shared< FE_P_disc<0> >(0);
			std::shared_ptr< FiniteElement<1> > fe1_ptr = std::make_shared< FE_P_disc<1> >(0);
			std::shared_ptr< FiniteElement<2> > fe2_ptr = std::make_shared< FE_P_disc<2> >(0);
			std::shared_ptr< FiniteElement<3> > fe3_ptr = std::make_shared< FE_P_disc<3> >(0);
			fe0 = new FESystem<0>(fe0_ptr, FEType::FEVector, 3);
			fe1 = new FESystem<1>(fe1_ptr, FEType::FEVector, 3);
			fe2 = new FESystem<2>(fe2_ptr, FEType::FEVector, 3);
			fe3 = new FESystem<3>(fe3_ptr, FEType::FEVector, 3);
			break;
		}
		case 9: { // tensor
			std::shared_ptr< FiniteElement<0> > fe0_ptr = std::make_shared< FE_P_disc<0> >(0);
			std::shared_ptr< FiniteElement<1> > fe1_ptr = std::make_shared< FE_P_disc<1> >(0);
			std::shared_ptr< FiniteElement<2> > fe2_ptr = std::make_shared< FE_P_disc<2> >(0);
			std::shared_ptr< FiniteElement<3> > fe3_ptr = std::make_shared< FE_P_disc<3> >(0);
			fe0 = new FESystem<0>(fe0_ptr, FEType::FETensor, 9);
			fe1 = new FESystem<1>(fe1_ptr, FEType::FETensor, 9);
			fe2 = new FESystem<2>(fe2_ptr, FEType::FETensor, 9);
			fe3 = new FESystem<3>(fe3_ptr, FEType::FETensor, 9);
			break;
		}
		default:
			ASSERT(false).error("Should not happen!\n");
	}

	// Prepare DOF handler
	dh = std::make_shared<DOFHandlerMultiDim>(mesh);
	std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( &mesh, fe0, fe1, fe2, fe3);
	dh->distribute_dofs(ds, true);

	// Prepare new empty VectorSeqDouble with corresponding size
	VectorSeqDouble *data_vec = new VectorSeqDouble();
	data_vec->resize( vec_seq.size() );

	// Construct FieldFE
	std::shared_ptr< FieldFE<spacedim, Value> > field_ptr = std::make_shared< FieldFE<spacedim, Value> >();
	field_ptr->set_fe_data(dh, &map1, &map2, &map3, data_vec);
	return field_ptr;
}

/**
 * Fill data to VecSeqDouble in order corresponding with element DOFs.
 *
 * Set data to data vector of field in correct order according to values of DOF handler indices.
 *
 * Temporary solution that will be remove after solving issue 995.
 */
template <int spacedim, class Value>
void fill_output_data(VectorSeqDouble & vec_seq, std::shared_ptr<FieldFE<spacedim, Value> > field_ptr)
{
	auto dh = field_ptr->get_dofhandler();
	unsigned int ndofs = dh->max_elem_dofs();
	unsigned int idof; // iterate over indices
	std::vector<LongIdx> indices(ndofs);

	// Fill DOF handler of FieldFE with correct permutation of data corresponding with DOFs.
	for (auto ele : dh->mesh()->elements_range()) {
		dh->get_dof_indices(ele, indices);
		for(idof=0; idof<ndofs; idof++) (*field_ptr->get_data_vec())[ indices[idof] ] = (*vec_seq.get_data_ptr())[ ndofs*ele.idx()+idof ];
	}
}


/**
 * Like VectorSeqDouble but for MPI PETSC vectors. Have acces to local part.
 */
class VectorMPI {
public:
    typedef typename std::vector<double> VectorData;
    typedef typename std::shared_ptr< VectorData > VectorDataPtr;

    VectorMPI(MPI_Comm comm = PETSC_COMM_WORLD)
    : communicator_(comm) {}

    /// Create shared pointer and PETSC vector with given size. COLLECTIVE.
    VectorMPI(unsigned int local_size, MPI_Comm comm = PETSC_COMM_WORLD)
    : communicator_(comm) {
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
        if (communicator_ == PETSC_COMM_SELF)
        	chkerr(VecCreateSeqWithArray(PETSC_COMM_SELF, 1, local_size, &((*data_ptr_)[0]), &data_petsc_));
        else
        	chkerr(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, local_size, PETSC_DECIDE, &((*data_ptr_)[0]), &data_petsc_));
        chkerr(VecZeroEntries( data_petsc_ ));
    }

    /// Return new vector with same parallel structure.
    void duplicate(VectorMPI other) {
    	ASSERT_EQ(this->communicator_, other.communicator_);
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
    	ASSERT_EQ(this->communicator_, other.communicator_);
        ASSERT_EQ(this->data_ptr_->size(), other.data_ptr_->size());
        uint size = this->data_ptr_->size();
        std::swap(this->data_ptr_, other.data_ptr_);
        chkerr(VecDestroy(&data_petsc_));
        chkerr(VecDestroy(&other.data_petsc_));
        chkerr(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, size, PETSC_DECIDE, &((*data_ptr_)[0]), &data_petsc_));
        chkerr(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, size, PETSC_DECIDE, &((*other.data_ptr_)[0]), &other.data_petsc_));
    }


    void copy(VectorMPI &other) {
    	ASSERT_EQ(this->communicator_, other.communicator_);
        ASSERT_EQ(this->data_ptr_->size(), other.data_ptr_->size());
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
    /// communicator
    MPI_Comm communicator_;
};



#endif /* VECTOR_SEQ_DOUBLE_HH_ */
