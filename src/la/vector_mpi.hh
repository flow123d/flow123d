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
 * @file    vector_mpi.hh
 * @brief
 */

#ifndef VECTOR_MPI_HH_
#define VECTOR_MPI_HH_

#include <vector>
#include <memory>
#include "system/system.hh"
#include "system/global_defs.h"
#include "mesh/long_idx.hh"

#include <petscvec.h>

/**
 * Auxiliary class for output elementwise concentration vectors
 * in convection transport, sorptions, dual porosity etc.
 *
 * Stores data in two formats:
 *  - shared pointer to std::vector of double
 *  - pointer to PETSC vector that use same data
 *
 * Allows the following functionalities:
 *  - access to local part
 *  - return shared pointer to std::vector of double
 *  - return pointer to PETSC vector
 */
class VectorMPI {
public:
    typedef typename std::vector<double> VectorData;
    typedef typename std::shared_ptr< VectorData > VectorDataPtr;

    VectorMPI(MPI_Comm comm = PETSC_COMM_SELF)
    : communicator_(comm) {}

    /// Create shared pointer and PETSC vector with given size. COLLECTIVE.
    VectorMPI(unsigned int local_size, MPI_Comm comm = PETSC_COMM_WORLD)
    : communicator_(comm) {
        resize(local_size);
    }
    
    /// Create PETSc vector with ghost values whose indices are specified in @p ghost_idx.
    VectorMPI(unsigned int local_size, std::vector<LongIdx> &ghost_idx)
    : communicator_(PETSC_COMM_WORLD) {
        resize(local_size, ghost_idx);
    }

    /**
     * Helper method creating VectorMPI of given size with serial Petsc communicator.
     *
     * Method is used for better readability of code.
     */
    static VectorMPI sequential(unsigned int size)
    {
    	return VectorMPI(size, PETSC_COMM_SELF);
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
    
    /**
     * Resize the vector to given local size with ghost values. Indices of ghost values are in ghost_idx.
     */
    void resize(unsigned int local_size, std::vector<LongIdx> &ghost_idx) {
        ASSERT_DBG(communicator_ == PETSC_COMM_WORLD).error("Cannot allocate ghost values in sequential vector.");
        if (data_ptr_.use_count() ==0) {
            data_ptr_ = std::make_shared< std::vector<double> >(local_size + ghost_idx.size());
        } else {
            ASSERT_DBG( data_ptr_.use_count() ==  1 ) ( data_ptr_.use_count() ).error("Object referenced by other pointer. Can not resize.");
            chkerr(VecDestroy(&data_petsc_));
            data_ptr_->resize(local_size + ghost_idx.size());
        }
        chkerr(VecCreateGhostWithArray(PETSC_COMM_WORLD, local_size, PETSC_DECIDE, ghost_idx.size(), ghost_idx.data(), data_ptr_->data(), &data_petsc_));
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
    
    const VectorData &data() const
    {
        ASSERT_DBG(data_ptr_);
        return *data_ptr_;
    }

    void swap(VectorMPI &other) {
    	ASSERT_EQ(this->communicator_, other.communicator_);
        ASSERT_EQ(this->data_ptr_->size(), other.data_ptr_->size());
        uint size = this->data_ptr_->size();
        std::swap(this->data_ptr_, other.data_ptr_);
        chkerr(VecDestroy(&data_petsc_));
        chkerr(VecDestroy(&other.data_petsc_));
        if (communicator_ == PETSC_COMM_SELF) {
        	chkerr(VecCreateSeqWithArray(PETSC_COMM_SELF, 1, size, &((*data_ptr_)[0]), &data_petsc_));
        	chkerr(VecCreateSeqWithArray(PETSC_COMM_SELF, 1, size, &((*other.data_ptr_)[0]), &other.data_petsc_));
        } else {
            chkerr(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, size, PETSC_DECIDE, &((*data_ptr_)[0]), &data_petsc_));
            chkerr(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, size, PETSC_DECIDE, &((*other.data_ptr_)[0]), &other.data_petsc_));
        }
    }


    void copy(VectorMPI &other) {
    	ASSERT_EQ(this->communicator_, other.communicator_);
        ASSERT_EQ(this->data_ptr_->size(), other.data_ptr_->size());
        chkerr(VecCopy(other.data_petsc_, data_petsc_));
    }
    

    /// local_to_ghost_{begin,end} updates the ghost values on neighbouring processors from local values
    void local_to_ghost_begin() {
        VecGhostUpdateBegin(data_petsc_, INSERT_VALUES, SCATTER_FORWARD);
    }
    
    /// local_to_ghost_{begin,end} updates the ghost values on neighbouring processors from local values
    void local_to_ghost_end() {
        VecGhostUpdateEnd(data_petsc_, INSERT_VALUES, SCATTER_FORWARD);
    }
    
    /// ghost_to_local_{begin,end} updates the local values by adding ghost values from neighbouring processors
    void ghost_to_local_begin() {
        VecGhostUpdateBegin(data_petsc_, ADD_VALUES, SCATTER_REVERSE);
    }
    
    /// ghost_to_local_{begin,end} updates the local values by adding ghost values from neighbouring processors
    void ghost_to_local_end() {
        VecGhostUpdateEnd(data_petsc_, ADD_VALUES, SCATTER_REVERSE);
    }

	/// Return size of output data.
	unsigned int size() const
	{
		ASSERT_PTR(data_ptr_).error("Uninitialized data vector.\n");
		return data_ptr_->size();
	}


    /// Destructor.
    ~VectorMPI()
    {
        if (data_ptr_.use_count() == 1)
            if (data_petsc_) chkerr(VecDestroy(&data_petsc_));
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
    
    /**
     * Access to the vector element on index @p idx (const version).
     */
    inline double &operator[](unsigned int idx) const
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


#endif /* VECTOR_MPI_HH_ */
