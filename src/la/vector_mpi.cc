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
 * @file    vector_mpi.cc
 * @brief
 */

#include "la/vector_mpi.hh"
#include <vector>
#include <memory>
#include "system/system.hh"
#include "system/global_defs.h"
#include "system/index_types.hh"

#include <petscvec.h>
#include <armadillo>


VectorMPI::VectorMPI(MPI_Comm comm)
: communicator_(comm)
{}


VectorMPI::VectorMPI(unsigned int local_size, MPI_Comm comm)
: communicator_(comm) {
    resize(local_size);
}

VectorMPI::VectorMPI(unsigned int local_size, std::vector<LongIdx> &ghost_idx)
: communicator_(PETSC_COMM_WORLD) {
    resize(local_size, ghost_idx);
}

VectorMPI VectorMPI::sequential(unsigned int size)
{
    return VectorMPI(size, PETSC_COMM_SELF);
}


void VectorMPI::resize(unsigned int local_size) {
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


void VectorMPI::resize(unsigned int local_size, std::vector<LongIdx> &ghost_idx) {
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


void VectorMPI::duplicate_from(VectorMPI other) {
    ASSERT_EQ(this->communicator_, other.communicator_);
    this->resize(other.data().size());
}


void VectorMPI::swap(VectorMPI &other) {
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


void VectorMPI::copy_from(VectorMPI &other) {
    ASSERT_EQ(this->communicator_, other.communicator_);
    ASSERT_EQ(this->data_ptr_->size(), other.data_ptr_->size());
    chkerr(VecCopy(other.data_petsc_, data_petsc_));
}


arma::vec VectorMPI::get_subvec(const LocDofVec& loc_indices)
{
    ASSERT_DBG(data_ptr_);
    arma::vec vec(loc_indices.n_elem);
    
    for(unsigned int i=0; i<loc_indices.n_elem; i++){
        unsigned int idx = loc_indices(i);
        ASSERT_DBG(idx < data_ptr_->size()) (idx) (data_ptr_->size());
        vec(i) = (*data_ptr_)[idx];
    }
    return vec;
}


arma::vec VectorMPI::get_subvec(const LocDofVec& loc_indices) const
{
    ASSERT_DBG(data_ptr_);
    arma::vec vec(loc_indices.n_elem);
    
    for(unsigned int i=0; i<loc_indices.n_elem; i++){
        unsigned int idx = loc_indices(i);
        ASSERT_DBG(idx < data_ptr_->size()) (idx) (data_ptr_->size());
        vec(i) = (*data_ptr_)[idx];
    }
    return vec;
}


void VectorMPI::set_subvec(const LocDofVec& loc_indices, const arma::vec& values)
{
    ASSERT_DBG(data_ptr_);
    ASSERT_EQ_DBG(loc_indices.n_elem, values.n_elem);
    
    for(unsigned int i=0; i<loc_indices.n_elem; i++){
        unsigned int idx = loc_indices(i);
        ASSERT_DBG(idx < data_ptr_->size()) (idx) (data_ptr_->size());
        (*data_ptr_)[idx] = values(i);
    }
}

/// Destructor.
VectorMPI::~VectorMPI()
{
    if (data_ptr_.use_count() == 1)
        if (data_petsc_) chkerr(VecDestroy(&data_petsc_));
}