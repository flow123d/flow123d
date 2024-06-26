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
 * @file    distribution.cc
 * @ingroup system
 * @brief   Objects and functions for mesh partitioning.
 */

#include "system/global_defs.h"
#include "system/system.hh"
#include "distribution.hh"

/**
 * create a Distribution from local sizes (dim = np )
 * (collective context)
 */
Distribution::Distribution(const unsigned int size, MPI_Comm comm)
:communicator(comm),
lsizes(NULL)
{
    chkerr(MPI_Comm_rank(communicator, &(my_proc)));
    chkerr(MPI_Comm_size(communicator, &(num_of_procs)));
// TODO: zavest odchytavani vyjimek a pouzivat new a delete

    // communicate global sizes array
    starts= new unsigned int [np()+1];
    unsigned int lsize=size; // since size is const
    MPI_Allgather(&lsize,1,MPI_INT,starts+1,1,MPI_INT,communicator);
    // count starts
    starts[0]=0;
    for( unsigned int i=1 ; i<=np(); i++) starts[i]+=starts[i-1];
}

/**
 * create a Distribution from local sizes (dim = np )
 * (local context)
 */
Distribution::Distribution(const unsigned int * const sizes, MPI_Comm comm)
:communicator(comm),
 lsizes(NULL)
{
    chkerr(MPI_Comm_rank(communicator, &(my_proc)));
    chkerr(MPI_Comm_size(communicator, &(num_of_procs)));
// TODO: zavest odchytavani vyjimek a pouzivat new a delete
    starts= new unsigned int [np()+1];
    starts[0]=0;
    for(unsigned  int i=0 ; i<np(); i++) starts[i+1]=starts[i]+sizes[i];
}

/**
 * constructor from existing PETSC vector
 * (collective context)
 */
Distribution::Distribution(const Vec &petsc_vector)
:communicator(PETSC_COMM_WORLD),
 lsizes(NULL)
{
    chkerr(MPI_Comm_rank(communicator, &(my_proc)));
    chkerr(MPI_Comm_size(communicator, &(num_of_procs)));

    const PetscInt *petsc_starts;
    chkerr(VecGetOwnershipRanges(petsc_vector,&petsc_starts));

    starts= new unsigned int [np()+1];
    for(unsigned  int i=0 ; i<=np(); i++) starts[i]=petsc_starts[i];
}

/**
 * construct from given global size
 * (collective context)
 */
Distribution::Distribution(const DistributionType &type, unsigned int global_size, MPI_Comm comm)
:communicator(comm),
 lsizes(NULL)
{
    chkerr(MPI_Comm_rank(communicator, &(my_proc)));
    chkerr(MPI_Comm_size(communicator, &(num_of_procs)));
    ASSERT_PERMANENT_GT( num_of_procs, 0).error("MPI size is not positive, possibly broken MPI communicator.\n");

    if (type.type_ == Block) {
        unsigned int reminder, per_proc;

        reminder=global_size % np(); per_proc=global_size / np();
        // set perproc rows to each proc, but for first "reminder" procs set one row more
        starts= new unsigned int [np()+1];
        starts[0]=0;
        for(unsigned int i=0; i<np(); i++)
            starts[i+1]=starts[i]+per_proc+(i<reminder?1:0);

    } else if (type.type_ == Localized) {

        starts= new unsigned int [np()+1];
        starts[0]=0;
        for(unsigned int i=1; i<=np(); i++) starts[i]=global_size;
    }
    else {
    	ASSERT_PERMANENT(0).error("Cyclic distribution is not yet implemented.\n");
    }
 }
/**
 * copy constructor
 *
 */
Distribution::Distribution(const Distribution &distr)
: communicator(distr.communicator)
{
    num_of_procs=distr.num_of_procs;
    my_proc=distr.my_proc;
    starts= new unsigned int [np()+1];
    memcpy(starts,distr.starts,(np()+1) * sizeof(unsigned int));
    lsizes=NULL;
}


/**
 * find the proc to which belongs index "idx" in the distribution
 * use simple linear search, better binary search could be implemented
 * (local context)
 */
unsigned int Distribution::get_proc(unsigned  int idx) const
{
	ASSERT_PTR( starts ).error("Distribution is not initialized.\n");
	ASSERT_LT(idx, size()).error("Index is greater than distribution size.\n");

    for(unsigned int i=0; i<np(); i++) {
        if (is_on_proc(idx,i)) return (i);
    }
    ASSERT_PERMANENT(0)(idx).error("Can not find owner of given index.\n");
    return (-1);
}

const unsigned int * Distribution::get_lsizes_array()
{
    if ( lsizes == NULL ) {
        lsizes= new unsigned int [np()];
        for(unsigned int i=0;i<np();i++) lsizes[i]=lsize(i);
    }

    return lsizes;
}



const unsigned int * Distribution::get_starts_array() const {
    return starts;
}



void Distribution::view(std::ostream &stream) const
{
    stream << "[" <<myp() << "]" << "Distribution size: " << size() << " lsize: " << lsize() << " offset: " << begin() << " mpi_size: " << np() << endl;
    for(unsigned int i=0; i<np();++i)
        stream << "[" <<myp() << "]" << "proc: " << i << " offset: " << begin(i) << " lsize: " << lsize(i) << endl;
}

/**
 * Destructor.
 */
Distribution::~Distribution()
{
// TODO: zavest odchytavani vyjimek a pouzivat new a delete
    delete [] starts;
    if (lsizes) delete [] lsizes;
}

