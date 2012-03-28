/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @ingroup system
 * @brief  Objects and functions for mesh partitioning.
 *
 */

#include "system/system.hh"
#include "distribution.hh"

/**
 * create a Distribution from local sizes (dim = np )
 * (collective context)
 */
Distribution::Distribution(const unsigned int size)
:communicator(PETSC_COMM_WORLD),
lsizes(NULL)
{
    F_ENTRY;

    int ierr;
    ierr=MPI_Comm_rank(communicator, &(my_proc));
    ASSERT( ! ierr , "Can not get MPI rank.\n" );
    ierr=MPI_Comm_size(communicator, &(num_of_procs));
    ASSERT( ! ierr  , "Can not get MPI size.\n" );
// TODO: zavest odchytavani vyjimek a pouzivat new a delete

    // communicate global sizes array
    starts=(int *)xmalloc((np()+1)*sizeof(int));
    int lsize=size; // since size is const
    MPI_Allgather(&lsize,1,MPI_INT,starts+1,1,MPI_INT,communicator);
    // count starts
    starts[0]=0;
    for( int i=1 ; i<=np(); i++) starts[i]+=starts[i-1];
}

/**
 * create a Distribution from local sizes (dim = np )
 * (local context)
 */
Distribution::Distribution(const unsigned int * const sizes)
:communicator(PETSC_COMM_WORLD),
 lsizes(NULL)
{
    F_ENTRY;

    int ierr;
    ierr=MPI_Comm_rank(communicator, &(my_proc));
    ASSERT( ! ierr , "Can not get MPI rank.\n" );
    ierr=MPI_Comm_size(communicator, &(num_of_procs));
    ASSERT( ! ierr  , "Can not get MPI size.\n" );
// TODO: zavest odchytavani vyjimek a pouzivat new a delete
    starts=(int *)xmalloc((np()+1)*sizeof(int));
    starts[0]=0;
    for( int i=0 ; i<np(); i++) starts[i+1]=starts[i]+sizes[i];
}

/**
 * constructor from existing PETSC vector
 * (collective context)
 */
Distribution::Distribution(const Vec &petsc_vector)
:communicator(PETSC_COMM_WORLD),
 lsizes(NULL)
{
    F_ENTRY;

    int ierr;
    ierr=MPI_Comm_rank(communicator, &(my_proc));
    ASSERT( ! ierr , "Can not get MPI rank.\n" );
    ierr=MPI_Comm_size(communicator, &(num_of_procs));
    ASSERT( ! ierr  , "Can not get MPI size.\n" );

    const PetscInt *petsc_starts;
    VecGetOwnershipRanges(petsc_vector,&petsc_starts);
    ASSERT( ! ierr , "Can not get vector ownership range.\n" );

    starts=(int *)xmalloc((np()+1)*sizeof(int));
    for( int i=0 ; i<=np(); i++) starts[i]=petsc_starts[i];
}

/**
 * construct from given global size
 * (collective context)
 */
Distribution::Distribution(const SpecialDistribution type, unsigned int global_size)
:communicator(PETSC_COMM_WORLD),
 lsizes(NULL)
{
    F_ENTRY;
    
    int ierr;
    ierr=MPI_Comm_rank(communicator, &(my_proc));
    ASSERT( ! ierr , "Can not get MPI rank.\n" );
    ierr=MPI_Comm_size(communicator, &(num_of_procs));
    ASSERT( ! ierr  , "Can not get MPI size.\n" );

    if (type == Block) {
        int reminder, per_proc;

        reminder=global_size % np(); per_proc=global_size / np();
        // set perproc rows to each proc, but for first "reminder" procs set one row more
        starts=(int *)xmalloc((np()+1)*sizeof(int));
        starts[0]=0;
        for(int i=0; i<np(); i++)
            starts[i+1]=starts[i]+per_proc+(i<reminder?1:0);

    } else if (type == Localized) {

        starts=(int *)xmalloc((np()+1)*sizeof(int));
        starts[0]=0;
        for(int i=1; i<=np(); i++) starts[i]=global_size;
    }
    else {
        ASSERT( 0 , "Cyclic distribution is not yet implemented.\n");
    }
 }
/**
 * copy constructor
 *
 */
Distribution::Distribution(const Distribution &distr)
: communicator(distr.communicator)
{
    DBGMSG("coping distribution\n");
    num_of_procs=distr.num_of_procs;
    my_proc=distr.my_proc;
    starts=(int *)xmalloc((np()+1)*sizeof(int));
    memcpy(starts,distr.starts,(np()+1) * sizeof(int));
    lsizes=NULL;
}


/**
 * find the proc to which belongs index "idx" in the distribution
 * use simple linear search, better binary search could be implemented
 * (local context)
 */
int Distribution::get_proc(int idx) const
{
    ASSERT(NONULL(starts),"Distribution is not initialized.\n");
    ASSERT(idx < size(), "Index %d greater then distribution size %d.\n", idx, size());

    for(int i=0; i<np(); i++) {
        if (is_on_proc(idx,i)) return (i);
    }
    ASSERT( 0 , "Can not find owner of index %d. \n", idx);
    return (-1);
}

const int * Distribution::get_lsizes_array()
{
    if ( lsizes == NULL ) {
        lsizes=(int *) xmalloc(np()*sizeof(int));
        for(int i=0;i<np();i++) lsizes[i]=lsize(i);
    }

    return lsizes;
}

void Distribution::view()
{
    xprintf(Msg,"Distribution:\n");
    xprintf(Msg,"size: %d lsize: %d n. proc: %d\n", size(), lsize(), np());
    for(int i=0; i<np();++i)
        xprintf(Msg,"proc: %d start: %d: lsize: %d\n",i,begin(i),lsize(i));
}

/**
 * Destructor.
 */
Distribution::~Distribution()
{
// TODO: zavest odchytavani vyjimek a pouzivat new a delete
    xfree(starts);
    if (lsizes) xfree(lsizes);
}

