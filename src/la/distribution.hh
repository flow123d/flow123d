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
 * @brief  Support classes for parallel programing.
 *
 * TODO:
 *
 * * need better resolution of constructors
 */


#ifndef MESH_PARTITION_HH_
#define MESH_PARTITION_HH_

#include <mpi.h>
#include <ostream>
#include <petscvec.h>

class DistributionType {
public:
explicit DistributionType(int type) : type_(type) {}
int type_;
};

class DistributionBlock : public DistributionType {
public: DistributionBlock() : DistributionType(-1) {}
};

class DistributionLocalized : public DistributionType {
public: DistributionLocalized() : DistributionType(-2) {}
};

/**
class DistributionCyclic : public DistributionType {
};
*/

/**
 * Continuous distribution of an 1D array of indexes.
 * Default is PETSC_COMM_WORLD communicator, but you can provide another.
 */

class Distribution {
public:
    /**
     * Type of distribution automatically created only from global number of indices.
     */
    typedef enum  {
        /** Distribute indices evenly. */
        Block=-1,
        /** Put all on the first processor. */
        Localized=-2
    } SpecialDistribution;

    /**
     * Constructor. It makes distribution from given local number of indices on every processor.
     *
     * COLLECTIVE
     * @param size Local size on calling processor.
     */
    //@param comm (optional) MPI Communicator. Default PETSC_COMM_WORLD.
    Distribution(const unsigned int size, MPI_Comm comm);

    /**
     * Constructor. It makes distribution from given array of sizes of processors.
     *
     * NOT COLLECTIVE, but still use MPI to provide information about processors.
     *
     * @param sizes Int array with sizes.
     */
    //@param comm (optional) MPI Communicator. Default PETSC_COMM_WORLD.
    Distribution(const unsigned int * const sizes, MPI_Comm comm);

    /**
     * Constructor. It makes distribution from distribution of a PETSC vector.
     *
     * NOT COLLECTIVE, but still use MPI to provide information about processors.
     * @param petsc_vector
     */
    //@param comm (optional) MPI Communicator. Default PETSC_COMM_WORLD.
    Distribution(const Vec &petsc_vector);

    /**
     * Constructor. It makes distribution from global number of indices according to given scheme.
     *
     * NOT COLLECTIVE, but still use MPI to provide information about processors.
     * @param type Either Block or Localized distribution of indices.
     * @param global_size Total number of indices to distribute.
     */
    //@param comm (optional) MPI Communicator. Default PETSC_COMM_WORLD.
    Distribution(const DistributionType &type, unsigned int global_size, MPI_Comm comm);

    /**
     * Copy Constructor.
     */
    Distribution(const Distribution &distr);

    /// get num of processors
    inline unsigned int np() const {return num_of_procs;}
    /// get my processor
    inline unsigned int myp() const {return my_proc;}
    /// get starting local index
    inline unsigned int begin(int proc) const {return (starts[proc]);}
    inline unsigned int begin() const {return ( begin(myp()) );}
    /// get last local index +1
    inline unsigned int end(int proc) const {return (starts[proc+1]);}
    inline unsigned int end()  const {return ( end(myp()) );}
    /// get local size
    inline unsigned int lsize(int proc) const {return (end(proc)-begin(proc));}
    inline unsigned int lsize() const {return ( lsize(myp()) );}
    /// get global size
    inline unsigned int size() const {return (starts[np()]);}
    /// identify local index
    inline bool is_local(unsigned int idx) const {return ( begin()<=(idx) && (idx)<end() );}
    inline bool is_on_proc(unsigned int idx, unsigned int proc) const {return ( begin(proc)<=(idx) && (idx)<end(proc) );}
    /// get processor of the given index
    unsigned int get_proc(unsigned int idx) const;
    /// get local sizes array
    const unsigned int * get_lsizes_array();
    /// get local starts array
    const unsigned int * get_starts_array() const;
    /// Returns communicator.
    inline MPI_Comm get_comm() {return communicator;}
    /// distribution view
    void view(std::ostream &stream) const;
    ~Distribution();
private:
    /// communicator
    const MPI_Comm communicator;
    /// number of procs
    int num_of_procs;
    /// my proc number
    int my_proc;
    /// starts[i] index of the first index on the proc i; starts[n_procs]=size of whole array
    unsigned int *starts;
    /// local sizes
    unsigned int *lsizes;
};

inline std::ostream & operator <<(std::ostream &stream, const Distribution &distr)
    { distr.view(stream); return stream; }



#endif // MESH_PARTITION_HH_
