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
 * $Id: la_linsys.hh 1299 2011-08-23 21:42:50Z jakub.sistek $
 * $Revision: 1299 $
 * $LastChangedBy: jakub.sistek $
 * $LastChangedDate: 2011-08-23 23:42:50 +0200 (Tue, 23 Aug 2011) $
 *
 * @file
 * @brief   Solver based on Multilevel BDDC - using corresponding class of OpenFTL package
 * @author  Jakub Sistek
 *
 *
 */

#include <mpi.h>

// need OpenFTL BDDC wrapper
#include "la/SystemSolveBddc.hpp"

#include "la/linsys_BDDC.hh"


LinSys_BDDC::LinSys_BDDC( const unsigned lsize,
                          const unsigned numDofsSub,
                          const MPI_Comm comm, 
                          const int matrixTypeInt,
                          const int  numSubLoc )
        : LinSys( lsize, comm )
{
    // set type
    type = LinSys::BDDC;

    la::SystemSolveBddc::MatrixType matrixType;
    switch ( matrixTypeInt ) {
        case 0:
            matrixType = la::SystemSolveBddc::GENERAL;
            break;
        case 1:
            matrixType = la::SystemSolveBddc::SPD;
            break;
        case 2:
            matrixType = la::SystemSolveBddc::SYMMETRICGENERAL;
            break;
        case 3:
            matrixType = la::SystemSolveBddc::SPD_VIA_SYMMETRICGENERAL;
            break;
        default:
            ASSERT( true, "Unknown matrix type %d", matrixTypeInt );
    }

    solver_ = new Solver_( size_,
                           numDofsSub,
                           matrixType,
                           comm, 
                           numSubLoc );

    // set type
    type = LinSys::BDDC;
}

void LinSys_BDDC::load_mesh( const int nDim, const int numNodes, const int numDofs,
                             const std::vector<int> & inet, 
                             const std::vector<int> & nnet, 
                             const std::vector<int> & nndf, 
                             const std::vector<int> & isegn, 
                             const std::vector<int> & isngn, 
                             const std::vector<int> & isvgvn,
                             const std::vector<double> & xyz,
                             const int meshDim ) 
{
    // simply pass the data to BDDCML solver
    isngn_.resize(isngn.size());
    std::copy( isngn.begin(), isngn.end(), isngn_.begin() );
    ASSERT( numDofs == size_, "Global problem size mismatch!" );

    solver_ -> loadRawMesh( nDim, numNodes, numDofs, inet, nnet, nndf, isegn, isngn, isvgvn, xyz, meshDim );
}

void LinSys_BDDC::mat_set_values( int nrow, int *rows, int ncol, int *cols, double *vals )
{
    namespace ublas = boost::numeric::ublas;

    std::vector< unsigned >  myRows( nrow ); 
    std::vector< unsigned >  myCols( ncol ); 
    ublas::matrix< double >  mat( nrow, ncol ); 

    std::copy( &(rows[0]), &(rows[nrow]), myRows.begin() );
    std::copy( &(cols[0]), &(cols[ncol]), myCols.begin() );

    for ( unsigned i = 0; i < nrow; i++ ) {
        for ( unsigned j = 0; j < ncol; j++ ) {
            mat( i, j ) = vals[i*ncol + j];
        }
    }

    solver_ -> insertToMatrix( mat, myRows, myCols );
} 

void LinSys_BDDC::rhs_set_values( int nrow, int *rows, double *vals)
{
    namespace ublas = boost::numeric::ublas;

    std::vector< unsigned >  myRows( nrow ); 
    ublas::vector< double >  vec( nrow ); 

    std::copy( &(rows[0]), &(rows[nrow]), myRows.begin() );

    for ( unsigned i = 0; i < nrow; i++ ) {
        vec( i ) = vals[i];
    }

    solver_ -> insertToRhs( vec, myRows );
}

void LinSys_BDDC::finish_assembly( )
{
    solver_ -> finishMatAssembly( );
}

void LinSys_BDDC::apply_constrains( double scalar )
{
    solver_ -> applyConstraints( constraints_, 1., scalar );
}

int LinSys_BDDC::solve( )
{
    double              tol            = 1.e-7; //!< tolerance on relative residual ||res||/||rhs||
    int                 numLevels      = 2;     //!< number of levels
    std::vector<int> *  numSubAtLevels = NULL;  //!< number of subdomains at levels
    int                 verboseLevel   = 1;     //!< level of verbosity of BDDCML library 
                                                //!< ( 0 - only fatal errors reported, 
                                                //!<   1 - mild output, 
                                                //!<   2 - detailed output )
    int                 maxIt          = 1000;  //!< maximum number of iterations
    int                 ndecrMax       = 30;    //!< maximum number of iterations with non-decreasing residual 
                                                //!< ( used to stop diverging process )
    bool                use_adaptive   = false; //!< should adaptive BDDC be used?

    solver_ -> solveSystem( tol, numLevels, numSubAtLevels, verboseLevel, maxIt, ndecrMax, use_adaptive );

    DBGMSG("BDDCML converged reason: %d ( 0 means OK ) \n", solver_ -> giveConvergedReason() );
    DBGMSG("BDDCML converged in %d iterations. \n", solver_ -> giveNumIterations() );
    DBGMSG("BDDCML estimated condition number is %f \n", solver_ -> giveCondNumber() );

    return solver_ -> giveConvergedReason();
}

void LinSys_BDDC::get_whole_solution( std::vector<double> & globalSolution )
{
    this -> gatherSolution_( );
    globalSolution.resize( globalSolution_.size( ) );
    std::copy( globalSolution_.begin(), globalSolution_.end(), globalSolution.begin() );
}

void LinSys_BDDC::set_whole_solution( std::vector<double> & globalSolution )
{
    globalSolution_.resize( globalSolution.size( ) );
    std::copy( globalSolution.begin(), globalSolution.end(), globalSolution_.begin() );
}

LinSys_BDDC::~LinSys_BDDC()
{ 
    isngn_.clear(); 
    delete solver_; 
};

// construct global solution
void LinSys_BDDC::gatherSolution_( )
{
    // download local solution
    std::vector<double> locSolution( isngn_.size() );
    solver_ -> giveSolution( isngn_, locSolution ); 
    int ierr;

    // merge solution on root
    int rank;
    MPI_Comm_rank( comm_, &rank );
    int nProc;
    MPI_Comm_size( comm_, &nProc );

    globalSolution_.resize( size_ );
    if ( rank == 0 ) {
        // merge my own data
        for ( int i = 0; i < isngn_.size(); i++ ) {
            int ind = isngn_[i];
            globalSolution_[ind] = locSolution[i];
        }
        for ( int iProc = 1; iProc < nProc; iProc++ ) {
            // receive length
            int length;
            MPI_Status status;
            ierr = MPI_Recv( &length, 1, MPI_INT, iProc, iProc, comm_, &status ); 

            // receive indices
            std::vector<int> inds( length );
            ierr = MPI_Recv( &(inds[0]), length, MPI_INT, iProc, iProc, comm_, &status ); 

            // receive local solution
            locSolution.resize( length );
            ierr = MPI_Recv( &(locSolution[0]), length, MPI_DOUBLE, iProc, iProc, comm_, &status ); 

            // merge data
            for ( int i = 0; i < length; i++ ) {
                int ind = inds[i];
                globalSolution_[ind] = locSolution[i];
            }
        }
    }
    else {
        // send my solution to root
        int length = isngn_.size();
        ierr = MPI_Send( &length,                1, MPI_INT,    0, rank, comm_ ); 
        ierr = MPI_Send( &(isngn_[0]),      length, MPI_INT,    0, rank, comm_ ); 
        ierr = MPI_Send( &(locSolution[0]), length, MPI_DOUBLE, 0, rank, comm_ ); 
    }
    // broadcast global solution from root
    ierr = MPI_Bcast( &(globalSolution_[0]), globalSolution_.size(), MPI_DOUBLE, 0, comm_ );
}


