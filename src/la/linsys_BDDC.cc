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
 * @file    linsys_BDDC.cc
 * @brief   Solver based on Multilevel BDDC - using corresponding class of OpenFTL package
 * @author  Jakub Sistek
 */

#include <mpi.h>
#include "config.h"

// need BDDCML wrapper
#ifdef FLOW123D_HAVE_BDDCML
  #include <map>
  #include "la/bddcml_wrapper.hpp"
#endif // FLOW123D_HAVE_BDDCML

#include "la/linsys.hh"
#include "la/linsys_BDDC.hh"
#include "system/sys_profiler.hh"



namespace it = Input::Type;


const it::Record & LinSys_BDDC::get_input_type() {
	return it::Record("Bddc", "Solver setting.")
		.derive_from(LinSys::get_input_type())
		.declare_key("r_tol", it::Double(0.0, 1.0), it::Default("1.0e-7"),
					"Relative residual tolerance (to initial error).")
		.declare_key("max_it", it::Integer(0), it::Default("10000"),
					"Maximum number of outer iterations of the linear solver.")
		.declare_key("max_nondecr_it", it::Integer(0), it::Default("30"),
					 "Maximum number of iterations of the linear solver with non-decreasing residual.")
		.declare_key("number_of_levels", it::Integer(0), it::Default("2"),
					 "Number of levels in the multilevel method (=2 for the standard BDDC).")
		.declare_key("use_adaptive_bddc", it::Bool(), it::Default("false"),
					 "Use adaptive selection of constraints in BDDCML.")
		.declare_key("bddcml_verbosity_level", it::Integer(0,2), it::Default("0"),
					 "Level of verbosity of the BDDCML library:\n\n - 0 - no output\n - 1 - mild output\n - 2 - detailed output.")
		.close();
}


const int LinSys_BDDC::registrar = LinSys_BDDC::get_input_type().size();



LinSys_BDDC::LinSys_BDDC( const unsigned numDofsSub,
                          const Distribution * rows_ds,
                          const int matrixTypeInt,
                          const int  numSubLoc,
                          const bool swap_sign )
        : LinSys( rows_ds ),
          swap_sign_(swap_sign)
{
#ifdef FLOW123D_HAVE_BDDCML
    // from the point of view of assembly, BDDC linsys is in the ADD state
    status_ = LinSys::ADD;

    la::BddcmlWrapper::MatrixType matrixType;
    switch ( matrixTypeInt ) {
        case 0:
            matrixType = la::BddcmlWrapper::GENERAL;
            break;
        case 1:
            matrixType = la::BddcmlWrapper::SPD;
            break;
        case 2:
            matrixType = la::BddcmlWrapper::SYMMETRICGENERAL;
            break;
        case 3:
            matrixType = la::BddcmlWrapper::SPD_VIA_SYMMETRICGENERAL;
            break;
        default:
            matrixType = la::BddcmlWrapper::GENERAL;
            ASSERT( true, "Unknown matrix type %d", matrixTypeInt );
    }

    bddcml_ = new Bddcml_( size_,
                           numDofsSub,
                           matrixType,
                           rows_ds->get_comm(),
                           numSubLoc );

    // prepare space for local solution
    locSolution_.resize( numDofsSub );

    // identify it with PETSc vector
    PetscErrorCode ierr;
    PetscInt numDofsSubInt = static_cast<PetscInt>( numDofsSub );
    //ierr = VecCreateSeq( PETSC_COMM_SELF, numDofsSubInt, &locSolVec_ ); 
    ierr = VecCreateSeqWithArray( PETSC_COMM_SELF, 1, numDofsSub, &(locSolution_[0]), &locSolVec_ ); 
    CHKERRV( ierr );
#else
    xprintf(UsrErr, "Compiled without support for BDDCML solver.\n");  
#endif // FLOW123D_HAVE_BDDCML
}

void LinSys_BDDC::load_mesh( const int nDim, const int numNodes, const int numDofs,
                             const std::vector<int> & inet, 
                             const std::vector<int> & nnet, 
                             const std::vector<int> & nndf, 
                             const std::vector<int> & isegn, 
                             const std::vector<int> & isngn, 
                             const std::vector<int> & isvgvn,
                             const std::vector<double> & xyz,
                             const std::vector<double> & element_permeability,
                             const int meshDim ) 
{
#ifdef FLOW123D_HAVE_BDDCML
    // simply pass the data to BDDCML solver
    isngn_.resize(isngn.size());
    std::copy( isngn.begin(), isngn.end(), isngn_.begin() );
    ASSERT( numDofs == size_, "Global problem size mismatch!" );

    bddcml_ -> loadRawMesh( nDim, numNodes, numDofs, inet, nnet, nndf, isegn, isngn, isvgvn, xyz, element_permeability, meshDim );

    // create a map for BDDCML to PETSc vector
    PetscErrorCode ierr;

    // local index set
    PetscInt numDofsSubInt = static_cast<int>( isngn_.size( ) );
    std::vector<PetscInt> idx( isngn_ );

    ISLocalToGlobalMapping subdomainMapping;
    ierr = ISLocalToGlobalMappingCreate( comm_, numDofsSubInt, &(idx[0]), PETSC_COPY_VALUES, &subdomainMapping ); CHKERRV( ierr );
    
    IS subdomainIndexSet;
    IS from;
    ierr = ISCreateStride( PETSC_COMM_SELF, numDofsSubInt, 0, 1, &subdomainIndexSet ); 
    ierr = ISLocalToGlobalMappingApplyIS( subdomainMapping, subdomainIndexSet, &from ); 


    // create scatter from parallel PETSc vector to local indices of subdomain
    ierr = VecScatterCreate( solution_, from, locSolVec_, subdomainIndexSet, &VSpetscToSubScatter_ ); CHKERRV( ierr );
    ierr = ISDestroy( &subdomainIndexSet ); CHKERRV( ierr );
    ierr = ISDestroy( &from ); CHKERRV( ierr );

    //double * locSolVecArray;
    //ierr = VecGetArray( locSolVec_, &locSolVecArray ); CHKERRV( ierr );
    //std::copy( locSolution_.begin(), locSolution_.end(), locSolVecArray );
    //ierr = VecRestoreArray( locSolVec_, &locSolVecArray ); CHKERRV( ierr );

    // scatter local solutions back to global one
    VecScatterBegin( VSpetscToSubScatter_, locSolVec_, solution_, INSERT_VALUES, SCATTER_REVERSE ); 
    VecScatterEnd(   VSpetscToSubScatter_, locSolVec_, solution_, INSERT_VALUES, SCATTER_REVERSE );
#endif // FLOW123D_HAVE_BDDCML
}

void LinSys_BDDC::mat_set_values( int nrow, int *rows, int ncol, int *cols, double *vals )
{
#ifdef FLOW123D_HAVE_BDDCML
	namespace ublas = boost::numeric::ublas;

    std::vector< unsigned >  myRows( nrow ); 
    std::vector< unsigned >  myCols( ncol ); 
    ublas::matrix< double >  mat( nrow, ncol ); 

    std::copy( &(rows[0]), &(rows[nrow]), myRows.begin() );
    std::copy( &(cols[0]), &(cols[ncol]), myCols.begin() );

    for ( int i = 0; i < nrow; i++ ) {
        for ( int j = 0; j < ncol; j++ ) {
            mat( i, j ) = vals[i*ncol + j];
        }
    }
    if (swap_sign_) {
       mat = -mat;
    }

    bddcml_ -> insertToMatrix( mat, myRows, myCols );
#endif // FLOW123D_HAVE_BDDCML
} 

void LinSys_BDDC::rhs_set_values( int nrow, int *rows, double *vals)
{
#ifdef FLOW123D_HAVE_BDDCML
    namespace ublas = boost::numeric::ublas;

    std::vector< unsigned >  myRows( nrow ); 
    ublas::vector< double >  vec( nrow ); 

    std::copy( &(rows[0]), &(rows[nrow]), myRows.begin() );

    for ( int i = 0; i < nrow; i++ ) {
        vec( i ) = vals[i];
    }
    if (swap_sign_) {
       vec = -vec;
    }

    bddcml_ -> insertToRhs( vec, myRows );
#endif // FLOW123D_HAVE_BDDCML
}

void LinSys_BDDC::diagonal_weights_set_value( int global_index, double value )
{
#ifdef FLOW123D_HAVE_BDDCML
    bddcml_ -> insertToDiagonalWeights( global_index, value );
#endif // FLOW123D_HAVE_BDDCML
}

PetscErrorCode LinSys_BDDC::mat_zero_entries()
{
#ifdef FLOW123D_HAVE_BDDCML
    bddcml_ -> clearMatrix( );
#endif // FLOW123D_HAVE_BDDCML
    return 0;
}

PetscErrorCode LinSys_BDDC::rhs_zero_entries()
{
#ifdef FLOW123D_HAVE_BDDCML
    bddcml_ -> clearRhs( );
#endif // FLOW123D_HAVE_BDDCML
    return 0;
}

void LinSys_BDDC::finish_assembly( )
{
#ifdef FLOW123D_HAVE_BDDCML
    bddcml_ -> finishMatAssembly( );
#endif // FLOW123D_HAVE_BDDCML
}

void LinSys_BDDC::apply_constrains( double scalar )
{
#ifdef FLOW123D_HAVE_BDDCML
    bddcml_ -> applyConstraints( constraints_, 1., scalar );
#endif // FLOW123D_HAVE_BDDCML
}

int LinSys_BDDC::solve()    // ! params are not currently used
{
#ifdef FLOW123D_HAVE_BDDCML
    std::vector<int> *  numSubAtLevels = NULL;  //!< number of subdomains at levels

    {
		START_TIMER("BDDC linear solver");
		START_TIMER("BDDC linear iteration");
		bddcml_ -> solveSystem( r_tol_, number_of_levels_, numSubAtLevels, bddcml_verbosity_level_, max_it_, max_nondecr_it_, use_adaptive_bddc_ );


		xprintf(MsgLog,"BDDCML converged reason: %d ( 0 means OK ) \n", bddcml_ -> giveConvergedReason() );
		xprintf(MsgLog,"BDDCML converged in %d iterations. \n", bddcml_ -> giveNumIterations() );
		xprintf(MsgLog,"BDDCML estimated condition number is %f \n", bddcml_ -> giveCondNumber() );
		ADD_CALLS(bddcml_ -> giveNumIterations());
    }

    // download local solution
    bddcml_ -> giveSolution( isngn_, locSolution_ ); 


    double * locSolVecArray;
    VecGetArray( locSolVec_, &locSolVecArray ); 
    std::copy( locSolution_.begin(), locSolution_.end(), locSolVecArray );
    VecRestoreArray( locSolVec_, &locSolVecArray ); 

    // scatter local solutions back to global one
    VecScatterBegin( VSpetscToSubScatter_, locSolVec_, solution_, INSERT_VALUES, SCATTER_REVERSE ); 
    VecScatterEnd(   VSpetscToSubScatter_, locSolVec_, solution_, INSERT_VALUES, SCATTER_REVERSE );

    // upper bound on the residual error
    residual_norm_ = r_tol_ * bddcml_->normRhs( ) ;

    return bddcml_ -> giveConvergedReason();
#else
	return 0;
#endif // FLOW123D_HAVE_BDDCML

}

void LinSys_BDDC::set_from_input(const Input::Record in_rec)
{
    // common values
    r_tol_  = in_rec.val<double>("r_tol");   
    max_it_ = in_rec.val<int>("max_it");   

    // BDDCML specific parameters
    max_nondecr_it_         = in_rec.val<int>("max_nondecr_it");   
    number_of_levels_       = in_rec.val<int>("number_of_levels");   
    use_adaptive_bddc_      = in_rec.val<bool>("use_adaptive_bddc");
    bddcml_verbosity_level_ = in_rec.val<int>("bddcml_verbosity_level");
}

LinSys_BDDC::~LinSys_BDDC()
{ 
#ifdef FLOW123D_HAVE_BDDCML
    isngn_.clear();
    locSolution_.clear(); 

    PetscErrorCode ierr;
    ierr = VecDestroy( &locSolVec_ ); CHKERRV( ierr );

    ierr = VecScatterDestroy( &VSpetscToSubScatter_ ); CHKERRV( ierr );

    delete bddcml_;
#endif // FLOW123D_HAVE_BDDCML
};

// construct global solution
//void LinSys_BDDC::gatherSolution_( )
//{
//#ifdef FLOW123D_HAVE_BDDCML
//    int ierr;
//
//    // merge solution on root
//    int rank;
//    MPI_Comm_rank( comm_, &rank );
//    int nProc;
//    MPI_Comm_size( comm_, &nProc );
//
//    globalSolution_.resize( size_ );
//    std::vector<double> locSolutionNeib;
//    if ( rank == 0 ) {
//        // merge my own data
//        for ( unsigned int i = 0; i < isngn_.size(); i++ ) {
//            int ind = isngn_[i];
//            globalSolution_[ind] = locSolution_[i];
//        }
//        for ( int iProc = 1; iProc < nProc; iProc++ ) {
//            // receive length
//            int length;
//            MPI_Status status;
//            ierr = MPI_Recv( &length, 1, MPI_INT, iProc, iProc, comm_, &status ); 
//
//            // receive indices
//            std::vector<int> inds( length );
//            ierr = MPI_Recv( &(inds[0]), length, MPI_INT, iProc, iProc, comm_, &status ); 
//
//            // receive local solution
//            locSolutionNeib.resize( length );
//            ierr = MPI_Recv( &(locSolutionNeib[0]), length, MPI_DOUBLE, iProc, iProc, comm_, &status ); 
//
//            // merge data
//            for ( int i = 0; i < length; i++ ) {
//                int ind = inds[i];
//                globalSolution_[ind] = locSolutionNeib[i];
//            }
//        }
//    }
//    else {
//        // send my solution to root
//        int length = isngn_.size();
//        ierr = MPI_Send( &length,                1, MPI_INT,    0, rank, comm_ ); 
//        ierr = MPI_Send( &(isngn_[0]),      length, MPI_INT,    0, rank, comm_ ); 
//        ierr = MPI_Send( &(locSolution_[0]), length, MPI_DOUBLE, 0, rank, comm_ ); 
//    }
//    // broadcast global solution from root
//    ierr = MPI_Bcast( &(globalSolution_[0]), globalSolution_.size(), MPI_DOUBLE, 0, comm_ );
//#endif // FLOW123D_HAVE_BDDCML
//}

double LinSys_BDDC::get_solution_precision()
{
	double bnorm=0.0;
	VecNorm(locSolVec_, NORM_2, &bnorm);

	return max(a_tol_, r_tol_*bnorm);
}


