//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2010.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//! @author Jakub Sistek, Thomas Rueberg
//! @date   3/2011

//------------------------------------------------------------------------------
/** Routine for preparing MATRIX assembly. It reserves memory in triplet.
 */
void corlib::SystemSolvePetsc::prepareMatAssembly( unsigned numElements, unsigned elMatSize )
{
    unsigned length = numElements * elMatSize;
    triplet_.prepareAssembly( length );

    return;
}
//------------------------------------------------------------------------------
/** Insert element stiffness matrix into global system matrix
 */
void corlib::SystemSolvePetsc::insertToMatrix( const SubMat_  & subMat, 
                                               const VecUint_ & rowIndices,
                                               const VecUint_ & colIndices )
{
       
    // copy dense matrix into triplet
    for ( unsigned i = 0; i < rowIndices.size(); i++) {
        for ( unsigned j = 0; j < colIndices.size(); j++) {

            // insert value
            triplet_.insert( rowIndices[i], colIndices[j], subMat(i,j) );
        }
    }

    // set flag
    isMatAssembled_ = false;

    return;
}

//------------------------------------------------------------------------------
/** Routine for assembly of the underlying triplet.
 *  This routine can save some memory if the sparse matrix is put into triplet.
 *  By calling this routine, the triplet is 
 *  sorted and assembled, but it still allows to add new entries at the end.
 *  This might be useful especially if multiple sweeps on elements are used.
 */
void corlib::SystemSolvePetsc::partialMatAssembly( )
{
    // do assembly of triplet
    triplet_.finishAssembly( );
}

//------------------------------------------------------------------------------
/** Routine for Petsc's communication.
 *  To be called right after the finish of the MATRIX assembly and before
 *  the RHS assembly start.
 */
void corlib::SystemSolvePetsc::finishMatAssembly( )
{
    // do nothing if already called
    if ( isMatAssembled_ ) return;

    // finish assembly of triplet
    triplet_.finishAssembly( );

    if ( sysMat_ == NULL ) {
        // create MPI matrix
        ierr_ = MatCreate( comm_, &sysMat_ ); CHKERRV( ierr_ ); 
        ierr_ = MatSetType( sysMat_, MATMPIAIJ ); CHKERRV( ierr_ ); 
        if ( numLocalDofs_ >= 0 ) {
            // specify number of local rows if this is given
            ierr_ = MatSetSizes( sysMat_, numLocalDofs_, numLocalDofs_, 
                                 numTotalDofs_, numTotalDofs_ ); 
            CHKERRV( ierr_ ); 
        }
        else {
            // let PETSc decide about number of local matrix rows
            ierr_ = MatSetSizes( sysMat_, PETSC_DECIDE, PETSC_DECIDE, 
                                 numTotalDofs_, numTotalDofs_ ); 
            CHKERRV( ierr_ ); 
        }

        // get index range [matLow,matHigh)
        ierr_ = MatGetOwnershipRange(sysMat_, &matLow_, &matHigh_); CHKERRV( ierr_ ); 

        // let Triplet find data for preallocation
        unsigned locBlock = matHigh_ - matLow_ ; 
        std::vector<PetscInt> nnzRowsIn(  locBlock, 0 );
        std::vector<PetscInt> nnzRowsOut( locBlock, 0 );

        // find number of nonzeros in [ matLow_, matHigh_ ) interval and outside on a row block corresponding
        // to the processor
        triplet_.getIntervalCounts( matLow_ ,  matHigh_ , 
                                    nnzRowsIn, nnzRowsOut );

        // check data locality
        std::vector<PetscInt>::const_iterator iter;
        iter = nnzRowsIn.begin();
        unsigned sum1 = 0;
        for ( ; iter != nnzRowsIn.end(); ++iter ) {
            sum1 += *iter;
        }
        iter = nnzRowsOut.begin();
        unsigned sum2 = 0;
        for ( ; iter != nnzRowsOut.end(); ++iter ) {
            sum2 += *iter;
        }

        if ( sum2 > 0.5*sum1 ) {
            std::cout << "Warning: More than 50 percent of entries are nonlocal to process " 
                      << rank_ << " : "
                      << float(sum2)/float(sum1) * 100 << " percent" << std::endl;
            std::cout.flush();
        }

#ifdef PROFILE
        std::cout << "Number of local/non-local entries on rank " << rank_ << " is "
                  << float(sum1) << " / " << float(sum2) << ",  "
                  << "nonlocal entries in percent: " 
                  << float(sum2) / float(sum1) * 100
                  << std::endl;
        std::cout.flush();
#endif

        // preallocate the MPI matrix
        ierr_ = MatMPIAIJSetPreallocation( sysMat_, 0, &(nnzRowsIn[0]), 0, &(nnzRowsOut[0]) );
        CHKERRV( ierr_ ); 

    }
    else {
        this -> clearMatrix( );
    }

    // how many entries to copy?
    PetscInt row, col;
    PetscReal val;
    unsigned numEntries = triplet_.nnz( );
    // copy entries from triplet to PETSc matrix
    for ( unsigned i = 0; i < numEntries; i++ ) {

        // extract entry from triplet
        triplet_.getEntry( i, row, col, val );

        if ( not filterLocalRows_ ) {
            // insert entry to existing system matrix for possible addition
            ierr_ = MatSetValue( sysMat_, row, col, val, ADD_VALUES ); 
            CHKERRV( ierr_ ); 
        }
        else {
            if ( row >= matLow_  and row <  matHigh_ ) {
                // if filtering is active, insert final entry to existing system matrix - no addition to be performed
                ierr_ = MatSetValue( sysMat_, row, col, val, INSERT_VALUES ); 
                CHKERRV( ierr_ ); 
            }
        }
    }

    // clear triplet after it has been copied to PETSc
    diagScalar_ = triplet_.getDiagScalar();
    triplet_.clear();

    // do the assembly-begin and -end 
#ifdef PROFILE
    double time1 = MPI_Wtime();
#endif
    ierr_ = MatAssemblyBegin( sysMat_, MAT_FINAL_ASSEMBLY ); CHKERRV( ierr_ ); 
    ierr_ = MatAssemblyEnd(   sysMat_, MAT_FINAL_ASSEMBLY ); CHKERRV( ierr_ ); 
#ifdef PROFILE
    double time2 = MPI_Wtime();
#endif

    MatInfo info;
    double  mal, nz_a, nz_u, nz_un;

    MatGetInfo( sysMat_, MAT_GLOBAL_SUM, &info);
    mal   = info.mallocs;
    nz_a  = info.nz_allocated;
    nz_u  = info.nz_used;
    nz_un = info.nz_unneeded;

#ifdef PROFILE
    if ( rank_ == 0 ) {
        std::cout << "Time of assembly in Petsc [s]: " 
                  << time2 - time1
                  << std::endl;
        std::cout.flush();
        std::cout << "Number of mallocs: " << mal << std::endl
                  << "Number of entries allocated: " << nz_a << std::endl
                  << "Number of entries used: " << nz_u << std::endl
                  << "Number of entries unneeded: " << nz_un << std::endl;
        std::cout.flush();
    }
#endif

    isMatAssembled_ = true;

    return;
}
//------------------------------------------------------------------------------
/** Routine for preparing RHS assembly.
 *  If called after the finish of the MATRIX assembly,
 *  structure of RHS is deduced from the assembled matrix.
 */
void corlib::SystemSolvePetsc::prepareRhsAssembly( )
{
    if ( rhsVec_ == NULL ) {
        // create RHS

        // check matrix exists
        if ( sysMat_ != NULL ) {
            // if matrix is created, deduce vector structure from the matrix
            ierr_ = MatGetVecs( sysMat_, PETSC_NULL, &rhsVec_ ); CHKERRV( ierr_ ); 
        }
        else {
            // if matrix is not created, create rhs vector
            if ( numLocalDofs_ >= 0 ) {
                // number of local dofs specified
                ierr_ = VecCreateMPI( comm_, numLocalDofs_, numTotalDofs_, &rhsVec_ ); CHKERRV( ierr_ ); 
            }
            else {
                // number of local dofs obtained by PETSc
                ierr_ = VecCreateMPI( comm_, PETSC_DECIDE,  numTotalDofs_, &rhsVec_ ); CHKERRV( ierr_ ); 
            }
        }
        // get index range [rhsLow,rhsHigh)
        ierr_ = VecGetOwnershipRange(rhsVec_, &rhsLow_, &rhsHigh_); CHKERRV( ierr_ ); 
    }
    else {
        // zero RHS
        this -> clearRhs( );
    }

    return;
}
//------------------------------------------------------------------------------
/** Insert subvector to the system's RHS vector
 */
void corlib::SystemSolvePetsc::insertToRhs( const SubVec_  & subVec, 
                                            const VecUint_ & dofIndices )
{
    unsigned nDofs = dofIndices.size();

    if ( rhsVec_ == NULL ) {
        this -> prepareRhsAssembly();
    }

    // simply loop over dofs
    for ( unsigned iDof = 0; iDof < nDofs; iDof++ ) {
        this -> insertValueToRhs( subVec (iDof), dofIndices[iDof] );
    }

    return;
}
//------------------------------------------------------------------------------
/** Insert a scalar to RHS
 */
void corlib::SystemSolvePetsc::insertValueToRhs( const double & entry, 
                                                 const unsigned & dofIndex )
{
    if ( rhsVec_ == NULL ) {
        this -> prepareRhsAssembly();
    }

    PetscInt    index = static_cast<PetscInt>( dofIndex );
    if ( not filterLocalRows_ or
         ( filterLocalRows_ and 
           index >= rhsLow_ and 
           index <  rhsHigh_ ) )
        {

        PetscScalar value = static_cast<PetscScalar>( entry );

        ierr_ = VecSetValues( rhsVec_, 1, &index, &value, ADD_VALUES );
        CHKERRV( ierr_ );
        // set flag
        isRhsAssembled_ = false;
    }


    return;
}

//------------------------------------------------------------------------------
/** Routine for finishing RHS assembly.
 *  To be called right after the finish of the RHS assembly and before
 *  the constraint application.
 */
void corlib::SystemSolvePetsc::finishRhsAssembly( )
{
    // do nothing if already called
    if ( isRhsAssembled_ ) return;

    // if the RHS vector was not created by inserting values, create and empty vector now
    if ( rhsVec_ == NULL ) {
        this -> prepareRhsAssembly();
    }

    ierr_ = VecAssemblyBegin( rhsVec_ ); CHKERRV( ierr_ ); 
    ierr_ = VecAssemblyEnd(   rhsVec_ ); CHKERRV( ierr_ ); 
    isRhsAssembled_ = true;

    return;
}

//------------------------------------------------------------------------------
/** Routine for finishing assemblies of both matrix and RHS.
 *  This is provided for compatibility. It is generally better to finish matrix assembly 
 *  and then derive the parallel layout of the RHS vector in prepareRhsAssembly, 
 *  i.e. follow this call sequence:
 *  prepareMatAssembly
 *  insertToMatrix
 *  finishMatAssembly
 *  prepareRhsAssembly
 *  insertToRhs
 *  finishRhsAssembly
 */
void corlib::SystemSolvePetsc::finishAssembly( )
{
    this -> finishMatAssembly();
    this -> finishRhsAssembly();

    return;
}

//------------------------------------------------------------------------------
/** Apply constraints to the system by replacing row 'i' with vector 'alpha*e_i',
 *  i.e., a zero vector with value 'alpha' at its i-th position, if the dof 'i'
 *  is prescribed. The corresponding value of the prescribed dof is set as the
 *  'i'-th value in the RHS-vector multiplied by a given factor.
 *  Note that Petsc requires the function 'MatZeroRows' to be called by all
 *  processes with all global indices even if not owned by some processes.
 */
void corlib::SystemSolvePetsc::applyConstraints( ConstraintVec_ & constraints,
                                                 const double factor , 
                                                 double scalar )
{
    // check that system matrix is assembled
    if ( not isMatAssembled_ ) {
        std::cout << " applyConstraints: The matrix has not been assembled and will be assembled now." << std::endl;
        this -> finishMatAssembly();
    }
    // check that RHS is assembled
    if ( not isRhsAssembled_ ) {
        std::cout << " applyConstraints: The right-hand side has not been assembled and will be assembled now." << std::endl;
        this -> finishRhsAssembly();
    }

    // number of constraints
    PetscInt numConstraints = static_cast<PetscInt>( constraints.size() );

    // Additional multiplier for numerical reasons (criterion to be established)
    const PetscScalar diagScalar = scalar;
    const PetscScalar rhsScalar = diagScalar * static_cast<PetscScalar>( factor );

    std::vector<PetscInt> globalDofs;
    std::vector<PetscScalar>  values;

    // Constraint iterators
    ConstraintVec_::const_iterator cIter = constraints.begin( );
    ConstraintVec_::const_iterator cEnd  = constraints.end( );
    // collect global dof indices and the correpsonding values
    for ( ; cIter != cEnd; ++cIter ) {
        globalDofs.push_back( static_cast<PetscInt>( cIter -> first ) );
        values.push_back( static_cast<PetscScalar>( cIter -> second ) * rhsScalar );
    }

    // prepare pointers to be passed to PETSc
    PetscInt * globalDofPtr = corlib::detail_::makePetscPointer( globalDofs );
    PetscScalar * valuePtr  = corlib::detail_::makePetscPointer( values );

    // set matrix rows to zero 
    ierr_ = MatZeroRows( sysMat_, numConstraints, globalDofPtr, diagScalar );
    CHKERRV( ierr_ ); 

    // set RHS entries to values (crashes if called with NULL pointers)
    if ( numConstraints ) {
        ierr_ = VecSetValues( rhsVec_, numConstraints, globalDofPtr, valuePtr, INSERT_VALUES );
        CHKERRV( ierr_ ); 
    }

    // perform communication in the rhs vector
    ierr_ = VecAssemblyBegin( rhsVec_ ); CHKERRV( ierr_ ); 
    ierr_ = VecAssemblyEnd(   rhsVec_ ); CHKERRV( ierr_ ); 

    return;
}

//------------------------------------------------------------------------------
//! Simply set the dof with index 'index' to zero
void corlib::SystemSolvePetsc::fixDOF( const unsigned index, const double scalar )
{
    std::vector< std::pair<unsigned,double> > thisConstraint;
    thisConstraint.push_back( std::make_pair( index, 0. ) );
    this -> applyConstraints( thisConstraint, 1., scalar );
    thisConstraint.clear();
    return;
}

//------------------------------------------------------------------------------
/** Let Petsc solve the system  */
void corlib::SystemSolvePetsc::solveSystem( bool reuseOperators, bool removeConstantNullSpace )
{
    // if this is the first call to this routine
    if ( solver_ == NULL ) {
        // solver object
        ierr_ = KSPCreate( comm_, &solver_ ); CHKERRV( ierr_ ); 

        // insert extra options if they were specified
        // Warning: PETSs inserts options globally at this point, so other instances of this SystemSolvePetsc get them as well
        //          to avoid unexpected behaviour, one should either specify desired options for any solver, or specify
        //          options for all instances of the solver
        if ( not myPetscOptions_.empty( ) ) {
            ierr_ = PetscOptionsInsertString( myPetscOptions_.c_str( ) ); CHKERRV( ierr_ ); 
        }

        // reasonable default if preconditioner was not set on command line
        // this is needed for sequential runs with MPIAIJ matrix
        PetscTruth flg;
        ierr_ = PetscOptionsHasName( PETSC_NULL, "-pc_type", &flg ); CHKERRV( ierr_ ); 
        if ( flg == PETSC_FALSE ) {
            ierr_ = PetscOptionsInsertString("-pc_type bjacobi"); CHKERRV( ierr_ ); 
        }
        ierr_ = KSPSetFromOptions( solver_ ); CHKERRV( ierr_ ); 

        // avoid illeagal settings
        FTL_VERIFY_DESCRIPTIVE( not reuseOperators, "Nothing is loaded in solver, cannot reuse it. \n");  

        // create solution vector
        ierr_ = VecDuplicate( rhsVec_, &sol_ ); CHKERRV( ierr_ ); 

        // create a large local solution vector for gathering 
        ierr_ = VecCreateSeq( PETSC_COMM_SELF, numTotalDofs_, &solGathered_ ); CHKERRV( ierr_ ); 

        // prepare grathering scatter
        ierr_ = VecScatterCreate( sol_, PETSC_NULL, solGathered_, PETSC_NULL, &VSdistToLarge_ ); CHKERRV( ierr_ ); 
    }
    else {
        // clear solution
        this -> clearSol( );
    }

    // set operators - this is needed each time matrix changes to correctly update preconditioner
    // SAME_NONZERO_PATTERN may not always be the case - other options are DIFFERENT_NONZERO_PATTERN and SAME_PRECONDITIONER
    if ( not reuseOperators ) {
        ierr_ = KSPSetOperators( solver_, sysMat_, sysMat_, SAME_NONZERO_PATTERN ); CHKERRV( ierr_ ); 
    }

    // remove one-dimensional nullspace of constants if wanted
    if ( removeConstantNullSpace ) {
        MatNullSpace nullSpace;
        ierr_ = MatNullSpaceCreate( comm_, PETSC_TRUE, 0, PETSC_NULL, &nullSpace ); CHKERRV( ierr_ ); 
        ierr_ = KSPSetNullSpace( solver_, nullSpace ); CHKERRV( ierr_ ); 
        MatNullSpaceDestroy( nullSpace );
    }

    // solve
    ierr_ = KSPSolve( solver_, rhsVec_, sol_ ); CHKERRV( ierr_ ); 

    // analyze solution process
    ierr_ = KSPGetIterationNumber( solver_, &numIterations_ ); CHKERRV( ierr_ ); 
    ierr_ = KSPGetConvergedReason( solver_, &convergedReason_ ); CHKERRV( ierr_ ); 
    ierr_ = KSPGetType( solver_, &kspType_ ); CHKERRV( ierr_ ); 
    PC precond;
    ierr_ = KSPGetPC( solver_, &precond ); CHKERRV( ierr_ ); 
    ierr_ = PCGetType( precond, &pcType_ ); CHKERRV( ierr_ );

    // Compute residual res = A*u - b and store it into rhsVec_
    ierr_ = VecScale( rhsVec_, -1.0 ); CHKERRV( ierr_ );
    ierr_ = MatMultAdd( sysMat_, sol_, rhsVec_, rhsVec_ ); CHKERRV( ierr_ );

    // get global solution vector at each process
    ierr_ = VecScatterBegin( VSdistToLarge_, sol_, solGathered_, INSERT_VALUES, SCATTER_FORWARD ); CHKERRV( ierr_ ); 
    ierr_ = VecScatterEnd(   VSdistToLarge_, sol_, solGathered_, INSERT_VALUES, SCATTER_FORWARD ); CHKERRV( ierr_ ); 

    return;
}

//------------------------------------------------------------------------------
/** Assign values of the solution vector to an iterator using iterators which
 *  indicated the desired range of global DOF-indices.
 */
template< typename IT1, typename IT2 >
void corlib::SystemSolvePetsc::giveDofSolution( IT1 dofIter, IT1 dofEnd, 
                                                IT2 valIter ) const
{
    // PETSc's error code
    PetscErrorCode ierr;

    // Number of dof-indices
    PetscInt nDofs = static_cast<PetscInt>( std::distance( dofIter, dofEnd ) );
    std::vector<PetscInt>    effectiveDofs( dofIter, dofEnd );
    std::vector<PetscScalar> dofSolutions( nDofs );

    // simply get the dof solutions by using the global indices
    ierr = VecGetValues( solGathered_, nDofs, &(effectiveDofs[0]), &(dofSolutions[0]) ); CHKERRV( ierr );

    std::transform( dofSolutions.begin(), dofSolutions.end(), 
                    valIter , corlib::detail_::PetscScalar2Double() );

    return;
}

//------------------------------------------------------------------------------
void corlib::SystemSolvePetsc::giveSolution( const VecUint_ & dofIndices,
                                             SubVec_ & result ) const
{
    this->giveDofSolution( dofIndices.begin( ), dofIndices.end( ),
                           result.begin( ) );
}

//------------------------------------------------------------------------------
void corlib::SystemSolvePetsc::ownConvergenceTest( const std::vector<unsigned> & dofIndices, 
                                                   double & residual1,  
                                                   double & residual2 ) 
{
    // Number of dof-indices
    std::vector<PetscInt>   dofIndPetsc;
    std::vector<unsigned>::const_iterator itDofB = dofIndices.begin();
    for ( ; itDofB != dofIndices.end(); ++itDofB ) { 
        PetscInt index = static_cast<PetscInt>( *itDofB );
        if ( index >= rhsLow_ and index < rhsHigh_) {
            dofIndPetsc.push_back( index );
        }
    }
    // debug print
    //std::vector<PetscInt>::const_iterator itB = dofIndPetsc.begin();
    //std::cout << "list of pressure dofs on process " << rank_ << std::endl;
    //for ( ; itB != dofIndPetsc.end(); ++itB ) { 
    //    std::cout << *itB << std::endl;
    //}
    //std::cout.flush();

    PetscInt nDofs = ( dofIndPetsc.size() );

    // Create index set of selected dofs
    IS subsetIs1;
    IS subsetIs2;
    ierr_ = ISCreateGeneral(comm_, nDofs, &(dofIndPetsc[0]), &subsetIs1); CHKERRV( ierr_ );

    ierr_ = ISComplement( subsetIs1, rhsLow_, rhsHigh_, &subsetIs2 ); CHKERRV( ierr_ );

    Vec resPart;
    ierr_ = VecDuplicate( rhsVec_, &resPart ); CHKERRV( ierr_ );

    // filter first block
    VecScatter VSdistToSubvector1;
    ierr_ = VecScatterCreate( rhsVec_, subsetIs1, resPart, subsetIs1, &VSdistToSubvector1 ); CHKERRV( ierr_ );
    // filter second block
    VecScatter VSdistToSubvector2;
    ierr_ = VecScatterCreate( rhsVec_, subsetIs2, resPart, subsetIs2, &VSdistToSubvector2 ); CHKERRV( ierr_ );

    // get first subvector at each process
    ierr_ = VecSet( resPart, 0.0 ); CHKERRV( ierr_ );
    ierr_ = VecScatterBegin( VSdistToSubvector1, rhsVec_, resPart, INSERT_VALUES, SCATTER_FORWARD ); CHKERRV( ierr_ ); 
    ierr_ = VecScatterEnd(   VSdistToSubvector1, rhsVec_, resPart, INSERT_VALUES, SCATTER_FORWARD ); CHKERRV( ierr_ ); 

    PetscReal normVal;

    // and measure its norm
    ierr_ = VecNorm( resPart, NORM_2, &normVal ); CHKERRV( ierr_ ); 

    // || res_1 ||_2
    residual1 = static_cast<double> ( normVal );

    // get second subvector at each process
    ierr_ = VecSet( resPart, 0.0 ); CHKERRV( ierr_ );
    ierr_ = VecScatterBegin( VSdistToSubvector2, rhsVec_, resPart, INSERT_VALUES, SCATTER_FORWARD ); CHKERRV( ierr_ ); 
    ierr_ = VecScatterEnd(   VSdistToSubvector2, rhsVec_, resPart, INSERT_VALUES, SCATTER_FORWARD ); CHKERRV( ierr_ ); 

    // and measure its norm
    ierr_ = VecNorm( resPart, NORM_2, &normVal ); CHKERRV( ierr_ ); 

    // || res_2 ||_2
    residual2 = static_cast<double> ( normVal );


    ierr_ = VecScatterDestroy( VSdistToSubvector2 ); CHKERRV( ierr_ );
    ierr_ = VecScatterDestroy( VSdistToSubvector1 ); CHKERRV( ierr_ );
    VecDestroy( resPart );
    ISDestroy( subsetIs1 );
    ISDestroy( subsetIs2 );
}
