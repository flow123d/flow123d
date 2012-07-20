//! @author Jakub Sistek
//! @date   5/2011

//------------------------------------------------------------------------------
/** Constructor for the BDDCML solver object
 */
la::SystemSolveBddc::SystemSolveBddc( const unsigned numDofs,
                                      const unsigned numDofsSub,
                                      const MatrixType matrixType, 
                                      const MPI_Comm comm,
                                      int numSubLoc )
  : numDofs_( static_cast<int>( numDofs ) ),
    numDofsSub_( static_cast<int>( numDofsSub ) ),
    matrixType_( matrixType ),
    comm_ (comm),
    numSubLoc_ ( numSubLoc )
{

    // Orient in the communicator
    MPI_Comm_size( comm_, &nProc_ );
    MPI_Comm_rank( comm_, &rank_ );

    // prepare RHS vector
    rhsVec_.resize( numDofsSub_, 0. );

    // prepare arrays to mark boundary conditions
    ifix_.resize( numDofsSub_, 0 );
    fixv_.resize( numDofsSub_, 0. );

    // initialize subdomain solution vector
    sol_.resize( numDofsSub_, 0. );
    normSol_ = 0.;
    numIter_ = 0;
    convergedReason_ = 0;
    condNumber_ = 0.;

    // initialize state of mesh and matrix
    meshLoaded_     = false;
    isMatAssembled_ = false;
}
//------------------------------------------------------------------------------
/** Destructor frees the BDDC solver structures
 */
la::SystemSolveBddc::~SystemSolveBddc()
{
    // clear RHS
    rhsVec_.clear();

    // clear BC
    ifix_.clear();
    fixv_.clear();

    // clear solution
    sol_.clear();

    // clear mesh data
    inet_.clear();            
    nnet_.clear();            
    nndf_.clear();            
    xyz_.clear();             
    isngn_.clear();         
    isvgvn_.clear();         
    isegn_.clear();           

    // clear global2local map
    global2LocalDofMap_.clear();

    // clear triplet
    triplet_.clear();

}

//------------------------------------------------------------------------------
/** Routine for loading raw mesh data into the solver - for cases of strange meshes, 
 * where these are not Mesh or Grid objects, user can create own raw description
 * to exploit the flexibility of mesh format underlaying BDDCML.
 */
void la::SystemSolveBddc::loadRawMesh( const int nDim, const int numNodes, const int numDofs,
                                           const std::vector<int> & inet, 
                                           const std::vector<int> & nnet, 
                                           const std::vector<int> & nndf, 
                                           const std::vector<int> & isegn, 
                                           const std::vector<int> & isngn, 
                                           const std::vector<int> & isvgvn,
                                           const std::vector<double> & xyz,
                                           const int meshDim )
{
    // space dimension
    nDim_ = nDim;

    // set topological dimension of mesh
    if ( meshDim == 0 ) { // if not specified, use space dimension
       meshDim_ = nDim_;
    }
    else {
       meshDim_ = static_cast<int>( meshDim ); // if specified, use it
    }

    // fill internal variables about mesh
    numElemSub_  = isegn.size( );
    numNodesSub_ = isngn.size( );

    // global number of nodes
    numNodes_ = numNodes;

    std::vector<int> procElemStarts( nProc_ + 1 );
    int numElemSubInt = static_cast<int>( numElemSub_ );
    MPI_Allgather( &numElemSubInt,       1, MPI_INT,
                   &(procElemStarts[0]), 1, MPI_INT, comm_ );
    // now procElemStarts_ contains counts of local elements at each process
    // change it to starts
    for ( unsigned i = 1; i < static_cast<unsigned>( nProc_ ); i++ )
        procElemStarts[i] = procElemStarts[i-1] + procElemStarts[i];
    // shift it one back
    for ( unsigned i = nProc_; i > 0; i-- )
            procElemStarts[i] = procElemStarts[i-1];
    // start from zero
    procElemStarts[0] = 0;

    // find global number of elements - assumes nonoverlapping division
    numElem_ = procElemStarts[ nProc_ ];

    // check sizes of arrays
    ASSERT( std::accumulate( nnet.begin(), nnet.end(), 0 ) == inet.size(),
                            "array inet size mismatch \n " );
    ASSERT( std::accumulate( nndf.begin(), nndf.end(), 0 ) == numDofsSub_,
                            "array nndf content mismatch: %d %d \n ",  std::accumulate( nndf.begin(), nndf.end(), 0 ), numDofsSub_ );
    ASSERT( isvgvn.size() == numDofsSub_,
                            "array isvgvn size mismatch \n " );
    ASSERT( nnet.size() == numElemSub_,
                            "array nnet size mismatch \n " );
    ASSERT( nndf.size() == numNodesSub_,
                            "array nndf size mismatch \n " );
    ASSERT( xyz.size() == numNodesSub_ * nDim_,
                            "array xyz size mismatch: %d %d \n ", xyz.size(), numNodesSub_ * nDim_ );

    // Simply copy input data into the private object equivalents
    inet_.resize( inet.size() );
    std::copy( inet.begin(), inet.end(), inet_.begin() );
    nnet_.resize( nnet.size() );
    std::copy( nnet.begin(), nnet.end(), nnet_.begin() );
    nndf_.resize( nndf.size() );
    std::copy( nndf.begin(), nndf.end(), nndf_.begin() );
    isegn_.resize( isegn.size() );
    std::copy( isegn.begin(), isegn.end(), isegn_.begin() );
    isngn_.resize( isngn.size() );
    std::copy( isngn.begin(), isngn.end(), isngn_.begin() );
    isvgvn_.resize( isvgvn.size() );
    std::copy( isvgvn.begin(), isvgvn.end(), isvgvn_.begin() );
    xyz_.resize( xyz.size() );
    std::copy( xyz.begin(), xyz.end(), xyz_.begin() );

    // prepare auxiliary map for renumbering dofs in the global array
    for ( unsigned ind = 0; ind < isvgvn_.size(); ++ind ) {
        global2LocalDofMap_.insert( std::make_pair( static_cast<unsigned>( isvgvn_[ind] ), ind ) );
    }

    // change the state
    meshLoaded_     = true;

    return;
}

//------------------------------------------------------------------------------
/** Routine for preparing MATRIX assembly. It reserves memory in triplet.
 */
void la::SystemSolveBddc::prepareMatAssembly( unsigned numElements, unsigned elMatSize )
{
    unsigned length = numElements * elMatSize;
    triplet_.prepareAssembly( length );

    return;
}
//------------------------------------------------------------------------------
/** Insert element stiffness matrix into global system matrix
 */
void la::SystemSolveBddc::insertToMatrix( const SubMat_  & subMat, 
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
/** Routine for triplet finalization.
 */
void la::SystemSolveBddc::finishMatAssembly( )
{
    // do nothing if already called
    if ( isMatAssembled_ ) return;

    // finish assembly of triplet
    triplet_.finishAssembly( );

    // set flag
    isMatAssembled_ = true;

    return;
}
//------------------------------------------------------------------------------
/** Insert subvector to the system's RHS vector
 */
void la::SystemSolveBddc::insertToRhs( const SubVec_  & subVec, 
                                           const VecUint_ & dofIndices )
{
    for ( unsigned i = 0; i < dofIndices.size( ); i ++ ) {

        // map global dof index to subdomain dof
        Global2LocalMap_::iterator pos = global2LocalDofMap_.find( dofIndices[i] );
        ASSERT( pos != global2LocalDofMap_.end(),
                                "Cannot remap index %d to local indices in right-hand side," 
                                " perhaps missing call to loadMesh() member function. \n ", dofIndices[i] );
        const unsigned indLoc = pos -> second;

        rhsVec_[ indLoc ] += subVec[i];
    }
    return;
}

//------------------------------------------------------------------------------
/** Apply boundary conditions
 *  Since BDDCML handles constraints, just fill proper BC arrays.
 *  Scalar is not used, it is kept for compatibility with constraint functors.
 */
void la::SystemSolveBddc::applyConstraints( ConstraintVec_ & constraints,
                                                const double factor, const double scalar ) 
{
    // Constraint iterators
    ConstraintVec_::const_iterator cIter = constraints.begin( );
    ConstraintVec_::const_iterator cEnd  = constraints.end( );
    // collect global dof indices and the correpsonding values
    for ( ; cIter != cEnd; ++cIter ) {
        unsigned globalDof =  cIter -> first;
        double   value    = (cIter -> second );

        // map global dof index to subdomain dof
        Global2LocalMap_::iterator pos = global2LocalDofMap_.find( globalDof );
        if ( pos != global2LocalDofMap_.end() ) {
            // I contain the dof
            unsigned indLoc = pos -> second;

            ifix_[ indLoc ] = 1;
            fixv_[ indLoc ] = factor * value;
        }
    }

    return;
}
//------------------------------------------------------------------------------
/** Simply set the dof with index 'index' to zero
 */
void la::SystemSolveBddc::fixDOF( const unsigned index, const double scalar )
{
    std::vector< std::pair<unsigned,double> > thisConstraint;
    thisConstraint.push_back( std::make_pair( index, 0. ) );
    this -> applyConstraints( thisConstraint, 1., scalar );
    thisConstraint.clear();
    return;
}
//------------------------------------------------------------------------------
/** Solve the system by BDDCML
 */
void la::SystemSolveBddc::solveSystem( double tol, int  numLevels, std::vector<int> *  numSubAtLevels, 
                                           int verboseLevel, int  maxIt, int  ndecrMax, bool use_adaptive )
{

    // check that mesh was loaded 
    ASSERT( meshLoaded_ ,
                            "Subdomain mesh was not loaded, perhaps missing call to loadMesh() member function. \n " );
   
    //============= initialize BDDC
    int numSub;                  //!< number of subdomains on first level   
    // Number of subdomains at first level - if not given, deduce it from MPI communicator
    MPI_Allreduce ( &numSubLoc_, &numSub, 1, MPI_INT, MPI_SUM, comm_ );

    // Number of subdomains at levels
    std::vector<int> numSubLev( numLevels ); //!< numbers of subdomains on levels ( should be monotonically decreasing, e.g. for 3 levels 1024, 32, 1 )
    if ( numSubAtLevels != NULL ) {
        // numbers are given 
        
        ASSERT( numLevels == (*numSubAtLevels).size(), 
                                "Inconsistent size of numbers of subdomains at levels: %d %d \n",
                                numLevels, (*numSubAtLevels).size() );
        ASSERT( (*numSubAtLevels)[0] == numSub, 
                                "Inconsistent number of subdomains at first level: %d %d \n",
                                (*numSubAtLevels)[0], numSub );

        std::copy( (*numSubAtLevels).begin(), (*numSubAtLevels).end(),
                    numSubLev.begin() );
    }
    else {
        // assume 2-level BDDC
        // generate default vector for two levels
        if ( numLevels == 2 ) {
            numSubLev[0] = numSub;
            numSubLev[1] = 1;
        }
        else {
            ASSERT( false, "Missing numbers of subdomains at levels. \n" );
        }
    }

    int lnumSubLev   = numSubLev.size();
    int commInt      = MPI_Comm_c2f( comm_ ); // translate C MPI communicator hanlder to Fortran integer
    int numBase      = 0;                     // indexing starts with 0 for C++

    bddcml_init( &numLevels, &(numSubLev[0]), &lnumSubLev, &numSubLoc_, &commInt, &verboseLevel, &numBase );


    //============= uploading data for subdomain
    // at the moment, number of subdomains equals number of processors
    int iSub = rank_;

    int linet   = inet_.size();
    int lnnet   = nnet_.size();
    int lnndf   = nndf_.size();
    int lisngn  = isngn_.size();
    int lisvgvn = isvgvn_.size();
    int lisegn  = isegn_.size();

    int lxyz1   = numNodesSub_;
    int lxyz2   = nDim_;

    int lifix   = ifix_.size();
    int lfixv   = fixv_.size();

    int lrhsVec = rhsVec_.size();
    int isRhsCompleteInt = 0; // since only local subassembly of right-hand side was performed, it does not contain 
                              // repeated entries

    std::vector<int>    i_sparse;
    std::vector<int>    j_sparse;
    std::vector<double> a_sparse;

    // specify if matrix is symmetric, if yes, use only upper triangle
    bool onlyUpperTriangle;
    if ( matrixType_ == GENERAL ) {
        onlyUpperTriangle = false;
    }
    else {
        onlyUpperTriangle = true;
    }

    // copy out arrays from triplet
    this -> finishMatAssembly();
    triplet_.extractArrays( i_sparse, j_sparse, a_sparse, onlyUpperTriangle );
    triplet_.clear();

    // map global row and column indices to local subdomain indices
    int indRow    = -1;
    int indCol    = -1;
    int indRowLoc = -1;
    int indColLoc = -1;
    for ( unsigned inz = 0; inz < a_sparse.size(); inz++ ) {
        if ( i_sparse[inz] != indRow ) {
           indRow = i_sparse[inz];
           Global2LocalMap_::iterator pos = global2LocalDofMap_.find( static_cast<unsigned> ( indRow ) );
           ASSERT( pos != global2LocalDofMap_.end(),
                                   "Cannot remap index %d to local indices. \n ", indRow );
           indRowLoc = static_cast<int> ( pos -> second );
        }
        if ( j_sparse[inz] != indCol ) {
           indCol = j_sparse[inz];
           Global2LocalMap_::iterator pos = global2LocalDofMap_.find( static_cast<unsigned> ( indCol ) );
           ASSERT( pos != global2LocalDofMap_.end(),
                                   "Cannot remap index %d to local indices. \n ", indCol );
           indColLoc = static_cast<int> ( pos -> second );
        }

        i_sparse[ inz ] = indRowLoc;
        j_sparse[ inz ] = indColLoc;

        //std::cout << " row: " << i_sparse[ inz ] << " col: " << j_sparse[ inz ] << std::endl;
    }
    // due to remap, it is not guaranteed that the triplet is still correctly ordered
    int isMatAssembledInt = 0;


    int la = a_sparse.size();

    // remove const attribute
    int numDofsInt    = static_cast<int> ( numDofs_ );
    int numDofsSubInt = static_cast<int> ( numDofsSub_ );

    int matrixTypeInt;
    if      ( matrixType_ == GENERAL )                  matrixTypeInt = 0 ;
    else if ( matrixType_ == SPD )                      matrixTypeInt = 1 ;
    else if ( matrixType_ == SYMMETRICGENERAL )         matrixTypeInt = 2 ;
    else if ( matrixType_ == SPD_VIA_SYMMETRICGENERAL ) matrixTypeInt = 2 ; 
    else    ASSERT( false, "Illegal matrixType \n " );

    // download local solution
    int lsol = sol_.size();

    // upload local data to the solver
    bddcml_upload_subdomain_data( &numElem_, &numNodes_, &numDofsInt, &nDim_, &meshDim_, 
                                  &iSub, &numElemSub_, &numNodesSub_, &numDofsSubInt, 
                                  &(inet_[0]),&linet, &(nnet_[0]),&lnnet, &(nndf_[0]),&lnndf, 
                                  &(isngn_[0]),&lisngn, &(isvgvn_[0]),&lisvgvn, &(isegn_[0]),&lisegn, 
                                  &(xyz_[0]),&lxyz1,&lxyz2, 
                                  &(ifix_[0]),&lifix, &(fixv_[0]),&lfixv, 
                                  &(rhsVec_[0]),&lrhsVec, &isRhsCompleteInt,
                                  &(sol_[0]), &lsol,
                                  &matrixTypeInt, &(i_sparse[0]), &(j_sparse[0]), &(a_sparse[0]), &la, &isMatAssembledInt );
    i_sparse.clear();
    j_sparse.clear();
    a_sparse.clear();

    
    //============= setup of BDDC preconditioner
    
    int use_defaults_int          = 0;
    int parallel_division_int     = 0;
    int use_arithmetic_int        = 1;
    int use_adaptive_int;
    if ( use_adaptive ) {
        use_adaptive_int = 1;
    }
    else {
        use_adaptive_int = 0;
    }

    // setup multilevel BDDC preconditioner
    bddcml_setup_preconditioner( & matrixTypeInt, & use_defaults_int, 
                                 & parallel_division_int, & use_arithmetic_int, & use_adaptive_int );

    //============= compute the norm on arrays with overlapping entries - apply weights
    double normSquaredLoc = 0.;
    bddcml_dotprod_subdomain( &iSub, &(rhsVec_[0]), &lrhsVec, &(rhsVec_[0]), &lrhsVec, &normSquaredLoc );

    // sum over processes
    double normSquared = 0.;
    MPI_Allreduce ( &normSquaredLoc, &normSquared, 1, MPI_DOUBLE, MPI_SUM, comm_ );
    normRhs_ = std::sqrt( normSquared );

    //============= solution of the problem
    int method;
    if      ( matrixType_ == SPD )
        method = 0; // PCG - set for SPD matrix
    else if ( matrixType_ == SPD_VIA_SYMMETRICGENERAL)
        method = 0; // PCG - set for benign matrix
    else
        method = 1; // BICGSTAB - set for any other matrix

    // call iterative solver
    bddcml_solve( & commInt, & method, & tol, & maxIt, & ndecrMax, &numIter_, &convergedReason_, &condNumber_);


    //============= downloading solution from BDDCML solver
    
    // download local solution
    lsol = sol_.size();

    bddcml_download_local_solution( &iSub, &(sol_[0]), &lsol );

    //============= compute the norm on arrays with overlapping entries - apply weights
    normSquaredLoc = 0.;
    bddcml_dotprod_subdomain( &iSub, &(sol_[0]), &lsol, &(sol_[0]), &lsol, &normSquaredLoc );

    // sum over processes
    normSquared = 0.;
    MPI_Allreduce ( &normSquaredLoc, &normSquared, 1, MPI_DOUBLE, MPI_SUM, comm_ );
    normSol_ = std::sqrt( normSquared );

    //============= erase memory used by BDDCML solver
    bddcml_finalize();

    return;
}

//------------------------------------------------------------------------------
template<typename VEC1,typename VEC2>
void la::SystemSolveBddc::giveSolution( const VEC1 & dofIndices,
                                        VEC2 & result ) const
{
    // simply get the dof solutions by using the local indices
    typename VEC1::const_iterator dofIter = dofIndices.begin();
    typename VEC1::const_iterator dofEnd  = dofIndices.end();
    typename VEC2::iterator valIter = result.begin();
    for ( ; dofIter != dofEnd; ++dofIter ) {

        // map it to local dof
        Global2LocalMap_::const_iterator pos = global2LocalDofMap_.find( *dofIter );
        ASSERT( pos != global2LocalDofMap_.end(),
                                "Cannot remap index %d to local indices in solution distribution. \n ", *dofIter );
        unsigned indLoc = pos -> second;

        // give solution
        *valIter++ = sol_[ indLoc ];
    }

}
