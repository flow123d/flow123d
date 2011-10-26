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

#ifndef corlib_systemsolvepetsc_h
#define corlib_systemsolvepetsc_h
//------------------------------------------------------------------------------
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <ostream>
#include <vector>
#include <set>
#include <cmath>
#include <petscksp.h>
#include <corlib/Triplet.hpp>

//------------------------------------------------------------------------------
namespace corlib{
    class SystemSolvePetsc;
    namespace ublas = boost::numeric::ublas;

    namespace detail_{
        //! \cond SKIPDOX
        // unsigned int to PetscInt casting functor
        struct Uint2PetscInt : public std::unary_function< unsigned, PetscInt >
        {
            PetscInt operator()( unsigned int arg ) 
            {
                return static_cast<PetscInt>( arg );
            }
        };
        // PetscScalar to double casting functor
        struct PetscScalar2Double : public std::unary_function< PetscScalar, double >
        {
            double operator()( PetscScalar arg ) 
            {
                return static_cast<double>( arg );
            }
        };    
        // make a pointer to the data array out of a std::vector
        template<typename T> 
        T *  makePetscPointer( std::vector<T> & array )
        {
            if ( array.size() ) return &(array[0]);
            return PETSC_NULL;
        }
        //! \endcond
    }
}

//------------------------------------------------------------------------------
/** \brief PETSc based linear system solver
 *  \details This class implements solution of a system of linear equations 
 *  by PETSc solvers. Particular solution method and preconditioner 
 *  can be inserted to the executable at command line by -ksp_type and -pc_type 
 *  together with related options, e.g. for Additive Schwarz preconditioner for GMRES
 *  and overlap 10 elements:
 *  -ksp_type gmres -pc_type asm -pc_asm_overlap 10 
 *  or for direct solution by MUMPS:
 *  -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps
 *
 *  The class uses preallocation of the PETSc matrix - this is performed via
 *  the triplet format implemented in Triplet.hpp file.
 *  The matrix is first constructed into the triplet and then converted to PETSc 
 *  MPIAIJ matrix.
 *  Until revision 2934, the version was without preallocation.
 *
 *  Filtering of local values can be switched on. This leads to insertion only values strictly local
 *  to process and thus avoiding later assembly of PETSc matrix. An extreme of using this is when 
 *  whole mesh is available to each process. A more sophisticated use of this functionality 
 *  employs overlapping partitioning, which for each process contains all elements contributing 
 *  to rows of this process.
 *
 *  The class is used as follows:
 *
 *  prepareMatAssembly - (optional) - reserves memory for sparse triplet
 *  insertToMatrix     - called many times, inserts values into triplet
 *  insertToRhs        - called many times, inserts values into RHS.
 *                       If matrix is created and assembled, structure of the vector is deduced from it. 
 *                       Otherwise, it is created as empty.
 *  finishAssembly     - assembles matrix (if not already done for RHS creation) and assembles RHS
 *  applyConstraints   - enforces boundary conditions to matrix and RHS
 *  solveSystem        - solves the system
 *  clearMatrix        - blanks matrix values, keeps structure
 *  clearRhs           - blanks RHS values, keeps structure
 */
class corlib::SystemSolvePetsc
{
private:
    typedef ublas::matrix< double >                         SubMat_;
    typedef ublas::vector< double >                         SubVec_;
    typedef std::vector< unsigned >                         VecUint_;
    typedef std::vector< std::pair<unsigned,double> >       ConstraintVec_;

    // produce initial state of object
    void initMyValues () {
        sysMat_ = NULL;
        rhsVec_ = NULL;
        solver_ = NULL;
        sol_ = NULL;
        solGathered_ = NULL;
        VSdistToLarge_ = NULL;

        isMatAssembled_ = false;
        isRhsAssembled_ = false;

        // orient in the communicator
        MPI_Comm_size( comm_, &nProc_ );
        MPI_Comm_rank( comm_, &rank_ );

        myPetscOptions_.clear( );
    }

public:
    //! Constructor with global dof numbers and optional communicator ( PETSC_COMM_WORLD is the default )
    SystemSolvePetsc( const unsigned numTotalDofs ) 
      : numTotalDofs_( static_cast<PetscInt>( numTotalDofs ) ),
        numLocalDofs_( -1 ),
        comm_ ( PETSC_COMM_WORLD ),
        filterLocalRows_ ( false )
    {
        this -> initMyValues();
    }

    //! Constructor with 
    //! global dof number, 
    //! local dof number, 
    //! optional filtering ( false is the default )
    //! and optional communicator ( PETSC_COMM_WORLD is the default )
    SystemSolvePetsc( const unsigned numTotalDofs, const unsigned numLocalDofs, 
                      bool filterLocalRows = false ,
                      const MPI_Comm comm = PETSC_COMM_WORLD )
      : numTotalDofs_( static_cast<PetscInt>( numTotalDofs ) ),
        numLocalDofs_( static_cast<PetscInt>( numLocalDofs ) ),
        comm_ (comm),
        filterLocalRows_ ( filterLocalRows )
    {
        this -> initMyValues();

        procDofsStarts_.resize( nProc_ + 1 );

        int numLocalDofsInt = static_cast<int>( numLocalDofs_ );
        MPI_Allgather( &numLocalDofsInt,      1, MPI_INT,
                       &(procDofsStarts_[0]), 1, MPI_INT, comm_ );
        // now procDofsStarts_ contains counts of local dofs at each process
        // change it to starts
        for ( unsigned i = 1; i < static_cast<unsigned>( nProc_ ); i++ )
            procDofsStarts_[i] = procDofsStarts_[i-1] + procDofsStarts_[i];
        // shift it one back
        for ( unsigned i = nProc_; i > 0; i-- )
            procDofsStarts_[i] = procDofsStarts_[i-1];
        // start from zero
        procDofsStarts_[0] = 0;

    }


    //! Destructor frees the Petsc structures
    ~SystemSolvePetsc()
    {
        this -> destroyMatrix();
        this -> destroyRhs();
        if( solver_ != NULL ){
            ierr_ = KSPDestroy( solver_ ); CHKERRV( ierr_ );
        }
        if( sol_ != NULL ){
            ierr_ = VecDestroy( sol_ ); CHKERRV( ierr_ );
        }
        if( solGathered_ != NULL ){
            ierr_ = VecDestroy( solGathered_ ); CHKERRV( ierr_ );
        }
        if( VSdistToLarge_ != NULL ){
            ierr_ = VecScatterDestroy( VSdistToLarge_ ); CHKERRV( ierr_ );
        }
    }

    //! Prepare assembly of matrix - reserve space in underlying triplet object
    void prepareMatAssembly( unsigned numElements, unsigned elMatSize );

    //! insert submatrix to system matrix - into triplet object
    void insertToMatrix( const SubMat_ & eStiff, 
                         const VecUint_ & rowIndices,
                         const VecUint_ & colIndices );

    //! Assemble the triplet at the background, but do not convert it into PETSc matrix yet.
    void partialMatAssembly( );

    //! Finalize assembly of matrix - finish assembly of triplet, preallocate PETSc matrix, and copy triplet into it.
    void finishMatAssembly( );

    //! Prepare assembly of RHS - derive parallel layout of RHS from the assembled matrix
    void prepareRhsAssembly( );

    //! insert subvector to system rhs vector
    void insertToRhs( const SubVec_ & eForce, const VecUint_ & dofIndices );

    //! insert a scalar to RHS
    void insertValueToRhs( const double & entry, const unsigned & dofIndex );

    //! Finalize assembly of RHS
    void finishRhsAssembly( );

    //! Compatibility function - calls finishMatAssembly and finishRhsAssembly
    void finishAssembly( );

    //! Outputs
    void writeMatrix(   std::ostream & out ) { MatView(sysMat_,PETSC_VIEWER_STDOUT_WORLD); }
    void writeRhs(      std::ostream & out ) { VecView(rhsVec_,PETSC_VIEWER_STDOUT_WORLD); }
    void writeSolution( std::ostream & out ) { VecView(sol_,PETSC_VIEWER_STDOUT_WORLD); }
    void writeMatrixMatlab( const std::string filename ) { 

        PetscViewer myViewer;
        PetscViewerASCIIOpen(comm_,filename.c_str(),&myViewer);
        PetscViewerSetFormat(myViewer,PETSC_VIEWER_ASCII_MATLAB);
        MatView( sysMat_, myViewer );
        PetscViewerDestroy(myViewer);
    }
    void writeRhsMatlab( const std::string filename ) { 

        PetscViewer myViewer;
        PetscViewerASCIIOpen(comm_,filename.c_str(),&myViewer);
        PetscViewerSetFormat(myViewer,PETSC_VIEWER_ASCII_MATLAB);
        VecView( rhsVec_, myViewer );
        PetscViewerDestroy(myViewer);
    }
    void writeSolutionMatlab( const std::string filename ) { 

        PetscViewer myViewer;
        PetscViewerASCIIOpen(comm_,filename.c_str(),&myViewer);
        PetscViewerSetFormat(myViewer,PETSC_VIEWER_ASCII_MATLAB);
        VecView( sol_, myViewer );
        PetscViewerDestroy(myViewer);
    }

    //! Compute norm of rhs vector
    double normRhs( ) 
    { 
        FTL_VERIFY_DESCRIPTIVE( isRhsAssembled_,
                                "Cannot compute norm of unassembled vectors. Call finishRhsAssembly first.\n " );
        double result;
        ierr_ = VecNorm( rhsVec_, NORM_2, &result ); CHKERRQ( ierr_ );
        return result;
    }

    //! Compute 1-norm of rhs vector
    double norm1Rhs( ) 
    { 
        FTL_VERIFY_DESCRIPTIVE( isRhsAssembled_,
                                "Cannot compute norm of unassembled vectors. Call finishRhsAssembly first.\n " );
        double result;
        ierr_ = VecNorm( rhsVec_, NORM_1, &result ); CHKERRQ( ierr_ );
        return result;
    }


    //! Compute norm of the solution vector
    double normSol( ) 
    { 
        double result;
        ierr_ = VecNorm( sol_, NORM_2, &result ); CHKERRQ( ierr_ );
        return result;
    }

    //! Apply constraints
    void applyConstraints( ConstraintVec_ & constraints, const double factor, const double scalar );

    //! Set any DOF to zero
    void fixDOF( const unsigned index, const double scalar );

    //! Insert particular PETSc options to solver
    void insertOptions( const std::string additionalOptions ) {
        myPetscOptions_.append( " " + additionalOptions + " " );
    }

    //! Solve the system
    void solveSystem( bool reuseOperators = false, bool removeConstantNullSpace = false );

    //! give number of iteration
    unsigned    giveNumIterations()   const { return static_cast<unsigned>( numIterations_ ); }
    //! give reason of convergence
    int         giveConvergedReason() const { return static_cast<int>( convergedReason_ ); }
    //! give preconditioner type
    std::string givePCType()          const { return static_cast<std::string>( pcType_ ); }
    //! give iterative solver type
    std::string giveKSPType()         const { return static_cast<std::string>( kspType_ ); }

    //! Get scalar suitable for diagonal while fixing BC
    double getDiagScalar( ) { return diagScalar_; }
    
    //! Fill matrix with zeros
    void clearMatrix( ) { 
        if( sysMat_ != NULL ){ 
            ierr_ = MatZeroEntries(sysMat_); CHKERRV( ierr_ ); 
        }
    }
    //! Fill RHS with zeros
    void clearRhs( )    { 
        if( rhsVec_ != NULL ){
            ierr_ = VecSet( rhsVec_, 0.);  CHKERRV( ierr_ ); 
        }
    }
    //! Fill solution with zeros
    void clearSol( )    { 
        if( sol_ != NULL ){
            ierr_ = VecSet( sol_, 0.);  CHKERRV( ierr_ ); 
        }
    }

    //! Erase matrix completely from memory - this can be faster than zeroing
    void destroyMatrix( ) { 
        if( sysMat_ != NULL ){
            ierr_ = MatDestroy( sysMat_ ); CHKERRV( ierr_ );
            sysMat_ = NULL;
            isMatAssembled_ = false;
        }
    }
    //! Erase RHS completely from memory 
    void destroyRhs( ) { 
        if( rhsVec_ != NULL ){
            ierr_ = VecDestroy( rhsVec_ ); CHKERRV( ierr_ );
            rhsVec_ = NULL;
            isRhsAssembled_ = false;
        }
    }

    //! Give access to the solution vector
    template< typename IT1, typename IT2 >
    void giveDofSolution( IT1 db, IT1 de, IT2 vb ) const;

    //! Give access to the solution vector
    void giveSolution( const VecUint_ & dofIndices, SubVec_ & result ) const;

    //! Compute two residuals based on set of indices 
    //  - initially meant to distinguish velocity and pressure dofs in Taylor-Hood computations
    void ownConvergenceTest( const std::vector<unsigned> & dofIndices, 
                             double & residual1,  
                             double & residual2 );

private:
    const PetscInt  numTotalDofs_;       //!< number of total dofs in the system
    const PetscInt  numLocalDofs_;       //!< number of local dofs, i.e. rows of the MPI matrix at the process
    const MPI_Comm   comm_;              //!< communicator
    const bool  filterLocalRows_;        //!< if active, PETSc will only insert values local to process
    int   nProc_;                        //!< number of processors in communicator
    int   rank_;                         //!< ID of processors in communicator
    corlib::Triplet<PetscInt,PetscReal>  triplet_; //!< matrix in triplet format
    double     diagScalar_;              //!< recommended value for to be put on diagonal in fixing BC
    Mat        sysMat_;                  //!< system matrix
    PetscInt   matLow_, matHigh_;        //!< [matLow,matHigh) is the range of dof indices in matrix
    std::vector<int> procDofsStarts_;    //!< where local Dofs starts determined outside PETSc - initialized only
                                         //!< for constructor with number of local dofs given
    bool       isMatAssembled_;          //!< true if matrix is assembled
    Vec        rhsVec_;                  //!< vector with RHS values
    PetscInt   rhsLow_, rhsHigh_;        //!< [rhsLow,rhsHigh) is the range of dof indices in right-hand side
    bool       isRhsAssembled_;          //!< true if RHS is assembled
    KSP        solver_;                  //!< solver context
    std::string myPetscOptions_;         //!< extra options - to be able to specify solver type from applications
    Vec        sol_;                     //!< distributed solution vector
    Vec        solGathered_;             //!< large local solution vector
    VecScatter VSdistToLarge_;           //!< scatter for gathering local parts into large vector
    PetscInt   numIterations_;           //!< number of iterations
    KSPConvergedReason convergedReason_; //!< how convergence was obtained
    const PCType  pcType_;               //!< Petsc's preconditioner type - string with name
    const KSPType kspType_;              //!< Petsc's solver type - string with name
    PetscErrorCode ierr_;                //!< Petsc's error handler
};

//-----------------------------------------------------------------------------

#include "SystemSolvePetsc.ipp"

//-----------------------------------------------------------------------------

#endif
