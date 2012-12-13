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

//! @author Jakub Sistek
//! @date   5/2011

#ifndef corlib_systemsolvebddc_h
#define corlib_systemsolvebddc_h
//------------------------------------------------------------------------------
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <ostream>
#include <vector>
#include <set>
#include <cmath>
#include <corlib/Triplet.hpp>
#include <corlib/Mesh.hpp>

extern "C" { 
    #include <bddcml_interface_c.h>
}

//------------------------------------------------------------------------------
namespace corlib{
    class SystemSolveBddc;
    namespace ublas = boost::numeric::ublas;

}

//------------------------------------------------------------------------------
/** \brief Multilevel BDDC based linear system solver
 *  \details This class implements solution of a system of linear equations 
 *  by multilevel BDDC implementation in BDDCML solver. 
 *
 *  The class uses assembly of the matrix into 
 *  the triplet format implemented in Triplet.hpp file.
 *
 *  prepareMatAssembly - (optional) - reserves memory for sparse triplet
 *  insertToMatrix     - called many times, inserts values into triplet
 *  insertToRhs        - called many times, inserts values into RHS.
 *                       If matrix is created and assembled, structure of the vector is deduced from it. 
 *                       Otherwise, it is created as empty.
 *  applyConstraints   - currently only passes list of boundary conditions into the BDDCML solver, 
 *                       rather then modifying the matrix and rhs
 *                     - as a result, no nonsymmetry is forced to SPD and symmetric indefinite problems by constraints, 
 *                       and these can be accelerated and be more memory efficient by specifying proper matrixType
 *  solveSystem        - solves the system
 */
class corlib::SystemSolveBddc
{
public:
    enum matrixTypeEnum { 
        GENERAL,                   //!< general (default), 
        SPD,                       //!< symmetric positive definite, 
        SYMMETRICGENERAL,          //!< general symmetric ( e.g. saddle point ),
        SPD_VIA_SYMMETRICGENERAL   //!< (advanced) interface problem is symm. pos. def., 
                                   //!< but subdomain problems are general symmetric ( e.g. saddle point )
    };
    typedef enum matrixTypeEnum                             MatrixType;

private:
    typedef ublas::matrix< double >                         SubMat_;
    typedef ublas::vector< double >                         SubVec_;
    typedef std::vector< unsigned >                         VecUint_;
    typedef std::vector< std::pair<unsigned,double> >       ConstraintVec_;

public:
    //! Constructor with 
    //! total number of dofs, 
    //! subdomain number of dofs, 
    //! (optional) type of matrix,
    //! (optional) communicator ( MPI_COMM_WORLD is the default ),
    //! (optional) number of subdomains ( if not given, use number of processes ),
    SystemSolveBddc( const unsigned numDofs,
                     const unsigned numDofsSub,
                     const MatrixType matrixType = GENERAL,
                     const MPI_Comm comm = MPI_COMM_WORLD, 
                     int numSubLoc = 1 );

    //! Destructor frees the BDDCML structures
    ~SystemSolveBddc();

    //! Load data about mesh
    //! number of space dimensions
    //! mesh of subdomain
    //! (optional) topological dimension of mesh, i.e. 2 for shells in 3D, default is nDim
    template< typename MESH >
    void loadMesh( const MESH     & meshPart,
                   const unsigned meshDim = 0 );

    //! Load data about grid
    //! number of space dimensions
    //! sub-grid 
    //! (optional) topological dimension of mesh, i.e. 2 for shells in 3D, default is nDim
    template< typename GRID >
    void loadGrid( const GRID     & meshPart,
                   const unsigned meshDim = 0 );

    //! Load raw data about mesh
    void loadRawMesh( const int nDim, const int numNodes, const int numDofs,
                      const std::vector<int> & inet, 
                      const std::vector<int> & nnet, 
                      const std::vector<int> & nndf, 
                      const std::vector<int> & isegn, 
                      const std::vector<int> & isngn, 
                      const std::vector<int> & isvgvn,
                      const std::vector<double> & xyz,
                      const int meshDim = 0 );

    //! Prepare assembly of matrix - reserve space in underlying triplet object
    void prepareMatAssembly( unsigned numElements, unsigned elMatSize );

    //! insert submatrix to system matrix - into triplet object
    void insertToMatrix( const SubMat_ & eStiff, 
                         const VecUint_ & rowIndices,
                         const VecUint_ & colIndices );

    //! Finalize assembly of matrix triplet
    void finishMatAssembly( );

    //! insert subvector to system rhs vector
    void insertToRhs( const SubVec_ & eForce, const VecUint_ & dofIndices );

    //! insert a scalar to RHS
    void insertValueToRhs( const double & entry, const unsigned & dofIndex );

    //! Outputs
    void writeMatrix(   std::ostream & out ) { triplet_.write( out ); }

    //! Apply constraints
    //! Currently, scalar is not used in the function and it is kept only for compatibility with constraint functors.
    //! Proper scalar is enforced within the BDDC solver.
    void applyConstraints( ConstraintVec_ & constraints, const double factor, const double scalar );

    //! Set any DOF to zero
    void fixDOF( const unsigned index, const double scalar );

    //! Solve the system
    void solveSystem( double tol = 1.e-7,                        //!< tolerance on relative residual ||res||/||rhs||
                      int  numLevels = 2,                        //!< number of levels
                      std::vector<int> *  numSubAtLevels = NULL, //!< number of subdomains at levels
                      int verboseLevel = 0,                      //!< level of verbosity of BDDCML library 
                                                                 //!< ( 0 - only fatal errors reported, 
                                                                 //!<   1 - mild output, 
                                                                 //!<   2 - detailed output )
                      int  maxIt = 1000,                         //!< maximum number of iterations
                      int  ndecrMax = 30                         //!< maximum number of iterations with non-decreasing residual 
                                                                 //!< ( used to stop diverging process )
                     );

    //! Get norm of right-hand side
    double normRhs( ) { return normRhs_; }

    //! Get norm of solution
    double normSol( ) { return normSol_; }

    //! Give number of iteration
    int giveNumIterations()   { return numIter_ ; }

    //! Give reason of convergence
    int giveConvergedReason() { return convergedReason_ ; }

    //! Give estimate of condition number
    double giveCondNumber() { return condNumber_ ; }

    //! Fill matrix with zeros
    void clearMatrix( ) { 
        // remember length
        unsigned length = triplet_.nnz( );
        triplet_.clear( ); 
        triplet_.prepareAssembly( length );
    }
    //! Fill RHS with zeros
    void clearRhs( )    { 
        std::fill( rhsVec_.begin( ), rhsVec_.end( ), 0. ); 
        normRhs_ = 0.; 
    }  

    //!  Blank arrays for boundary conditions
    void clearBC( ) { 
        std::fill( ifix_.begin(), ifix_.end(), 0  );
        std::fill( fixv_.begin(), fixv_.end(), 0. );
    }

    //! Fill solution with zeros
    void clearSol( )    
    {   std::fill( sol_.begin( ), sol_.end( ), 0. ); 
        normSol_ = 0.; 
    }

    //! Destroy matrix
    void destroyMatrix( ) { 
        triplet_.clear( ); 
    }
    //! Destroy RHS
    void destroyRhs( )    { rhsVec_.clear( ); }
    //! Destroy boundary conditions
    void destroyBC( )    { 
        ifix_.clear( ); 
        fixv_.clear( ); 
    }
    //! Destroy solution
    void destroySol( )    { sol_.clear( ); }

    //! Give access to the solution vector
    template< typename IT1, typename IT2 >
    void giveDofSolution( IT1 db, IT1 de, IT2 vb ) const;

    //! Give access to the solution vector
    void giveSolution( const VecUint_ & dofIndices, SubVec_ & result ) const;

private:
    const int  numDofs_;                   //!< number of total dofs
    const int  numDofsSub_;                //!< number of subdomain dofs
    const MatrixType matrixType_;          //!< type of matrix
    const MPI_Comm comm_;                  //!< communicator
    int        numSubLoc_;                 //!< number of subdomains on process

    int        nProc_;                     //!< number of processors in communicator
    int        rank_;                      //!< ID of processors in communicator

    int        nDim_;                      //!< space dimension
    int        meshDim_;                   //!< mesh dimension, may differ from space dimension - e.g. 2 for shells in 3D

    int        numElemSub_;                //!< number of elements in subdomain
    int        numElem_;                   //!< total number of elements
    int        numNodesSub_;               //!< number of nodes in subdomain
    int        numNodes_;                  //!< total number of nodes 

    bool       meshLoaded_;                //! safety state variable to warn if solver is called without loading subdomain mesh
    std::vector<int>     inet_;            //!< indices of nodes on elements (linearized connectivity) - numbering w.r.t. local subdomain numbering
    std::vector<int>     nnet_;            //!< numbers of nodes on elements 
                                           //!< (determines chunks of nodes for each element in inet), 
                                           //!< e.g. [ 8, 8, 8, ... ] for linear hexahedra
    std::vector<int>     nndf_;            //!< number of nodal degrees of freedom ( entry for each node ), 
                                           //!< e.g. [ 3, 3, 3, ... ] for 3D solids
    std::vector<double>  xyz_;             //!< linearized array of coordinates - all x values, all y values, all z values ( if present )
                                           //!< i.e. [ x_0 x_1 x_2 ... x_(numNodes - 1) y_0 y_1 y_2 ... y_(numNodes - 1) z_0 z_1 z_2 ... z_(numNodes - 1) ]
                                           
    std::vector<int>     isngn_;           //!< Indices of Subdomain Nodes in Global Numbering
    std::vector<int>     isvgvn_;          //!< Indices of Subdomain Variables in Global Variable Numbering ( i.e. dof mapping )
    std::vector<int>     isegn_;           //!< Indices of Subdomain Elements in Global Numbering

    typedef std::map<unsigned,unsigned> Global2LocalMap_; //! type for storage of global to local map
    Global2LocalMap_ global2LocalDofMap_;     //!< map from global dof indices to subdomain local indices


    std::vector<int>     ifix_;            //!< indices of fixed boundary conditions - ( 0 for free variable, 1 for constrained dof )
    std::vector<double>  fixv_;            //!< values of fixed boundary conditions at places of ifix_ - values outside fixed dofs are ignored

    corlib::Triplet<int,double>  triplet_; //!< matrix in triplet format
    bool                 isMatAssembled_;  //!< true if matrix is assembled

    std::vector<double>  rhsVec_;          //!< vector with RHS values restricted to subdomain, i.e. values are repeated at shared nodes
    double               normRhs_;         //!< norm of the global right-hand side

    std::vector<double>  sol_;             //!< distributed solution vector
    double               normSol_;         //!< norm of the global solution

    int                  numIter_;         //!< required number of iterations
    int                  convergedReason_; //!< reason of convergence:
                                           //!<   0 - converged relative residual, desirable state
                                           //!<  -1 - reached maximal number of iterations
                                           //!<  -2 - reached maximal number of iterations with non-decreasing residual
    double               condNumber_;      //!< estimated condition number of the preconditioned operator
                                           //!<  only provided for SPD problems solved by PCG ( -1.0 if not meaningful )

};

//-----------------------------------------------------------------------------

#include "SystemSolveBddc.ipp"

//-----------------------------------------------------------------------------

#endif
