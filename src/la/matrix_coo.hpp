#ifndef la_matrix_coo_h
#define la_matrix_coo_h

#define ALLOWED_MEM_OVERHEAD 2.e+6 // entries, this corresponds to 32 MB of 
                                   //  overhead in <int, int, double> matrixCoo
//------------------------------------------------------------------------------
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <numeric>

//------------------------------------------------------------------------------
namespace la{
    template< typename INDT, typename VALT >
    class MatrixCoo;
}

//------------------------------------------------------------------------------
/** \brief Temporary storage for a sparse matrix 
 *  \details This object stores and handles matrix in coordinate format, i.e., 
 *  <row_index,col_index,value>-combinations for each matrix entry of a sparse matrix.
 *  Functionalities are:
 *  - preparing assembly of elements - allocating proper space to all elements 
 *    in unassembled form
 *  - inserting an element - it is user's responsibility not to insert zeros
 *  - finishing assembly of elements - this is the core routine which 
 *    - sorts triples, merges entries with the same pair of indices, 
 *      truncates the array to assembled size
 *  - 'zeroing' a row and setting its diagonal entry to a given value
 *  - query of number of entries - this counts also zeros if they were inserted
 *  - get inside/outide-a-range counts of elements - for PETSc preallocation
 *  - filling the arrays for a row-compressed data storage for external
 *    solver libraries
 *  - computing inf norm
 *  - cleanup and writing
 */
template< typename INDT, typename VALT >
class la::MatrixCoo
{
public:
    typedef INDT         IndexType;
    typedef VALT         ValueType;

public:
    MatrixCoo( unsigned nRows = 0, unsigned nCols = 0 ) 
      : nRows_( nRows ),
        nCols_( nCols ),
        sorted_ ( false )
    { 
        sortedSize_ = 0; 
    } 

    ~MatrixCoo( ) { 
        matrixCooVec_.clear();
    }

    //! Return number of non-zeros
    IndexType nnz() const { return matrixCooVec_.size(); }

    //! Prepare matrix for assembly - allocate buffer for UNASSEMBLED entries
    //! serial copies of element matrices to memory
    //! length = numElements * size( elementMatrix )
    void prepareAssembly( const IndexType length )
    {
        matrixCooVec_.reserve( length );

        return;
    }

    //! Insert a value to unassembled vector 
    void insert( const IndexType i, const IndexType j, const ValueType value )
    {
        // combine index pair
        IndexPair_ ip( i, j );

        // combine triples
        Triple_ trip( ip, value );

        // insert it to buffer in any case
        matrixCooVec_.push_back( trip );

        sorted_ = false;

        // resort if number of new unsorted entries excess a critical value
        if ( (matrixCooVec_.size() - sortedSize_) > ALLOWED_MEM_OVERHEAD ) {
            this -> finishAssembly();
        }

        return;
    }
    //! Finalize matrix assembly - build the assembled input and truncate matrix
    void finishAssembly( )
    {
        // sort entries to make sure
        std::sort( matrixCooVec_.begin( ), matrixCooVec_.end( ), TripleLessThan_( ) );

        // combine the sorted matrix into one vector
        typename MatrixCooVec_::const_iterator itGlobale = matrixCooVec_.end();
    
        typename MatrixCooVec_::const_iterator itSelectb = matrixCooVec_.begin();
        typename MatrixCooVec_::const_iterator itSelecte = itSelectb;

        Triple_ tmpTriple;

        typename MatrixCooVec_::iterator itStore = matrixCooVec_.begin();

        // for each pair of indices
        while (itSelectb != itGlobale) {

            // find any change in indices
            itSelecte = std::upper_bound( itSelectb, itGlobale, *itSelectb, TripleLessThan_( ) );

            // sum values in the range
            tmpTriple = combineTriples_ ( itSelectb, itSelecte );

            // store new  
            *itStore = tmpTriple;
            itStore++;

            // shift ranges
            itSelectb = itSelecte;
        }

        // determine number of nonzeros in assembled matrix
        unsigned nnz = std::distance( matrixCooVec_.begin( ),  itStore);

        // truncate the vector from the end
        matrixCooVec_.resize( nnz );

        // mark the amount of sorted entries
        sortedSize_ = nnz;
        
        sorted_ = true;

        return;
    }

    //! Determine number of nonzeros within a block [low,high) - used for PETSc preallocation
    void getIntervalCounts( IndexType low, IndexType high, 
                            std::vector<IndexType> & nnzRowsIn, std::vector<IndexType> & nnzRowsOut) 
    {
        // sort entries to make sure
        if ( not sorted_ ) {
            this -> finishAssembly( );
        }

        // run through all entries
        typename MatrixCooVec_::const_iterator iter  = matrixCooVec_.begin();
        typename MatrixCooVec_::const_iterator itere = matrixCooVec_.end();

        for ( ; iter != itere; ++iter ) {

            IndexType row = (*iter).first.first;
            IndexType col = (*iter).first.second;

            // check that it belongs to the interval
            if ( row >= low and row < high ) {

                IndexType rowLoc = row - low;
                // it is on my row, now where to put it
                
                if ( col >= low and col < high ) {
                    // it is inside the diagonal block
                    nnzRowsIn[rowLoc]++;
                }
                else {
                    // it is outside the diagonal block
                    nnzRowsOut[rowLoc]++;
                }
            }
        }
    }

    //! Get specific entry
    void getEntry( IndexType index, IndexType & row, IndexType & col, ValueType & val )
    {
        row = matrixCooVec_[ index ].first.first;
        col = matrixCooVec_[ index ].first.second;
        val = matrixCooVec_[ index ].second;
    }

    //! Zero a row and set diagonal to scalar 
    void zeroRow( const IndexType row, const ValueType factor )
    {

        // sort entries to make sure
        if ( not sorted_ ) {
            this -> finishAssembly( );
        }

        // erasing iterator
        typename MatrixCooVec_::iterator  tripEraseIter;
        typename MatrixCooVec_::iterator  tripEraseIterE;

        // find constrained row in matrix and set iterator there
        tripEraseIter = std::lower_bound ( matrixCooVec_.begin(), matrixCooVec_.end(), 
                                           row, TripleRowIndexLessThan_( ) );

        // find last entry at the row
        tripEraseIterE = std::upper_bound( tripEraseIter, matrixCooVec_.end(), 
                                           *tripEraseIter,  TripleRowLessThan_( ) );

        // find diagonal
        typename MatrixCooVec_::iterator tripEraseIterDiag;
        tripEraseIterDiag = std::find_if( tripEraseIter, tripEraseIterE, TripleDiagEntry_( ) ); 

        // set row entries to zero
        for ( ; tripEraseIter != tripEraseIterE; ++tripEraseIter ) {
            (*tripEraseIter).second = static_cast< ValueType >( 0.0 );
        }

        // put value from constraint on diagonal
        Triple_ tmpTrip = std::make_pair( std::make_pair ( row, row ), factor );

        if ( tripEraseIterDiag != tripEraseIterE ) {
            // diagonal entry already present
            *tripEraseIterDiag = tmpTrip;
        }
        else {
            // nothing on diagonal, add it
            matrixCooVec_.push_back( tmpTrip );
            sorted_ = false;
        }
    }

    //! Fill arrays with column indices, row pointers and values for
    //! row-compressed data storage
    void setCompressedRowArrays( std::vector<IndexType> & colIndices, 
                                 std::vector<IndexType> & rowPtrs,
                                 std::vector<ValueType> & values ) 
    {

        // sort entries to make sure
        if ( not sorted_ ) {
            this -> finishAssembly( );
        }

        // iterate over all elements and fill the nonzero, 
        // column index, and row pointer arrys
        typename MatrixCooVec_::const_iterator tvIter = matrixCooVec_.begin();
        typename MatrixCooVec_::const_iterator tvEnd  = matrixCooVec_.end();

        // reserve space 
        colIndices.reserve(matrixCooVec_.size());
        values.reserve(matrixCooVec_.size());

        unsigned nr = this -> deduceNumRows_( );
        rowPtrs.reserve( nr + 1 );

        int previousRow = -1;
        IndexType index  = 0;

        for ( ; tvIter != tvEnd; ++tvIter ) {
            const IndexType row =  (*tvIter).first.first;
            const IndexType col =  (*tvIter).first.second;
            
            values.push_back( tvIter -> second );
            colIndices.push_back( col );
            if ( row == previousRow ) index++;
            else rowPtrs.push_back( index++ );

            previousRow = row;
        }

        rowPtrs.push_back( values.size() );

        return;       
    }

    //! Fill independent arrays with column indices, row pointers and values
    void extractArrays( std::vector<IndexType> & rowIndices,
                        std::vector<IndexType> & colIndices,
                        std::vector<ValueType> & values,
                        bool onlyUpperTriangle = false ) 
    {

        // sort entries to make sure
        if ( not sorted_ ) {
            this -> finishAssembly( );
        }

        // iterate over all elements and fill the nonzero, 
        // column index, and row index arrys
        typename MatrixCooVec_::const_iterator tvIter = matrixCooVec_.begin();
        typename MatrixCooVec_::const_iterator tvEnd  = matrixCooVec_.end();

        // reserve space 
        colIndices.reserve(matrixCooVec_.size());
        rowIndices.reserve(matrixCooVec_.size());
        values.reserve(matrixCooVec_.size());

        for ( ; tvIter != tvEnd; ++tvIter ) {
            const IndexType row = (*tvIter).first.first  ;
            const IndexType col = (*tvIter).first.second ;
            
            if ( not onlyUpperTriangle or col >= row ) {
                values.push_back( (*tvIter).second );
                rowIndices.push_back( row );
                colIndices.push_back( col );
            }
        }
        return;       
    }

    //! destroy all contents
    void clear() 
    {
        matrixCooVec_.clear();
        // force deallocation of the vector
        MatrixCooVec_().swap( matrixCooVec_ ); 
        return;
    }

    //! erase all contents
    void erase() 
    {
        // erasing iterator
        typename MatrixCooVec_::iterator  tripEraseIter  = matrixCooVec_.begin( );
        typename MatrixCooVec_::iterator  tripEraseIterE = matrixCooVec_.end( );

        // loop to erase values
        for ( ; tripEraseIter != tripEraseIterE; ++tripEraseIter ) {
            (*tripEraseIter).second = static_cast< ValueType >(0.0);
        }

        return;
    }

    //! debug-output
    std::ostream & write( std::ostream & out )
    {
        // sort entries to make sure
        if ( not sorted_ ) {
            this -> finishAssembly( );
        }
        typename MatrixCooVec_::const_iterator tvIter = matrixCooVec_.begin();
        typename MatrixCooVec_::const_iterator tvEnd  = matrixCooVec_.end();
        for ( ; tvIter != tvEnd; ++tvIter ) {
            const IndexType row = (*tvIter).first.first;
            const IndexType col = (*tvIter).first.second;

            out << row << "  " << col << "  " 
                << std::setprecision(12)
                << (*tvIter).second << std::endl;

        }
        return out;
    }

    //----------------------------------------------------------------------------
    //! Compute maximal row sum (infinity-norm) of a matrix in coordinate format. 
    double infNorm() 
    {
        // sort entries to make sure
        if ( not sorted_ ) {
            this -> finishAssembly( );
        }
        unsigned nr = this -> deduceNumRows_( );
        std::vector<double> rowSums( nr, 0.0 );
        typename MatrixCooVec_::const_iterator  iter = matrixCooVec_.begin();
        typename MatrixCooVec_::const_iterator  last = matrixCooVec_.end();  
        for ( ; iter != last; ++iter ) {
            const IndexType iRow = (*iter).first.first;
            rowSums[ iRow ] += std::fabs( (*iter).second );
        }
        
        return *std::max_element( rowSums.begin(), rowSums.end() );
    }

    //--------------------------------------------------------------------------
    //! Return Frobenius norm of matrix \f$||a_{ij}||_F=\sqrt{\sum_i\sum_j |a_{ij}|^2}\f$
    ValueType frobeniusNorm() 
    {
        // sort entries to make sure
        if ( not sorted_ ) {
            this -> finishAssembly( );
        }

        ValueType sumOfSquares = 0.0;

        typename MatrixCooVec_::const_iterator  iter = matrixCooVec_.begin();
        typename MatrixCooVec_::const_iterator  last = matrixCooVec_.end();  
        for ( ; iter != last; ++iter ) {

            ValueType val = (*iter).second;
            sumOfSquares += val * val;
        }

        return std::sqrt( sumOfSquares );
    }

    //--------------------------------------------------------------------------
    //! Get scalar suitable to be put on diagonal while fixing BC
    //!  This is the mean value of extremal entries on diagonal - the aim of this guess is to 
    //!  be within the spectra of the original matrix - this value immediately corresponds to an 
    //!  eigenvalue 
    double getDiagScalar( )
    {

        ValueType val;
        // initialize maximum and minimum
        if ( matrixCooVec_.size() > 0 ) {
            val = std::fabs( matrixCooVec_[0].second ); 
        }
        else {
            val = 0.0;
        }
        ValueType maxVal = val;
        ValueType minVal = val;
        IndexType row, col;
        for ( IndexType i = 1; i < static_cast<IndexType>( matrixCooVec_.size() ); i++ ) {

            row = matrixCooVec_[i].first.first;
            col = matrixCooVec_[i].first.second;

            if ( row == col ) {
                val = std::fabs( matrixCooVec_[i].second ); 

                // update maximum
                if ( val > maxVal ) {
                    maxVal = val;
                }

                // update minimum
                if ( val < minVal ) {
                    minVal = val;
                }
            }
        }

        // recommend value in between extremal absolute values along diagonal
        ValueType scalar = ( minVal + maxVal ) / 2.0;

        return scalar;

    }

private:
    typedef std::pair<IndexType, IndexType>    IndexPair_;   //!< pair of indices
    typedef std::pair<IndexPair_,ValueType>    Triple_;      //!< triple
    typedef std::vector<Triple_>               MatrixCooVec_;//!< vector of triples
    MatrixCooVec_                              matrixCooVec_;//!< matrix in coordinate format <i,j, values>

    IndexType                                  nRows_;       //!< number of rows of matrix
    IndexType                                  nCols_;       //!< number of columns of matrix

    unsigned                                   sortedSize_;  //!< number of sorted entries

    bool                                       sorted_;      //!< state variable - is matrix sorted?

    //--------------------------------------------------------------------------
    // private functions
    
    //! deduce number of rows - if it was given by user, use it, otherwise, 
    //!                         deduce it from the index of the last entry
    IndexType deduceNumRows_( )
    {
        IndexType nr;

        if ( nRows_ > 0 ) { 
            // number of rows was given by user
            nr = nRows_;
        }
        else { 
            // deduce number of rows from  sorted matrix
            nr = (*matrixCooVec_.rbegin()).first.first + 1;
        }

        return nr;
    }

    
    //! returns true if the index pair of the first entry is smaller than of the second entry
    struct TripleLessThan_ : std::binary_function< Triple_, Triple_, bool > 
    {
        bool operator() ( Triple_ firstTriple, Triple_ secondTriple) {
        
            IndexPair_ firstTripleIndxPair  = firstTriple.first;
            IndexPair_ secondTripleIndxPair = secondTriple.first;

            return (firstTripleIndxPair < secondTripleIndxPair);
        }
    };

    //! returns true if the first entry has smaller ROW index than the second one
    struct TripleRowLessThan_ : std::binary_function< Triple_, Triple_, bool > 
    {
        bool operator() ( Triple_ firstTriple, Triple_ secondTriple) {
        
            IndexType firstRowIndex  = firstTriple.first.first;
            IndexType secondRowIndex = secondTriple.first.first;

            return (firstRowIndex < secondRowIndex);
        }
    };

    //! returns true if the row of an entry is smaller than a prescribed index
    struct TripleRowIndexLessThan_ 
    {
        bool operator() ( Triple_ triple, IndexType row ) {
        
            IndexType tripleRowIndex  = triple.first.first;

            return (tripleRowIndex < row);
        }
    };

    //! returns true if the entry is on diagonal, i.e. row == column
    struct TripleDiagEntry_ 
    {
        bool operator() ( Triple_ triple ) {
        
            IndexType tripleRowIndex  = triple.first.first;
            IndexType tripleColIndex  = triple.first.second;

            return (tripleRowIndex == tripleColIndex);
        }
    };

    //! merge range of triples by addition
    Triple_ combineTriples_ ( typename MatrixCooVec_::const_iterator itb, typename MatrixCooVec_::const_iterator ite ) {

        IndexPair_ tmpPair = (*itb).first;
        ValueType val = 0.0;

        for ( ; itb != ite; ++itb) {
            val += (*itb).second;
        }

        Triple_ combinedTriple = std::make_pair(tmpPair, val); 
    
        return combinedTriple;

    }
};

#endif
