//! @author Jakub Sistek, Thomas Rueberg
//! @date   6/2011

#ifndef la_triplet_h
#define la_triplet_h

#define ALLOWED_MEM_OVERHEAD 2.e+6 // entries, this corresponds to 32 MB of 
                                   //  overhead in <int, int, double> triplet
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
    class Triplet;
}

//------------------------------------------------------------------------------
/** \brief Temporary storage for a sparse matrix 
 *  \details This object stores and handles the so-called triplets, i.e., 
 *  <row_index,col_index,value>-combinations for each matrix entry of a sparse matrix.
 *  Functionalities are:
 *  - preparing assembly of elements - allocating proper space to all elements 
 *    in unassembled form
 *  - inserting an element - it is user's responsibility not to insert zeros
 *  - finishing assembly of elements - this is the core routine which 
 *    - sorts triplets, merges entries with the same pair of indices, 
 *      truncates the array to assembled size
 *  - 'zeroing' a row and setting its diagonal entry to a given value
 *  - query of number of entries - this counts also zeros if they were inserted
 *  - get inside/outide-a-range counts of elements - for PETSc preallocation
 *  - filling the arrays for a row-compressed data storage for external
 *    solver libraries
 *  - computing inf norm
 *  - cleanup and writing
 *
 *  Triplet first allocates memory for unassembled matrix, which is then 
 *  sorted once and triplets with same coordinates are merged.
 *  This is generally faster then inserting elements to the proper places 
 *  each time, but can require considerably more memory,
 *  needed at one time of this computation. 
 *  This memory is released after the assembly. 
 *  If memory is not a concern, this object should be preferred. Should this code 
 *  crash on allocation, the other approach does the job instead.
 *  The approach with putting entries directly to proper positions was implemented 
 *  in an eariler version of Triplet.hpp in this repository ( until revision 2877 ).
 *
 *  The object supports the interface of previous version. However, for optimal 
 *  performance, one should use the prepareAssembly() and finishAssembly() routines
 *  around insert().
 */
template< typename INDT, typename VALT >
class la::Triplet
{
public:
    typedef INDT         IndexType;
    typedef VALT         ValueType;

public:
    Triplet( unsigned nRows = 0, unsigned nCols = 0 ) 
      : nRows_( nRows ),
        nCols_( nCols ),
        sorted_ ( false )
    { 
        sortedSize_ = 0; 
    } 

    ~Triplet( ) { 
        tripletVec_.clear();
    }

    //! Return number of non-zeros
    IndexType nnz() const { return tripletVec_.size(); }

    //! Prepare triplet for assembly - allocate buffer for UNASSEMBLED entries
    //! serial copies of element matrices to memory
    //! length = numElements * size( elementMatrix )
    void prepareAssembly( const IndexType length )
    {
        tripletVec_.reserve( length );

        return;
    }

    //! Insert a value to unassembled vector 
    void insert( const IndexType i, const IndexType j, const ValueType value )
    {
        // combine index pair
        IndexPair_ ip( i, j );

        // combine triplet
        Triplet_ trip( ip, value );

        // insert it to buffer in any case
        tripletVec_.push_back( trip );

        sorted_ = false;

        // resort if number of new unsorted entries excess a critical value
        if ( (tripletVec_.size() - sortedSize_) > ALLOWED_MEM_OVERHEAD ) {
            //std::cout << " Performing resorting due to the critical value ..." << std::endl;
            this -> finishAssembly();
        }

        return;
    }
    //! Finalize triplet assembly - build the assembled input and truncate triplet
    void finishAssembly( )
    {
        // sort entries to make sure
        std::sort( tripletVec_.begin( ), tripletVec_.end( ), TripletLessThan_( ) );

        // combine the sorted triplet into one vector
        typename TripletVec_::const_iterator itGlobale = tripletVec_.end();
    
        typename TripletVec_::const_iterator itSelectb = tripletVec_.begin();
        typename TripletVec_::const_iterator itSelecte = itSelectb;

        Triplet_ tmpTriplet;

        typename TripletVec_::iterator itStore = tripletVec_.begin();

        // for each pair of indices
        while (itSelectb != itGlobale) {

            // find any change in indices
            itSelecte = std::upper_bound( itSelectb, itGlobale, *itSelectb, TripletLessThan_( ) );

            // sum values in the range
            tmpTriplet = combineTriplets_ ( itSelectb, itSelecte );

            // store new  
            *itStore = tmpTriplet;
            itStore++;

            // shift ranges
            itSelectb = itSelecte;
        }

        // determine number of nonzeros in assembled matrix
        unsigned nnz = std::distance( tripletVec_.begin( ),  itStore);

        // truncate the vector from the end
        tripletVec_.resize( nnz );

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
        typename TripletVec_::const_iterator iter  = tripletVec_.begin();
        typename TripletVec_::const_iterator itere = tripletVec_.end();

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
        row = tripletVec_[ index ].first.first;
        col = tripletVec_[ index ].first.second;
        val = tripletVec_[ index ].second;
    }

    //! Zero a row and set diagonal to scalar 
    void zeroRow( const IndexType row, const ValueType factor )
    {

        // sort entries to make sure
        if ( not sorted_ ) {
            this -> finishAssembly( );
        }

        // erasing iterator
        typename TripletVec_::iterator  tripEraseIter;
        typename TripletVec_::iterator  tripEraseIterE;

        // find constrained row in triplets and set iterator there
        tripEraseIter = std::lower_bound ( tripletVec_.begin(), tripletVec_.end(), 
                                           row, TripletRowIndexLessThan_( ) );

        // find last entry at the row
        tripEraseIterE = std::upper_bound( tripEraseIter, tripletVec_.end(), 
                                           *tripEraseIter,  TripletRowLessThan_( ) );

        // find diagonal
        typename TripletVec_::iterator tripEraseIterDiag;
        tripEraseIterDiag = std::find_if( tripEraseIter, tripEraseIterE, TripletDiagEntry_( ) ); 

        // set row entries to zero
        for ( ; tripEraseIter != tripEraseIterE; ++tripEraseIter ) {
            (*tripEraseIter).second = static_cast< ValueType >( 0.0 );
        }

        // put value from constraint on diagonal
        Triplet_ tmpTrip = std::make_pair( std::make_pair ( row, row ), factor );

        if ( tripEraseIterDiag != tripEraseIterE ) {
            // diagonal entry already present
            *tripEraseIterDiag = tmpTrip;
        }
        else {
            // nothing on diagonal, add it
            tripletVec_.push_back( tmpTrip );
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
        typename TripletVec_::const_iterator tvIter = tripletVec_.begin();
        typename TripletVec_::const_iterator tvEnd  = tripletVec_.end();

        // reserve space 
        colIndices.reserve(tripletVec_.size());
        values.reserve(tripletVec_.size());

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
        typename TripletVec_::const_iterator tvIter = tripletVec_.begin();
        typename TripletVec_::const_iterator tvEnd  = tripletVec_.end();

        // reserve space 
        colIndices.reserve(tripletVec_.size());
        rowIndices.reserve(tripletVec_.size());
        values.reserve(tripletVec_.size());

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
        tripletVec_.clear();
        // force deallocation of the vector
        TripletVec_().swap( tripletVec_ ); 
        return;
    }

    //! erase all contents
    void erase() 
    {
        // erasing iterator
        typename TripletVec_::iterator  tripEraseIter  = tripletVec_.begin( );
        typename TripletVec_::iterator  tripEraseIterE = tripletVec_.end( );

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
        typename TripletVec_::const_iterator tvIter = tripletVec_.begin();
        typename TripletVec_::const_iterator tvEnd  = tripletVec_.end();
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
    //! Compute maximal row sum (infinity-norm) of a matrix in triplet format. 
    double infNorm() 
    {
        // sort entries to make sure
        if ( not sorted_ ) {
            this -> finishAssembly( );
        }
        unsigned nr = this -> deduceNumRows_( );
        std::vector<double> rowSums( nr, 0.0 );
        typename TripletVec_::const_iterator  iter = tripletVec_.begin();
        typename TripletVec_::const_iterator  last = tripletVec_.end();  
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

        typename TripletVec_::const_iterator  iter = tripletVec_.begin();
        typename TripletVec_::const_iterator  last = tripletVec_.end();  
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
        if ( tripletVec_.size() > 0 ) {
            val = std::fabs( tripletVec_[0].second ); 
        }
        else {
            val = 0.0;
        }
        ValueType maxVal = val;
        ValueType minVal = val;
        IndexType row, col;
        for ( IndexType i = 1; i < static_cast<IndexType>( tripletVec_.size() ); i++ ) {

            row = tripletVec_[i].first.first;
            col = tripletVec_[i].first.second;

            if ( row == col ) {
                val = std::fabs( tripletVec_[i].second ); 

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
    typedef std::pair<IndexPair_,ValueType>    Triplet_;     //!< triplet
    typedef std::vector<Triplet_>              TripletVec_;  //!< vector of triplets
    TripletVec_                                tripletVec_;  //!< matrix in triplet <i,j, values>

    IndexType                                  nRows_;       //!< number of rows of matrix
    IndexType                                  nCols_;       //!< number of columns of matrix

    unsigned                                   sortedSize_;  //!< number of sorted entries

    bool                                       sorted_;      //!< state variable - is triplet sorted?

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
            nr = (*tripletVec_.rbegin()).first.first + 1;
        }

        return nr;
    }

    
    //! returns true if the index pair of the first entry is smaller than of the second entry
    struct TripletLessThan_ : std::binary_function< Triplet_, Triplet_, bool > 
    {
        bool operator() ( Triplet_ firstTriplet, Triplet_ secondTriplet) {
        
            IndexPair_ firstTripletIndxPair  = firstTriplet.first;
            IndexPair_ secondTripletIndxPair = secondTriplet.first;

            return (firstTripletIndxPair < secondTripletIndxPair);
        }
    };

    //! returns true if the first entry has smaller ROW index than the second one
    struct TripletRowLessThan_ : std::binary_function< Triplet_, Triplet_, bool > 
    {
        bool operator() ( Triplet_ firstTriplet, Triplet_ secondTriplet) {
        
            IndexType firstRowIndex  = firstTriplet.first.first;
            IndexType secondRowIndex = secondTriplet.first.first;

            return (firstRowIndex < secondRowIndex);
        }
    };

    //! returns true if the row of an entry is smaller than a prescribed index
    struct TripletRowIndexLessThan_ 
    {
        bool operator() ( Triplet_ triplet, IndexType row ) {
        
            IndexType tripletRowIndex  = triplet.first.first;

            return (tripletRowIndex < row);
        }
    };

    //! returns true if the entry is on diagonal, i.e. row == column
    struct TripletDiagEntry_ 
    {
        bool operator() ( Triplet_ triplet ) {
        
            IndexType tripletRowIndex  = triplet.first.first;
            IndexType tripletColIndex  = triplet.first.second;

            return (tripletRowIndex == tripletColIndex);
        }
    };

    //! merge range of triplets by addition
    Triplet_ combineTriplets_ ( typename TripletVec_::const_iterator itb, typename TripletVec_::const_iterator ite ) {

        IndexPair_ tmpPair = (*itb).first;
        ValueType val = 0.0;

        for ( ; itb != ite; ++itb) {
            val += (*itb).second;
        }

        Triplet_ combinedTriplet = std::make_pair(tmpPair, val); 
    
        return combinedTriplet;

    }
};

#endif
