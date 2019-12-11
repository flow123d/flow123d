#ifndef ARMOR_HH
#define ARMOR_HH

//#define ARMA_DONT_USE_WRAPPER
//#define ARMA_NO_DEBUG
#include <armadillo>
#include <array>
#include "system/asserts.hh"

namespace Armor {

template <class Type, uint nRows, uint nCols>
class Mat
{
private:
    Type * const data;
    typedef typename arma::Mat<Type>::template fixed<nRows,nCols> ArmaType;
public:
    Mat()
    : data(new Type[nRows * nCols])
    {}

    Mat(Type * const mem)
    : data(mem)
    {}

    inline Mat(const ArmaType & arma) : data(new Type[nRows * nCols]) {
        for (uint i = 0; i < nRows * nCols; ++i) {
            data[i] = arma[i];
        }
    }

    inline Mat(std::initializer_list<std::initializer_list<Type>> list) : data(new Type[nRows * nCols]) {
        const auto * listIt = list.begin();
        const Type * it;
        for (uint i = 0; i < nRows; ++i) {
            it = (listIt + i)->begin();
            for (uint j = 0; j < nCols; ++j) {
                data[i+j*nRows] = *(it + j);
            }
        }
    }

    inline Mat(std::initializer_list<Type> list) : data(new Type[nRows * nCols]) {
        const Type * it = list.begin();
        for (uint i = 0; i < nRows * nCols; ++i) {
            data[i] = *(it + i);
        }
    }

    inline Mat(const Armor::Mat<Type, nRows, nCols> & other) 
    : data(other.data)
    {
        for (uint i = 0; i < nRows * nCols; ++i) {
            data[i] = other[i];
        }
    }

    inline const Type * begin() const {
        return data;
    }
    inline Type * begin() {
        return data;
    }
    inline const Type * end() const {
        return begin() + nCols * nRows;
    }
    inline Type * end() {
        return begin() + nCols * nRows;
    }
    inline uint size() const {
        return nRows * nCols;
    }
    inline Type * memptr() {
        return begin();
    }
    inline const Type & operator[](uint index) const {
        return data[index];
    }
    inline Type & operator[](uint index) {
        return data[index];
    }
    inline const Type & operator()(uint row, uint col) const {
        return data[row+col*nRows];
    }
    inline Type & operator()(uint row, uint col) {
        return data[row+col*nRows];
    }
    inline ArmaType arma() const {
        return ArmaType(begin());
    }
    inline const Mat<Type, nRows, nCols> & operator=(const Mat<Type, nRows, nCols> & other) {
        for (uint i = 0; i < nRows * nCols; ++i) {
            data[i] = other[i];
        }
        return *this;
    }
    inline const Mat<Type, nRows, nCols> & operator=(const ArmaType & other) {        
        for (uint i = 0; i < nRows * nCols; ++i) {
            data[i] = other[i];
        }
        return *this;
    }
    inline const Mat<Type, nRows, nCols> & operator=(std::initializer_list<std::initializer_list<Type>> list) {
        const auto * listIt = list.begin();
        const Type * it;
        for (uint i = 0; i < nRows; ++i) {
            it = (listIt + i)->begin();
            for (uint j = 0; j < nCols; ++j) {
                data[i+j*nRows] = *(it + j);
            }
        }
        return *this;
    }
    inline const Mat<Type, nRows, nCols> & operator=(std::initializer_list<Type> list) {
        const Type * it = list.begin();
        for (uint i = 0; i < nRows * nCols; ++i) {
            data[i] = *(it + i);
        }
        return *this;
    }
    inline bool operator==(const ArmaType & other) {
        const Type * first1 = begin();
        const Type * last1 = end();
        const Type * first2 = other.begin();
        for (; first1 != last1; ++first1, ++first2) {
            if (*first1 != *first2) {
                return false;
            }
        }
        return true;
    }
    inline bool operator==(const Mat<Type, nRows, nCols> & other) {
        const Type * first1 = begin();
        const Type * last1 = end();
        const Type * first2 = other.begin();
        for (; first1 != last1; ++first1, ++first2) {
            if (*first1 != *first2) {
                return false;
            }
        }
        return true;
    }
};


/**
 * Array of Armor::Mat with given shape. Provides contiguous storage for the data and access to the array elements.
 * The shape of the matrices is specified at run time, so the class Array is independent of additional template parameters.
 * However, to access the array elements, one must use the templated method get().
 */
template<class Type>
class Array {
public:
    /**
     * Construct array of Armor matrices.
     * @param nv    Number of matrices in the array.
     * @param nr    Number of rows in each matrix.
     * @param nc    Number of columns in each matrix.
     */
    Array(uint nv, uint nr, uint nc = 1)
    : nRows(nr),
      nCols(nc),
      nVals(nv),
      data(nVals*nRows*nCols, 0)
    {}
    
    /**
     * Change number of elements in the array, while keeping the shape of arrays.
     * @param size  New size of array.
     */
    inline void resize(uint nv)
    {
    	nVals = nv;
        data.resize(nVals*nRows*nCols);
    }
    
    /**
     * Return number of matrices.
     */
    inline uint n_vals() const
    {
        return nVals;
    }

    inline uint n_rows() const
    {
        return nRows;
    }

    inline uint n_cols() const
    {
        return nCols;
    }

    /**
     * Insert new matrix to the end of the array.
     * @param p  Vector of values for the new matrix.
     */
    inline void push_back(const std::vector<Type> &p)
    {
        ASSERT_DBG( p.size() == nRows*nCols );
        nVals++;
        for (auto i : p) data.push_back( i );
    }
    
    /**
     * Return matrix at given position in array. The returned object is a Armor::Mat
     * pointing to the respective data block in the Array's storage.
     * @param i  Index of matrix.
     *
     * TODO: Should be renamed to item(), but we have compilation problem in Field::loc_point_value
     */
    template<uint nr, uint nc = 1>
    inline Mat<Type,nr,nc> get(uint i) const
    {
        ASSERT_DBG( (nr == nRows) && (nc == nCols) );
        return Mat<Type,nr,nc>( (Type *)data.data() + i*nRows*nCols );
    }
    
    /**
     * Return matrix at given position in array. The returned object is a Armor::Mat
     * pointing to the respective data block in the Array's storage.
     * @param i  Index of matrix.
     *
     * TODO: Should be renamed to item(), but we have compilation problem in Field::loc_point_value
     */
    template<uint nr, uint nc = 1>
    inline Mat<Type,nr,nc> get(uint i)
    {
        ASSERT_DBG( (nr == nRows) && (nc == nCols) );
        return Mat<Type,nr,nc>( (Type*)(data.data()) + i*nRows*nCols );
    }

    /**
     * Return armadillo matrix at given position in array.
     * @param i  Index of matrix.
     */
    inline arma::mat arma_mat(uint i) const
    {
    	return arma::vec( (Type*)(data.data()) + i*nRows*nCols, nRows*nCols );
    }

    /**
     * Return armadillo vector at given position in array.
     * Warning! Method can be used only if nCols == 1.
     * @param i  Index of matrix.
     */
    inline arma::vec arma_vec(uint i) const
    {
        ASSERT_EQ_DBG(nCols, 1);
    	return arma::vec( (Type*)(data.data()) + i*nRows, nRows );
    }

private:
    uint nRows;
    uint nCols;
    uint nVals;
    std::vector<Type> data;
};


template <class Type, uint nRows, uint nCols>
inline Type dot(const Mat<Type, nRows, nCols> & a, const Mat<Type, nRows, nCols> & b) {
    return arma::dot(a.arma(), b.arma());
}

template <class Type, uint nRows, uint nCols>
inline typename arma::Mat<Type>::template fixed<nRows,nCols> operator+(const Mat<Type, nRows, nCols> & a, const Mat<Type, nRows, nCols> & b) {
    return a.arma() + b.arma();
}

template <class Type, uint nRows, uint nCols>
inline typename arma::Mat<Type>::template fixed<nRows,nCols> operator-(const Mat<Type, nRows, nCols> & a, const Mat<Type, nRows, nCols> & b) {
    return a.arma() - b.arma();
}

template <class Type, uint resRows, uint commonDimension, uint resCols>
inline typename arma::Mat<Type>::template fixed<resRows,resCols> operator*(const Mat<Type, resRows, commonDimension> & a, const Mat<Type, commonDimension, resCols> & b) {
    return a.arma() * b.arma();
}

template <class Type, uint resRows, uint commonDimension, uint resCols>
inline typename arma::Mat<Type>::template fixed<resRows,resCols> operator*(const Mat<Type, resRows, commonDimension> & a,
                                                                           const typename arma::Mat<Type>::template fixed<commonDimension, resCols> & b) {
    return a.arma() * b;
}

template <class Type, uint nRows, uint nCols>
inline typename arma::Mat<Type>::template fixed<nRows,nCols> operator%(const Mat<Type, nRows, nCols> & a, const Mat<Type, nRows, nCols> & b) {
    return a.arma() % b.arma();
}

template <class Type, uint nRows, uint nCols>
inline typename arma::Mat<Type>::template fixed<nRows,nCols> operator*(Type number, const Mat<Type, nRows, nCols> & a) {
    return number * a.arma();
}

template <class Type, uint nRows, uint nCols>
inline typename arma::Mat<Type>::template fixed<nRows,nCols> operator/(const Mat<Type, nRows, nCols> & a, Type number) {
    return a.arma() / number;
}

template <class Type, uint nRows, uint nCols>
inline typename arma::Mat<Type>::template fixed<nCols,nRows> trans(const Mat<Type, nRows, nCols> & mat) {
    return trans(mat.arma());
}

template <uint N>
using vec = Mat<double, N, 1>;

template <uint N, uint M>
using mat = Mat<double, N, M>;

using array = Array<double>;

}

#endif
