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
public:
    typedef typename arma::Mat<Type>::template fixed<nRows, nCols> ArmaType;

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

    inline Mat(ArmaType arma_mat)
    : data(arma_mat.memptr())
    {}



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

    inline operator ArmaType() const {return arma();}

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
 * Check for close vectors or matrices.
 * |a - b| < a_tol  or  |a-b| < r_tol * |a|
 *
 * Note: Similar function arma::approx_equal.
 */
template <class Type, uint nRows, uint nCols>
inline bool is_close(Mat<Type, nRows, nCols> a, Mat<Type, nRows, nCols> b, double a_tol = 1e-10, double r_tol = 1e-6) {
    double a_diff = arma::norm(a - b, 1);
    if (a_diff < a_tol) return true;
    if (a_diff < r_tol * arma::norm(a, 1)) return true;
    return false;
}







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
    Array(uint nr, uint nc = 1, uint size = 0)
    : n_rows_(nr),
      n_cols_(nc),
      size_(0),
      reserved_(size),
      data_(new Type[nr * nc * size])
    {}
    
    /**
     * Drop all data and allocate new space of given size.
     * Change number of elements in the array, while keeping the shape of arrays.
     * @param size  New size of array.
     */
    void reinit(uint size) {
        delete data_;
        reserved_ = size;
        size_ = 0;
        data_ = new Type[n_rows_ * n_cols_ * reserved_];
    


    /**
     * Resize active part of the allocated space.
     */
    void resize(uint size) {
        ASSERT_LE_DBG(size, reserved_);
        size_ = size;
    }

    /**
     * Get size of active space.
     */
    inline unsigned int size() const {
        return size_;
    }

    /**
     * Increase active space by 1 and store given Mat value to the end of the active space.
     */
    template<uint nr, uint nc = 1>
    inline void append(Mat<Type,nr,nc> item)
    {
        ASSERT_LE_DBG(size_, reserved_);
        get<nr,nc>(size_) = item;
        size_ += 1;
    }
    
    /**
     * Return matrix at given position in array. The returned object is a Armor::Mat
     * pointing to the respective data_ block in the Array's storage.
     * One can assign to the Armor::Mat which performs postponed evaluation and storing the result to the array.
     *
     * @param i  Index of the matrix.
     * TODO: Should be renamed to item(), but we have compilation problem in Field::loc_point_value
     */
    template<uint nr, uint nc = 1>
    inline Mat<Type,nr,nc> get(uint mat_index) const
    {
        ASSERT_DBG( (nr == n_rows_) && (nc == n_cols_) );
        return Mat<Type,nr,nc>( data_ + mat_index * n_rows_ * n_cols_ );
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
    inline uint space_() { return n_rows_ * n_cols_ * reserved_; }
    uint n_rows_;
    uint n_cols_;
    uint size_;
    uint reserved_;
    Type * data_;
};


template <class Type, uint nRows, uint nCols>
inline Type dot(const Mat<Type, nRows, nCols> & a, const Mat<Type, nRows, nCols> & b) {
    return arma::dot(a.arma(), b.arma());
}

template <class Type, uint nRows, uint nCols>
inline typename Mat<Type, nRows, nCols>::ArmaType operator+(const Mat<Type, nRows, nCols> & a, const Mat<Type, nRows, nCols> & b) {
    return a.arma() + b.arma();
}

template <class Type, uint nRows, uint nCols>
inline typename Mat<Type, nRows, nCols>::ArmaType operator-(const Mat<Type, nRows, nCols> & a, const Mat<Type, nRows, nCols> & b) {
    return a.arma() - b.arma();
}

template <class Type, uint nRows, uint nCols>
inline typename Mat<Type, nRows, nCols>::ArmaType operator*(const Mat<Type, nRows, nCols> & a, const Mat<Type, nRows, nCols> & b) {
    return a.arma() * b.arma();
}

template <class Type, uint nRows, uint nCols>
inline typename Mat<Type, nRows, nCols>::ArmaType operator%(const Mat<Type, nRows, nCols> & a, const Mat<Type, nRows, nCols> & b) {
    return a.arma() % b.arma();
}

template <class Type, uint nRows, uint nCols>
inline typename Mat<Type, nRows, nCols>::ArmaType operator*(Type number, const Mat<Type, nRows, nCols> & a) {
    return number * a.arma();
}

template <class Type, uint nRows, uint nCols>
inline typename Mat<Type, nRows, nCols>::ArmaType operator/(const Mat<Type, nRows, nCols> & a, Type number) {
    return a.arma() / number;
}


template <uint N>
using vec = Mat<double, N, 1>;

template <uint N, uint M>
using mat = Mat<double, N, M>;



using array = Array<double>;

}

#endif
