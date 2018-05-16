#ifndef ARMOR_H
#define ARMOR_H

//#define ARMA_DONT_USE_WRAPPER
//#define ARMA_NO_DEBUG
#include <armadillo>
#include <array>

namespace armor {

template <class Type, uint nRows, uint nCols>
class Mat
{
private:
    Type data[nCols][nRows];
public:
    typedef typename arma::Mat<Type>::template fixed<nRows,nCols> ArmaType;
    Mat() {
        for (uint i = 0; i < nRows * nCols; ++i) {
            *(*data + i) = 0;
        }
    }
    Mat(std::initializer_list<std::initializer_list<Type>> list) {
        auto listIt = list.begin();
        auto it = listIt->begin();
        for (uint i = 0; i < nRows; ++i, ++listIt) {
            it = listIt->begin();
            for (uint j = 0; j < nCols; ++j, ++it) {
                data[j][i] = *it;
            }
        }
    }
    Mat(std::initializer_list<Type> list) {
        auto it = list.begin();
        for (uint i = 0; i < nRows * nCols; ++i, ++it) {
            *(*data + i) = *it;
        }
    }
    inline Mat(const Mat<Type, nRows, nCols> & other) {
        for (uint i = 0; i < nRows * nCols; ++i) {
            *(*data + i) = other[i];
        }
    }
    inline Mat(const ArmaType & other) {
        for (uint i = 0; i < nRows * nCols; ++i) {
            *(*data + i) = other[i];
        }
    }
    /*
    inline const Type * begin() const {
        return *data;
    }
    inline Type * begin() {
        return *data;
    }
    inline const Type * end() const {
        return begin() + nCols * nRows;
    }
    inline Type * end() {
        return begin() + nCols * nRows;
    }
    */
    inline uint size() const {
        return nRows * nCols;
    }
    inline const Type * memptr() const {
        return *data;
    }
    inline Type * memptr() {
        return *data;
    }
    inline const Type operator[](uint index) const {
        return *(*data + index);
    }
    inline Type & operator[](uint index) {
        return *(*data + index);
    }
    inline const Type operator()(uint index) const {
        return *(*data + index);
    }
    inline Type & operator()(uint index) {
        return *(*data + index);
    }
    inline const Type operator()(uint row, uint col) const {
        return data[col][row];
    }
    inline Type & operator()(uint row, uint col) {
        return data[col][row];
    }
    inline ArmaType arma() const {
        return ArmaType(memptr());
    }
    inline const Mat<Type, nRows, nCols> & operator=(const Mat<Type, nRows, nCols> & other) {
        for (uint i = 0; i < nRows * nCols; ++i) {
            *(*data + i) = other[i];
        }
        return *this;
    }
    inline const Mat<Type, nRows, nCols> & operator=(const ArmaType & other) {
        for (uint i = 0; i < nRows * nCols; ++i) {
            *(*data + i) = other[i];
        }
        return *this;
    }
    inline const Mat<Type, nRows, nCols> & operator=(std::initializer_list<std::initializer_list<Type>> list) {
        const auto listIt = list.begin();
        const auto it = listIt->begin();
        for (uint i = 0; i < nRows; ++i, ++listIt) {
            it = listIt->begin();
            for (uint j = 0; j < nCols; ++j, ++it) {
                data[j][i] = *it;
            }
        }
        return *this;
    }
    inline const Mat<Type, nRows, nCols> & operator=(std::initializer_list<Type> list) {
        const auto it = list.begin();
        for (uint i = 0; i < nRows * nCols; ++i, ++it) {
            *(*data + i) = *it;
        }
        return *this;
    }
    /*
    inline bool operator==(const ArmaType & other) {
        const Type * first1 = memptr();
        const Type * last1 = memptr() + nRows * nCols;
        const Type * first2 = other.memptr();
        for (; first1 != last1; ++first1, ++first2) {
            if (*first1 != *first2) {
                return false;
            }
        }
        return true;
    }
    inline bool operator==(const Mat<Type, nRows, nCols> & other) {
        const Type * first1 = memptr();
        const Type * last1 = memptr() + nRows * nCols;
        const Type * first2 = other.memptr();
        for (; first1 != last1; ++first1, ++first2) {
            if (*first1 != *first2) {
                return false;
            }
        }
        return true;
    }
    */
    inline void zeros() {
        for (uint i = 0; i < nRows * nCols; ++i) {
            *(*data + i) = 0;
        }
    }
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
inline bool approx_equal(const Mat<Type, nRows, nCols> & a, const Mat<Type, nRows, nCols> & b, const char * method, Type tol) {
//    return arma::approx_equal(a.arma(), b.arma(), method, tol);
}

template <class Type, uint nRows, uint nCols>
inline bool approx_equal(const Mat<Type, nRows, nCols> & a, const typename arma::Mat<Type>::template fixed<nRows,nCols> & b, const char * method, Type tol) {
//    return arma::approx_equal(a.arma(), b, method, tol);
}

template <class Type, uint nRows, uint nCols>
inline bool approx_equal(const typename arma::Mat<Type>::template fixed<nRows,nCols> & a, const Mat<Type, nRows, nCols> & b, const char * method, Type tol) {
//    return arma::approx_equal(a, b.arma(), method, tol);
}

template <class Type, uint nRows, uint nCols>
inline const Mat<Type, nRows, nCols> & operator+=(Mat<Type, nRows, nCols> & a, const Mat<Type, nRows, nCols> & b) {
    for (uint i = 0; i < nRows * nCols; ++i) {
        a[i] += b[i];
    }
    return a;
}

template <class Type, uint nRows, uint nCols>
inline const Mat<Type, nRows, nCols> & operator/=(Mat<Type, nRows, nCols> & a, Type b) {
    for (uint i = 0; i < nRows * nCols; ++i) {
        a[i] /= b;
    }
    return a;
}

template <uint N>
using vec = Mat<double, N, 1>;

template <uint N, uint M>
using mat = Mat<double, N, M>;

}

#endif // ARMOR_H
