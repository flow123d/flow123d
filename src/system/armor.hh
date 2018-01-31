//#define ARMA_DONT_USE_WRAPPER
//#define ARMA_NO_DEBUG
#include <armadillo>
#include <array>

namespace Armor {

template <class Type, uint nRows, uint nCols>
class Mat
{
private:
    std::array<std::array<Type, nRows>, nCols> data;
    typedef typename arma::Mat<Type>::template fixed<nRows,nCols> ArmaType;
public:
    Mat() {
        for (uint i = 0; i < nRows * nCols; ++i) {
            data[0][i] = 0;
        }
    }
    Mat(std::initializer_list<std::initializer_list<Type>> list) {
        const auto * listIt = list.begin();
        const Type * it;
        for (uint i = 0; i < nRows; ++i) {
            it = (listIt + i)->begin();
            for (uint j = 0; j < nCols; ++j) {
                data[j][i] = *(it + j);
            }
        }
    }
    Mat(std::initializer_list<Type> list) {
        const Type * it = list.begin();
        for (uint i = 0; i < nRows * nCols; ++i) {
            data[0][i] = *(it + i);
        }
    }
    inline Mat(const Armor::Mat<Type, nRows, nCols> & other) {
        for (uint i = 0; i < nRows * nCols; ++i) {
            data[0][i] = other[i];
        }
    }
    inline Mat(const ArmaType & other) {
        for (uint i = 0; i < nRows * nCols; ++i) {
            data[0][i] = other[i];
        }
    }
    inline const Type * begin() const {
        return data.data()->data();
    }
    inline Type * begin() {
        return data.data()->data();
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
        return data[0][index];
    }
    inline Type & operator[](uint index) {
        return data[0][index];
    }
    inline const Type & operator()(uint row, uint col) const {
        return data[col][row];
    }
    inline Type & operator()(uint row, uint col) {
        return data[col][row];
    }
    inline ArmaType arma() const {
        return ArmaType(begin());
    }
    inline const Mat<Type, nRows, nCols> & operator=(const Mat<Type, nRows, nCols> & other) {
        for (uint i = 0; i < nRows * nCols; ++i) {
            data[0][i] = other[i];
        }
        return *this;
    }
    inline const Mat<Type, nRows, nCols> & operator=(const ArmaType & other) {        
        for (uint i = 0; i < nRows * nCols; ++i) {
            data[0][i] = other[i];
        }
        return *this;
    }
    inline const Mat<Type, nRows, nCols> & operator=(std::initializer_list<std::initializer_list<Type>> list) {
        const auto * listIt = list.begin();
        const Type * it;
        for (uint i = 0; i < nRows; ++i) {
            it = (listIt + i)->begin();
            for (uint j = 0; j < nCols; ++j) {
                data[j][i] = *(it + j);
            }
        }
        return *this;
    }
    inline const Mat<Type, nRows, nCols> & operator=(std::initializer_list<Type> list) {
        const Type * it = list.begin();
        for (uint i = 0; i < nRows * nCols; ++i) {
            data[0][i] = *(it + i);
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

template <uint N>
using vec = Mat<double, N, 1>;

template <uint N, uint M>
using mat = Mat<double, N, M>;

}
