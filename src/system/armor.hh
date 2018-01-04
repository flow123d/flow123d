#define ARMA_DONT_USE_WRAPPER
#define ARMA_NO_DEBUG
#include <armadillo>
#include <algorithm>

namespace Armor {

template <class Type, uint nRows, uint nCols>
class Mat
{
private:
    Type data[nCols][nRows];
    typedef typename arma::Mat<Type>::template fixed<nRows,nCols> ArmaType;
public:
    Mat() {
        std::fill(begin(), end(), 0);
    }
    Mat(std::initializer_list<std::initializer_list<Type>> list) {
        uint i{0}, j{0};
        for (auto & row : list) {
            for (auto & item : row) {
                data[i][j] = item;
                ++i;
            }
            i = 0;
            ++j;
        }
    }
    Mat(std::initializer_list<Type> list) {
        std::copy(list.begin(), list.end(), begin());
    }
    inline const Type * begin() const {
        return *data;
    }
    inline Type * begin() {
        return *data;
    }
    inline const Type * end() const {
        return *data + nCols * nRows;
    }
    inline Type * end() {
        return *data + nCols * nRows;
    }
    inline Mat(const ArmaType & a) {
        std::copy(a.begin(), a.end(), begin());
    }
    inline const Type operator [](uint index) const {
        return *(*data + index);
    }
    inline Type & operator[](uint index) {
        return *(*data + index);
    }
    inline const Type operator()(uint row, uint col) const {
        return data[col][row];
    }
    inline Type & operator()(uint row, uint col) {
        return data[col][row];
    }
    inline ArmaType arma() const {
        return ArmaType(begin());
    }
    inline const Mat<Type, nRows, nCols> & operator=(const Mat<Type, nRows, nCols> & other) {
        std::copy(other.begin(), other.end(), begin());
        return *this;
    }
    inline const Mat<Type, nRows, nCols> & operator=(const ArmaType & other) {
        std::copy(other.begin(), other.end(), begin());
        return *this;
    }
    inline const Mat<Type, nRows, nCols> & operator=(std::initializer_list<std::initializer_list<Type>> list) {
        uint i{0}, j{0};
        for (auto & row : list) {
            for (auto & item : row) {
                data[i][j] = item;
                ++i;
            }
            i = 0;
            ++j;
        }
        return *this;
    }
    inline const Mat<Type, nRows, nCols> & operator=(std::initializer_list<Type> list) {
        std::copy(list.begin(), list.end(), begin());
        return *this;
    }
    inline bool operator==(const ArmaType & other) {
        uint size = nRows * nCols;
        for (uint i{0}; i < size; ++i) {
            if (*(*data + i) != other[i]) {
                return false;
            }
        }
        return true;
    }
    inline bool operator==(const Mat<Type, nRows, nCols> & other) {
        return std::equal(begin(), end(), other.begin());
    }
    inline uint size() const {
        return nRows * nCols;
    }
    inline Type * memptr() {
        return begin();
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
using Vec = Mat<double, N, 1>;

/*
template <uint N, uint M>
class Mat : Mat<double, N, M> {

};
*/

}

