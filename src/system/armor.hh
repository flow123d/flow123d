//#define ARMA_DONT_USE_WRAPPER
//#define ARMA_NO_DEBUG
#include <armadillo>

namespace Armor {

template <class Type, uint nRows, uint nCols>
class Mat
{
private:
    Type data[nCols][nRows];
    typedef typename arma::Mat<Type>::template fixed<nRows,nCols> ArmaType;
public:
    Mat() {
        uint size{nCols * nRows};
        for (uint i{0}; i < size; ++i) {
            *(*data + i) = 0;
        }
    }
    Mat(std::initializer_list<std::initializer_list<Type>> list) {
        uint i{(uint) -1};
        uint j{0};
        for (auto & row : list) {
            for (auto & item : row) {
                data[++i][j] = item;
            }
            i = (uint) -1;
            ++j;
        }
    }
    Mat(std::initializer_list<Type> list) {
        uint i{(uint) -1};
        for (auto & item : list) {
            *(*data + ++i) = item;
        }
    }
    inline Mat(const ArmaType & a) {
        uint size{nCols * nRows};
        for (uint i{0}; i < size; ++i) {
            *(*data + i) = a[i];
        }
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
        return ArmaType(data[0]);
    }
    inline operator ArmaType() const {
        return this->arma();
    }
    inline const Mat<Type, nRows, nCols> & operator=(const Mat<Type, nRows, nCols> & other) {
        uint size{nCols * nRows};
        for (uint i{0}; i < size; ++i) {
            *(*data + i) = other[i];
        }
        return *this;
    }
    inline const Mat<Type, nRows, nCols> & operator=(const ArmaType & other) {
        uint size{nCols * nRows};
        for (uint i{0}; i < size; ++i) {
            *(*data + i) = other[i];
        }
        return *this;
    }
    inline const Mat<Type, nRows, nCols> & operator=(std::initializer_list<std::initializer_list<Type>> list) {
        uint i{(uint) -1};
        uint j{0};
        for (auto & row : list) {
            for (auto & item : row) {
                data[++i][j] = item;
            }
            i = (uint) -1;
            ++j;
        }
        return *this;
    }
    inline const Mat<Type, nRows, nCols> & operator=(std::initializer_list<Type> list) {
        uint i{(uint) -1};
        for (auto & item : list) {
            *(*data + ++i) = item;
        }
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
        uint size = nRows * nCols;
        for (uint i{0}; i < size; ++i) {
            if (*(*data + i) != other[i]) {
                return false;
            }
        }
        return true;
    }
    uint size() const {
        return nRows * nCols;
    }
    double * memptr() {
        return *data;
    }
};

template <class Type, uint nRows, uint nCols>
inline double dot(const Mat<Type, nRows, nCols> & a, const Mat<Type, nRows, nCols> & b) {
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

}

