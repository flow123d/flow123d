#ifndef ARMOR_HH
#define ARMOR_HH

/**
 * The Mat class template is used for small vectors and matricies.
 * It is just a wrapper class for an actual storage, for calculations the armadillo library
 * is used, in combination with inlining and expression templates all auxiliary constructors
 * and data copy are usually optimized out leaving just optimal code for the desired calculation.
 *
 * When deducing template parameters for the function templates the compiler consider only
 * exact type match with possible const conversion and conversion to the parent. In particular
 * implicit conversions are not considered. This makes transition between Armor and armadillo
 * a bit complicated. In order to make everythink as smooth as possible we use several
 * tricks:
 * 1) defining a friend function inside of the class template creates an explicit function
 * instance. As the instance is already created the implicit conversion is considered.
 * See: https://stackoverflow.com/questions/9787593/implicit-type-conversion-with-template
 * At least with GCC it seems, that this approach works for operators. Argument dependent lookup (ADL)
 * finds the operator or function if at least one of its parameters is Mat template. However
 * for operators the implicit conversion of the other argument is applied, while for the function it is not.
 *
 *
 * 2) In order to consturct Mat<Type, nr, 1> and Mat<Type, 1, 1> also from arma::Col and from
 * the Type variable, we use specialization that is derived from the generic case with nr>1, nc>1.
 *
 * 3) In order to prevent ambiguous method resolution, we can have just single
 * TODO: try to transpose storage format used in Armor (prefered row first) and Armadillo (column first).
 */

//#define ARMA_DONT_USE_WRAPPER
//#define ARMA_NO_DEBUG
#include <armadillo>
#include <array>
#include "system/asserts.hh"
#include "system/logger.hh"


//
//template<typename T1, typename T2>
//inline bool is_close(
//        const arma::BaseCube<typename T1::elem_type,T1>& a,
//        const arma::BaseCube<typename T1::elem_type,T2>& b,
//        const typename T1::pod_type abs_tol=1e-10, const typename T1::pod_type rel_tol=1e-6)
//{
//    double a_diff = arma::norm(a.get_ref() - b.get_ref(), 1);
//    if (a_diff < abs_tol) return true;
//    if (a_diff < rel_tol * arma::norm(a.get_ref(), 1)) return true;
//    return false;
//}
//
//template<typename T1, typename T2>
//inline bool is_close(
//        const arma::Base<typename T1::elem_type,T1>& a,
//        const arma::Base<typename T1::elem_type,T2>& b,
//        const typename T1::pod_type abs_tol=1e-10, const typename T1::pod_type rel_tol=1e-6)
//{
//    double a_diff = arma::norm(a.get_ref() - b.get_ref(), 1);
//    if (a_diff < abs_tol) return true;
//    if (a_diff < rel_tol * arma::norm(a.get_ref(), 1)) return true;
//    return false;
//}


namespace Armor {

//template <class type, uint nr, uint nc>
//struct _MatSimpleType {
//    typedef typename arma::Mat<type>::template fixed<nr, nc> AType;
//    inline static AType convert(const type *begin) {return AType(begin);}
//};
//
//template <class type, uint nr>
//struct _MatSimpleType<type, nr, 1> {
//    typedef typename arma::Col<type>::template fixed<nr> AType;
//    inline static AType convert(const type *begin) {return AType(begin);}
//};
//
//template <class type>
//struct _MatSimpleType<type, 1, 1> {
//    typedef type AType;
//    inline static AType convert(const type *begin) {return *begin;}
//};


//
//
//    template <class Type, uint nRows, uint nCols>
//    class _Mat
//    {
//    protected:
//        Type * const data;
//
//    public:
//        static const uint n_rows = nRows;
//        static const uint n_cols = nCols;
//
//        typedef typename arma::Mat<Type>::template fixed<nRows, nCols> Arma;
//    //    typedef typename arma::_Mat<Type>::template fixed<nRows, nCols> BaseAr_Matype;
//    //    typedef typename __MatSimpleType<Type, nRows, nCols>::AType Ar_Matype;
//    //    typedef typename arma::Col<Type>::template fixed<nRows> VecType;
//
//        _Mat(Type * const mem)
//        : data(mem)
//        {}
//
//        inline _Mat(const Armor::_Mat<Type, nRows, nCols> & other)
//        : data(other.data)
//        {}
//
//        inline _Mat(const Arma &arma_mat)
//        : data(const_cast<Type *>(arma_mat.memptr()))
//        // TODO: Here possible problem with creating non-constant Mat interface to
//        // a temporary Arma object.
//        {}
//
//
//        inline operator Arma() const {return arma();}
//
//    //    inline operator Arma() {return arma();}
//
//    //    inline const Arma &arma() const {
//    //        return Arma(begin());
//    //    }
//
//        inline Arma arma() const {
//            // See proper way how to deal with l-value wrapper and r-value wrapper distinction:
//            // https://stackoverflow.com/questions/20928044/best-way-to-write-a-value-wrapper-class
//            // Currently we provide only r-value convertion to Arma types.
//            // Use Index access and assignement to modify values.
//            return Arma(begin());
//        }
//
//
//        inline const Type * begin() const {
//            return data;
//        }
//
//        inline Type * begin() {
//            return data;
//        }
//
//        inline const Type * end() const {
//            return begin() + nCols * nRows;
//        }
//
//        inline Type * end() {
//            return begin() + nCols * nRows;
//        }
//
//        inline uint size() const {
//            return nRows * nCols;
//        }
//
//        inline Type * memptr() {
//            return begin();
//        }
//
//        inline const Type & operator[](uint index) const {
//            return data[index];
//        }
//
//        inline Type & operator[](uint index) {
//            return data[index];
//        }
//
//        inline const Type & operator()(uint row, uint col) const {
//            return data[row+col*nRows];
//        }
//
//        inline Type & operator()(uint row, uint col) {
//            return data[row+col*nRows];
//        }
//
//        inline const Type & at(uint row, uint col) const {
//            return data[row+col*nRows];
//        }
//
//        inline Type & at(uint row, uint col) {
//            return data[row+col*nRows];
//        }
//
//
//        inline const _Mat<Type, nRows, nCols> & operator=(const _Mat<Type, nRows, nCols> & other) {
//            _data_copy(other.data);
//            return *this;
//        }
//
//    //    inline const _Mat<Type, nRows, nCols> & operator=(const Arma & other) {
//    //        _data_copy(other.memptr());
//    //        return *this;
//    //    }
//
//    private:
//
//        inline void _data_copy(Type *other_data) {
//    //        DebugOut() << "data: " << data << "\n";
//    //        DebugOut() << "other: " << other_data << "\n";
//            //DebugOut() << "dc size : " << n_rows * n_cols << "\n";
//            for (uint i = 0; i < nRows * nCols; ++i) {
//                data[i] = other_data[i];
//            }
//        }
//
//    public:
//    //    inline const _Mat<Type, nRows, nCols> & operator=(const VecType & other) {
//    //        for (uint i = 0; i < nRows * nCols; ++i) {
//    //            data[i] = other[i];
//    //        }
//    //        return *this;
//    //    }
//
//
//    //    inline const _Mat<Type, nRows, nCols> & operator=(std::initializer_list<std::initializer_list<Type>> list) {
//    //        *this = (Arma(list));
//    //
//    ////        const auto * listIt = list.begin();
//    ////        const Type * it;
//    ////        for (uint i = 0; i < nRows; ++i) {
//    ////            it = (listIt + i)->begin();
//    ////            for (uint j = 0; j < nCols; ++j) {
//    ////                data[i+j*nRows] = *(it + j);
//    ////            }
//    ////        }
//    ////        return *this;
//    //    }
//    //
//    //    inline const _Mat<Type, nRows, nCols> & operator=(std::initializer_list<Type> list) {
//    //        *this = (Arma(list));
//    //
//    ////        const Type * it = list.begin();
//    ////        for (uint i = 0; i < nRows * nCols; ++i) {
//    ////            data[i] = *(it + i);
//    ////        }
//    ////        return *this;
//    //    }
//
//
//
//    //        inline const _Mat<Type, nRows, nCols> & init(std::initializer_list<Type> list) {
//    //            *this = (Arma(list));
//    //        }
//
//
//    //    inline bool operator==(const Arma & other) {
//    //        const Type * first1 = begin();
//    //        const Type * last1 = end();
//    //        const Type * first2 = other.begin();
//    //        for (; first1 != last1; ++first1, ++first2) {
//    //            if (*first1 != *first2) {
//    //                return false;
//    //            }
//    //        }
//    //        return true;
//    //    }
//    //
//    //    inline bool operator==(const _Mat<Type, nRows, nCols> & other) {
//    //        const Type * first1 = begin();
//    //        const Type * last1 = end();
//    //        const Type * first2 = other.begin();
//    //        for (; first1 != last1; ++first1, ++first2) {
//    //            if (*first1 != *first2) {
//    //                return false;
//    //            }
//    //        }
//    //        return true;
//    //    }
//
//    //
//    //    inline friend bool operator==(const _Mat& a, const _Mat& b) {
//    //        return a.arma() == b.arma();
//    //    }
//
//
//
//        inline friend Arma operator+=(_Mat& a, const _Mat& b) {
//            return (a = (a.arma() + b.arma()));
//        }
//
//        inline friend Arma operator-(const _Mat& a, const _Mat& b) {
//            return a.arma() - b.arma();
//        }
//
//        inline friend Arma operator-=(_Mat& a, const _Mat& b) {
//            return (a = (a.arma() - b.arma()));
//        }
//
//
//    //    // Matrix multiplication, matrix-vector multiplication.
//    //    inline friend Arma operator*=(_Mat& a, const _Mat& b) {
//    //        return (a = (a.arma() * b.arma()));
//    //    }
//
//        // Elementwise multiplication.
//        inline friend Arma operator%(const _Mat& a, const _Mat& b) {
//            return a.arma() % b.arma();
//        }
//
//        // Elementwise multiplication.
//        inline friend Arma operator%=(_Mat& a, const _Mat& b) {
//            return (a = (a.arma() % b.arma()));
//        }
//
//        // Elementwise division.
//        inline friend Arma operator/(const _Mat& a, const _Mat& b) {
//            return a.arma() / b.arma();
//        }
//
//        // Elementwise division.
//        inline friend Arma operator/=(_Mat& a, const _Mat& b) {
//            return (a = (a.arma() / b.arma()));
//        }
//
//
//
//        inline friend Type dot(const _Mat& a, const _Mat& b) {
//            return arma::dot(a.arma(), b.arma());
//        }
//
//        inline friend bool approx_equal(const _Mat &a, const _Mat &b,
//                const char* method,
//                double a_tol = 1e-10, double r_tol = 1e-6)
//        {
//            return arma::approx_equal(a.arma(), b.arma(), method, a_tol, r_tol);
//        }
//
//
//        //template <class Type, uint nRows, uint nCols>
//        inline friend bool is_close(const _Mat &a, const _Mat &b, double a_tol = 1e-10, double r_tol = 1e-6) {
//            is_close(a.arma(), b.arma(), a_tol, r_tol);
//        }
//
//
//
//
//    };
//
//
//
//

//    // Matrix multiplication, matrix-vector multiplication.
//    template<class Type, uint nrA, uint ncA, uint nrB, uint ncB>
//    inline auto operator*(
//            const _Mat<Type, nrA, ncA>& a,
//            const _Mat<Type, nrB, ncB>& b)
//            -> decltype(a.arma() * b.arma())
//    {
//        return a.arma() * b.arma();
//    }
//
//    template<class Type, uint nrA, uint ncA, uint nrB, uint ncB>
//    inline auto operator*(
//            const typename _Mat<Type, nrA, ncA>::Arma& a,
//            const _Mat<Type, nrB, ncB>& b)
//            -> decltype(a * b.arma())
//    {
//        return a * b.arma();
//    }
//
//    template<class Type, uint nrA, uint ncA, uint nrB, uint ncB>
//    inline auto operator*(
//            const _Mat<Type, nrA, ncA>& a,
//            const typename _Mat<Type, nrB, ncB>::Arma& b)
//            -> decltype(a.arma() * b)
//    {
//        return a.arma() * b;
//    }
//
//    /**
//     * Check for close vectors or matrices.
//     * |a - b| < a_tol  or  |a-b| < r_tol * |a|
//     *
//     * Note: Similar function arma::approx_equal.
//     */
//
//
//    /** dot ********************************/
//
//
//    /** operator + ********************************/
//    template <class Type, uint nRows, uint nCols>
//    inline auto operator+(const _Mat<Type, nRows, nCols> & a, const _Mat<Type, nRows, nCols> & b)
//                ->decltype(a.arma() + b.arma())
//    {
//        return a.arma() + b.arma();
//    }
//
//    template <class Type, uint nRows, uint nCols, class TB>
//    inline auto operator+(const _Mat<Type, nRows, nCols> & a, const TB & b)
//                ->decltype(a.arma() + b)
//    {
//        return a.arma() + b;
//    }
//
//    template <class Type, uint nRows, uint nCols, class TB>
//    inline auto operator+(const TB& a, const _Mat<Type, nRows, nCols> & b)
//                ->decltype(a + b.arma())
//    {
//        return a + b.arma();
//    }
//
//
//
//
//    /** operator - ********************************/
//    //
//    //template <class Type, uint nRows, uint nCols>
//    //inline typename Mat<Type, nRows, nCols>::ArmaType operator-(const Mat<Type, nRows, nCols> & a, const Mat<Type, nRows, nCols> & b) {
//    //    return a.arma() - b.arma();
//    //}
//    //
//    //template <class Type, uint nRows, uint nCols>
//    //inline typename Mat<Type, nRows, nCols>::ArmaType operator-(const Mat<Type, nRows, nCols> & a, const typename Mat<Type, nRows, nCols>::ArmaType & b) {
//    //    return a.arma() - b;
//    //}
//    //
//    //template <class Type, uint nRows, uint nCols>
//    //inline typename Mat<Type, nRows, nCols>::ArmaType operator-(const typename Mat<Type, nRows, nCols>::ArmaType & a, const Mat<Type, nRows, nCols> & b) {
//    //    return a - b.arma();
//    //}
//    //
//    //template <class Type, uint nRows, uint nCols>
//    //inline typename Mat<Type, nRows, nCols>::ArmaType operator-=(typename Mat<Type, nRows, nCols>::ArmaType & a, const Mat<Type, nRows, nCols> & b) {
//    //    return (a = (a - b.arma()));
//    //}
//
//    /** operator *  ********************************/
//    //template <class Type, uint nRows, uint nCols>
//    //inline typename Mat<Type, nRows, nCols>::ArmaType operator*(const Mat<Type, nRows, nCols> & a, const Mat<Type, nRows, nCols> & b) {
//    //    return a.arma() * b.arma();
//    //}
//    //
//    //template <class Type, uint nRows, uint nCols>
//    //inline typename Mat<Type, nRows, nCols>::ArmaType operator*(const Mat<Type, nRows, nCols> & a, const typename Mat<Type, nRows, nCols>::ArmaType & b) {
//    //    return a.arma() * b;
//    //}
//    //
//    //template <class Type, uint nRows, uint nCols>
//    //inline typename Mat<Type, nRows, nCols>::ArmaType operator*(const typename Mat<Type, nRows, nCols>::ArmaType & a, const Mat<Type, nRows, nCols> & b) {
//    //    return a * b.arma();
//    //}
//
//
//    /*
//
//    template <class Type, uint nRows, uint nCols>
//    inline typename Mat<Type, nRows, nCols>::ArmaType operator%(const Mat<Type, nRows, nCols> & a, const Mat<Type, nRows, nCols> & b) {
//        return a.arma() % b.arma();
//    }
//
//    template <class Type, uint nRows, uint nCols>
//    inline typename Mat<Type, nRows, nCols>::ArmaType operator%(const Mat<Type, nRows, nCols> & a, const typename Mat<Type, nRows, nCols>::ArmaType & b) {
//        return a.arma() % b;
//    }
//
//    template <class Type, uint nRows, uint nCols>
//    inline typename Mat<Type, nRows, nCols>::ArmaType operator%(const typename Mat<Type, nRows, nCols>::ArmaType & a, const Mat<Type, nRows, nCols> & b) {
//        return a % b.arma();
//    }
//    */
//
//
//
//    /** scalar * and / **/
//    //template <class Type, uint nRows, uint nCols>
//    //inline typename Mat<Type, nRows, nCols>::ArmaType operator*(Type number, const Mat<Type, nRows, nCols> & a) {
//    //    return number * a.arma();
//    //}
//    //
//    //template <class Type, uint nRows, uint nCols>
//    //inline typename Mat<Type, nRows, nCols>::ArmaType operator/(const Mat<Type, nRows, nCols> & a, Type number) {
//    //    return a.arma() / number;
//    //}
//
//
//
//    template <uint N>
//    using vec = Mat<double, N, 1>;
//
//    template <uint N, uint M>
//    using mat = Mat<double, N, M>;
//
//
//    template<class Type, int nr, int nc>
//    Mat<Type, nr, nc> mat_type(typename arma::Mat<Type>::template fixed<nr, nc> &x) {return Mat<Type, nr, nc>(x.memptr());}
//
//    template<class Type>
//    Mat<Type, 1, 1> mat_type(double &x) {return Mat<Type, 1, 1>(&x);}
//
//    template<class eT, int nr>
//    Mat<eT, nr, 1> mat_type(typename arma::Col<eT>::template fixed<nr> &x) {return Mat<eT, nr, 1>(x.memptr());}
//
//
//
//    //template<Type, int nr>
//    //class Mat<Type, nr, 1> : public Mat<Type, nr, 1>{
//    //
//    //}
//

template <class Type, uint nr, uint nc>
using ArmaMat = typename arma::Mat<Type>::template fixed<nr, nc>;

template <class Type, uint nr>
using ArmaVec = typename arma::Col<Type>::template fixed<nr>;

//// assignment wrapper
//
//template <class Type, uint nRows, uint nCols>
//class Mat
//{
//    Type *ptr;
//public:
//    typedef typename ArmaMat<Type, nRows, nCols> Arma;
//
//    // inherit constructors
//    //using ArmaMat<Type, nRows, nCols>::ArmaMat;
//
//    Mat(Type *ptr)
//    : ptr(ptr)
//    {}
//
//    const Arma &operator=(const Arma &other)
//    {
//        for (uint i = 0; i < nRows * nCols; ++i) {
//            *(ptr + i) = *(other.mem + i);
//        }
//        return *this;
//    }
//
//};
//
//
//template <class Type, uint nRows>
//class Mat<Type, nRows, 1> : public ArmaMat<Type, nRows, 1>
//{
//public:
//    typedef typename arma::Col<Type>::template fixed<nRows> ArmaVec;
//    // inherit constructors
//    using ArmaMat<Type, nRows, nCols>::ArmaMat;
//
//    // Add Col constructor
//    Mat<Type, nRows, 1>(const ArmaVec &other)
//    : ArmaMat<Type, nRows, 1>(other.memptr())
//    {}
//
////    const Mat<Type, nRows, 1> &operator=(const ArmaVec &other)
////    {
////        this->_data_copy(other.memptr());
////        return *this;
////    }
//};
//
//
//template <class Type>
//class Mat<Type, 1, 1> : public _Mat<Type, 1, 1>
//{
//public:
//    typedef typename _Mat<Type, 1, 1>::Arma Arma;
//    typedef typename arma::Col<Type>::template fixed<1> ArmaVec;
//    typedef Type Scalar;
//
//    // inherit constructors
//    using _Mat<Type, 1, 1>::_Mat;
//
//    // Add Col and Scalar constructor
//    Mat<Type, 1, 1>(const ArmaVec &other)
//    : _Mat<Type, 1, 1>(other.memptr())
//    {}
//
//    Mat<Type, 1, 1>(const Scalar &other)
//    : _Mat<Type, 1, 1>(&other)
//    {}
//
////    const Mat<Type, 1, 1> &operator=(const ArmaVec &other)
////    {
////        this->_data_copy(other.memptr());
////        return *this;
////    }
////
////    const Mat<Type, 1, 1> &operator=(const Scalar &other)
////    {
////        this->_data_copy(&other);
////        return *this;
////    }
//
//};
//


/**
 * Array of Armor::Mat with given shape. Provides contiguous storage for the data and access to the array elements.
 * The shape of the matrices is specified at run time, so the class Array is independent of additional template parameters.
 * However, to access the array elements, one must use the templated method get().
 */
template<class Type>
class Array {
public:
    class ArrayMatSet {
        Type * ptr_;
        uint n_rows_, n_cols_;
    public:
        inline ArrayMatSet(Type *ptr,  uint n_rows, uint n_cols)
        : ptr_(ptr), n_rows_(n_rows), n_cols_(n_cols)  {}

//        template<class T>
//        ArrayMatSet &operator=(const typename arma::Base<Type, T>& arma_x)
//        {
//            const T &derived = arma_x.get_ref();
//            ASSERT_EQ(n_rows_, T::n_rows);
//            ASSERT_EQ(n_cols_, T::n_cols);
//            copy<T::n_rows, T::n_cols>(derived.memptr());
//        }

        template<long long unsigned int nr, long long unsigned int nc>
        ArrayMatSet &operator=(const ArmaMat<Type, nr, nc>& arma_x)
        {
            ASSERT_EQ(n_rows_, nr);
            ASSERT_EQ(n_cols_, nc);
            copy<nr, nc>(arma_x.memptr());
            return *this;
        }

        template<long long unsigned int nr>
        ArrayMatSet &operator=(const ArmaVec<Type, nr>& arma_x)
        {
            ASSERT_EQ(n_rows_, nr);
            ASSERT_EQ(n_cols_, 1);
            copy<nr, 1>(arma_x.memptr());
            return *this;
        }


        template <uint nr, uint nc>
        void copy(const Type *other_ptr) {
            for (uint i = 0; i < nr * nc; ++i) {
                *(ptr_ + i) = *(other_ptr + i);
            }
        }
    };

public:
    /**
     * Construct array of Armor matrices.
     * @param nv    Number of matrices in the array.
     * @param nr    Number of rows in each matrix.
     * @param nc    Number of columns in each matrix.
     */
    Array(uint nr, uint nc = 1, uint size = 0)
    : data_(new Type[nr * nc * size]),
      n_rows_(nr),
      n_cols_(nc),
      size_(size),
      reserved_(size)
    {
    }
    
    Array(const Array &other)
    : Array(other.n_rows_, other.n_cols_, other.size_)
    {
        for(uint i = 0; i < n_rows_ * n_cols_ * size(); i++) {
            data_[i] = other.data_[i];
        }
    }

    ~Array() {
        delete [] data_;
        data_ = nullptr;
    }

    Array &operator=(const Array &other)
    {
        ASSERT_DBG( (n_rows_ == other.n_rows_) && (n_cols_ == other.n_cols_) );
        reinit(other.size());
        resize(other.size());

        for(uint i = 0; i < n_rows_ * n_cols_ * size(); i++) {
            data_[i] = other.data_[i];
        }
        return *this;
    }

    /**
     * Drop all data and allocate new space of given size.
     * Change number of elements in the array, while keeping the shape of arrays.
     * @param size  New size of array.
     */
    void reinit(uint size) {
        delete [] data_;
        data_ = nullptr;
        reserved_ = size;
        size_ = 0;
        data_ = new Type[n_rows_ * n_cols_ * reserved_];
    }


    /**
     * Resize active part of the allocated space.
     */
    void resize(uint size) {
        ASSERT_LE_DBG(size, reserved_);
        size_ = size;
    }

    inline uint n_rows() const
    {
        return n_rows_;
    }

    inline uint n_cols() const
    {
        return n_cols_;
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
    template<unsigned long long int nr, unsigned long long int nc = 1>
    inline void append(const ArmaMat<Type,nr,nc> &item)
    {
        ASSERT_LE_DBG(size_, reserved_);
        size_ += 1;
        set(size_ - 1) = item;

    }

    template<unsigned long long int nr>
    inline void append(const ArmaVec<Type,nr> &item)
    {
        append<nr, 1>(ArmaMat<Type,nr,1>(item));
    }

    /**
     * Return matrix at given position in array. The returned object is a Armor::Mat
     * pointing to the respective data_ block in the Array's storage.
     * One can assign to the Armor::Mat which performs postponed evaluation and storing the result to the array.
     *
     * @param i  Index of the matrix.
     * TODO: Should be renamed to item(), but we have compilation problem in Field::loc_point_value
     */
//    template<uint nr, uint nc = 1>
//    inline Mat<Type,nr,nc> get(uint mat_index) const
//    {
//        ASSERT_DBG( (nr == n_rows_) && (nc == n_cols_) );
//        return Mat<Type,nr,nc>( data_ + mat_index * n_rows_ * n_cols_ );
//    }
//
//    template<uint nr, uint nc = 1>
//    inline Mat<Type,nr,nc> get(uint mat_index)
//    {
//        ASSERT_DBG( (nr == n_rows_) && (nc == n_cols_) );
//        return Mat<Type,nr,nc>( data_ + mat_index * n_rows_ * n_cols_ );
//    }


//    template<uint nr, uint nc = 1>
//    inline Mat<Type,nr,nc> get(uint mat_index) const
//    {
//        ASSERT_DBG( (nr == n_rows_) && (nc == n_cols_) );
//        return Mat<Type,nr,nc>( data_ + mat_index * n_rows_ * n_cols_ );
//    }

    template<uint nr, uint nc = 1>
    inline ArmaMat<Type,nr,nc> mat(uint mat_index) const
    {
        ASSERT_DBG( (nr == n_rows_) && (nc == n_cols_) );
//        ArmaMat<Type,nr,nc> m;
//        double ** ptr = const_cast<double **>(&(m.mem));
//        *ptr = data_ + mat_index * n_rows_ * n_cols_;
//        return m;
        return ArmaMat<Type,nr,nc>(data_ + mat_index * n_rows_ * n_cols_);
    }

//    template<long long unsigned int nr, long long unsigned int nc = 1>
//    inline void set(uint mat_index, const ArmaMat<Type,nr,nc> &mat)
//    {
//        ASSERT_DBG( (nr == n_rows_) && (nc == n_cols_) );
//
//        for (uint i = 0; i < nr * nc; ++i) {
//            *(data_ + i) = *(mat.mem + i);
//        }
//
//
////        ArmaMat<Type, nr, nc> m;
////        double ** ptr = const_cast<double **>(&(m.mem));
////        *ptr = data_ + mat_index * n_rows_ * n_cols_;
////        m = mat;
//    }



//    template<uint nr, uint nc = 1>
//    inline Vec<Type,nr,nc> get_vec(uint mat_index) const
//    {
//        ASSERT_DBG( (nr == n_rows_) && (nc == n_cols_) );
//        return Vec<Type,nr,nc>( data_ + mat_index * n_rows_ * n_cols_ );
//    }

    template<uint nr>
    inline ArmaVec<Type, nr> vec(uint mat_index) const
    {
        ASSERT_DBG( (nr == n_rows_) && (1 == n_cols_) )(n_rows_)(n_cols_);
        ASSERT_LT_DBG(mat_index, size());
        return ArmaVec<Type, nr>( data_ + mat_index * n_rows_ * n_cols_ );
    }

    inline Type scalar(uint mat_index) const
    {
        ASSERT_DBG( (1 == n_rows_) && (1 == n_cols_) )(n_rows_)(n_cols_);
        ASSERT_LT_DBG(mat_index, size());
        return ArmaMat<Type,1,1>( data_ + mat_index * n_rows_ * n_cols_ )(0);
    }

    inline ArrayMatSet set(uint index) {
        ASSERT_LT_DBG(index, size());
        return ArrayMatSet(data_ + index * n_rows_ * n_cols_, n_rows_, n_cols_);
    }


    /**
     * Return armadillo matrix at given position in array.
     * @param i  Index of matrix.
     */
    inline arma::mat arma_mat(uint i) const
    {
        ASSERT_LT_DBG(i, size());
   	    return arma::mat( data_ + i*n_rows_*n_cols_, n_rows_, n_cols_ );
    }

    /**
     * Return armadillo vector at given position in array.
     * Warning! Method can be used only if nCols == 1.
     * @param i  Index of matrix.
     */
    inline arma::vec arma_vec(uint i) const
    {
        ASSERT_LT_DBG(i, size());
        ASSERT_EQ_DBG(n_cols_, 1);
   	    return arma::vec( data_ + i*n_rows_, n_rows_ );
    }

    Type * data_;

private:
    inline uint space_() { return n_rows_ * n_cols_ * reserved_; }
    uint n_rows_;
    uint n_cols_;
    uint size_;
    uint reserved_;
};


template <uint N>
using vec = ArmaVec<double, N>;

template <uint N, uint M>
using mat = ArmaMat<double, N, M>;

using array = Array<double>;
}

#endif
