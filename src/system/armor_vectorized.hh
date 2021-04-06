#ifndef ARMOR_VECTORIZED_HH
#define ARMOR_VECTORIZED_HH

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

#include "system/armor.hh"
#include "include/processor.hh"     // bparser


namespace Armor {

using namespace bparser;

class ArrayVectorized : public Array<double> {
public:
	ArrayVectorized(uint nr, uint nc = 1, uint size = 0)
    : Array<double>(nr, nc), subset_(nullptr)
    {
	    change_reserved(size);
	    change_size(this->reserved_);
    }

	ArrayVectorized(const ArrayVectorized &other)
    : ArrayVectorized(other.n_rows_, other.n_cols_, other.size_)
    {
        for(uint i = 0; i < n_rows_ * n_cols_ * size(); i++) {
            data_[i] = other.data_[i];
        }
    }

    ~ArrayVectorized() {
        delete [] data_;
        delete [] subset_;
    }

    ArrayVectorized &operator=(const ArrayVectorized &other)
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
        change_reserved(size);
        change_size(0);
    }


    /**
     * Resize active part of the allocated space.
     */
    void resize(uint size) {
        while (size%4 != 0) size++;
        ASSERT_LE_DBG(size, reserved_);
        change_size(size);
    }

    double4 *data_;
    uint *subset_;

private:
    inline void change_reserved(uint size) {
        while (size%4 != 0) size++;
        data_ = new double4[this->n_rows_ * this->n_cols_ * size];
        this->reserved_ = size;
    }

    inline void change_size(uint size) {
        ASSERT_LE_DBG(size%4, 0);
        if (subset_ != nullptr) delete [] subset_;
        this->size_ = size;
        subset_ = new uint[this->n_rows_ * this->n_cols_ * size / 4];
        uint i_subset=0;
        for (uint i=0; i<this->n_rows_ * this->n_cols_; ++i) {
        	uint c_shift = i * this->reserved_ / 4;
            for (uint j=0; j<size/4; ++j, ++i_subset) {
                subset_[i_subset] = c_shift+j;
            }
        }
    }
};

}

#endif // ARMOR_VECTORIZED_HH
