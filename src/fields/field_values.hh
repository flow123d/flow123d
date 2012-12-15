/*
 * field_values.hh
 *
 *  Created on: Dec 6, 2012
 *      Author: jb
 */

#ifndef FIELD_VALUES_HH_
#define FIELD_VALUES_HH_

#include <armadillo>
#include <boost/type_traits.hpp>
#include <boost/format.hpp>


std::string type_name_(double)
{ return "Real"; }

std::string type_name_(int)
{ return "Int"; }


/**
 * Template for class representing all possible Filed values. It is just common interface to
 * scalar (double, int) values, vector and tensor values with fixed size and vector/tensor values with variable size.
 *
 * ET is type of elements, n_cols and n_rows gives fixed dimensions of the tensor value (nx1 is vector, 1x1 is scalar,
 * 0x1 is variable size vector, 0x0 is variable size tensor (not implemented yet) )
 *
 * TODO: change (i,j) access to .at(i,j) access ( latter is without bound check)
 */
template <int NRows, int NCols, class ET>
class FieldValue_ {
public:
    typedef ET element_type;
    typedef typename arma::Mat<ET>::template fixed<NRows, NCols> return_type;
    static std::string type_name() { return boost::str(boost::format("%s[%d,%d]") % type_name_( ET() ) % NRows % NCols); }

    FieldValue_(return_type &val) : value_(val) {}
    void set_n_comp(unsigned int) {};
    inline unsigned int n_cols() const
        { return NCols; }
    inline unsigned int n_rows() const
        { return NRows; }
    inline ET &operator() ( unsigned int i, unsigned int j)
        { return value_(i,j); }
    inline operator return_type()
        { return value_;}

private:
    return_type &value_;
};




/// Specialization for scalars
template <class ET>
class FieldValue_<1,1,ET> {
public:
    typedef ET element_type;
    typedef ET return_type;
    static std::string type_name() { return boost::str(boost::format("%s") % type_name_( ET() ) ); }

    FieldValue_(return_type &val) : value_(val) {}
    void set_n_comp(unsigned int) {};
    inline unsigned int n_cols() const
        { return 1; }
    inline unsigned int n_rows() const
        { return 1; }
    inline ET &operator() ( unsigned int, unsigned int )
        { return value_; }
    inline operator return_type()
        { return value_;}

private:
    return_type &value_;
};

/// Specialization for variable size vectors
template <class ET>
class FieldValue_<0,1,ET> {
public:
    typedef ET element_type;
    typedef arma::Col<ET> return_type;
    static std::string type_name() { return boost::str(boost::format("%s[n]") % type_name_( ET() ) ); }

    FieldValue_(return_type &val) : value_(val) {}
    void set_n_comp(unsigned int n_comp) { value_ = return_type(n_comp); };
    inline unsigned int n_cols() const
        { return 1; }
    inline unsigned int n_rows() const
        { return value_.n_rows; }
    inline ET &operator() ( unsigned int i, unsigned int )
        { return value_[i]; }
    inline operator return_type()
        { return value_;}

private:
    return_type &value_;
};

/// Specialization for fixed size vectors
template <int NRows, class ET>
class FieldValue_<NRows,1,ET> {
public:
    typedef ET element_type;
    typedef typename arma::Col<ET>::template fixed<NRows> return_type;
    static std::string type_name() { return boost::str(boost::format("%s[%d]") % type_name_( ET() ) % NRows ); }

    FieldValue_(return_type &val) : value_(val) {}
    void set_n_comp(unsigned int) {};
    inline unsigned int n_cols() const
        { return 1; }
    inline unsigned int n_rows() const
        { return NRows; }
    inline ET &operator() ( unsigned int i, unsigned int )
        { return value_[i]; }
    inline operator return_type()
        { return value_;}

private:
    return_type &value_;
};






/**
 * Class that provides common template-less interface to all templated Field structes.
 *
 * Maybe we can delete this since common interface will be given by "Quantity".
 */

template <int spacedim>
struct FieldValue {
    // typedefs for possible field values
    typedef FieldValue_<1,1,int>     Discrete;
    typedef FieldValue_<1,1,double>   Scalar;
    typedef FieldValue_<spacedim,1,double> VectorFixed;
    typedef FieldValue_<0,1,double> Vector;
    typedef FieldValue_<spacedim,spacedim,double> TensorFixed;
};










#endif /* FIELD_VALUES_HH_ */
