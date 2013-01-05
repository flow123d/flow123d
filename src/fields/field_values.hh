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

#include "input/input_type.hh"
#include "input/accessors.hh"
namespace IT=Input::Type;


TYPEDEF_ERR_INFO( EI_InputMsg, const string );
DECLARE_INPUT_EXCEPTION( ExcFV_Input, << "Wrong field value input: " << EI_InputMsg::val );

namespace internal {

// Helper functions to get scalar type name
std::string type_name_(double);
std::string type_name_(int);


template <class ET>
struct InputType { typedef Input::Type::Double type; };

template <>
struct InputType<int> { typedef Input::Type::Integer type; };

// for FieldFormula
template <>
struct InputType<std::string> { typedef Input::Type::String type; };


} // namespace internal



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
    typedef Input::Array InputType;
    static std::string type_name() { return boost::str(boost::format("%s[%d,%d]") % internal::type_name_( ET() ) % NRows % NCols); }
    static Input::Type::Array get_input_type() {
        if (NRows == NCols)
            // for square tensors allow initialization by diagonal vector, etc.
            return IT::Array( IT::Array( typename internal::InputType<ET>::type(), 1), 1 );
        else
            return IT::Array( IT::Array( typename internal::InputType<ET>::type(), NCols, NCols), NRows, NRows );
    }

    inline FieldValue_(return_type &val) : value_(val) {}

    void init_from_input( InputType rec ) {
        Input::Iterator<Input::Array> it = rec.begin<Input::Array>();
        if (NRows == NCols) {
            // square tensor
            // input = 3  expands  to [ [ 3 ] ]; init to  3 times identity matrix
            // input = [1, 2, 3] expands to [[1], [2], [3]]; init to diag. matrix
            // input = [1, 2, 3, .. , (N+1)*N/2], ....     ; init to symmetric matrix [ [1 ,2 ,3], [2, 4, 5], [ 3, 5, 6] ]

            if (it->size() == 1) {
                if (rec.size() == 1)  {// scalar times identity
                    value_.zeros();
                    ET scalar=*(it->begin<ET>());
                    for(unsigned int i=0; i< NRows; i++) value_.at(i,i)=scalar;
                } else if (rec.size() == NRows) { // diagonal vector
                    value_.zeros();
                    for(unsigned int i=0; i< NRows; i++, ++it) value_.at(i,i)=*(it->begin<ET>());
                } else if (rec.size() == (NRows+1)*NRows/2) { // symmetric part
                    for( unsigned int row=0; row<NRows; row++)
                        for( unsigned int col=0; col<NCols; col++)
                            if (row <= col) {
                                value_.at(row,col) = *(it->begin<ET>());
                                ++it;
                            } else value_.at(row,col) = value_.at(col,row);
                } else {
                    THROW( ExcFV_Input() << EI_InputMsg(
                            boost::str(boost::format("Initializing symmetric matrix %dx%d by vector of wrong size %d, should be 1, %d, or %d.")
                                % NRows % NCols % rec.size() % NRows % ((NRows+1)*NRows/2))
                         ));
                }
            }
        } else {
            // accept only full tensor
            if (rec.size() == NRows && it->size() == NCols) {

                for (unsigned int row = 0; row < NRows; row++, ++it) {
                    if (it->size() != NCols)
                        THROW( ExcFV_Input());
                    Input::Iterator<ET> col_it = it->begin<ET>();
                    for (unsigned int col = 0; col < NCols; col++, ++col_it)
                        value_.at(row, col) = *col_it;
                }
            } else {
                THROW( ExcFV_Input() << EI_InputMsg(
                        boost::str(boost::format("Initializing matrix %dx%d by matrix of wrong size %dx%d.")
                            % NRows % NCols % rec.size() % it->size() )
                     ));
            }
        }
    }

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
    typedef ET InputType;

    static std::string type_name() { return boost::str(boost::format("%s") % internal::type_name_( ET() ) ); }
    static typename internal::InputType<ET>::type get_input_type() { return typename internal::InputType<ET>::type(); }

    inline FieldValue_(return_type &val) : value_(val) {}
    void init_from_input( InputType val ) { value_ = val; }

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
    typedef Input::Array InputType;

    static std::string type_name() { return boost::str(boost::format("%s[n]") % internal::type_name_( ET() ) ); }
    static Input::Type::Array get_input_type() {
        return Input::Type::Array( typename internal::InputType<ET>::type(), 1);
    }

    inline FieldValue_(return_type &val) : value_(val) {}
    void init_from_input( InputType rec ) {
        Input::Iterator<ET> it = rec.begin<ET>();

        if ( rec.size() == 1 ) value_.fill( *it );
        else if ( rec.size() == n_rows() ) {
            for(unsigned int i=0; i< n_rows(); i++, ++it)
                value_.at(i)=*it;
        } else {
            THROW( ExcFV_Input() << EI_InputMsg(
                    boost::str(boost::format("Initializing vector of size %d by vector of size %d.")
                        % n_rows() % rec.size() )
                 ));
        }
    }

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
    typedef Input::Array InputType;

    static std::string type_name() { return boost::str(boost::format("%s[%d]") % internal::type_name_( ET() ) % NRows ); }
    static IT::Array get_input_type() {
        return IT::Array( typename internal::InputType<ET>::type(), 1, NRows);
    }

    inline FieldValue_(return_type &val) : value_(val) {}

    void init_from_input( InputType rec ) {
        Input::Iterator<ET> it = rec.begin<ET>();

        if ( rec.size() == 1 ) value_.fill( *it );
        else if ( rec.size() == NRows ) {
            for(unsigned int i=0; i< NRows; i++, ++it)
                value_.at(i)=*it;
        } else {
            THROW( ExcFV_Input() << EI_InputMsg(
                    boost::str(boost::format("Initializing fixed vector of size %d by vector of size %d.")
                        % n_rows() % rec.size() )
                 ));
        }
    }

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
