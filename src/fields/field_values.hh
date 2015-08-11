/*
 * field_values.hh
 *
 *  Created on: Dec 6, 2012
 *      Author: jb
 *
 *
 */

#ifndef FIELD_VALUES_HH_
#define FIELD_VALUES_HH_

#include <armadillo>
#include <boost/format.hpp>
#include <system/exceptions.hh>
#include "input/input_type.hh"
#include "input/accessors.hh"
#include <ostream>

namespace IT=Input::Type;

/**
 * @file
 *
 * This file contains various dispatch classes to simplify implementation of Fields. Essential is class template
 * @p FieldValues_  which provides unified access and initialization to scalar, vector and matrix type object in Armadillo library
 */

TYPEDEF_ERR_INFO( EI_InputMsg, const string );
DECLARE_INPUT_EXCEPTION( ExcFV_Input, << "Wrong field value input: " << EI_InputMsg::val );


/**
 * Mimics arma::mat<std::string>.
 */
class StringTensor {
public:
    StringTensor( unsigned int n_rows, unsigned int n_cols )
    : n_rows(n_rows),n_cols(n_cols),n_elem(n_rows*n_cols), values_(n_elem) {}

    StringTensor(const std::string &value)
    : n_rows(1), n_cols(1), n_elem(1), values_(1, value) {}

    std::string & at(unsigned int row) { return at(row,0); }
    std::string & at(unsigned int row, unsigned int col) { return values_[col*n_rows+row]; }
    void zeros() {
        for( auto &elem: values_) elem = "0.0";
    }
    unsigned int n_rows;
    unsigned int n_cols;
    unsigned int n_elem;
    operator std::string() {
        ASSERT_EQUAL( n_elem, 1);
        return values_[0];
    }
    const std::string * memptr() {
    	return &(values_[0]);
    }
private:
    std::vector<std::string>  values_;

};


typedef unsigned int FieldEnum;



namespace internal {

// Helper functions to get scalar type name
std::string type_name_(double);
std::string type_name_(int);
std::string type_name_(std::string);
std::string type_name_(FieldEnum);


inline double &scalar_value_conversion(double &ref) {return ref;}
inline int &scalar_value_conversion(int &ref) {return ref;}
inline FieldEnum &scalar_value_conversion(FieldEnum &ref) {return ref;}
inline std::string &scalar_value_conversion(StringTensor &ref) {
    ASSERT( ref.n_rows==1 , "Converting StringTensor(n,m) too std::string with m!=1 or n!=1.");
    return ref.at(0,0);
}


/**
 * InputType dispatch from elementary type @p ET of FieldValues_ to elementary Input::Type, i.e. descendant of Input::Type::Scalar.
 */
template <class ET>
struct InputType { typedef Input::Type::Double type; };

template <>
struct InputType<int> { typedef Input::Type::Integer type; };

template <>
struct InputType<std::string> { typedef Input::Type::String type; };    // for FieldFormula

template <>
struct InputType<FieldEnum> { typedef Input::Type::Selection type; };


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// resolution of Value::return_type

// general element type
template<int NRows, int NCols, class ET>
struct ReturnType { typedef typename arma::Mat<ET>::template fixed<NRows, NCols> return_type; };

template<class ET>
struct ReturnType<1,1,ET> { typedef ET return_type; };

template <class ET>
struct ReturnType<0,1,ET> { typedef arma::Col<ET>  return_type; };

template <int NRows, class ET>
struct ReturnType<NRows,1,ET> { typedef typename arma::Col<ET>::template fixed<NRows> return_type; };


// string element type (for FieldFormula)
template<int NRows, int NCols>
struct ReturnType<NRows, NCols, std::string> { typedef StringTensor return_type; };

template<>
struct ReturnType<1,1, std::string> { typedef StringTensor return_type; };

template <>
struct ReturnType<0,1, std::string> { typedef StringTensor return_type; };

template <int NRows>
struct ReturnType<NRows,1, std::string> { typedef StringTensor return_type; };


// FiledEnum element type - this just returns types with ET=unsigned int, however input should be different
template<int NRows, int NCols>
struct ReturnType<NRows, NCols, FieldEnum> { typedef typename arma::Mat<unsigned int>::template fixed<NRows, NCols> return_type; };

template<>
struct ReturnType<1,1, FieldEnum> { typedef unsigned int return_type; };

template <>
struct ReturnType<0,1, FieldEnum> { typedef arma::Col<unsigned int> return_type; };

template <int NRows>
struct ReturnType<NRows,1, FieldEnum> { typedef typename arma::Col<unsigned int>::template fixed<NRows> return_type; };


// Resolution of helper functions for raw constructor
template <class RT> inline RT & set_raw_scalar(RT &val, double *raw_data) { return *raw_data;}
template <class RT> inline RT & set_raw_scalar(RT &val, int *raw_data) { return *raw_data;}
template <class RT> inline RT & set_raw_scalar(RT &val, string *raw_data) { return val;}
template <class RT> inline RT & set_raw_scalar(RT &val, FieldEnum *raw_data) { return *raw_data;}

template <class RT> inline RT & set_raw_vec(RT &val, double *raw_data) { arma::access::rw(val.mem) = raw_data; return val;}
template <class RT> inline RT & set_raw_vec(RT &val, int *raw_data) { arma::access::rw(val.mem) = raw_data; return val;}
template <class RT> inline RT & set_raw_vec(RT &val, string *raw_data) { return val;}
template <class RT> inline RT & set_raw_vec(RT &val, FieldEnum *raw_data) { arma::access::rw(val.mem) = raw_data; return val;}

template <class RT> inline RT & set_raw_fix(RT &val, double *raw_data) {  val = RT(raw_data); return val;}
template <class RT> inline RT & set_raw_fix(RT &val, int *raw_data) { val = RT(raw_data); return val;}
template <class RT> inline RT & set_raw_fix(RT &val, string *raw_data) { return val;}
template <class RT> inline RT & set_raw_fix(RT &val, FieldEnum *raw_data) { val = RT(raw_data); return val;}

} // namespace internal



/**
 * Template for class representing all possible Filed values. It is just common interface to
 * scalar (double, int) values, vector and tensor values with fixed size and vector/tensor values with variable size.
 *
 * ET is type of elements, n_cols and n_rows gives fixed dimensions of the tensor value (nx1 is vector, 1x1 is scalar,
 * 0x1 is variable size vector, 0x0 is variable size tensor (not implemented yet) )
 *
 * TODO:
 * This wrapper serves at least to several different things:
 * - Unified reading of input values (for FieldConstant, FieldFormula, etc.)
 *    provided by init_from_input
 * - Unified InputType objects, provided by type_name(), get_input_type()
 *
 * - For unified matrix-like access even to scalar and vector values, without compromising performance.
 *    provided by operator(); n_cols, n_rows, from_raw, ...
 *
 * Maybe it could be better to split these two functions into two distinguish but related classes.
 *
 */
template <int NRows, int NCols, class ET>
class FieldValue_ {
public:
    typedef ET element_type;
    typedef typename internal::ReturnType<NRows, NCols, ET>::return_type return_type;
    typedef typename internal::InputType<ET>::type ElementInputType;
    typedef Input::Array AccessType;
    const static int NRows_ = NRows;
    const static int NCols_ = NCols;

    static std::string type_name() { return boost::str(boost::format("%s[%d,%d]") % internal::type_name_( ET() ) % NRows % NCols); }
    static IT::Array get_input_type() {
		if (NRows == NCols)
			// for square tensors allow initialization by diagonal vector, etc.
			return IT::Array( IT::Array( IT::Parameter("element_input_type"), 1), 1 );
		else
			return IT::Array( IT::Array( IT::Parameter("element_input_type"), NCols, NCols), NRows, NRows );

    }


    inline FieldValue_(return_type &val) : value_(val) {}
    inline static const return_type &from_raw(return_type &val, ET *raw_data) {return internal::set_raw_fix(val, raw_data);}
    const ET * mem_ptr() { return value_.memptr(); }


    void init_from_input( AccessType rec ) {
        Input::Iterator<Input::Array> it = rec.begin<Input::Array>();
        if (it->size() == 1 && NRows == NCols) {
            // square tensor
            // input = 3  expands  to [ [ 3 ] ]; init to  3 * (identity matrix)
            // input = [1, 2, 3] expands to [[1], [2], [3]]; init to diag. matrix
            // input = [1, 2, 3, .. , (N+1)*N/2], ....     ; init to symmetric matrix [ [1 ,2 ,3], [2, 4, 5], [ 3, 5, 6] ]


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
                    THROW( ExcFV_Input()
                    		<< EI_InputMsg(
                    				boost::str(boost::format("Initializing symmetric matrix %dx%d by vector of wrong size %d, should be 1, %d, or %d.")
                                	% NRows % NCols % rec.size() % NRows % ((NRows+1)*NRows/2)))
                    		<< rec.ei_address()

                         );
                }

        } else {
            // accept only full tensor
            if (rec.size() == NRows && it->size() == NCols) {

                for (unsigned int row = 0; row < NRows; row++, ++it) {
                    if (it->size() != NCols)
                        THROW( ExcFV_Input() << EI_InputMsg("Wrong number of columns.")
                        		             << rec.ei_address());
                    Input::Iterator<ET> col_it = it->begin<ET>();
                    for (unsigned int col = 0; col < NCols; col++, ++col_it)
                        value_.at(row, col) = *col_it;
                }
            } else {
                THROW( ExcFV_Input()
                		<< EI_InputMsg(
                				boost::str(boost::format("Initializing matrix %dx%d by matrix of wrong size %dx%d.")
                					% NRows % NCols % rec.size() % it->size() ))
                		<< rec.ei_address()
                     );
            }
        }
    }

    void set_n_comp(unsigned int) {};
    inline unsigned int n_cols() const
        { return NCols; }
    inline unsigned int n_rows() const
        { return NRows; }
    inline ET &operator() ( unsigned int i, unsigned int j)
        { return value_.at(i,j); }
    inline ET operator() ( unsigned int i, unsigned int j) const
        { return value_.at(i,j); }
    inline operator return_type() const
        { return value_;}

private:
    return_type &value_;
};

template <class ET>
struct AccessTypeDispatch { typedef ET type;};
template <>
struct AccessTypeDispatch<unsigned int> { typedef Input::Enum type; };



/// **********************************************************************
/// Specialization for scalars
template <class ET>
class FieldValue_<1,1,ET> {
public:
    typedef ET element_type;
    typedef typename internal::ReturnType<1, 1, ET>::return_type return_type;
    typedef typename internal::InputType<ET>::type ElementInputType;
    typedef typename AccessTypeDispatch<ET>::type AccessType;
    const static int NRows_ = 1;
    const static int NCols_ = 1;

    static std::string type_name() { return boost::str(boost::format("%s") % internal::type_name_( ET() ) ); }
    static IT::Parameter get_input_type()
    {
        return IT::Parameter("element_input_type");
    }

    inline FieldValue_(return_type &val) : value_(val) {}

    /**
     * Returns reference to the return_type (i.e. double, or arma::vec or arma::mat); with data provided by the parameter @p raw_data.
     * A reference to a work space @p val has to be provided for efficient work with vector and matrix values.
     */
    inline static const return_type &from_raw(return_type &val, ET *raw_data) {return internal::set_raw_scalar(val, raw_data);}
    const ET * mem_ptr() { return &(internal::scalar_value_conversion(value_)); }

    void init_from_input( AccessType val ) { value_ = return_type(val); }

    void set_n_comp(unsigned int) {};
    inline unsigned int n_cols() const
        { return 1; }
    inline unsigned int n_rows() const
        { return 1; }
    inline ET &operator() ( unsigned int, unsigned int )
        { return internal::scalar_value_conversion(value_); }
    inline ET operator() ( unsigned int i, unsigned int j) const
        { return internal::scalar_value_conversion(value_); }

    inline operator return_type() const
        { return value_;}

private:
    return_type &value_;
};




/// **********************************************************************
/// Specialization for variable size vectors
template <class ET>
class FieldValue_<0,1,ET> {
public:
    typedef ET element_type;
    typedef typename internal::ReturnType<0, 1, ET>::return_type return_type;
    typedef typename internal::InputType<ET>::type ElementInputType;
    typedef Input::Array AccessType;
    const static int NRows_ = 0;
    const static int NCols_ = 1;


    static std::string type_name() { return boost::str(boost::format("%s[n]") % internal::type_name_( ET() ) ); }
    static IT::Array get_input_type() {
        return IT::Array( IT::Parameter("element_input_type"), 1);
    }
    inline static const return_type &from_raw(return_type &val, ET *raw_data) {return internal::set_raw_vec(val, raw_data);}
    const ET * mem_ptr() { return value_.memptr(); }

    inline FieldValue_(return_type &val) : value_(val) {}


    void init_from_input( AccessType rec ) {
        typedef typename AccessTypeDispatch<ET>::type InnerType;
        Input::Iterator<InnerType> it = rec.begin<InnerType>();

        if ( rec.size() == 1 ) {
            for(unsigned int i=0; i< n_rows(); i++)
                value_.at(i)=ET(*it);
        } else if ( rec.size() == n_rows() ) {
            for(unsigned int i=0; i< n_rows(); i++, ++it) {
                value_.at(i)=ET(*it);
            }
        } else {
            THROW( ExcFV_Input()
            		<< EI_InputMsg(
            				boost::str(boost::format("Initializing vector of size %d by vector of size %d.")
                        		% n_rows() % rec.size() ))
                    << rec.ei_address()
                 );
        }
    }

    void set_n_comp(unsigned int n_comp) { value_ = return_type(n_comp,1); };
    inline unsigned int n_cols() const
        { return 1; }
    inline unsigned int n_rows() const
        { return value_.n_rows; }
    inline ET &operator() ( unsigned int i, unsigned int )
        { return value_.at(i); }
    inline ET operator() ( unsigned int i, unsigned int j) const
        { return value_.at(i); }

    inline operator return_type() const
        { return value_;}

private:
    return_type &value_;
};

/// **********************************************************************
/// Specialization for fixed size vectors
template <int NRows, class ET>
class FieldValue_<NRows,1,ET> {
public:
    typedef ET element_type;
    typedef typename internal::ReturnType<NRows, 1, ET>::return_type return_type;
    typedef typename internal::InputType<ET>::type ElementInputType;
    typedef Input::Array AccessType;
    const static int NRows_ = NRows;
    const static int NCols_ = 1;


    static std::string type_name() { return boost::str(boost::format("%s[%d]") % internal::type_name_( ET() ) % NRows ); }
    static IT::Array get_input_type() {
        return IT::Array( IT::Parameter("element_input_type"), 1, NRows);
    }

    inline FieldValue_(return_type &val) : value_(val) {}
    inline static const return_type &from_raw(return_type &val, ET *raw_data) {return internal::set_raw_fix(val, raw_data);}
    const ET * mem_ptr() { return value_.memptr(); }

    void init_from_input( AccessType rec ) {
        Input::Iterator<ET> it = rec.begin<ET>();

        if ( rec.size() == 1 ) {
            for(unsigned int i=0; i< n_rows(); i++)
                value_.at(i)=*it;
        } else if ( rec.size() == NRows ) {
            for(unsigned int i=0; i< NRows; i++, ++it)
                value_.at(i)=*it;
        } else {
            THROW( ExcFV_Input()
            		<< EI_InputMsg(
            				boost::str(boost::format("Initializing fixed vector of size %d by vector of size %d.")
                        		% n_rows() % rec.size() ))
                    << rec.ei_address()
                 );
        }
    }

    void set_n_comp(unsigned int) {};
    inline unsigned int n_cols() const
        { return 1; }
    inline unsigned int n_rows() const
        { return NRows; }
    inline ET &operator() ( unsigned int i, unsigned int )
        { return value_.at(i); }
    inline ET operator() ( unsigned int i, unsigned int j) const
        { return value_.at(i); }

    inline operator return_type() const
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
    typedef FieldValue_<1,1,int>            Integer;
    typedef FieldValue_<1,1, FieldEnum>     Enum;
    typedef FieldValue_<0,1, FieldEnum>     EnumVector;
    typedef FieldValue_<1,1,double>         Scalar;
    typedef FieldValue_<spacedim,1,double>  VectorFixed;
    typedef FieldValue_<0,1,double>         Vector;
    typedef FieldValue_<0,1,int>            IntVector;
    typedef FieldValue_<spacedim,spacedim,double> TensorFixed;
};










#endif /* FIELD_VALUES_HH_ */
