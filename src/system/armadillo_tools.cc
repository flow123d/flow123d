/*
 * armadillo_setup.cc
 *
 *  Created on: Jul 4, 2016
 *      Author: jb
 */





#include "armadillo_tools.hh"

//#include <ostream>
//#include <sstream>
#include <type_traits>
#include "system/exc_common.hh"
#include "system/logger.hh"

#include <armadillo>

using namespace std;

/**
 * Auxiliary exception in order to change standard frames_to_cut_.
 */
class ExcArmadillo : public ExceptionBase {
public:
    ExcArmadillo()
    {
        this->frames_to_cut_ = {"arma_stop"};
    }

    void print_info(std::ostringstream &out) const override {
        ::internal::ExcStream estream(out, *this);
        estream << "Armadillo Message: " << EI_Message::val << "\n";
    }
};

/**
 * Auxiliary output buffer to catch the error messages and report stacktrace.
 */
class ArmaStreamBuf : public std::stringbuf {
protected:
    /**
     * Override sync to throw on error message.
     */
    int sync() override;
};

int ArmaStreamBuf::sync() {
    if (this->str().find("error") !=  string::npos) {
        // Can not throw here, since any exception is cached somewhere between armadillo call point and this method.
        // Must write out the stack here.
        auto e = ExcArmadillo();
        e << EI_Message(this->str());
        _LOG( Logger::MsgType::message ) << e.what();

    }
    std::cout << this->str() << std::endl;
    this->str().clear();
    return 1;
}


void armadillo_setup()
{
    static ArmaStreamBuf stream_buf;
    static std::ostream arma_stream(&stream_buf);
    arma::set_stream_err1(arma_stream);
    arma::set_stream_err2(arma_stream);
}




// internal implementation template
template <class T>
string field_value_to_yaml_matrix(const T &mat, unsigned int prec) {
    stringstream ss;
    ss.precision(prec);
    ss << "[ ";
    for(unsigned int i_row=0; i_row < mat.n_rows; i_row ++ ) {
        if (i_row != 0) ss << " , ";
        ss << "[ ";
        for(unsigned int i_col=0; i_col < mat.n_cols; i_col++) {
            if (i_col != 0) ss << " , ";
            ss << mat.at(i_row, i_col);
        }
        ss << " ]";
    }
    ss << " ]";
    return ss.str();
}

// internal implementation template
template <class T>
string field_value_to_yaml_vector(const T &vec, unsigned int prec) {
    stringstream ss;
    ss.precision(prec);
    ss <<  "[ ";
    for(unsigned int i=0; i < vec.n_elem; i++) {
        if (i != 0) ss << " , ";
        ss << vec(i);
    }
    ss << " ]";
    return ss.str();
}

template <class IsScalar>
struct field_value_scalar_resolution;

template <>
struct field_value_scalar_resolution<std::true_type> {
    template <class T>
    inline static string print(const  T &mat, unsigned int prec) {
        stringstream ss;
        ss.precision(prec);
        ss << mat;
        return ss.str();
    }
};

template <>
struct field_value_scalar_resolution<std::false_type> {
    template <class T>
    inline static string print(const  T &mat, unsigned int prec) {
        if (mat.n_cols == 1 || mat.n_rows == 1) {
            return field_value_to_yaml_vector(mat, prec);
        } else {
            return field_value_to_yaml_matrix(mat, prec);
        }
    }
};


template <class T>
inline string field_value_to_yaml(const T &mat, unsigned int prec) {
    return field_value_scalar_resolution< typename std::is_scalar<T>::type >::print(mat, prec);
}



#define FIELD_VALUE_TO_YAML_DIM(dim) \
    template string field_value_to_yaml(arma::Mat<double>::fixed<dim,dim> const  &mat, unsigned int prec);  \
    template string field_value_to_yaml(const arma::Mat<int>::fixed<dim,dim> &mat, unsigned int prec); \
    template string field_value_to_yaml(const arma::Mat<unsigned int>::fixed<dim,dim> &mat, unsigned int prec); \
    template string field_value_to_yaml(const arma::Col<double>::fixed<dim> &mat, unsigned int prec);  \
    template string field_value_to_yaml(const arma::Col<int>::fixed<dim> &mat, unsigned int prec);   \
    template string field_value_to_yaml(const arma::Col<unsigned int>::fixed<dim> &mat, unsigned int prec);

FIELD_VALUE_TO_YAML_DIM(2)
FIELD_VALUE_TO_YAML_DIM(3)

template string field_value_to_yaml(const arma::Col<double> &mat, unsigned int prec);
template string field_value_to_yaml(const arma::Col<unsigned int> &mat, unsigned int prec);

template string field_value_to_yaml(const double &mat, unsigned int prec);
template string field_value_to_yaml(const int &mat, unsigned int prec);
template string field_value_to_yaml(const unsigned int &mat, unsigned int prec);




