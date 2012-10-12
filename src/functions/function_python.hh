/*
 * function_python.hh
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */
// TODO: make FunctionPython dummy class if we do not have python, so that we
// need not optional code elsewhere

#include "functions/functions_all.hh"

#ifndef FUNCTION_PYTHON_HH_
#define FUNCTION_PYTHON_HH_


#include "system/system.hh"
#include "system/python_loader.hh"
#include "functions/function_base.hh"

#include <string>
using namespace std;

/**
 *
 * This class assumes function python function with @p dim arguments containing coordinates of the given point.
 * The function should return  a tuple representing a vector value (possibly of size one for scalar functions)
 *
 * TODO:
 * - use rather only one argument - tuple representing the whole point
 * - time functions
 * - set some parameter in the python module
 * - for functions with one component allow python functions returning the value directly
 *
 */
template <int dim>
class FunctionPython : public FunctionBase<dim>
{
public:
    typedef typename FunctionBase<dim>::Point Point;

    FunctionPython(const unsigned int n_components=1, const double init_time=0.0);

    static Input::Type::Record &get_input_type();

    void init_from_input( Input::Record rec);

    /**
     * Set the file and function to be called.
     */
    void set_python_function_from_file( const string &file_name, const string &func_name);

    /**
     * Set the source in a string and name of the function to be called.
     */
    void set_python_function_from_string( const string &python_source, const string &func_name);

    /**
     * Returns one scalar value in one given point.
     */
    virtual double value(const Point &p, const unsigned int  component = 0) const;
    /**
     * Returns one vector value in one given point.
     */
    virtual void   vector_value(const Point &p, std::vector<double>     &value) const;

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void   value_list (const std::vector< Point >  &point_list,
                       std::vector<double>         &value_list,
                       const unsigned int  component = 0) const;

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void   vector_value_list (const std::vector< Point >    &point_list,
                              std::vector< std::vector<double> >      &value_list) const;

    virtual ~FunctionPython();

private:
    /**
     * Common part of set_python_function_from_* methods
     */
    void set_func(const string &func_name);
#ifdef HAVE_PYTHON
    PyObject *p_func_;
    PyObject *p_module_;
    mutable PyObject *p_args_;
    mutable PyObject *p_value_;
#endif // HAVE_PYTHON

};



#endif /* FUNCTION_PYTHON_HH_ */
