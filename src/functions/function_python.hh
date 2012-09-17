/*
 * function_python.hh
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */

#ifndef FUNCTION_PYTHON_HH_
#define FUNCTION_PYTHON_HH_

#ifdef HAVE_PYTHON

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

    FunctionPython(const unsigned int n_components=1, const double init_time=0.0)
    : FunctionBase<dim>(n_components, init_time),
      p_func_(NULL), p_module_(NULL), p_args_(NULL), p_value_(NULL)
    {}

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

    PyObject *p_func_;
    PyObject *p_module_;
    mutable PyObject *p_args_;
    mutable PyObject *p_value_;
};


/// Implementation.

template <int dim>
void FunctionPython<dim>::set_python_function_from_string(const string &python_source, const string &func_name)
{
    p_module_ = PythonLoader::load_module_from_string("python_function_"+func_name, python_source);
    set_func(func_name);
}



template <int dim>
void FunctionPython<dim>::set_python_function_from_file(const string &file_name, const string &func_name)
{
    p_module_ = PythonLoader::load_module_from_file(file_name);
    set_func(func_name);
}




template <int dim>
void FunctionPython<dim>::set_func(const string &func_name)
{
    p_func_ = PyObject_GetAttrString(p_module_, func_name.c_str() );
    if (! p_func_) {
        if (PyErr_Occurred()) PyErr_Print();
        xprintf(UsrErr, "Function '%s' not found in the python module: %s\n", func_name.c_str(), PyModule_GetName(p_module_) );
        Py_XDECREF(p_func_);
        Py_XDECREF(p_module_);

    }
    if (! PyCallable_Check(p_func_)) {
        xprintf(UsrErr, "Function '%s' from the python module: %s is not callable.\n", func_name.c_str(), PyModule_GetName(p_module_) );
        Py_XDECREF(p_func_);
        Py_XDECREF(p_module_);

    }

    p_args_ = PyTuple_New( dim );

    // try function call
    for(unsigned int i = 0; i < dim; i++) {
        p_value_ = PyFloat_FromDouble( double(i) );
        PyTuple_SetItem(p_args_, i, p_value_);
    }
    p_value_ = PyObject_CallObject(p_func_, p_args_);

    if (p_value_ == NULL) {
        PyErr_Print();
        // TODO: use cout and output also function arguments
        xprintf(Err,"Failed to call function '%s' from the python module: %s\n", func_name.c_str(), PyModule_GetName(p_module_) );
    }

    if ( ! PyTuple_Check( p_value_)) {
        xprintf(UsrErr, "Function '%s' from the python module: %s doesn't return Tuple.\n", func_name.c_str(), PyModule_GetName(p_module_) );
    }

    unsigned int size = PyTuple_Size( p_value_);
    if ( size != this->n_components_) {
        xprintf(UsrErr, "Function '%s' from the python module: %s returns %d components but should return %d components.\n"
                ,func_name.c_str(), PyModule_GetName(p_module_), size, this->n_components_ );
    }


}



template <int dim>
double FunctionPython<dim>::value(const Point &p, const unsigned int  component) const
{
    for(unsigned int i = 0; i < dim; i++) {
        p_value_ = PyFloat_FromDouble( p[i] );
        PyTuple_SetItem(p_args_, i, p_value_);
    }
    p_value_ = PyObject_CallObject(p_func_, p_args_);

    if (p_value_ == NULL) {
        PyErr_Print();
        xprintf(Err,"FunctionPython call failed\n");
    }

    return PyFloat_AsDouble( PyTuple_GET_ITEM( p_value_, component) );
}

/**
* Returns one vector value in one given point.
*/
template <int dim>
void   FunctionPython<dim>::vector_value(const Point &p, std::vector<double>     &value) const
{
    ASSERT_SIZES( this->n_components_, value.size() );
    for(unsigned int i = 0; i < dim; i++) {
        p_value_ = PyFloat_FromDouble( p[i] );
        PyTuple_SetItem(p_args_, i, p_value_);
    }
    p_value_ = PyObject_CallObject(p_func_, p_args_);

    if (p_value_ == NULL) {
        PyErr_Print();
        xprintf(Err,"FunctionPython call failed\n");
    }

    for(unsigned int i = 0; i < this->n_components_; i++) {
        value[i] = PyFloat_AsDouble( PyTuple_GetItem( p_value_, i ) );
    }
}


/**
* Returns std::vector of scalar values in several points at once.
*/
template <int dim>
void   FunctionPython<dim>::value_list (const std::vector< Point >  &point_list,
                  std::vector<double>         &value_list,
                  const unsigned int  component) const
{
    ASSERT_SIZES( point_list.size(), value_list.size() );
    for(unsigned int i=0; i< point_list.size(); i++)
        value_list[i] = value(point_list[i], component);
}

/**
* Returns std::vector of scalar values in several points at once.
*/
template <int dim>
void   FunctionPython<dim>::vector_value_list (const std::vector< Point >    &point_list,
                         std::vector< std::vector<double> >      &value_list) const
{
    ASSERT_SIZES( point_list.size(), value_list.size() );
    for(unsigned int i=0; i< point_list.size(); i++)
        vector_value( point_list[i], value_list[i]);
}


template <int dim>
FunctionPython<dim>::~FunctionPython() {
    Py_CLEAR(p_module_);
    Py_CLEAR(p_func_);
    Py_CLEAR(p_value_);
    Py_CLEAR(p_args_);
}


#endif // HAVE_PYTHON

#endif /* FUNCTION_PYTHON_HH_ */
