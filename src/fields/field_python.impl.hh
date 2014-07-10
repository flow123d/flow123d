/*
 * function_python_impl.hh
 *
 *  Created on: Oct 10, 2012
 *      Author: jb
 */

#ifndef FIELD_PYTHON_IMPL_HH_
#define FIELD_PYTHON_IMPL_HH_


#include <boost/type_traits.hpp>
#include "fields/field_python.hh"

/// Implementation.

namespace it = Input::Type;

template <int spacedim, class Value>
it::Record FieldPython<spacedim, Value>::input_type= get_input_type( FieldAlgorithmBase<spacedim, Value>::input_type, NULL);



template <int spacedim, class Value>
Input::Type::Record FieldPython<spacedim, Value>::get_input_type(
        Input::Type::AbstractRecord &a_type, const typename Value::ElementInputType *eit
        )
{
    it::Record type
    = it::Record("FieldPython", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field given by a Python script.")
    .derive_from(a_type)
    .declare_key("script_string", it::String(), it::Default::read_time("Obligatory if 'script_file' is not given."),
            "Python script given as in place string")
    .declare_key("script_file", it::FileName::input(), it::Default::read_time("Obligatory if 'script_striong' is not given."),
            "Python script given as external file")
    .declare_key("function", it::String(), it::Default::obligatory(),
            "Function in the given script that returns tuple containing components of the return type.\n"
            "For NxM tensor values: tensor(row,col) = tuple( M*row + col ).");

    return type;
}



template <int spacedim, class Value>
FieldPython<spacedim, Value>::FieldPython(unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>( n_comp)
{
#ifdef HAVE_PYTHON
    p_func_=NULL;
    p_module_=NULL;
    p_args_=NULL;
    p_value_=NULL;
#else
    xprintf(UsrErr, "Flow123d compiled without support for Python, FieldPython can not be used.\n");
#endif // HAVE_PYTHON
}



template <int spacedim, class Value>
void FieldPython<spacedim, Value>::set_python_field_from_string(const string &python_source, const string &func_name)
{
#ifdef HAVE_PYTHON
    p_module_ = PythonLoader::load_module_from_string("python_field_"+func_name, python_source);
    set_func(func_name);
#endif // HAVE_PYTHON
}





template <int spacedim, class Value>
void FieldPython<spacedim, Value>::init_from_input(const Input::Record &rec) {
    Input::Iterator<string> it = rec.find<string>("script_string");
    if (it) {
        set_python_field_from_string( *it, rec.val<string>("function") );
    } else {
        Input::Iterator<FilePath> it = rec.find<FilePath>("script_file");
        if (! it) xprintf(UsrErr, "Either 'script_string' or 'script_file' has to be specified in PythonField initialization.");
        set_python_field_from_file( *it, rec.val<string>("function") );
    }
}



template <int spacedim, class Value>
void FieldPython<spacedim, Value>::set_python_field_from_file(const FilePath &file_name, const string &func_name)
{
#ifdef HAVE_PYTHON
    p_module_ = PythonLoader::load_module_from_file( string(file_name) );
    set_func(func_name);
#endif // HAVE_PYTHON
}




template <int spacedim, class Value>
void FieldPython<spacedim, Value>::set_func(const string &func_name)
{
#ifdef HAVE_PYTHON
    char func_char[func_name.size()+2];
    std::strcpy(func_char, func_name.c_str());
    p_func_ = PyObject_GetAttrString(p_module_, func_char );
    if (! p_func_) {
        if (PyErr_Occurred()) PyErr_Print();
        xprintf(UsrErr, "Field '%s' not found in the python module: %s\n", func_name.c_str(), PyModule_GetName(p_module_) );
        Py_XDECREF(p_func_);
        Py_XDECREF(p_module_);

    }
    if (! PyCallable_Check(p_func_)) {
        xprintf(UsrErr, "Field '%s' from the python module: %s is not callable.\n", func_name.c_str(), PyModule_GetName(p_module_) );
        Py_XDECREF(p_func_);
        Py_XDECREF(p_module_);

    }

    p_args_ = PyTuple_New( spacedim );

    // try field call
    for(unsigned int i = 0; i < spacedim; i++) {
        p_value_ = PyFloat_FromDouble( double(i) );
        PyTuple_SetItem(p_args_, i, p_value_);
    }
    p_value_ = PyObject_CallObject(p_func_, p_args_);

    if (p_value_ == NULL) {
        PyErr_Print();
        // TODO: use cout to output also field arguments
        xprintf(Err,"Failed to call field '%s' from the python module: %s\n", func_name.c_str(), PyModule_GetName(p_module_) );
    }

    if ( ! PyTuple_Check( p_value_)) {
        xprintf(UsrErr, "Field '%s' from the python module: %s doesn't return Tuple.\n", func_name.c_str(), PyModule_GetName(p_module_) );
    }

    unsigned int size = PyTuple_Size( p_value_);

    unsigned int value_size=this->value_.n_rows() * this->value_.n_cols();
    if ( size !=  value_size) {
        xprintf(UsrErr, "Field '%s' from the python module: %s returns %d components but should return %d components.\n"
                ,func_name.c_str(), PyModule_GetName(p_module_), size, value_size);
    }

#endif // HAVE_PYTHON

}

/**
 * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
 */
template <int spacedim, class Value>
typename Value::return_type const & FieldPython<spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
    set_value(p,elm, this->value_);
    return this->r_value_;
}


/**
 * Returns std::vector of scalar values in several points at once.
 */
template <int spacedim, class Value>
void FieldPython<spacedim, Value>::value_list (const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{
    ASSERT_EQUAL( point_list.size(), value_list.size() );
    for(unsigned int i=0; i< point_list.size(); i++) {
        Value envelope(value_list[i]);
        set_value(point_list[i], elm, envelope );

    }
}

/**
* Returns one vector value in one given point.
*/
template <int spacedim, class Value>
void FieldPython<spacedim, Value>::set_value(const Point &p, const ElementAccessor<spacedim> &elm, Value &value)
{
#ifdef HAVE_PYTHON
    for(unsigned int i = 0; i < spacedim; i++) {
        p_value_ = PyFloat_FromDouble( p[i] );
        PyTuple_SetItem(p_args_, i, p_value_);
    }
    p_value_ = PyObject_CallObject(p_func_, p_args_);

    if (p_value_ == NULL) {
        PyErr_Print();
        xprintf(Err,"FieldPython call failed\n");
    }

    unsigned int pos =0;
    for(unsigned int row=0; row < value.n_rows(); row++)
        for(unsigned int col=0; col < value.n_cols(); col++, pos++)
            if ( boost::is_integral< typename Value::element_type >::value ) value(row,col) = PyLong_AsLong( PyTuple_GetItem( p_value_, pos ) );
            else value(row,col) = PyFloat_AsDouble( PyTuple_GetItem( p_value_, pos ) );

#endif // HAVE_PYTHON
}




template <int spacedim, class Value>
FieldPython<spacedim, Value>::~FieldPython() {
#ifdef HAVE_PYTHON
    Py_CLEAR(p_module_);
    Py_CLEAR(p_func_);
    Py_CLEAR(p_value_);
    Py_CLEAR(p_args_);
#endif // HAVE_PYTHON
}


#endif /* FIELD_PYTHON_IMPL_HH_ */
