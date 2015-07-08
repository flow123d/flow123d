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

//FLOW123D_FORCE_LINK_IN_CHILD(field_python)



template <int spacedim, class Value>
const Input::Type::Record & FieldPython<spacedim, Value>::get_input_type()
{
    return it::Record("FieldPython", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field given by a Python script.")
		.derive_from(FieldAlgorithmBase<spacedim, Value>::get_abstract_input_type())
		.declare_key("script_string", it::String(), it::Default::read_time("Obligatory if 'script_file' is not given."),
				"Python script given as in place string")
		.declare_key("script_file", it::FileName::input(), it::Default::read_time("Obligatory if 'script_striong' is not given."),
				"Python script given as external file")
		.declare_key("function", it::String(), it::Default::obligatory(),
				"Function in the given script that returns tuple containing components of the return type.\n"
				"For NxM tensor values: tensor(row,col) = tuple( M*row + col ).")
		.close();
}


template <int spacedim, class Value>
const int FieldPython<spacedim, Value>::registrar =
		Input::register_class< FieldPython<spacedim, Value>, unsigned int >("FieldPython") +
		FieldPython<spacedim, Value>::get_input_type().size();



template <int spacedim, class Value>
FieldPython<spacedim, Value>::FieldPython(unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>( n_comp)
{
#ifdef FLOW123D_HAVE_PYTHON
    p_func_=NULL;
    p_module_=NULL;
    p_args_=NULL;
    p_value_=NULL;
#else
    xprintf(UsrErr, "Flow123d compiled without support for Python, FieldPython can not be used.\n");
#endif // FLOW123D_HAVE_PYTHON
}



template <int spacedim, class Value>
void FieldPython<spacedim, Value>::set_python_field_from_string(const string &python_source, const string &func_name)
{
#ifdef FLOW123D_HAVE_PYTHON
    p_module_ = PythonLoader::load_module_from_string("python_field_"+func_name, python_source);
    set_func(func_name);
#endif // FLOW123D_HAVE_PYTHON
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
#ifdef FLOW123D_HAVE_PYTHON
    p_module_ = PythonLoader::load_module_from_file( string(file_name) );
    set_func(func_name);
#endif // FLOW123D_HAVE_PYTHON
}




template <int spacedim, class Value>
void FieldPython<spacedim, Value>::set_func(const string &func_name)
{
#ifdef FLOW123D_HAVE_PYTHON
	p_func_ = PythonLoader::get_callable(p_module_, func_name);

    p_args_ = PyTuple_New( spacedim );

    // try field call
    for(unsigned int i = 0; i < spacedim; i++) {
        p_value_ = PyFloat_FromDouble( double(i) );
        PyTuple_SetItem(p_args_, i, p_value_);
    }
    p_value_ = PyObject_CallObject(p_func_, p_args_);
    PythonLoader::check_error();

    if ( ! PyTuple_Check( p_value_)) {
    	stringstream ss;
    	ss << "Field '" << func_name << "' from the python module: " << PyModule_GetName(p_module_) << " doesn't return Tuple." << endl;
        THROW( ExcMessage() << EI_Message( ss.str() ));
    }

    unsigned int size = PyTuple_Size( p_value_);

    unsigned int value_size=this->value_.n_rows() * this->value_.n_cols();
    if ( size !=  value_size) {
        xprintf(UsrErr, "Field '%s' from the python module: %s returns %d components but should return %d components.\n"
                ,func_name.c_str(), PyModule_GetName(p_module_), size, value_size);
    }

#endif // FLOW123D_HAVE_PYTHON

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
        ASSERT( envelope.n_rows()==this->value_.n_rows(),
                "value_list[%d] has wrong number of rows: %d; should match number of components: %d\n",
                i, envelope.n_rows(),this->value_.n_rows());
        set_value(point_list[i], elm, envelope );
    }
}

/**
* Returns one vector value in one given point.
*/
template <int spacedim, class Value>
void FieldPython<spacedim, Value>::set_value(const Point &p, const ElementAccessor<spacedim> &elm, Value &value)
{
#ifdef FLOW123D_HAVE_PYTHON
    for(unsigned int i = 0; i < spacedim; i++) {
        p_value_ = PyFloat_FromDouble( p[i] );
        PyTuple_SetItem(p_args_, i, p_value_);
    }
    p_value_ = PyObject_CallObject(p_func_, p_args_);
    PythonLoader::check_error();

    unsigned int pos =0;
    for(unsigned int row=0; row < value.n_rows(); row++)
        for(unsigned int col=0; col < value.n_cols(); col++, pos++)
            if ( boost::is_integral< typename Value::element_type >::value ) value(row,col) = PyLong_AsLong( PyTuple_GetItem( p_value_, pos ) );
            else value(row,col) = PyFloat_AsDouble( PyTuple_GetItem( p_value_, pos ) );

#endif // FLOW123D_HAVE_PYTHON
}




template <int spacedim, class Value>
FieldPython<spacedim, Value>::~FieldPython() {
#ifdef FLOW123D_HAVE_PYTHON
    Py_CLEAR(p_module_);
    Py_CLEAR(p_func_);
    Py_CLEAR(p_value_);
    Py_CLEAR(p_args_);
#endif // FLOW123D_HAVE_PYTHON
}


#endif /* FIELD_PYTHON_IMPL_HH_ */
