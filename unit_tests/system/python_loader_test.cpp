/*
 * python_loader_test.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */



#include <flow_gtest.hh>

#include "system/global_defs.h"


#ifdef FLOW123D_HAVE_PYTHON

#include <string>
#include <cmath>

#include "system/python_loader.hh"
#include "system/file_path.hh"

using namespace std;


string python_function = R"CODE(
import math

def func_xyz(x,y,z):
    return ( x*y*z+ , )     # one value tuple

)CODE";

string python_print = R"CODE(
import math

def func_xyz(x,y,z):
    return ( x*y*z+a , )     # one value tuple

print func_xyz(1, 2, 3)

)CODE";

TEST(PythonLoader, print_error) {
	EXPECT_THROW_WHAT( { PythonLoader::load_module_from_string("func_xyz", python_print); }, PythonLoader::ExcPythonError,
        "Python Error: global name 'a' is not defined");
}


TEST(PythonLoader, function_error) {
	EXPECT_THROW( { PythonLoader::load_module_from_string("func_xyz", python_function); }, PythonLoader::ExcPythonError);
}


TEST(PythonLoader, file_error) {
    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    EXPECT_THROW( { PythonLoader::load_module_from_file(FilePath::get_absolute_working_dir() + "/python_loader_script.py"); },
    		PythonLoader::ExcPythonError);
}

TEST(PythonLoader, module_retrieve) {
  PyObject * python_module;

  // grab module and function by importing module profiler_formatter_module.py
  EXPECT_THROW ( { PythonLoader::load_module_by_name ("profiler.non_existing_module"); }, PythonLoader::ExcPythonError);
  python_module = PythonLoader::load_module_by_name ("profiler.profiler_formatter_module");
  EXPECT_THROW ( { PythonLoader::get_callable (python_module, "non_existing_callable" ); }, PythonLoader::ExcPythonError);
  PythonLoader::get_callable (python_module, "convert" );
  PythonLoader::check_error();
}

TEST(PythonLoader, profiler_broken_json) {
  PyObject * python_module;
  PyObject * convert_method;
  PyObject * arguments;
  PyObject * tmp;
  int argument_index = 0;

  // grab module and function by importing module profiler_formatter_module.py
  python_module = PythonLoader::load_module_by_name ("profiler.profiler_formatter_module");
  convert_method  = PythonLoader::get_callable (python_module, "convert" );

  // generate argument tuple
  arguments = PyTuple_New (3);
  // set json path location as first argument
  tmp = PyString_FromString ("python_profiler_files/broken.json");
  PyTuple_SetItem (arguments, argument_index++, tmp);
  // set output path location as second argument
  tmp = PyString_FromString ("non-existent-file.txt");
  PyTuple_SetItem (arguments, argument_index++, tmp);
  // set Formatter class as third value
  tmp = PyString_FromString ("SimpleTableFormatter");
  PyTuple_SetItem (arguments, argument_index++, tmp);

  // execute method with arguments
  cerr << "Error traceback will appear below" << endl;
  PyObject_CallObject (convert_method, arguments);
  EXPECT_THROW ( { PythonLoader::check_error(); }, PythonLoader::ExcPythonError);
}

TEST(PythonLoader, profiler_missing_json) {
  PyObject * python_module;
  PyObject * convert_method;
  PyObject * arguments;
  PyObject * tmp;
  int argument_index = 0;

  // grab module and function by importing module profiler_formatter_module.py
  python_module = PythonLoader::load_module_by_name ("profiler.profiler_formatter_module");
  convert_method  = PythonLoader::get_callable (python_module, "convert" );

  // generate argument tuple
  arguments = PyTuple_New (3);
  // set json path location as first argument
  tmp = PyString_FromString ("python_profiler_files/missing.json");
  PyTuple_SetItem (arguments, argument_index++, tmp);
  // set output path location as second argument
  tmp = PyString_FromString ("non-existent-file.txt");
  PyTuple_SetItem (arguments, argument_index++, tmp);
  // set Formatter class as third value
  tmp = PyString_FromString ("SimpleTableFormatter");
  PyTuple_SetItem (arguments, argument_index++, tmp);

  // execute method with arguments
  cerr << "Error traceback will appear below" << endl;
  PyObject_CallObject (convert_method, arguments);
  EXPECT_THROW ( { PythonLoader::check_error(); }, PythonLoader::ExcPythonError);
}

TEST(PythonLoader, profiler_empty_json) {
  PyObject * python_module;
  PyObject * convert_method;
  PyObject * arguments;
  PyObject * tmp;
  int argument_index = 0;

  // grab module and function by importing module profiler_formatter_module.py
  python_module = PythonLoader::load_module_by_name ("profiler.profiler_formatter_module");
  convert_method  = PythonLoader::get_callable (python_module, "convert" );

  // generate argument tuple
  arguments = PyTuple_New (3);
  // set json path location as first argument
  tmp = PyString_FromString ("python_profiler_files/empty.json");
  PyTuple_SetItem (arguments, argument_index++, tmp);
  // set output path location as second argument
  tmp = PyString_FromString ("non-existent-file.txt");
  PyTuple_SetItem (arguments, argument_index++, tmp);
  // set Formatter class as third value
  tmp = PyString_FromString ("SimpleTableFormatter");
  PyTuple_SetItem (arguments, argument_index++, tmp);

  // execute method with arguments
  cerr << "Error traceback will appear below" << endl;
  PyObject_CallObject (convert_method, arguments);
  EXPECT_THROW ( { PythonLoader::check_error(); }, PythonLoader::ExcPythonError);
}

TEST(PythonLoader, profiler_ok_json) {
  PyObject * python_module;
  PyObject * convert_method;
  PyObject * arguments;
  PyObject * return_value;
  PyObject * tmp;
  int argument_index = 0;

  // grab module and function by importing module profiler_formatter_module.py
  python_module = PythonLoader::load_module_by_name ("profiler.profiler_formatter_module");
  convert_method  = PythonLoader::get_callable (python_module, "convert" );

  // generate argument tuple
  arguments = PyTuple_New (3);
  // set json path location as first argument
  tmp = PyString_FromString ("python_profiler_files/ok.json");
  PyTuple_SetItem (arguments, argument_index++, tmp);
  // set output path location as second argument
  tmp = PyString_FromString ("python_profiler_files/ok.txt");
  PyTuple_SetItem (arguments, argument_index++, tmp);
  // set Formatter class as third value
  tmp = PyString_FromString ("SimpleTableFormatter");
  PyTuple_SetItem (arguments, argument_index++, tmp);

  // execute method with arguments
  return_value = PyObject_CallObject (convert_method, arguments);
  EXPECT_EQ (return_value, Py_True);
}

#endif // FLOW123D_HAVE_PYTHON
