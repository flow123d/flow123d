/*
 * python_loader_test.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */


#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>

#include "system/global_defs.h"


#include <string>
#include <cmath>

#include "system/python_loader.hh"
#include "system/file_path.hh"

using namespace std;
namespace py = pybind11;


string python_function = R"CODE(
import math

def func_xyz(x,y,z):
    return ( x*y*z+ , )     # one value tuple

)CODE";

string python_print = R"CODE(
import math

def func_xyz(x,y,z):
    return ( x*y*z+a , )     # one value tuple

print (func_xyz(1, 2, 3))

)CODE";


string invalid_code = R"CODE(
import math

def func_xyz(x,y,z):
    return ( x*y*z+a , )     # one value tuple

print func_xyz(1, 2, 3)

)CODE";

string invalid_code2 = R"CODE(
this is invalid python code
)CODE";


string produce_error = R"CODE(
def func_xyz():
    return a()
    
def a():
    b()
    
def b():
    return division_by_zero_origin()
    
def division_by_zero_origin():
    return 1/0

func_xyz()
)CODE";


/**
 * We are testing that load_module_from_string call will fail because
 * variable is not defined in the code
 */
TEST(PythonLoader, print_error) {
    EXPECT_THROW_WHAT(
        { PythonLoader::load_module_from_string("test1", "func_xyz", python_print); },
        PythonLoader::ExcPythonError,
        "name 'a' is not defined"
    );
}


/**
 * We are testing that compilation here will fail, since code itself is invalid
 * after Py_CompileString call check_error will react and raise Error
 */
TEST(PythonLoader, compilation_error) {
    EXPECT_THROW_WHAT(
        { PythonLoader::load_module_from_string("test2", "func_xyz", invalid_code); },
        PythonLoader::ExcPythonError,
        "invalid syntax"
    );
    EXPECT_THROW_WHAT(
        { PythonLoader::load_module_from_string("test3", "func_xyz", invalid_code2); },
        PythonLoader::ExcPythonError,
        "invalid syntax"
    );
}


// only test embedded python if we actually copied out Python
// this tests only checks if embedded python is loading modules from correct
// location. This cannot be tested if python was not copied out.
//#ifdef FLOW123D_PYTHON_COPY
//TEST(PythonLoader, test_embedded_python) {
//    FilePath::set_io_dirs(".", UNIT_TESTS_SRC_DIR, "", ".");
//    PythonLoader::initialize();
//
//    // string which must be present in the output
//    string embedded_path = "build_tree/lib";
//
//    // get callable object from file
//    PyObject * arguments = PyTuple_New (0);
//    PyObject * module = PythonLoader::load_module_from_file(string(UNIT_TESTS_SRC_DIR) + "/system/python_embedded.py");
//    PyObject * callable  = PythonLoader::get_callable (module, "test");
//    PyObject * result = PyObject_CallObject (callable, arguments);
//    PythonLoader::check_error();
//
//    // check whether result from python call was indeed string
//    if (PyString_Check(result)) {
//        string result_string = string(PyString_AsString(result));
//
//        stringstream lines(result_string);
//        string line;
//        while(std::getline(lines,line,'\n')) {
//            if (line.find(embedded_path) == string::npos) {
//                FAIL() << "Python is not using embedded library! Path must contain '" << embedded_path << "' part :" << line << endl;
//            } else {
//                cout << "OK Using embedded python library: " << line << endl;
//            }
//        }
//    } else {
//        FAIL() << "Returned value from module is not type of string. Embedded Python is not working properly!";
//    }
//}
//#endif // FLOW123D_PYTHON_COPY

TEST(PythonLoader, function_error) {
	EXPECT_THROW_WHAT( { PythonLoader::load_module_from_string("test5", "func_xyz", python_function); },
	        PythonLoader::ExcPythonError,
	        "invalid syntax");
    /**
     * We are testing that compilation here will succeed but execution of this code
     * will fail, causing traceback to be displayed
     */
    EXPECT_THROW_WHAT(
        { PythonLoader::load_module_from_string("test4", "func_xyz", produce_error); },
        PythonLoader::ExcPythonError,
        "division by zero"  // it should be "division_by_zero_origin"
    );
}


TEST(PythonLoader, from_file) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    py::module_ module = PythonLoader::load_module_from_file("fields/field_python_script.py");
    py::object result = module.attr("func_multi")(2, 3, 4);
    EXPECT_EQ(result.cast<int>(), 24);
}


TEST(PythonLoader, file_error) {
    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    EXPECT_THROW_WHAT( { PythonLoader::load_module_from_file("system/python_loader_script.py"); },
            PythonLoader::ExcPythonError, "invalid syntax");
    EXPECT_THROW_WHAT( { PythonLoader::load_module_from_file("system/field_python_script.py"); }, // file doesn't exist
            FilePath::ExcFileOpen, "Program Error: Can not open file");
}

