/*
 * python_loader_test.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */



#define FEAL_OVERRIDE_ASSERTS

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
	EXPECT_THROW_WHAT( { PythonLoader::load_module_from_string("func_xyz", python_print); },
	    PythonLoader::ExcPythonError,
        "Message: global name 'a' is not defined\nTraceback");
}

TEST(PythonLoader, test_embedded_python) {
    FilePath::set_io_dirs(".", UNIT_TESTS_SRC_DIR, "", ".");
    
    PythonLoader::initialize();
    
    string embedded_path = "build_tree/lib";
    char* pPath = getenv ("PYTHONPATH");
    PySys_SetPath(pPath);
    
    PyObject * arguments = PyTuple_New (0);
    PyObject * module = PythonLoader::load_module_from_file(string(UNIT_TESTS_SRC_DIR) + "/system/python_embedded.py");
    PyObject * callable  = PythonLoader::get_callable (module, "test");
    PyObject * result = PyObject_CallObject (callable, arguments);

    if (PyString_Check(result)) {
        string result_string = string(PyString_AsString(result));
        cout << "result: " << endl << result_string << endl;
        
        stringstream lines(result_string);
        string line;
        while(std::getline(lines,line,'\n')) {
            if (line.find(embedded_path) == string::npos) {
                FAIL() << "Python is not using embedded library! Path must contain '" << embedded_path << "' part :" << line << endl;
            } else {
                cout << "OK Using embedded python library: " << line << endl;
            }
        }
    } else {
        FAIL() << "Returned value from module is not type of string. Embedded Python now working propertly!";
    }
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


#endif // FLOW123D_HAVE_PYTHON
