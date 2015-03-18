/*
 * python_loader_test.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */



#include <flow_gtest.hh>
#include <string>
#include <cmath>


#ifdef FLOW123D_HAVE_PYTHON

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


#endif // FLOW123D_HAVE_PYTHON
