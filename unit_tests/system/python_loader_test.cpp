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

using namespace std;


string python_function = R"CODE(
import math

def func_xyz(x,y,z):
    return ( x*y*z+ , )     # one value tuple

)CODE";


TEST(PythonLoader, python_error) {
	EXPECT_THROW( { PythonLoader::load_module_from_string("func_xyz", python_function); }, PythonLoader::ExcPythonError);
    //"Program Error: Python Error: invalid syntax (flow123d_python_loader, line 5)"
}


#endif // FLOW123D_HAVE_PYTHON
