/*
 * python_function_test.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */



#include <gtest/gtest.h>
#include <string>
#include <cmath>


#ifdef HAVE_PYTHON

#include "system/python_loader.hh"
#include "functions/function_python.hh"

using namespace std;

string python_code = R"CODE(
def testFunc():
    print "spam!"

class testClass:
    def testMethod(self):
        print "eggs!"
)CODE";

string python_function = R"CODE(
import math

def func_xyz(x,y,z):
    return ( x*y*z , )     # one value tuple

def func_circle(r,phi):
    return ( r * math.cos(phi), r * math.sin(phi) )
)CODE";



TEST(PythonLoader, all) {
    PyObject * module = PythonLoader::load_module_from_string("my_module", python_code);
    PyObject * p_func = PyObject_GetAttrString(module, "testFunc" );
    PyObject * p_args = PyTuple_New( 0 );
    PyObject_CallObject(p_func, p_args); // this should print out 'spam!'
}

TEST(FunctionPython, two_args) {
    double pi = 4.0 * atan(1);

    Point<2> point_1, point_2;
    point_1(0)=1.0; point_1(1)= pi / 2.0;
    point_2(0)= sqrt(2.0); point_2(1)= 3.0 * pi / 4.0;

    FunctionPython<2> vec_func(2);
    vec_func.set_python_function_from_string(python_function, "func_circle");

    std::vector<double> result(2);
    vec_func.vector_value( point_1, result);
    EXPECT_DOUBLE_EQ( cos(pi /2.0 ) , result[0]); // should be 0.0
    EXPECT_DOUBLE_EQ( 1, result[1]);

    vec_func.vector_value( point_2, result);
    EXPECT_DOUBLE_EQ( -1, result[0]);
    EXPECT_DOUBLE_EQ( 1, result[1]);

}


TEST(FunctionPython, three_args) {
    Point<3> point_1, point_2;
    point_1(0)=1; point_1(1)=0; point_1(2)=0;
    point_2(0)=1; point_2(1)=2; point_2(2)=3;

    FunctionPython<3> scalar_func;
    scalar_func.set_python_function_from_string(python_function, "func_xyz");

    EXPECT_EQ( 0, scalar_func.value(point_1));
    EXPECT_EQ( 6, scalar_func.value(point_2));


}


#endif // HAVE_PYTHON

