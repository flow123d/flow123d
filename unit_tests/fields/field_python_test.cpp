/*
 * python_function_test.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: jb
 */



#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <string>
#include <cmath>



#include "system/global_defs.h"



#include "system/python_loader.hh"
#include "fields/field_python.hh"
#include "tools/unit_si.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"

#include <pybind11.h>
#include <embed.h> // everything needed for embedding

using namespace std;


string test_pybind = R"CODE(
test = testClass()
test.testMethod()
)CODE";

string python_function = R"CODE(
import math

def func_xyz(x,y,z):
    return ( x*y*z , )     # one value tuple

def func_circle(r,phi,n):
    return ( r * math.cos(phi), r * math.sin(phi), 1 )
)CODE";

string python_call_object_err = R"CODE(
import math

def func_xyz(x,y,z,a):
    return ( x*y*z+a , )     # one value tuple
)CODE";


string input = R"INPUT(
{   
   field_string={
       TYPE="FieldPython",
       function="func_circle",
       script_string="import math\ndef func_circle(r,phi,n): return ( r * math.cos(phi), r * math.sin(phi), 1 )"
   },
   field_string_unit_conversion={
       TYPE="FieldPython",
       function="func_circle",
       script_string="import math\ndef func_circle(r,phi,n): return ( r * math.cos(phi), r * math.sin(phi), 1 )",
       unit="cm"
   },
   field_file={
       TYPE="FieldPython",
       function="func_xyz",
       script_file="fields/field_python_script.py"
   }
}
)INPUT";


//TEST(FieldPython, vector_3D) {
//
//    double pi = 4.0 * atan(1);
//
//    Space<3>::Point point_1, point_2;
//    point_1(0)=1.0; point_1(1)= pi / 2.0; point_1(2)=1.0;
//    point_2(0)= sqrt(2.0); point_2(1)= 3.0 * pi / 4.0; point_2(2)= pi / 2.0;
//
//    FieldPython<3, FieldValue<3>::Vector > vec_func;
//    vec_func.set_python_field_from_string(python_function, "func_circle");
//    ElementAccessor<3> elm;
//
//    arma::vec3 result;
//    {
//    result = vec_func.value( point_1, elm);
//    EXPECT_DOUBLE_EQ( cos(pi /2.0 ) , result[0]); // should be 0.0
//    EXPECT_DOUBLE_EQ( 1, result[1]);
//    EXPECT_DOUBLE_EQ( 1, result[2]);
//    }
//
//    {
//    result = vec_func.value( point_2, elm);
//    EXPECT_DOUBLE_EQ( -1, result[0]);
//    EXPECT_DOUBLE_EQ( 1, result[1]);
//    EXPECT_DOUBLE_EQ( 1, result[2]);
//    }
//}


//TEST(FieldPython, double_3D) {
//    Space<3>::Point point_1, point_2;
//    point_1(0)=1; point_1(1)=0; point_1(2)=0;
//    point_2(0)=1; point_2(1)=2; point_2(2)=3;
//
//    ElementAccessor<3> elm;
//    FieldPython<3, FieldValue<3>::Scalar> scalar_func;
//    scalar_func.set_python_field_from_string(python_function, "func_xyz");
//
//    EXPECT_EQ( 0, scalar_func.value(point_1, elm));
//    EXPECT_EQ( 6, scalar_func.value(point_2, elm));
//
//
//}


// TODO Fix test
//TEST(FieldPython, read_from_input) {
//    typedef FieldAlgorithmBase<3, FieldValue<3>::Vector > VectorField;
//    typedef FieldAlgorithmBase<3, FieldValue<3>::Scalar > ScalarField;
//    double pi = 4.0 * atan(1);
//
//    // setup FilePath directories
//    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
//
//    Input::Type::Record rec_type = Input::Type::Record("FieldPythonTest","")
//        .declare_key("field_string", VectorField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
//        .declare_key("field_string_unit_conversion", VectorField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
//        .declare_key("field_file", ScalarField::get_input_type_instance(), Input::Type::Default::obligatory(), "" )
//        .close();
//
//    // read input string
//    Input::ReaderToStorage reader( input, rec_type, Input::FileFormat::format_JSON );
//    Input::Record in_rec=reader.get_root_interface<Input::Record>();
//    UnitSI unit = UnitSI().m();
//    FieldAlgoBaseInitData init_data("field_python", 3, unit);
//
//    auto flux=VectorField::function_factory(in_rec.val<Input::AbstractRecord>("field_string"), init_data);
//    {
//        Space<3>::Point point_1, point_2;
//        point_1(0)=1.0; point_1(1)= pi / 2.0; point_1(2)=1.0;
//        point_2(0)= sqrt(2.0); point_2(1)= 3.0 * pi / 4.0; point_2(2)= pi / 2.0;
//
//        ElementAccessor<3> elm;
//        arma::vec3 result;
//
//        result = flux->value( point_1, elm);
//        EXPECT_DOUBLE_EQ( cos(pi /2.0 ) , result[0]); // should be 0.0
//        EXPECT_DOUBLE_EQ( 1, result[1]);
//        EXPECT_DOUBLE_EQ( 1, result[2]);
//
//        result = flux->value( point_2, elm);
//        EXPECT_DOUBLE_EQ( -1, result[0]);
//        EXPECT_DOUBLE_EQ( 1, result[1]);
//        EXPECT_DOUBLE_EQ( 1, result[2]);
//    }
//
//    auto flux_unit_conv=VectorField::function_factory(in_rec.val<Input::AbstractRecord>("field_string_unit_conversion"), init_data);
//    {
//        Space<3>::Point point_1, point_2;
//        point_1(0)=1.0; point_1(1)= pi / 2.0; point_1(2)=1.0;
//        point_2(0)= sqrt(2.0); point_2(1)= 3.0 * pi / 4.0; point_2(2)= pi / 2.0;
//
//        ElementAccessor<3> elm;
//        arma::vec3 result;
//
//        result = flux_unit_conv->value( point_1, elm);
//        EXPECT_DOUBLE_EQ( 0.01*cos(pi /2.0 ) , result[0]); // should be 0.0
//        EXPECT_DOUBLE_EQ( 0.01, result[1]);
//        EXPECT_DOUBLE_EQ( 0.01, result[2]);
//
//        result = flux_unit_conv->value( point_2, elm);
//        EXPECT_DOUBLE_EQ( -0.01, result[0]);
//        EXPECT_DOUBLE_EQ( 0.01, result[1]);
//        EXPECT_DOUBLE_EQ( 0.01, result[2]);
//    }
//
//    auto conc=ScalarField::function_factory(in_rec.val<Input::AbstractRecord>("field_file"), init_data);
//    {
//        Space<3>::Point point_1, point_2;
//        point_1(0)=1; point_1(1)=0; point_1(2)=0;
//        point_2(0)=1; point_2(1)=2; point_2(2)=3;
//        ElementAccessor<3> elm;
//
//        EXPECT_EQ( 0, conc->value(point_1, elm));
//        EXPECT_EQ( 6, conc->value(point_2, elm));
//    }
//
//
//}

//TEST(FieldPython, python_exception) {
//    FieldPython<3, FieldValue<3>::Scalar> scalar_func;
//	EXPECT_THROW_WHAT( { scalar_func.set_python_field_from_string(python_function, "func_xxx"); }, PythonLoader::ExcPythonError,
//        "func_xxx");
//
//}


//TEST(FieldPython, call_object_error) {
//    FieldPython<3, FieldValue<3>::Scalar> scalar_func;
//	EXPECT_THROW_WHAT( { scalar_func.set_python_field_from_string(python_call_object_err, "func_xyz"); }, PythonLoader::ExcPythonError,
//        "missing 1 required positional argument: 'a'");
//
//}


TEST(FieldPython, new_assembly) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    PythonLoader::initialize();

    FieldPython<3, FieldValue<3>::Vector> vec_func;
    vec_func.set_python_field_from_class("fields/field_python_class.py", "PythonAsm");
}

