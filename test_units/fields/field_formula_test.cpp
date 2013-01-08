/*
 * field_formula_test.cpp
 *
 *  Created on: Jan 8, 2013
 *      Author: jb
 */


#include <gtest/gtest.h>


#include "fields/field_constant.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/json_to_storage.hh"


string input = R"INPUT(
{   
   init_conc={ # formula on 2d 
       TYPE="FieldFormula",
       value=["x", "x*y", "y+t"]
   },
   conductivity_3d={ #3x3 tensor
       TYPE="FieldFormula",
       value=["sin(x)+cos(y)","exp(x)+y^2", "base:=(x+y); base+base^2"]
   }
}
)INPUT";


TEST(FieldFormula, read_from_input) {
    typedef FieldBase<2, FieldValue<3>::TensorFixed > TensorField;
    typedef FieldBase<2, FieldValue<3>::Vector > VectorField;

    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Input::Type::Record  rec_type("FieldFormulaTest","");
    rec_type.declare_key("conductivity_3d", TensorField::input_type, Input::Type::Default::obligatory(),"" );
    rec_type.declare_key("init_conc", VectorField::input_type, Input::Type::Default::obligatory(), "" );
    rec_type.finish();

    // read input string
    std::stringstream ss(input);
    Input::JSONToStorage reader;
    reader.read_stream( ss, rec_type );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();

    Point<2> point_1, point_2;
    point_1(0)=1.0; point_1(1)=  2.0;
    point_2(0)= 2.0; point_2(1)= 4.0;
    ElementAccessor<2> elm;

    VectorField  *conc=VectorField::function_factory(in_rec.val<Input::AbstractRecord>("init_conc"), 0.0, 3);
    {
        arma::vec result;

        result = conc->value( point_1, elm);
        EXPECT_DOUBLE_EQ( point_1(0) ,              result[0]);
        EXPECT_DOUBLE_EQ( point_1(0)*point_1(1),    result[1]);
        EXPECT_DOUBLE_EQ( point_1(1),               result[2]);

        result = conc->value( point_2, elm);
        EXPECT_DOUBLE_EQ( point_2(0) ,              result[0]);
        EXPECT_DOUBLE_EQ( point_2(0)*point_2(1),    result[1]);
        EXPECT_DOUBLE_EQ( point_2(1),               result[2]);
    }

    TensorField  *cond=TensorField::function_factory(in_rec.val<Input::AbstractRecord>("conductivity_3d"), 0.0);
    {
        arma::mat::fixed<3,3> result;
        double x,y,base;

        result = cond->value( point_1, elm);
        x=point_1(0); y=point_1(1); base=(x+y);
        EXPECT_DOUBLE_EQ( sin(x)+cos(y),              result.at(0,0));
        EXPECT_DOUBLE_EQ( 0.0,                      result.at(0,1));
        EXPECT_DOUBLE_EQ( 0.0,                      result.at(0,2));

        EXPECT_DOUBLE_EQ( 0.0 ,                     result.at(1,0));
        EXPECT_DOUBLE_EQ( exp(x)+y*y,               result.at(1,1));
        EXPECT_DOUBLE_EQ( 0.0 ,                     result.at(1,2));

        EXPECT_DOUBLE_EQ( 0.0 ,                     result.at(2,0));
        EXPECT_DOUBLE_EQ( 0.0 ,                     result.at(2,1));
        EXPECT_DOUBLE_EQ( base+base*base,               result.at(2,2));

        result = cond->value( point_2, elm);
        x=point_2(0); y=point_2(1); base=(x+y);
        EXPECT_DOUBLE_EQ( sin(x)+cos(y),              result.at(0,0));
        EXPECT_DOUBLE_EQ( 0.0,                      result.at(0,1));
        EXPECT_DOUBLE_EQ( 0.0,                      result.at(0,2));

        EXPECT_DOUBLE_EQ( 0.0 ,                     result.at(1,0));
        EXPECT_DOUBLE_EQ( exp(x)+y*y,               result.at(1,1));
        EXPECT_DOUBLE_EQ( 0.0 ,                     result.at(1,2));

        EXPECT_DOUBLE_EQ( 0.0 ,                     result.at(2,0));
        EXPECT_DOUBLE_EQ( 0.0 ,                     result.at(2,1));
        EXPECT_DOUBLE_EQ( base+base*base,               result.at(2,2));

    }
}
