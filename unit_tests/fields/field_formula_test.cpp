/*
 * field_formula_test.cpp
 *
 *  Created on: Jan 8, 2013
 *      Author: jb
 */


#include <flow_gtest.hh>


#include "fields/field_constant.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/json_to_storage.hh"


string input = R"INPUT(
{   
   init_conc={ // formula on 2d 
       TYPE="FieldFormula",
       value=["x", "x*y", "y+t"]
   },
   conductivity_3d={ // 3x3 tensor
       TYPE="FieldFormula",
       value=["sin(x)+cos(y)","exp(x)+y^2", "base:=(x+y); base+base^2"]
   }
}
)INPUT";


TEST(FieldFormula, read_from_input) {
    typedef FieldBase<3, FieldValue<3>::TensorFixed > TensorField;
    typedef FieldBase<3, FieldValue<3>::Vector > VectorField;

    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Input::Type::Record  rec_type("FieldFormulaTest","");
    rec_type.declare_key("conductivity_3d", TensorField::input_type, Input::Type::Default::obligatory(),"" );
    rec_type.declare_key("init_conc", VectorField::input_type, Input::Type::Default::obligatory(), "" );
    rec_type.finish();

    // read input string
    Input::JSONToStorage reader( input, rec_type );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();

    Space<3>::Point point_1, point_2;
    point_1(0)=1.0; point_1(1)=2.0; point_1(2)=3.0;
    point_2(0)=2.0; point_2(1)=4.0; point_2(2)=6.0;
    ElementAccessor<3> elm;

    auto conc=VectorField::function_factory(in_rec.val<Input::AbstractRecord>("init_conc"), 3);
    {
        arma::vec result;

        conc->set_time(0.0);
        result = conc->value( point_1, elm);
        EXPECT_DOUBLE_EQ( point_1(0) ,              result[0]);
        EXPECT_DOUBLE_EQ( point_1(0)*point_1(1),    result[1]);
        EXPECT_DOUBLE_EQ( point_1(1),               result[2]);

        result = conc->value( point_2, elm);
        EXPECT_DOUBLE_EQ( point_2(0) ,              result[0]);
        EXPECT_DOUBLE_EQ( point_2(0)*point_2(1),    result[1]);
        EXPECT_DOUBLE_EQ( point_2(1),               result[2]);

        conc->set_time(1.0);
        result = conc->value( point_1, elm);
        EXPECT_DOUBLE_EQ( point_1(0) ,              result[0]);
        EXPECT_DOUBLE_EQ( point_1(0)*point_1(1),    result[1]);
        EXPECT_DOUBLE_EQ( point_1(1) +1.0,               result[2]);
    }

    auto cond=TensorField::function_factory(in_rec.val<Input::AbstractRecord>("conductivity_3d"));
    cond->set_time(0.0);
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


string set_time_input = R"INPUT(
[ 
      { TYPE="FieldFormula",  value=["x", "x*y", "y+t"] },
      { TYPE="FieldFormula",  value=["x", "x*y", "y"] },
      { TYPE="FieldFormula",  value=["x+t", "x*y+t", "y+t"] },
      { TYPE="FieldFormula",  value=["x", "x*y", "y"] }
]

)INPUT";


TEST(FieldFormula, set_time) {
    typedef FieldBase<2, FieldValue<3>::Vector > VectorField;

    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Input::Type::Array  input_type(VectorField::input_type);

    // read input string
    Input::JSONToStorage reader( set_time_input, input_type );
    Input::Array in_array=reader.get_root_interface<Input::Array>();

    auto it = in_array.begin<Input::AbstractRecord>();

    {
        auto field=VectorField::function_factory(*it, 3);
        EXPECT_TRUE( field->set_time(1.0) );
        EXPECT_TRUE( field->set_time(2.0) );
    }
    ++it;

    {
        auto field=VectorField::function_factory(*it, 3);
        EXPECT_TRUE( field->set_time(3.0) );
        EXPECT_FALSE( field->set_time(4.0) );
    }
    ++it;

    {
        auto field=VectorField::function_factory(*it, 3);
        EXPECT_TRUE( field->set_time(1.5) );
        EXPECT_TRUE( field->set_time(2.5) );
    }
    ++it;

    {
        auto field=VectorField::function_factory(*it, 3);
        EXPECT_TRUE( field->set_time(0.0) );
        EXPECT_FALSE( field->set_time(2.0) );
    }

}

