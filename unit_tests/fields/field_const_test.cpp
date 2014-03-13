/*
 * field_const_test.cpp
 *
 *  Created on: Dec 18, 2012
 *      Author: jb
 */



#include <flow_gtest.hh>


#include "fields/field_constant.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/json_to_storage.hh"



string input = R"INPUT(
{   
   conductivity_3d={
       TYPE="FieldConstant",
       value=[1,2, 3]
   },
   init_conc={
       TYPE="FieldConstant",
       value=[1.2, 2.3, 3.4]
   }
}
)INPUT";


TEST(FieldConst, read_from_input) {
    typedef FieldBase<3, FieldValue<3>::TensorFixed > TensorField;
    typedef FieldBase<3, FieldValue<3>::Vector > VectorField;

    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Input::Type::Record  rec_type("FieldConstTest","");
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

        result = conc->value( point_1, elm);
        EXPECT_DOUBLE_EQ( 1.2 , result[0]);
        EXPECT_DOUBLE_EQ( 2.3, result[1]);
        EXPECT_DOUBLE_EQ( 3.4, result[2]);

        result = conc->value( point_2, elm);
        EXPECT_DOUBLE_EQ( 1.2 , result[0]);
        EXPECT_DOUBLE_EQ( 2.3, result[1]);
        EXPECT_DOUBLE_EQ( 3.4, result[2]);
    }

    auto cond=TensorField::function_factory(in_rec.val<Input::AbstractRecord>("conductivity_3d"));
    {
        arma::mat::fixed<3,3> result;

        result = cond->value( point_2, elm);
        arma::umat match = ( arma::mat::fixed<3,3>("1 0 0; 0 2 0; 0 0 3") == result );
        EXPECT_TRUE( match.max());

        match = ( arma::mat::fixed<3,3>("1 0 0; 0 2 0; 0 0 3") == cond->value( point_1, elm) );
        EXPECT_TRUE( match.max());

    }
}
