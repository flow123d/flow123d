/*
 * field_const_test.cpp
 *
 *  Created on: Dec 18, 2012
 *      Author: jb
 */



#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <memory>


#include "fields/field_constant.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"


FLOW123D_FORCE_LINK_IN_PARENT(field_constant)


string input = R"INPUT(
{   
   tensor1=3.14,
   tensor2=[1,2, 3],
   tensor3=[1,2, 3,4,5,6],
   tensor4=[ [1,2, 3], [4,5,6], [7,8,9]],

   init_conc={
       TYPE="FieldConstant",
       value=[1.2, 2.3, 3.4]
   }
}
)INPUT";


typedef FieldAlgorithmBase<3, FieldValue<3>::TensorFixed > TensorField;
typedef FieldAlgorithmBase<3, FieldValue<3>::VectorFixed > VectorField;


void check_tensor_field(std::shared_ptr<TensorField> field,const  std::string& expected, vector<Space<3>::Point> points,ElementAccessor<3> elem) {
    arma::mat::fixed<3,3> result;
    arma::mat::fixed<3,3> exp(expected);

    for(auto point : points) {
    	result = field->value( point, elem);
    	arma::umat mat_match = (  exp == result );
    	bool mat_equal = arma::min( arma::min(mat_match));
    	if (! mat_equal) cout << "Expected: " << exp << endl << "Result: "<< result << endl;
    	EXPECT_TRUE( mat_equal );
    }
}


TEST(FieldConst, read_from_input) {

    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Input::Type::Record rec_type = Input::Type::Record("FieldConstTest","")
        .declare_key("tensor1", TensorField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
        .declare_key("tensor2", TensorField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
        .declare_key("tensor3", TensorField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
        .declare_key("tensor4", TensorField::get_input_type_instance(), Input::Type::Default::obligatory(),"" )
        //.declare_key("init_conc", VectorField::get_input_type_instance(), Input::Type::Default::obligatory(), "" )
        .close();

    // read input string
    Input::ReaderToStorage reader( input, rec_type, Input::FileFormat::format_JSON );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();

    Space<3>::Point point_1, point_2;
    point_1(0)=1.0; point_1(1)=2.0; point_1(2)=3.0;
    point_2(0)=2.0; point_2(1)=4.0; point_2(2)=6.0;
    ElementAccessor<3> elm;

    /*
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
    }*/

    auto tensor1=TensorField::function_factory(in_rec.val<Input::AbstractRecord>("tensor1"));
    check_tensor_field(tensor1, "3.14 0 0; 0 3.14 0; 0 0 3.14", {point_1, point_2}, elm);

    auto tensor2=TensorField::function_factory(in_rec.val<Input::AbstractRecord>("tensor2"));
    check_tensor_field(tensor2, "1 0 0; 0 2 0; 0 0 3", {point_1, point_2}, elm);

    auto tensor3=TensorField::function_factory(in_rec.val<Input::AbstractRecord>("tensor3"));
    check_tensor_field(tensor3, "1 2 3; 2 4 5; 3 5 6", {point_1, point_2}, elm);

    auto tensor4=TensorField::function_factory(in_rec.val<Input::AbstractRecord>("tensor4"));
    check_tensor_field(tensor4, "1 2 3; 4 5 6; 7 8 9", {point_1, point_2}, elm);
}

