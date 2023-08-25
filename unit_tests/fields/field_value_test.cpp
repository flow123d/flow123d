/*
 * field_value_test.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: jb
 */


#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <fields/field_values.hh>
#include <fields/field_constant.hh>
#include "arma_expect.hh"
#include <system/armor.hh>

#include <iostream>
using namespace std;

TEST(FieldValue_, all) {

    {
        typedef FieldValue_<2,2, double> T;
        T::return_type x_val;
        T val(x_val);
        val(0,0)=3.0; val(0,1)=4.0;
        EXPECT_EQ(3.0,  val(0,0) );
        EXPECT_EQ(4.0,  val(0,1) );

        T::return_type r_val = val;
        EXPECT_EQ(3.0,  r_val(0,0) );
        EXPECT_EQ(4.0,  r_val(0,1) );

    }

}


TEST(FieldValue_, construction_from_raw) {

    double raw_data[6]={1,2,3,4,5,6};
    FieldEnum ui_raw[6]={10, 20, 30 ,40, 50, 60};
    int i_raw[6]={10, 20, 30 ,40, 50, 60};
    // scalars
    {
        typedef FieldValue_<1,1,double> T; T::return_type x_val;
        x_val=0;
        const T::return_type & val = T::from_raw(x_val, raw_data);
        EXPECT_DOUBLE_EQ(1, double(val));
    }
    {
        typedef FieldValue_<1,1,int> T; T::return_type x_val;
        x_val=0;
        const T::return_type & val = T::from_raw(x_val, i_raw);
        EXPECT_DOUBLE_EQ(10, int(val));
    }
    {
        typedef FieldValue_<1,1,FieldEnum> T; T::return_type x_val;
        x_val=0;
        const T::return_type & val = T::from_raw(x_val, ui_raw);
        EXPECT_DOUBLE_EQ(10, FieldEnum(val));
    }

    // vectors
    {
        typedef FieldValue_<3,1,double> T; T::return_type x_val;
        x_val.zeros();
        const T::return_type & val = T::from_raw(x_val, raw_data);
        EXPECT_TRUE( arma::min(T::return_type("1 2 3") == T::return_type(val)) );
    }
    {
        typedef FieldValue_<0,1,double> T; T::return_type x_val(2);
        x_val.zeros();
        const T::return_type & val = T::from_raw(x_val, raw_data);
        EXPECT_TRUE( arma::min(T::return_type("1 2") == T::return_type(val)) );
    }
    {
        typedef FieldValue_<0,1,FieldEnum> T; T::return_type x_val(2);
        x_val.zeros();
        const T::return_type & val = T::from_raw(x_val, ui_raw);
        cout << T::return_type(val);
        EXPECT_TRUE( arma::min(T::return_type("10 20") == T::return_type(val)) );
    }
    {
        typedef FieldValue_<0,1,int> T; T::return_type x_val(2);
        x_val.zeros();
        const T::return_type & val = T::from_raw(x_val, i_raw);
        cout << T::return_type(val);
        EXPECT_TRUE( arma::min(T::return_type("10 20") == T::return_type(val)) );
    }

    // tensor
    {
        typedef FieldValue_<2,3,double> T; T::return_type x_val;
        x_val.zeros();
        const T::return_type & val = T::from_raw(x_val, raw_data);
        arma::umat match = (T::return_type("1 3 5; 2 4 6") == T::return_type(val));

        EXPECT_TRUE( match.min());
    }
}



string input = R"INPUT(
{   
double_scalar=1.3,

double_fix_vector_full=[1.2, 3.4, 5.6],
int_fix_vector_full=[1,2,3],
double_fix_vector_const=1.3,
int_fix_vector_const=23,

double_vector_full=[1.2,3.4],
int_vector_full=[1,2,3,4],
enum_vector_full=["zero", "one", "one"],
double_vector_const=1.2,
int_vector_const=23,

double_fix_tensor_full=[ [1.1, 1.2, 1.3], [2.1, 2.2, 2.3] ],
double_fix_tensor_symm=[ 1, 2, 3],
double_fix_tensor_diag=[1,2],
double_fix_tensor_cdiag=1.3
}
)INPUT";



#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"


template <class Value>
Input::Type::TypeBase::MakeInstanceReturnType get_instance(const Input::Type::Selection *sel = NULL) {
	std::vector<Input::Type::TypeBase::ParameterPair> param_vec;
	if (sel) {
		param_vec.push_back( std::make_pair("element_input_type", std::make_shared<Input::Type::Selection>(*sel)) );
	} else {
		param_vec.push_back( std::make_pair("element_input_type", std::make_shared<typename Value::ElementInputType>()) );
	}

	static auto value_type = FieldConstant<3, Value>::get_tensor_input_type();
	return value_type.make_instance(param_vec);
}


TEST(FieldValue_, init_from_input) {
    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Input::Type::Selection aux_sel = Input::Type::Selection("AuxSel")
    	.add_value(0,"zero","")
    	.add_value(1,"one","")
		.close();

    Input::Type::Record rec_type = Input::Type::Record("FieldValueTest","")
    	.declare_key("double_scalar",get_instance< FieldValue_<3,3,double> >().first, Input::Type::Default::obligatory(),"" )

    	.declare_key("double_fix_vector_full",get_instance< FieldValue_<3,3,double> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("double_fix_vector_const",get_instance< FieldValue_<3,3,double> >().first, Input::Type::Default::obligatory(),"" )

    	.declare_key("double_vector_full",get_instance< FieldValue_<3,3,double> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("double_vector_const",get_instance< FieldValue_<3,3,double> >().first, Input::Type::Default::obligatory(),"" )

    	.declare_key("double_fix_tensor_full",get_instance< FieldValue_<3,3,double> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("double_fix_tensor_symm",get_instance< FieldValue_<3,3,double> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("double_fix_tensor_diag",get_instance< FieldValue_<3,3,double> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("double_fix_tensor_cdiag",get_instance< FieldValue_<3,3,double> >().first, Input::Type::Default::obligatory(),"" )

    	.close();

    // read input string
    Input::ReaderToStorage reader( input, rec_type, Input::FileFormat::format_JSON );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();


    {
        typedef FieldValue_<1,1,double> T; T::return_type x_val; T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_scalar"));
        EXPECT_EQ(T::return_type(val), 1.3);
    }


    {
        typedef FieldValue_<3,1,double> T; T::return_type x_val; T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_fix_vector_full"));
        T::return_type expected("1.2 3.4 5.6");
        EXPECT_TRUE( arma::min(expected == T::return_type(val)) );
    }
    {
        typedef FieldValue_<3,1,double> T; T::return_type x_val; T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_fix_vector_const"));
        EXPECT_TRUE( arma::min(T::return_type("1.3 1.3 1.3") == T::return_type(val)) );
    }


    {
        typedef FieldValue_<0,1,double> T; T::return_type x_val(2); T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_vector_full"));
        EXPECT_TRUE( arma::min(T::return_type("1.2 3.4") == T::return_type(val)) );
    }
    {
        typedef FieldValue_<0,1,double> T; T::return_type x_val(3); T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_vector_const"));
        EXPECT_TRUE( arma::min(T::return_type("1.2 1.2 1.2") == T::return_type(val)) );
    }


    {
        typedef FieldValue_<2,3,double> T; T::return_type x_val; T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_fix_tensor_full"));
        arma::umat match = (T::return_type("1.1 1.2 1.3; 2.1 2.2 2.3") == T::return_type(val));
        EXPECT_TRUE( match.min());
    }
    {
        typedef FieldValue_<2,2,double> T; T::return_type x_val; T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_fix_tensor_symm"));
        arma::umat match = (T::return_type("1 2; 2 3") == T::return_type(val));
        EXPECT_TRUE( match.min());
    }
    {
        typedef FieldValue_<2,2,double> T; T::return_type x_val; T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_fix_tensor_diag"));
        arma::umat match = (T::return_type("1 0; 0 2") == T::return_type(val));
        EXPECT_TRUE( match.min());
    }
    {
        typedef FieldValue_<2,2,double> T; T::return_type x_val; T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_fix_tensor_cdiag"));
        arma::umat match = (T::return_type("1.3 0; 0 1.3") == T::return_type(val));
        EXPECT_TRUE( match.min());
    }
}


TEST(FieldValue_, get_from_array) {
    {  // scalar
        Armor::Array<double> arr(1, 1, 1);
        Armor::ArmaMat<double, 1, 1> m1{1.5}; // first item
        arr.set(0) = m1;

        typedef FieldValue_<1,1,double> T;
        EXPECT_EQ(T::get_from_array(arr, 0), 1.5);
    }
    {  // vector
        Armor::Array<double> arr(3, 1, 1);
        Armor::ArmaMat<double, 3, 1> m1{1, 2, 3}; // first item
        arr.set(0) = m1;

        typedef FieldValue_<3,1,double> T;
        arma::vec3 expected = {1, 2, 3};
        EXPECT_ARMA_EQ(T::get_from_array(arr, 0), expected);
    }
    {  // tensor
        Armor::Array<double> arr(3, 3, 1);
        Armor::ArmaMat<double, 3, 3> m1{1, 2, 3, 4, 5, 6, 7, 8, 9}; // first item
        arr.set(0) = m1;

        typedef FieldValue_<3,3,double> T;
        arma::mat33 expected = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        EXPECT_ARMA_EQ(T::get_from_array(arr, 0), expected);
    }
    {  // FieldEnum
        Armor::Array<unsigned int> arr(1, 1, 1);
        Armor::ArmaMat<unsigned int, 1, 1> m1{0}; // first item
        arr.set(0) = m1;

        typedef FieldValue_<1,1,FieldEnum> T;
        EXPECT_EQ(T::get_from_array(arr, 0), 0);
    }
    {  // int
        Armor::Array<int> arr(1, 1, 1);
        Armor::ArmaMat<int, 1, 1> m1{1}; // first item
        arr.set(0) = m1;

        typedef FieldValue_<1,1,int> T;
        EXPECT_EQ(T::get_from_array(arr, 0), 1);
    }
}
