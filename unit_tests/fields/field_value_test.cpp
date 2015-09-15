/*
 * field_value_test.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: jb
 */


#include <flow_gtest.hh>
#include <fields/field_values.hh>

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
int_scalar=23,
enum_scalar="one",

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
		param_vec.push_back( std::make_pair("element_input_type", boost::make_shared<Input::Type::Selection>(*sel)) );
	} else {
		param_vec.push_back( std::make_pair("element_input_type", boost::make_shared<typename Value::ElementInputType>()) );
	}

	static auto value_type = Value::get_input_type();
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
    	.declare_key("double_scalar",get_instance< FieldValue_<1,1,double> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("int_scalar",get_instance< FieldValue_<1,1,int> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("enum_scalar",aux_sel, Input::Type::Default::obligatory(),"" )

    	.declare_key("double_fix_vector_full",get_instance< FieldValue_<3,1,double> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("int_fix_vector_full",get_instance< FieldValue_<3,1,int> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("double_fix_vector_const",get_instance< FieldValue_<3,1,double> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("int_fix_vector_const",get_instance< FieldValue_<3,1,int> >().first, Input::Type::Default::obligatory(),"" )

    	.declare_key("double_vector_full",get_instance< FieldValue_<0,1,double> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("int_vector_full",get_instance< FieldValue_<0,1,int> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("enum_vector_full",get_instance< FieldValue_<0,1,FieldEnum> >(&aux_sel).first, Input::Type::Default::obligatory(),"" )
    	.declare_key("double_vector_const",get_instance< FieldValue_<0,1,double> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("int_vector_const",get_instance< FieldValue_<0,1,int> >().first, Input::Type::Default::obligatory(),"" )

    	.declare_key("double_fix_tensor_full",get_instance< FieldValue_<2,3,double> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("double_fix_tensor_symm",get_instance< FieldValue_<2,2,double> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("double_fix_tensor_diag",get_instance< FieldValue_<2,2,double> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("double_fix_tensor_cdiag",get_instance< FieldValue_<2,2,double> >().first, Input::Type::Default::obligatory(),"" )

    	.close();

    // read input string
    Input::ReaderToStorage reader( input, rec_type, Input::FileFormat::format_JSON );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();


    {
        typedef FieldValue_<1,1,double> T; T::return_type x_val; T val(x_val);
        val.init_from_input(in_rec.val<double>("double_scalar"));
        EXPECT_EQ(T::return_type(val), 1.3);
    }
    {
        typedef FieldValue_<1,1,int> T; T::return_type x_val; T val(x_val);
        val.init_from_input(in_rec.val<int>("int_scalar"));
        EXPECT_EQ(T::return_type(val), 23);
    }
    {
        typedef FieldValue_<1,1,FieldEnum> T; T::return_type x_val; T val(x_val);
        val.init_from_input(in_rec.val<Input::Enum>("enum_scalar"));
        EXPECT_EQ(T::return_type(val), 1);
    }


    {
        typedef FieldValue_<3,1,double> T; T::return_type x_val; T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_fix_vector_full"));
        T::return_type expected("1.2 3.4 5.6");
        EXPECT_TRUE( arma::min(expected == T::return_type(val)) );
    }
    {
        typedef FieldValue_<3,1,int> T; T::return_type x_val; T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("int_fix_vector_full"));
        EXPECT_TRUE( arma::min(T::return_type("1 2 3") == T::return_type(val)) );
    }
    {
        typedef FieldValue_<3,1,double> T; T::return_type x_val; T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_fix_vector_const"));
        EXPECT_TRUE( arma::min(T::return_type("1.3 1.3 1.3") == T::return_type(val)) );
    }
    {
        typedef FieldValue_<3,1,int> T; T::return_type x_val; T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("int_fix_vector_const"));
        EXPECT_TRUE( arma::min(T::return_type("23 23 23") == T::return_type(val)) );
    }


    {
        typedef FieldValue_<0,1,double> T; T::return_type x_val(2); T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_vector_full"));
        EXPECT_TRUE( arma::min(T::return_type("1.2 3.4") == T::return_type(val)) );
    }
    {
        typedef FieldValue_<0,1,int> T; T::return_type x_val(4); T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("int_vector_full"));
        EXPECT_TRUE( arma::min(T::return_type("1 2 3 4") == T::return_type(val)) );
    }
    {
        typedef FieldValue_<0,1,double> T; T::return_type x_val(3); T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_vector_const"));
        EXPECT_TRUE( arma::min(T::return_type("1.2 1.2 1.2") == T::return_type(val)) );
    }
    {
        typedef FieldValue_<0,1,FieldEnum> T; T::return_type x_val(3); T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("enum_vector_full"));
        EXPECT_TRUE( arma::min(T::return_type("0 1 1") == T::return_type(val)) );
    }
    {
        typedef FieldValue_<0,1,int> T; T::return_type x_val(3); T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("int_vector_const"));
        EXPECT_TRUE( arma::min(T::return_type("23 23 23") == T::return_type(val)) );
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

string formula_input = R"INPUT(
{   
double_scalar="x",

double_fix_vector_full=["x", "y", "z"],
double_fix_vector_const="x",

double_vector_full=["x","y"],
double_vector_const="x",

double_fix_tensor_full=[ ["x", "y", "z"], ["x*x", "y*y", "z*z"] ],
double_fix_tensor_symm=[ "x*x", "x*y", "y*y"],
double_fix_tensor_diag=[ "x*x", "y*y"],
double_fix_tensor_cdiag="x*y*z"
}
)INPUT";

TEST(FieldValue_, string_values_init_from_input) {
    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Input::Type::Record rec_type = Input::Type::Record("FieldValueTest","")
		.declare_key("double_scalar",get_instance< FieldValue_<1,1,std::string> >().first, Input::Type::Default::obligatory(),"" )

    	.declare_key("double_fix_vector_full",get_instance< FieldValue_<3,1,std::string> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("double_fix_vector_const",get_instance< FieldValue_<3,1,std::string> >().first, Input::Type::Default::obligatory(),"" )

    	.declare_key("double_vector_full",get_instance< FieldValue_<0,1,std::string> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("double_vector_const",get_instance< FieldValue_<0,1,std::string> >().first, Input::Type::Default::obligatory(),"" )

    	.declare_key("double_fix_tensor_full",get_instance< FieldValue_<2,3,std::string> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("double_fix_tensor_symm",get_instance< FieldValue_<2,2,std::string> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("double_fix_tensor_diag",get_instance< FieldValue_<2,2,std::string> >().first, Input::Type::Default::obligatory(),"" )
    	.declare_key("double_fix_tensor_cdiag",get_instance< FieldValue_<2,2,std::string> >().first, Input::Type::Default::obligatory(),"" )

    	.close();

    // read input string
    Input::ReaderToStorage reader( formula_input, rec_type, Input::FileFormat::format_JSON );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();


    {
        typedef FieldValue_<1,1,std::string> T; T::return_type x_val(1,1); T val(x_val);
        val.init_from_input(in_rec.val<std::string>("double_scalar"));
        EXPECT_EQ(x_val.at(0,0), "x");
    }


    {
        typedef FieldValue_<3,1,std::string> T; T::return_type x_val(3,1); T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_fix_vector_full"));
        EXPECT_EQ( "x",x_val.at(0));
        EXPECT_EQ( "y",x_val.at(1));
        EXPECT_EQ( "z",x_val.at(2));
    }
    {
        typedef FieldValue_<3,1,std::string> T; T::return_type x_val(3,1); T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_fix_vector_const"));
        EXPECT_EQ( "x",x_val.at(0));
        EXPECT_EQ( "x",x_val.at(1));
        EXPECT_EQ( "x",x_val.at(2));
    }


    {
        typedef FieldValue_<0,1,std::string> T; T::return_type x_val(2,1); T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_vector_full"));
        EXPECT_EQ( "x",x_val.at(0));
        EXPECT_EQ( "y",x_val.at(1));
    }
    {
        typedef FieldValue_<0,1,std::string> T; T::return_type x_val(3,1); T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_vector_const"));
        EXPECT_EQ( "x",x_val.at(0));
        EXPECT_EQ( "x",x_val.at(1));
        EXPECT_EQ( "x",x_val.at(2));
    }


    {
        typedef FieldValue_<2,3,std::string> T; T::return_type x_val(2,3); T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_fix_tensor_full"));
        EXPECT_EQ( "x",x_val.at(0,0));
        EXPECT_EQ( "y",x_val.at(0,1));
        EXPECT_EQ( "z",x_val.at(0,2));
        EXPECT_EQ( "x*x",x_val.at(1,0));
        EXPECT_EQ( "y*y",x_val.at(1,1));
        EXPECT_EQ( "z*z",x_val.at(1,2));
    }
    {
        typedef FieldValue_<2,2,std::string> T; T::return_type x_val(2,2); T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_fix_tensor_symm"));
        EXPECT_EQ( "x*x",x_val.at(0,0));
        EXPECT_EQ( "x*y",x_val.at(0,1));
        EXPECT_EQ( "x*y",x_val.at(1,0));
        EXPECT_EQ( "y*y",x_val.at(1,1));
    }
    {
        typedef FieldValue_<2,2,std::string> T; T::return_type x_val(2,2); T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_fix_tensor_diag"));
        EXPECT_EQ( "x*x",x_val.at(0,0));
        EXPECT_EQ( "0.0",x_val.at(0,1));
        EXPECT_EQ( "0.0",x_val.at(1,0));
        EXPECT_EQ( "y*y",x_val.at(1,1));
    }
    {
        typedef FieldValue_<2,2,std::string> T; T::return_type x_val(2,2); T val(x_val);
        val.init_from_input(in_rec.val<Input::Array>("double_fix_tensor_cdiag"));
        EXPECT_EQ( "x*y*z",x_val.at(0,0));
        EXPECT_EQ( "0.0",x_val.at(0,1));
        EXPECT_EQ( "0.0",x_val.at(1,0));
        EXPECT_EQ( "x*y*z",x_val.at(1,1));
    }
}

