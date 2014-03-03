/*
 * type_use_case_test.cpp
 *
 *  Created on: May 4, 2012
 *      Author: jb
 */



#include <flow_gtest.hh>

#include <fstream>
#include <input/input_type.hh>
#include <input/json_to_storage.hh>
#include <input/accessors.hh>


class Equation {
public:
    static Input::Type::AbstractRecord input_type;

};


class EquationA : public Equation {
public:
    static Input::Type::Record input_type;
    EquationA(Input::Record rec);
};


class EquationB : public Equation {
public:
    static Input::Type::Record input_type;
    EquationB(Input::Record rec);
};


/**
 * gtest fixture class
 */
class Application : public testing::Test{
public:
    static Input::Type::Record input_type;

    inline Input::Record input()
    { return  json_reader.get_root_interface<Input::Record>(); }
    virtual void SetUp();
    virtual void TearDown();
private:
    Input::JSONToStorage json_reader;

};

TEST_F(Application, init) {
    using namespace Input;

    FilePath::set_io_dirs("/root","/root","variant", "/output");

    Array eq_arr = input().val<Array>("equations");

    Iterator<AbstractRecord> it = eq_arr.begin<AbstractRecord>();

    Record x = *it;

    EquationA first_eq( x );

    ++it;

    EquationB second_eq( *it );

}

/******************************************************************************************
 * Implementation - typically in *.cc file.
 *
 */

namespace it = Input::Type;

it::Record Application::input_type = it::Record("Application", "Root record of the whole application.")
    // Array of equations with types given by method of class Equation
    .declare_key("equations", it::Array( Equation::input_type, 1, 10 ), it::Default::obligatory(), "");




it::Record EquationA::input_type = it::Record("EquationA", "For example explicit transport equation solver.")
    .derive_from( Equation::input_type )
    .declare_key("parameter_a", it::Double(), "");


it::Record EquationB::input_type = it::Record("EquationB", "For example implicit transport equation solver.")
    .derive_from( Equation::input_type )
    .declare_key("parameter_b", it::Integer(), it::Default("111"), "")
    .declare_key("default_str", it::String(), it::Default("str value"), "" )
    .declare_key("substances", it::Array( it::String() ), it::Default::obligatory(), "" );

it::AbstractRecord Equation::input_type = it::AbstractRecord("AbstractEquation","Abstract input Record type for any equation.")
	// keys that will be derived by every equation, but their type can be overridden
    .declare_key("mesh",it::FileName::input(),it::Default::obligatory(),"");




EquationA::EquationA(Input::Record rec) {
    using namespace Input;

    string mesh_file = rec.val<FilePath>("mesh");
    EXPECT_EQ("/root/some.msh", mesh_file);

    Iterator<double> it = rec.find<double>("parameter_a");
    EXPECT_TRUE(it);
    EXPECT_EQ(3.14, *it);
}






EquationB::EquationB(Input::Record rec) {
    using namespace Input;

    string mesh_file = rec.val<FilePath>("mesh");
    EXPECT_EQ("/root/some.msh", mesh_file);
    int param = rec.val<int>("parameter_b");
    EXPECT_EQ(314, param);

    EXPECT_EQ("str value", rec.val<string>("default_str"));
    EXPECT_TRUE( rec.find<string>("default_str"));

    Array array( rec.val<Array>("substances") );
}






void Application::SetUp() {
    std::string f_name = string(UNIT_TESTS_SRC_DIR) + "/input/type_use_case_test.con";
    std::ifstream in_stream(f_name.c_str());
    json_reader.read_stream(in_stream, input_type );
}



void Application::TearDown() {

}

