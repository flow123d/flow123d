/*
 * type_use_case_test.cpp
 *
 *  Created on: May 4, 2012
 *      Author: jb
 */



#include <gtest/gtest.h>

#include <input/input_type.hh>
#include <input/json_to_storage.hh>
#include <input/accessors.hh>


class Equation {
public:
    static Input::Type::AbstractRecord &get_input_type();

};


class EquationA : public Equation {
public:
    static Input::Type::Record &get_input_type();
    EquationA(Input::Record rec);
};


class EquationB : public Equation {
public:
    static Input::Type::Record &get_input_type();
    EquationB(Input::Record rec);
};


/**
 * gtest fixture class
 */
class Application : public testing::Test{
public:
    static Input::Type::Record &get_input_type();

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

Input::Type::AbstractRecord &Equation::get_input_type() {
    using namespace Input::Type;
    static AbstractRecord abstr_rec("AbstractEquation","Abstract input Record type for any equation.");
    if (! abstr_rec.is_finished() ) {
        // keys that will be derived by every equation, but their type can be overridden
        abstr_rec.declare_key("mesh",FileName::input(),Default::obligatory(),"");
        abstr_rec.finish(); // finish declaration of keys to allow equations to derive keys

        // declare descendants, they automatically register themself to the abstract record
        EquationA::get_input_type();
        EquationB::get_input_type();

        // finish registering of descendants
        abstr_rec.no_more_descendants();
    }
    return abstr_rec;
}



Input::Type::Record &EquationA::get_input_type() {
    using namespace Input::Type;
    static Record input_record("EquationA", "For example explicit transport equation solver.");

    if (! input_record.is_finished()) {
        input_record.derive_from( Equation::get_input_type() );
        input_record.declare_key("parameter_a", Double(), "");
        input_record.finish();
    }
    return input_record;
}



EquationA::EquationA(Input::Record rec) {
    using namespace Input;

    string mesh_file = rec.val<FilePath>("mesh");
    EXPECT_EQ("/root/some.msh", mesh_file);

    Iterator<double> it = rec.find<double>("parameter_a");
    EXPECT_TRUE(it);
    EXPECT_EQ(3.14, *it);
}


Input::Type::Record &EquationB::get_input_type() {
    using namespace Input::Type;
    static Record input_record("EquationB", "For example implicit transport equation solver.");

    if (! input_record.is_finished()) {
        input_record.derive_from( Equation::get_input_type() );
        input_record.declare_key("parameter_b", Integer(), Default("111"), "");
        input_record.declare_key("default_str", String(), Default("str value"), "" );
        input_record.declare_key("substances", Array( String() ), Default::obligatory(), "" );
        input_record.finish();
    }
    return input_record;
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



Input::Type::Record &Application::get_input_type() {
    using namespace Input::Type;
    static Record input_record("Application", "Root record of the whole application.");

    if (! input_record.is_finished()) {
        // Array of equations with types given by method of class Equation
        Array eq_array( Equation::get_input_type(), 1, 10 );
        input_record.declare_key("equations", eq_array, Default::obligatory(), "");
        input_record.finish();
    }
    return input_record;
}



/**
 * In place input file.
 */
const string flow_json = R"JSON(
{
global_mesh = "some.msh",

equations= 
[
   {
      TYPE="EquationA",
      mesh={REF:"/global_mesh"},
      parameter_a=3.14 
   },
   {
      TYPE="EquationB",
      mesh={REF:"/global_mesh"},
      parameter_b=314,
      substances = [ "Rn", "Cs", "I", "C" ]
   }   
]
}
)JSON";


void Application::SetUp() {
    std::stringstream in_stream(flow_json);
    json_reader.read_stream(in_stream, get_input_type() );
}



void Application::TearDown() {

}

