/*
 * type_use_case_test.cpp
 *
 *  Created on: May 4, 2012
 *      Author: jb
 */



#include <gtest/gtest.h>

#include <input/input_type.hh>
#include <input/json_to_storage.hh>
#include <input/input_interface.hh>


class Equation {
public:
    static Input::Type::AbstractRecord &get_input_type();

};


class EquationA : public Equation {
public:
    static Input::Type::Record &get_input_type();
    EquationA(Input::Interface::Record rec);
};


class EquationB : public Equation {
public:
    static Input::Type::Record &get_input_type();
    EquationB(Input::Interface::Record rec);
};


/**
 * gtest fixture class
 */
class Application : public testing::Test{
public:
    static Input::Type::Record &get_input_type();
    inline Input::Interface::Record input()
    { return *( json_reader.get_root_interface<Input::Interface::Record>() ); }
    virtual void SetUp();
    virtual void TearDown();
private:
    Input::JSONToStorage json_reader;

};

TEST_F(Application, init) {
    using namespace Input::Interface;

    Array eq_arr = input().key<Array>("equations");
    Iterator<Record> it = eq_arr.begin<Record>();
    EquationA first_eq( *it );

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
        abstr_rec.declare_key("mesh",FileName(input_file),DefaultValue(DefaultValue::obligatory),"");
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
        input_record.declare_key("parametr_a", Double(), "");
        input_record.finish();
    }
    return input_record;
}



EquationA::EquationA(Input::Interface::Record rec) {
    string mesh_file = rec.key<string>("mesh");
    double param = rec.key<double>("parameter_a");
}


Input::Type::Record &EquationB::get_input_type() {
    using namespace Input::Type;
    static Record input_record("EquationB", "For example implicit transport equation solver.");

    if (! input_record.is_finished()) {
        input_record.derive_from( Equation::get_input_type() );
        input_record.declare_key("parametr_b", Integer(), "");
        input_record.finish();
    }
    return input_record;
}



EquationB::EquationB(Input::Interface::Record rec) {
    string mesh_file = rec.key<string>("mesh");
    int param = rec.key<int>("parameter_a");
}



Input::Type::Record &Application::get_input_type() {
    using namespace Input::Type;
    static Record input_record("Application", "Root record of the whole application.");

    if (! input_record.is_finished()) {
        // Array of equations with types given by method of class Equation
        Array eq_array( Equation::get_input_type(), 1, 10 );
        input_record.declare_key("equations", eq_array, DefaultValue( DefaultValue::obligatory), "");
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

