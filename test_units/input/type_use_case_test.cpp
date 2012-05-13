/*
 * type_use_case_test.cpp
 *
 *  Created on: May 4, 2012
 *      Author: jb
 */



#include <gtest/gtest.h>

#include <input/input_type.hh>

class Application {
    static Input::Type::TypeBase get_input_type() {
        using namespace Input::Type;
        static Record input_record("Application", "Root record of the whole application.");
        Array eq_array( Equation::get_input_type(), 1, 10 );
        input_record.declare_key("equations", eq_array, DefaultValue( DefalutValue::obligatory), "");
    }
};


class Equation {

    static Input::Type::TypeBase get_input_type() {

    }

};


class EquationA : public Equation {

    static Input::Type::TypeBase get_input_type() {

    }

};

class EquationB : public Equation {

    static Input::Type::TypeBase get_input_type() {

    }

};
