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
#include <input/factory.hh>


class Equation {
public:
    static Input::Type::AbstractRecord & get_input_type();

};


class EquationA : public Equation {
public:
	typedef Equation FactoryBaseType;

    static Input::Type::Record & get_input_rec();
    EquationA(Input::Record rec);

private:
    /// Registrar of class to factory
    static const int registrar;
};


class EquationB : public Equation {
public:
	typedef Equation FactoryBaseType;

	static Input::Type::Record & get_input_rec();
    EquationB(Input::Record rec);

private:
    /// Registrar of class to factory
    static const int registrar;
};


/**
 * gtest fixture class
 */
class Application : public testing::Test{
public:
    static Input::Type::Record & get_input_type();

    inline Input::Record input()
    {
    	root_rec = json_reader->get_root_interface<Input::Record>();
    	return root_rec;
    }
    void SetUp() override {
    	std::string f_name = string(UNIT_TESTS_SRC_DIR) + "/input/type_use_case_test.con";
    	ifstream in(f_name);
    	json_reader = new Input::JSONToStorage(in, get_input_type() );
    }
    void TearDown() override {
    	delete json_reader;
    }
private:
    Input::JSONToStorage *json_reader;
    Input::Record root_rec;

};

TEST_F(Application, init) {
    using namespace Input;

    FilePath::set_io_dirs("./root","/root","variant", "./output");

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

it::Record & Application::get_input_type() {
	static it::Record type = it::Record("Application", "Root record of the whole application.")
    // Array of equations with types given by method of class Equation
    .declare_key("equations", it::Array( Equation::get_input_type(), 1, 10 ), it::Default::obligatory(), "");
	type.close();
	return type;
}




it::Record & EquationA::get_input_rec() {
	static it::Record type = it::Record("EquationA", "For example explicit transport equation solver.")
    .derive_from( Equation::get_input_type() )
	.declare_key("mesh",it::FileName::input(),it::Default::obligatory(),"")
    .declare_key("parameter_a", it::Double(), "");
	type.close();
	return type;
}


const int EquationA::registrar =
		Input::register_class< EquationA, Input::Record >("EquationA");



it::Record & EquationB::get_input_rec() {
	static it::Record type = it::Record("EquationB", "For example implicit transport equation solver.")
    .derive_from( Equation::get_input_type() )
	.declare_key("mesh",it::FileName::input(),it::Default::obligatory(),"")
    .declare_key("parameter_b", it::Integer(), it::Default("111"), "")
    .declare_key("default_str", it::String(), it::Default("str value"), "" )
    .declare_key("substances", it::Array( it::String() ), it::Default::obligatory(), "" );
	type.close();
	return type;
}


const int EquationB::registrar =
		Input::register_class< EquationB, Input::Record >("EquationB");



it::AbstractRecord & Equation::get_input_type() {
	static it::AbstractRecord type = it::AbstractRecord("AbstractEquation","Abstract input Record type for any equation.");
	type.close();
	return type;
}




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



