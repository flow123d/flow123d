/*
 * path_base_test.cpp
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

/**
 * TODO: test catching of errors in JSON file format.
 */

#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <fstream>


#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include "input/reader_internal_base.hh"
#include "input/path_base.hh"
#include "input/path_json.hh"
#include "input/path_yaml.hh"
#include "input/path_python.hh"

using namespace std;
using namespace Input;

namespace py = pybind11;


TEST(PathJSON, all) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    ifstream in_str((string(UNIT_TESTS_SRC_DIR) + "/input/reader_to_storage_test.con").c_str());
    PathJSON path(in_str);

    { ostringstream os;
    os << path;
    EXPECT_EQ("/",os.str());
    }

    path.down(6);
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6",os.str());
    }

    path.down("a");
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6/a",os.str());
    }
    EXPECT_EQ(1,path.find_ref_node()->get_int_value() );

    path.go_to_root();
    path.down(6);
    path.down("b");
    EXPECT_EQ("ctyri",path.find_ref_node()->get_string_value() );
}


TEST(PathJSON, errors) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";
	ifstream in_str((string(UNIT_TESTS_SRC_DIR) + "/input/reader_to_storage_test.con").c_str());

    PathJSON path(in_str);
    path.down(8);
    EXPECT_THROW_WHAT( { path.find_ref_node();}, PathBase::ExcRefOfWrongType,"has wrong type, should by string." );

    path.go_to_root();
    path.down(9); // "REF":"/5/10"
    EXPECT_THROW_WHAT( { path.find_ref_node();}, PathBase::ExcReferenceNotFound, "index out of size of Array" );

    path.go_to_root();
    path.down(10); // "REF":"/6/../.."
    EXPECT_THROW_WHAT( { path.find_ref_node();}, PathBase::ExcReferenceNotFound, "can not go up from root" );

    path.go_to_root();
    path.down(11); // "REF":"/key"
    EXPECT_THROW_WHAT( { path.find_ref_node();}, PathBase::ExcReferenceNotFound, "there should be Record" );

    path.go_to_root();
    path.down(12); // "REF":"/6/key"
    EXPECT_THROW_WHAT( { path.find_ref_node();}, PathBase::ExcReferenceNotFound, "key 'key' not found" );
}



TEST(PathYAML, all) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";

	ifstream in_str((string(UNIT_TESTS_SRC_DIR) + "/input/reader_to_storage_test.yaml").c_str());
	PathYAML path(in_str);

    { ostringstream os;
    os << path;
    EXPECT_EQ("/",os.str());
    }

    path.down(6);
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6",os.str());
    }

    path.down("a");
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6/a",os.str());
    }

    path.up();
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6",os.str());
    }

    path.down("b");
    { ostringstream os;
    os << path;
    EXPECT_EQ("/6/b",os.str());
    }

}


TEST(PathYAML, values) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";

	ifstream in_str((string(UNIT_TESTS_SRC_DIR) + "/input/reader_to_storage_test.yaml").c_str());
	PathYAML path(in_str);

	path.down(0); // bool value
	EXPECT_FALSE(path.get_bool_value());
	EXPECT_THROW( { path.get_int_value(); }, ReaderInternalBase::ExcInputError );
	EXPECT_THROW( { path.get_double_value(); }, ReaderInternalBase::ExcInputError );
	EXPECT_STREQ("false", path.get_string_value().c_str());
	path.up();

	path.down(1); // int value
	EXPECT_THROW( { path.get_bool_value(); }, ReaderInternalBase::ExcInputError );
	EXPECT_EQ(1, path.get_int_value());
	EXPECT_FLOAT_EQ(1.0, path.get_double_value());
	EXPECT_STREQ("1", path.get_string_value().c_str());
	path.up();

	path.down(3); // double value
	EXPECT_THROW( { path.get_bool_value(); }, ReaderInternalBase::ExcInputError );
	EXPECT_THROW( { path.get_int_value(); }, ReaderInternalBase::ExcInputError );
	EXPECT_FLOAT_EQ(3.3, path.get_double_value());
	EXPECT_STREQ("3.3", path.get_string_value().c_str());
	path.up();

	path.down(4); // string value
	EXPECT_THROW( { path.get_bool_value(); }, ReaderInternalBase::ExcInputError );
	EXPECT_THROW( { path.get_int_value(); }, ReaderInternalBase::ExcInputError );
	EXPECT_THROW( { path.get_double_value(); }, ReaderInternalBase::ExcInputError );
	EXPECT_STREQ("ctyri", path.get_string_value().c_str());
	path.up();

	path.down(9); // int64 value
	EXPECT_EQ(5000000000000, path.get_int_value());
	path.up();

	path.down(6); // record
	std::set<std::string> set;
	path.get_record_key_set(set);
	EXPECT_EQ(2, set.size());
	EXPECT_TRUE( set.find("a")!=set.end() );
	EXPECT_TRUE( set.find("b")!=set.end() );
	EXPECT_FALSE( set.find("c")!=set.end() );

	path.down("b"); // reference
	EXPECT_STREQ("ctyri", path.get_string_value().c_str());
}



py::dict create_input_dict() {
    py::dict input_dict;

    input_dict["flow_version"] = "4.0.0a01";
    input_dict["description"] = "Simple test - Steady flow with sources";
    input_dict["time_limit"] = 20.5;
    input_dict["pause_after_run"] = false;

    py::dict mesh;
    mesh["mesh"] = "../00_mesh/square_1x1_shift.msh";
    mesh["optimize"] = true;
    input_dict["mesh"] = mesh;

    py::dict flow_equation;

    py::dict solver;
    solver["type"] = "Petsc";
    solver["r_tol"] = 1.0e-5;
    solver["a_tol"] = 1.0e-5;
    flow_equation["solver"] = solver;

    py::list input_fields;
    py::dict plane_reg;
    plane_reg["region"] = "plane";
    plane_reg["anisotropy"] = 1;
    plane_reg["water_source_density"] = 0.8;
    input_fields.append(plane_reg);
    py::dict plane_bdr_reg;
    plane_bdr_reg["region"] = ".plane_boundary";
    plane_bdr_reg["bc_type"] = "dirichlet";
    plane_bdr_reg["bc_pressure"] = 0;
    input_fields.append(plane_bdr_reg);
    flow_equation["input_fields"] = input_fields;

    py::list output_fields;
    output_fields.append("pressure_p0");
    output_fields.append("velocity");
    flow_equation["output_fields"] = output_fields;

    input_dict["flow_equation"] = flow_equation;

    return input_dict;
}

// Register Python module
PYBIND11_MODULE(input_module, m) {
    m.def("create_dict", &create_input_dict, "Creates and fills Python dict");
}

TEST(PathPython, all) {
::testing::FLAGS_gtest_death_test_style = "threadsafe";

    py::scoped_interpreter guard{};

    py::dict input_dict = create_input_dict();
    std::cout << py::str(input_dict).cast<std::string>() << std::endl;
    PathPython path(input_dict);

    { ostringstream os;
    os << path;
    EXPECT_EQ("/",os.str());
    }

    path.down("time_limit");
    { ostringstream os;
    os << path;
    EXPECT_EQ("/time_limit",os.str());
    }

    path.up();
    path.down("flow_equation");
    path.down("solver");
    path.down("a_tol");
    { ostringstream os;
    os << path;
    EXPECT_EQ("/flow_equation/solver/a_tol",os.str());
    }

    path.up();
    path.up();
    path.down("input_fields");
    path.down(0);
    path.down("anisotropy");
    { ostringstream os;
    os << path;
    EXPECT_EQ("/flow_equation/input_fields/0/anisotropy",os.str());
    }

    path.up();
    path.up();
    path.down(1);
    path.down("bc_type");
    { ostringstream os;
    os << path;
    EXPECT_EQ("/flow_equation/input_fields/1/bc_type",os.str());
    }

}


