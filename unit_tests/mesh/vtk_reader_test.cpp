/*
 * vtk_reader_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */

#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>

#include <string>
#include <iostream>
#include <pugixml.hpp>

#include "mesh/msh_vtkreader.hh"


class VtkMeshReaderTest : public testing::Test, public VtkMeshReader {
protected:
	virtual void SetUp() {
	    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
	    FilePath vtu_file("output/test_output_vtk_ascii_ref.vtu", FilePath::input_file);

	    parse_result_ = doc_.load_file( ((std::string)vtu_file).c_str() );
	    read_nodes_elms_count();
	}

	virtual void TearDown() {}

	std::string get_parse_result() {
		std::stringstream ss; ss << this->parse_result_.description();
		return ss.str();
	}
};


// simple test of short xml input, check functionality of pugixml
TEST(PugiXml, read_simple_xml) {
	std::stringstream ss;
	ss << "<mesh name='sphere'>\n";
	ss << "<bounds>0 0 1 1</bounds>\n";
	ss << "</mesh>\n";

	pugi::xml_document doc;
    pugi::xml_parse_result parse_result = doc.load(ss);

    std::stringstream result_desc; result_desc << parse_result.description();
	std::stringstream mesh_name; mesh_name << doc.child("mesh").attribute("name").value();
	std::stringstream mesh_bounds; mesh_bounds << doc.child("mesh").child("bounds").child_value();

	EXPECT_EQ( "No error", result_desc.str() );
	EXPECT_EQ( "sphere",   mesh_name.str() );
	EXPECT_EQ( "0 0 1 1",  mesh_bounds.str() );
}


// test of reading of VTU file
TEST_F(VtkMeshReaderTest, read_xml_file) {

    EXPECT_EQ( "No error", this->get_parse_result() );

    // get inner tag with VTU data
    pugi::xml_node piece_node = doc_.child("VTKFile").child("UnstructuredGrid").child("Piece");

    {
    	// test of attributes
		EXPECT_EQ( 8, this->n_nodes() );
		EXPECT_EQ( 6, this->n_elements() );
    }

    {
    	// test of Points data section
    	pugi::xml_node points_data_node = piece_node.child("Points").child("DataArray");
    	std::stringstream type; type << points_data_node.attribute("type").value();
    	std::stringstream value; value << points_data_node.child_value();
    	EXPECT_EQ( "Float64", type.str() );
    	EXPECT_EQ( "\n1",     value.str().substr(0,2) );
    }

    {
    	// test of connectivity data array
    	pugi::xml_node cells_data_node = piece_node.child("Cells").find_child_by_attribute("DataArray", "Name", "connectivity");
    	std::stringstream value; value << cells_data_node.child_value();
    	EXPECT_EQ( "\n2 6 0 1 2 6 1 7 2 6 7 5 2 6 5 4 2 6 4 3 2 6 3 0 \n", value.str() );
    }
}
