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


class VtkMeshReaderTest : public VtkMeshReader {
public:
	VtkMeshReaderTest(const FilePath &file_name)
	: VtkMeshReader(file_name) {}

	VtkMeshReaderTest(std::istream &in)
	: VtkMeshReader(in) {}

	std::string get_parse_result() {
		std::stringstream ss; ss << this->parse_result_.description();
		return ss.str();
	}

	pugi::xml_document & get_xml_doc() {
		return this->doc_;
	}
};


// simple test of short xml input, check functionality of pugixml
TEST(VTKReader, read_simple_xml) {
	std::stringstream ss;
	ss << "<mesh name='sphere'>\n";
	ss << "<bounds>0 0 1 1</bounds>\n";
	ss << "</mesh>\n";

	VtkMeshReaderTest vtk_reader( ss );

	std::stringstream mesh_name; mesh_name << vtk_reader.get_xml_doc().child("mesh").attribute("name").value();
	std::stringstream mesh_bounds; mesh_bounds << vtk_reader.get_xml_doc().child("mesh").child("bounds").child_value();

	EXPECT_EQ( "No error", vtk_reader.get_parse_result() );
	EXPECT_EQ( "sphere",   mesh_name.str() );
	EXPECT_EQ( "0 0 1 1",  mesh_bounds.str() );
}


// test of reading of VTU file
TEST(VTKReader, read_xml_file) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    FilePath vtu_file("output/test_output_vtk_ascii_ref.vtu", FilePath::input_file);

    VtkMeshReaderTest vtk_reader( vtu_file );
    EXPECT_EQ( "No error", vtk_reader.get_parse_result() );

    // get inner tag with VTU data
    pugi::xml_node piece_node = vtk_reader.get_xml_doc().child("VTKFile").child("UnstructuredGrid").child("Piece");

    {
    	// test of attributes
		std::stringstream points; points << piece_node.attribute("NumberOfPoints").value();
		std::stringstream cells; cells << piece_node.attribute("NumberOfCells").value();
		EXPECT_EQ( "8",  points.str() );
		EXPECT_EQ( "6",  cells.str() );
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
