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

    {
    	// test of attributes
		EXPECT_EQ( 8, this->n_nodes() );
		EXPECT_EQ( 6, this->n_elements() );
    }

    {
    	// test of Points data section
    	auto data_attr = this->get_data_array_attr(DataSections::points);
    	EXPECT_EQ( DataType::float64, data_attr.type_ );
    	EXPECT_EQ( DataFormat::ascii, data_attr.format_ );
    	EXPECT_EQ( 3, data_attr.n_components_ );
    }

    {
    	// test of connectivity data array
    	auto data_attr = this->get_data_array_attr(DataSections::cells, "connectivity");
    	EXPECT_EQ( DataType::uint32, data_attr.type_ );
    	EXPECT_EQ( DataFormat::ascii, data_attr.format_ );
    	EXPECT_EQ( 1, data_attr.n_components_ );
    }
}
