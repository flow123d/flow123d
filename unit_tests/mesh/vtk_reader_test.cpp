/*
 * vtk_reader_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */

#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <mesh_constructor.hh>

#include <string>
#include <iostream>
#include <pugixml.hpp>

#include "mesh/msh_vtkreader.hh"


class VtkMeshReaderTest : public testing::Test, public VtkMeshReader {
protected:
	virtual void SetUp() {
	    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
	}

	virtual void TearDown() {}

	void read_input_file(const std::string &file_name) {
	    FilePath vtu_file(file_name, FilePath::input_file);
	    parse_result_ = doc_.load_file( ((std::string)vtu_file).c_str() );
	    this->read_base_vtk_attributes();
		// data of appended tag
		if (this->header_type_==DataType::undefined) { // no AppendedData tag
			appended_pos_ = 0;
		} else {
			this->set_appended_stream(vtu_file);
		}
	}

	std::string get_parse_result() {
		std::stringstream ss; ss << this->parse_result_.description();
		return ss.str();
	}

	void read_nodes(Mesh* mesh)
	{
		std::vector<double> nodes_data;
		auto data_attr = this->get_data_array_attr(DataSections::points);
		switch (data_format_) {
			case DataFormat::ascii: {
				nodes_data = parse_ascii_data<double>( this->n_nodes_*data_attr.n_components_ , data_attr.tag_value_ );
				break;
			}
			case DataFormat::binary_uncompressed: {
				ASSERT_PTR(appended_stream_).error();
				nodes_data = parse_binary_data<double>( appended_pos_+data_attr.offset_, data_attr.type_);
				break;
			}
			case DataFormat::binary_zlib: {
				ASSERT_PTR(appended_stream_).error();
				nodes_data = parse_compressed_data<double>( appended_pos_+data_attr.offset_, data_attr.type_);
				break;
			}
			default: {
				ASSERT(false).error(); // should not happen
				break;
			}
		}

		mesh->node_vector.reserve(this->n_nodes_);
		for (unsigned int i=0, ivec=0; i<this->n_nodes_; ++i) {
	        NodeFullIter node = mesh->node_vector.add_item(i);
	        node->point()(0)=nodes_data[ivec]; ++ivec;
	        node->point()(1)=nodes_data[ivec]; ++ivec;
	        node->point()(2)=nodes_data[ivec]; ++ivec;
		}
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
TEST_F(VtkMeshReaderTest, read_ascii_vtu) {
	read_input_file("output/test_output_vtk_ascii_ref.vtu");

    EXPECT_EQ( "No error", this->get_parse_result() );

    {
    	// test of attributes
		EXPECT_EQ( 8, this->n_nodes() );
		EXPECT_EQ( 6, this->n_elements() );
    }

    Mesh *mesh = mesh_constructor();

    {
    	// test of Points data section
        this->read_nodes(mesh);
        EXPECT_EQ(8, mesh->n_nodes());
    }

    {
    	// test of connectivity data array
    	auto data_attr = this->get_data_array_attr(DataSections::cells, "connectivity");
    	EXPECT_EQ( DataType::uint32, data_attr.type_ );
    	EXPECT_EQ( DataFormat::ascii, data_format_ );
    	EXPECT_EQ( 1, data_attr.n_components_ );
    }
}


TEST_F(VtkMeshReaderTest, read_binary_vtu) {
	read_input_file("output/test_output_vtk_binary_ref.vtu");

    {
    	// test of attributes
		EXPECT_EQ( 8, this->n_nodes() );
		EXPECT_EQ( 6, this->n_elements() );
    }

    Mesh *mesh = mesh_constructor();

    {
    	// test of Points data section
        this->read_nodes(mesh);
        EXPECT_EQ(8, mesh->n_nodes());
    }

    {
    	// test of connectivity data array
    	auto data_attr = this->get_data_array_attr(DataSections::cells, "connectivity");
    	EXPECT_EQ( DataType::uint32, data_attr.type_ );
    	EXPECT_EQ( DataFormat::binary_uncompressed, data_format_ );
    	EXPECT_EQ( 1, data_attr.n_components_ );
    }
}


TEST_F(VtkMeshReaderTest, read_compressed_vtu) {
	read_input_file("output/test_output_vtk_zlib_ref.vtu");

    {
    	// test of attributes
		EXPECT_EQ( 8, this->n_nodes() );
		EXPECT_EQ( 6, this->n_elements() );
    }

    Mesh *mesh = mesh_constructor();

    {
    	// test of Points data section
        this->read_nodes(mesh);
        EXPECT_EQ(8, mesh->n_nodes());
    }

    {
    	// test of connectivity data array
    	auto data_attr = this->get_data_array_attr(DataSections::cells, "connectivity");
    	EXPECT_EQ( DataType::uint32, data_attr.type_ );
    	EXPECT_EQ( DataFormat::binary_zlib, data_format_ );
    	EXPECT_EQ( 1, data_attr.n_components_ );
    }
}
