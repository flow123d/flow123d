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
#include "system/tokenizer.hh"


class VtkMeshReaderTest : public VtkMeshReader {
public:
	VtkMeshReaderTest(const FilePath &file_name)
	: VtkMeshReader(file_name) {
	}

	DataFormat get_data_format() {
		return data_format_;
	}

	void read_nodes(Mesh* mesh)
	{
		auto actual_header = this->find_header(0.0, "Points");

		// set new cache
	    delete current_cache_;
	    typename ElementDataCache<double>::CacheData data_cache
			= ElementDataCache<double>::create_data_cache(1, actual_header.n_components*actual_header.n_entities);
	    current_cache_ = new ElementDataCache<double>(actual_header.time, actual_header.field_name, data_cache);

		switch (data_format_) {
			case DataFormat::ascii: {
				parse_ascii_data<double>( 1, actual_header.n_components, actual_header.n_entities, actual_header.position );
				break;
			}
			case DataFormat::binary_uncompressed: {
				ASSERT_PTR(data_stream_).error();
				parse_binary_data<double>( 1, actual_header.n_components, actual_header.n_entities, actual_header.position,
						actual_header.type );
				break;
			}
			case DataFormat::binary_zlib: {
				ASSERT_PTR(data_stream_).error();
				parse_compressed_data<double>( 1, actual_header.n_components, actual_header.n_entities, actual_header.position,
						actual_header.type);
				break;
			}
			default: {
				ASSERT(false).error(); // should not happen
				break;
			}
		}

		mesh->node_vector.reserve(actual_header.n_entities);
		std::vector<double> &vect = *( static_cast< ElementDataCache<double> *>(current_cache_)->get_component_data(0).get() );
		for (unsigned int i=0, ivec=0; i<actual_header.n_entities; ++i) {
	        NodeFullIter node = mesh->node_vector.add_item(i);
	        node->point()(0)=vect[ivec]; ++ivec;
	        node->point()(1)=vect[ivec]; ++ivec;
	        node->point()(2)=vect[ivec]; ++ivec;
		}
	}

	MeshDataHeader & find_header(double time, std::string field_name) {
		return VtkMeshReader::find_header(time, field_name);
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
TEST(VtkReaderTest, read_ascii_vtu) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    FilePath vtu_file("output/test_output_vtk_ascii_ref.vtu", FilePath::input_file);
    VtkMeshReaderTest reader(vtu_file);

    Mesh *mesh = mesh_constructor();

    {
    	// test of Points data section
    	reader.read_nodes(mesh);
        EXPECT_EQ(8, mesh->n_nodes());
    }

    {
    	// test of connectivity data array
    	auto data_attr = reader.find_header(0.0, "connectivity");
    	EXPECT_EQ( DataType::uint32, data_attr.type );
    	EXPECT_EQ( VtkMeshReader::DataFormat::ascii, reader.get_data_format() );
    	EXPECT_EQ( 1, data_attr.n_components );
    	EXPECT_EQ( 6, data_attr.n_entities );
    }
}


TEST(VtkReaderTest, read_binary_vtu) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    FilePath vtu_file("output/test_output_vtk_binary_ref.vtu", FilePath::input_file);
    VtkMeshReaderTest reader(vtu_file);

    Mesh *mesh = mesh_constructor();

    {
    	// test of Points data section
    	reader.read_nodes(mesh);
        EXPECT_EQ(8, mesh->n_nodes());
    }

    {
    	// test of connectivity data array
    	auto data_attr = reader.find_header(0.0, "connectivity");
    	EXPECT_EQ( DataType::uint32, data_attr.type );
    	EXPECT_EQ( VtkMeshReader::DataFormat::binary_uncompressed, reader.get_data_format() );
    	EXPECT_EQ( 1, data_attr.n_components );
    	EXPECT_EQ( 6, data_attr.n_entities );
    }

    std::vector<int> el_ids;
    unsigned int i, j;
    for (i=0; i<6; ++i) el_ids.push_back(i);


    // read data by components for MultiField
    bool actual_data = false;
    for (i=0; i<3; ++i) {
        typename ElementDataCache<double>::ComponentDataPtr multifield_data =
        		reader.get_element_data<double>("vector_field", 0.0, 6, 1, actual_data, el_ids, i);
    	std::vector<double> &vec = *( multifield_data.get() );
    	EXPECT_EQ(6, vec.size());
    	for (j=0; j<vec.size(); j++) {
    		EXPECT_DOUBLE_EQ( 0.5*(i+1), vec[j] );
    	}
    }

    // read data to one vector for Field
    actual_data=false;
    {
    	typename ElementDataCache<double>::ComponentDataPtr field_data =
    			reader.get_element_data<double>("vector_field", 1.0, 6, 3, actual_data, el_ids, 0);
    	std::vector<double> &vec = *( field_data.get() );
    	EXPECT_EQ(18, vec.size());
    	for (j=0; j<vec.size(); j++) {
    		EXPECT_DOUBLE_EQ( 0.5*(j%3+1), vec[j] );
    	}
    }

}


TEST(VtkReaderTest, read_compressed_vtu) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    FilePath vtu_file("output/test_output_vtk_zlib_ref.vtu", FilePath::input_file);
    VtkMeshReaderTest reader(vtu_file);

    Mesh *mesh = mesh_constructor();

    {
    	// test of Points data section
    	reader.read_nodes(mesh);
        EXPECT_EQ(8, mesh->n_nodes());
    }

    {
    	// test of connectivity data array
    	auto data_attr = reader.find_header(0.0, "connectivity");
    	EXPECT_EQ( DataType::uint32, data_attr.type );
    	EXPECT_EQ( VtkMeshReader::DataFormat::binary_zlib, reader.get_data_format() );
    	EXPECT_EQ( 1, data_attr.n_components );
    	EXPECT_EQ( 6, data_attr.n_entities );
    }
}
