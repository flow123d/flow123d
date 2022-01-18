/*
 * vtk_reader_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>

#include <string>
#include <iostream>
#include <pugixml.hpp>

#include "io/msh_vtkreader.hh"
#include "io/msh_gmshreader.h"
#include "io/reader_cache.hh"
#include "system/tokenizer.hh"


class VtkMeshReaderTest : public VtkMeshReader {
public:
	/// Helper factory, creates shared pointer to instance of VtkMeshReaderTest
	static std::shared_ptr<VtkMeshReaderTest> test_factory(const std::string &input_str,
			Input::FileFormat format = Input::FileFormat::format_JSON) {
		istringstream is(input_str);
	    Input::ReaderToStorage reader;
	    IT::Record &in_rec = const_cast<IT::Record &>(Mesh::get_input_type());
	    in_rec.finish();
	    reader.read_stream(is, in_rec, format);

		return std::make_shared<VtkMeshReaderTest>(reader.get_root_interface<Input::Record>().val<FilePath>("mesh_file"));
	}

	VtkMeshReaderTest(const FilePath &file_name)
	: VtkMeshReader(file_name) {
	}

	DataFormat get_data_format() {
		return data_format_;
	}

	void read_nodes(Mesh* mesh)
	{
		HeaderQuery header_params("Points", 0.0, OutputTime::DiscreteSpace::MESH_DEFINITION);
		auto actual_header = this->find_header(header_params);

		bulk_elements_id_.resize(actual_header.n_entities);
		for (unsigned int i=0; i<bulk_elements_id_.size(); ++i) bulk_elements_id_[i]=i;

		// set new cache
	    ElementDataCacheBase *current_cache = new ElementDataCache<double>(actual_header.field_name, actual_header.time,
	    		actual_header.n_components*actual_header.n_entities);

		switch (data_format_) {
			case DataFormat::ascii: {
				parse_ascii_data( *current_cache, actual_header.n_components, actual_header.n_entities, actual_header.position );
				break;
			}
			case DataFormat::binary_uncompressed: {
				ASSERT_PERMANENT_PTR(data_stream_).error();
				parse_binary_data( *current_cache, actual_header.n_components, actual_header.n_entities, actual_header.position);
				break;
			}
			case DataFormat::binary_zlib: {
				ASSERT_PERMANENT_PTR(data_stream_).error();
				parse_compressed_data(* current_cache, actual_header.n_components, actual_header.n_entities, actual_header.position);
				break;
			}
			default: {
				ASSERT_PERMANENT(false).error(); // should not happen
				break;
			}
		}

		mesh->init_node_vector(actual_header.n_entities);
		std::vector<double> &vect = *( static_cast< ElementDataCache<double> *>(current_cache)->get_data().get() );
		arma::vec3 point;
		for (unsigned int i=0, ivec=0; i<actual_header.n_entities; ++i) {
	        point[0]=vect[ivec]; ++ivec;
	        point[1]=vect[ivec]; ++ivec;
	        point[2]=vect[ivec]; ++ivec;
	        mesh->add_node(i, point);
		}

		bulk_elements_id_.clear();
	    delete current_cache;
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
    std::string mesh_in_string = "{mesh_file=\"output/test_output_vtk_ascii_ref.vtu\"}";
    auto reader = VtkMeshReaderTest::test_factory(mesh_in_string);

	Mesh * mesh = mesh_constructor(mesh_in_string);

    {
    	// test of Points data section
    	reader->read_nodes(mesh);
        EXPECT_EQ(8, mesh->n_nodes());
    }

    {
    	// test of connectivity data array
    	BaseMeshReader::HeaderQuery header_params("connectivity", 0.0, OutputTime::DiscreteSpace::MESH_DEFINITION);
    	auto data_attr = reader->find_header(header_params);
    	EXPECT_EQ( DataType::uint32, data_attr.type );
    	EXPECT_EQ( VtkMeshReader::DataFormat::ascii, reader->get_data_format() );
    	EXPECT_EQ( 1, data_attr.n_components );
    	EXPECT_EQ( 6, data_attr.n_entities );
    	EXPECT_EQ( data_attr.discretization, OutputTime::DiscreteSpace::MESH_DEFINITION );
    }

    {
    	// test of connectivity data array
    	BaseMeshReader::HeaderQuery header_params("flow_data", 0.0, OutputTime::DiscreteSpace::MESH_DEFINITION, 6);
    	auto data_attr = reader->find_header(header_params);
    	EXPECT_EQ( DataType::float64, data_attr.type );
    	EXPECT_EQ( VtkMeshReader::DataFormat::ascii, reader->get_data_format() );
    	EXPECT_EQ( 1, data_attr.n_components );
    	EXPECT_EQ( 6, data_attr.n_entities );
    	EXPECT_EQ( data_attr.discretization, OutputTime::DiscreteSpace::NATIVE_DATA );
    	EXPECT_EQ( 6, data_attr.dof_handler_hash );
    }

    //delete mesh;
}


TEST(VtkReaderTest, read_binary_vtu) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    FilePath mesh_file("output/test_output_vtk_binary_ref.vtu", FilePath::input_file);

    {
    	std::string mesh_in_string = "{ mesh_file=\"fields/simplest_cube_3d.msh\", optimize_mesh=false }";
    	auto gmsh_reader = reader_constructor( mesh_in_string );
    	Mesh * source_mesh = mesh_constructor( mesh_in_string );
    	gmsh_reader->read_physical_names(source_mesh);
    	gmsh_reader->read_raw_mesh(source_mesh);
    	vector<LongIdx> bulk_elms_id, boundary_elms_id;
    	source_mesh->setup_topology();
    	source_mesh->check_and_finish();
    	ReaderCache::get_mesh(mesh_file)->check_compatible_mesh(const_cast<Mesh &>(*source_mesh));
    	ReaderCache::get_element_ids(mesh_file, *source_mesh);
        delete source_mesh;
    }

    {
    	// test of Points data section
        EXPECT_EQ(8, ReaderCache::get_mesh(mesh_file)->n_nodes());
    }

    {
    	// test of connectivity data array
    	BaseMeshReader::HeaderQuery header_params("connectivity", 0.0, OutputTime::DiscreteSpace::MESH_DEFINITION);
    	auto data_attr = ReaderCache::get_reader(mesh_file)->find_header(header_params);
    	EXPECT_EQ( DataType::uint32, data_attr.type );
    	//EXPECT_EQ( VtkMeshReader::DataFormat::binary_uncompressed, ReaderCache::get_reader(mesh_file)->get_data_format() );
    	EXPECT_EQ( 1, data_attr.n_components );
    	EXPECT_EQ( 6, data_attr.n_entities );
    }

    std::vector<int> el_ids;
    unsigned int i, j;
    for (i=0; i<6; ++i) el_ids.push_back(i);


    bool boundary_domain = false; // bulk data
    BaseMeshReader::HeaderQuery vector_header_params("vector_field", 0.0, OutputTime::DiscreteSpace::ELEM_DATA);
    auto header = ReaderCache::get_reader(mesh_file)->find_header(vector_header_params);
	{
        // test exception on wrong number of components
        uint wrong_number_of_components = 7;
        EXPECT_EQ( 3, header.n_components );
		EXPECT_THROW_WHAT( { ReaderCache::get_reader(mesh_file)->template get_element_data<double>(
            header, 6, wrong_number_of_components, boundary_domain); }, BaseMeshReader::ExcWrongComponentsCount,
            "Wrong number of components of field 'vector_field' at time 0 in the input file: ");
	}
	{
        // test reading vector field
        typename ElementDataCache<double>::CacheData multifield_data =
                ReaderCache::get_reader(mesh_file)->template get_element_data<double>(header, 6, 3, boundary_domain);
        std::vector<double> &vec = *( multifield_data.get() );
        EXPECT_EQ(18, vec.size());
        for (j=0; j<vec.size(); j++) {
            // DebugOut() << i << ": " << vec[j] << "\n";
            EXPECT_DOUBLE_EQ( 0.5*(j%3+1), vec[j] );
        }
	}

    // read data to one vector for Field
    BaseMeshReader::HeaderQuery tensor_header_params("tensor_field", 1.0, OutputTime::DiscreteSpace::ELEM_DATA);
    {
    	std::vector<double> ref_data = { 1, 4, 7, 2, 5, 8, 3, 6, 9 };
    	auto header = ReaderCache::get_reader(mesh_file)->find_header(tensor_header_params);
    	typename ElementDataCache<double>::CacheData field_data =
    	        ReaderCache::get_reader(mesh_file)->template get_element_data<double>(header, 6, 9, boundary_domain);
    	std::vector<double> &vec = *( field_data.get() );
    	EXPECT_EQ(54, vec.size());
    	for (j=0; j<vec.size(); j++) {
    		EXPECT_DOUBLE_EQ( ref_data[j%9], vec[j] );
    	}
    }

   // delete mesh;
}


TEST(VtkReaderTest, read_compressed_vtu) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    std::string mesh_in_string = "{mesh_file=\"output/test_output_vtk_zlib_ref.vtu\"}";
    auto reader = VtkMeshReaderTest::test_factory(mesh_in_string);

    Mesh * mesh = mesh_constructor(mesh_in_string);

    {
    	// test of Points data section
    	reader->read_nodes(mesh);
        EXPECT_EQ(8, mesh->n_nodes());
    }

    {
    	// test of connectivity data array
    	BaseMeshReader::HeaderQuery header_params("connectivity", 0.0, OutputTime::DiscreteSpace::MESH_DEFINITION);
    	auto data_attr = reader->find_header(header_params);
    	EXPECT_EQ( DataType::uint32, data_attr.type );
    	EXPECT_EQ( VtkMeshReader::DataFormat::binary_zlib, reader->get_data_format() );
    	EXPECT_EQ( 1, data_attr.n_components );
    	EXPECT_EQ( 6, data_attr.n_entities );
    }

    //delete mesh;
}


TEST(VtkReaderTest, read_mesh) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    std::string mesh_in_string = "{ mesh_file=\"output/test_output_vtk_ascii_ref.vtu\", optimize_mesh=false }";
    auto reader = reader_constructor(mesh_in_string);
    Mesh * mesh = mesh_constructor(mesh_in_string);
    reader->read_raw_mesh(mesh);
    mesh->setup_topology();
    mesh->check_and_finish();

    EXPECT_EQ(8, mesh->n_nodes());
    EXPECT_EQ(6, mesh->n_elements());
}
