/*
 * reader_cache_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <sstream>
#include <string>
#include <vector>
#include <mesh_constructor.hh>

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "io/reader_cache.hh"



Input::Record get_input_record(const std::string &input_str, Input::FileFormat format = Input::FileFormat::format_JSON) {

	istringstream is(input_str);
    Input::ReaderToStorage reader;
    IT::Record &in_rec = const_cast<IT::Record &>(Mesh::get_input_type());
    in_rec.finish();
    reader.read_stream(is, in_rec, format);

    return reader.get_root_interface<Input::Record>();
}


TEST(ReaderCache, get_bulk_element_) {
	Profiler::instance();
	unsigned int i, j;

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Input::Record i_rec = get_input_record("{mesh_file=\"fields/simplest_cube_data.msh\"}");
    FilePath file_name = i_rec.val<FilePath>("mesh_file");
    Mesh * mesh = new Mesh(i_rec);
    auto reader = ReaderCache::get_reader(file_name);
    reader->read_physical_names(mesh);
    reader->read_raw_mesh(mesh);
    mesh->setup_topology();
    ReaderCache::check_compatible_mesh(file_name, *mesh);

    // read  by components for MultiField
    BaseMeshReader::HeaderQuery header_params("vector_fixed", 0.0, OutputTime::DiscreteSpace::ELEM_DATA);
    for (i=0; i<3; ++i) {
    	ReaderCache::get_reader(file_name)->find_header(header_params);
        typename ElementDataCache<int>::ComponentDataPtr multifield_ =
        		ReaderCache::get_reader(file_name)->get_element_data<int>(9, 1, false, i);
    	std::vector<int> &vec = *( multifield_.get() );
    	EXPECT_EQ(9, vec.size());
    	for (j=0; j<mesh->n_elements(); j++) EXPECT_EQ( i+1, vec[j] );
    }

    // read  to one vector for Field
    {
    	BaseMeshReader::HeaderQuery header_params("vector_fixed", 1.0, OutputTime::DiscreteSpace::ELEM_DATA);
    	ReaderCache::get_reader(file_name)->find_header(header_params);
    	typename ElementDataCache<int>::ComponentDataPtr field_ =
    			ReaderCache::get_reader(file_name)->get_element_data<int>(9, 3, false, 0);
    	std::vector<int> &vec = *( field_.get() );
    	EXPECT_EQ(27, vec.size());
    	for (j=0; j<3*mesh->n_elements(); j++) EXPECT_EQ( 2+(j%3), vec[j] );
    }

    delete mesh;
}


TEST(ReaderCache, get_boundary_element_) {
	Profiler::instance();
	unsigned int i, j;

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Input::Record i_rec = get_input_record("{mesh_file=\"fields/simplest_cube_data.msh\"}");
    FilePath file_name = i_rec.val<FilePath>("mesh_file");
    Mesh * mesh = new Mesh(i_rec);
    auto reader = ReaderCache::get_reader(file_name);
    reader->read_physical_names(mesh);
    reader->read_raw_mesh(mesh);
    reader->check_compatible_mesh(*mesh);

    // read  by components for MultiField
    BaseMeshReader::HeaderQuery header_params("vector_fixed", 0.0, OutputTime::DiscreteSpace::ELEM_DATA);
    for (i=0; i<3; ++i) {
    	ReaderCache::get_reader(file_name)->find_header(header_params);
        typename ElementDataCache<int>::ComponentDataPtr multifield_ =
        		ReaderCache::get_reader(file_name)->get_element_data<int>(4, 1, true, i);
    	std::vector<int> &vec = *( multifield_.get() );
    	EXPECT_EQ(4, vec.size());
    	for (j=0; j<mesh->n_elements(true); j++) EXPECT_EQ( i+4, vec[j] );
    }

    // read  to one vector for Field
    {
    	BaseMeshReader::HeaderQuery header_params("vector_fixed", 1.0, OutputTime::DiscreteSpace::ELEM_DATA);
    	ReaderCache::get_reader(file_name)->find_header(header_params);
    	typename ElementDataCache<int>::ComponentDataPtr field_ =
    			ReaderCache::get_reader(file_name)->get_element_data<int>(4, 3, true, 0);
    	std::vector<int> &vec = *( field_.get() );
    	EXPECT_EQ(12, vec.size());
    	for (j=0; j<3*mesh->n_elements(true); j++) EXPECT_EQ( 5+(j%3), vec[j] );
    }

    delete mesh;
}


TEST(ReaderCache, find_header) {
    Profiler::instance();

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Input::Record i_rec = get_input_record("{mesh_file=\"fields/simplest_cube_data.msh\"}");
    FilePath file_name = i_rec.val<FilePath>("mesh_file");

    Mesh * mesh = new Mesh(i_rec);
    auto reader = ReaderCache::get_reader(file_name);
    reader->read_physical_names(mesh);
    reader->read_raw_mesh(mesh);
    mesh->setup_topology();
    ReaderCache::check_compatible_mesh(file_name, *mesh);
    delete mesh;

    BaseMeshReader::HeaderQuery header_params("vector_fixed", 0.0, OutputTime::DiscreteSpace::ELEM_DATA);

    {
    BaseMeshReader::MeshDataHeader &header = ReaderCache::get_reader(file_name)->find_header(header_params);
	EXPECT_EQ( header.field_name, "vector_fixed" );
	EXPECT_EQ( header.time, 0.0 );
    }

    {
    header_params.time = 0.1;
    BaseMeshReader::MeshDataHeader &header = ReaderCache::get_reader(file_name)->find_header(header_params);
	EXPECT_EQ( header.field_name, "vector_fixed" );
	EXPECT_EQ( header.time, 0.0 );
    }

    {
    header_params.time = 0.9;
    BaseMeshReader::MeshDataHeader &header = ReaderCache::get_reader(file_name)->find_header(header_params);
	EXPECT_EQ( header.field_name, "vector_fixed" );
	EXPECT_EQ( header.time, 0.0 );
    }

    {
    header_params.time = 1.0;
    BaseMeshReader::MeshDataHeader &header = ReaderCache::get_reader(file_name)->find_header(header_params);
    EXPECT_EQ( header.field_name, "vector_fixed" );
    EXPECT_EQ( header.time, 1.0 );
    }

    {
   	header_params.time = 1.1;
    BaseMeshReader::MeshDataHeader &header = ReaderCache::get_reader(file_name)->find_header(header_params);
    EXPECT_EQ( header.field_name, "vector_fixed" );
    EXPECT_EQ( header.time, 1.0 );
    }

    {
    header_params.time = 2.1;
    BaseMeshReader::MeshDataHeader &header = ReaderCache::get_reader(file_name)->find_header(header_params);
    EXPECT_EQ( header.field_name, "vector_fixed" );
    EXPECT_EQ( header.time, 1.0 );
    }

    {
    header_params.time = 200;
    BaseMeshReader::MeshDataHeader &header = ReaderCache::get_reader(file_name)->find_header(header_params);
    EXPECT_EQ( header.field_name, "vector_fixed" );
    EXPECT_EQ( header.time, 1.0 );
    }

}


TEST(ReaderCache, get_reader) {
	Profiler::instance();

	// has to introduce some flag for passing absolute path to 'test_units' in source tree
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	{
		FilePath mesh1("mesh/test_input.msh", FilePath::input_file);
		FilePath mesh2("mesh/simplest_cube.msh", FilePath::input_file);
		ReaderCache::get_reader(mesh1);
		ReaderCache::get_reader(mesh2);
	}

	{
	    Input::Record i_rec = get_input_record("{mesh_file=\"mesh/test_input.msh\"}");
	    Mesh * mesh = new Mesh(i_rec);
	    auto mesh_reader = ReaderCache::get_reader( i_rec.val<FilePath>("mesh_file") );
		mesh_reader->read_physical_names(mesh);
		mesh_reader->read_raw_mesh(mesh);

		EXPECT_EQ(118, mesh->n_nodes());

		delete mesh;
	}

	{
	    Input::Record i_rec = get_input_record("{mesh_file=\"mesh/simplest_cube.msh\"}");
	    Mesh * mesh = new Mesh(i_rec);
	    auto mesh_reader = ReaderCache::get_reader( i_rec.val<FilePath>("mesh_file") );
		mesh_reader->read_physical_names(mesh);
		mesh_reader->read_raw_mesh(mesh);

		EXPECT_EQ(8, mesh->n_nodes());

		delete mesh;
	}
}


TEST(ReaderCache, repeat_call) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Input::Record i_rec = get_input_record("{mesh_file=\"mesh/test_108_elem.msh\"}");
    for (unsigned int i=0; i<2; ++i) {
        auto mesh_reader = ReaderCache::get_reader( i_rec.val<FilePath>("mesh_file") );
        Mesh * mesh = new Mesh(i_rec);
        mesh_reader->read_physical_names(mesh);
        mesh_reader->read_raw_mesh(mesh);
        delete mesh;
    }
}
