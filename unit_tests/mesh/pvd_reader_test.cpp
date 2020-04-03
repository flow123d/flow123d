/*
 * mpi_reader_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest_mpi.hh>
#include <sstream>
#include <string>
#include <mesh_constructor.hh>

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "io/msh_pvdreader.hh"


class PvdMeshReaderTest : public PvdMeshReader {
public:
	PvdMeshReaderTest(const FilePath &file_name)
	: PvdMeshReader(file_name) {}

	void test_find_header(double time, double expected_time, std::string field_name, OutputTime::DiscreteSpace disc) {
		HeaderQuery header_params(field_name, time, disc);
		MeshDataHeader &header = this->find_header(header_params);
		EXPECT_DOUBLE_EQ(expected_time, header.time);
		EXPECT_EQ(field_name, header.field_name);
	}
};


TEST(PVDReader, read_pvd) {
    Profiler::instance();

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    FilePath mesh_file("mesh/pvd-test.pvd", FilePath::input_file);

    PvdMeshReaderTest reader(mesh_file);
    reader.test_find_header(0.0,  0.0, "offsets", OutputTime::DiscreteSpace::MESH_DEFINITION);
    reader.test_find_header(0.05, 0.0, "offsets", OutputTime::DiscreteSpace::MESH_DEFINITION);
    reader.test_find_header(0.11, 0.1, "offsets", OutputTime::DiscreteSpace::MESH_DEFINITION);
    reader.test_find_header(0.21, 0.2, "offsets", OutputTime::DiscreteSpace::MESH_DEFINITION);
}


TEST(PVDReader, get_element_data) {
    Profiler::instance();

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    FilePath mesh_file("mesh/pvd-test.pvd", FilePath::input_file);

    PvdMeshReaderTest reader(mesh_file);
    {
        std::string mesh_in_string = "{mesh_file=\"fields/simplest_cube_3d.msh\"}";
        auto gmsh_reader = reader_constructor( mesh_in_string );
        Mesh * source_mesh = mesh_constructor( mesh_in_string );
        gmsh_reader->read_physical_names(source_mesh);
        gmsh_reader->read_raw_mesh(source_mesh);
        source_mesh->setup_topology();
        source_mesh->check_and_finish();
        reader.check_compatible_mesh(*source_mesh);
        delete source_mesh;
    }

    for (unsigned int i=0; i<3; ++i) {
    	BaseMeshReader::HeaderQuery header_params("scalar_field", i*0.1, OutputTime::DiscreteSpace::ELEM_DATA);
    	reader.find_header(header_params);
        typename ElementDataCache<double>::ComponentDataPtr scalar_data = reader.get_element_data<double>(6, 1, false, 0);
        std::vector<double> &vec = *( scalar_data.get() );
        EXPECT_EQ(6, vec.size());
        for (unsigned int j=0; j<vec.size(); j++) {
        	EXPECT_DOUBLE_EQ( (i+j+1)*0.5, vec[j] );
        }
    }

    for (unsigned int i=0; i<3; ++i) {
    	BaseMeshReader::HeaderQuery header_params("vector_field", i*0.1, OutputTime::DiscreteSpace::ELEM_DATA);
    	reader.find_header(header_params);
        typename ElementDataCache<double>::ComponentDataPtr vector_data = reader.get_element_data<double>(6, 3, false, 0);
        std::vector<double> &vec = *( vector_data.get() );
        EXPECT_EQ(18, vec.size());
        for (unsigned int j=0; j<vec.size(); j++) {
        	EXPECT_DOUBLE_EQ( (i*2+1)*(j%3+1)*0.5, vec[j] );
        }
    }
}
