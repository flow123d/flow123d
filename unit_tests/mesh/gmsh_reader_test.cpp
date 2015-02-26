/*
 * gmsh_reader_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */


#include <flow_gtest.hh>
#include <sstream>
#include <string>
#include <vector>

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "mesh/reader_instances.hh"



TEST(GMSHReader, read_mesh_from_stream) {
    Profiler::initialize();
    
    string fname = string(UNIT_TESTS_SRC_DIR) + "/mesh/simplest_cube.msh";
    ifstream ifs( fname.c_str() );
    if (! ifs) cout << "Can not open file!" << endl;
    stringstream ss;
    ss << ifs.rdbuf();

    Mesh mesh;
    GmshMeshReader reader(ss);

    reader.read_mesh(&mesh);

    EXPECT_EQ(9, mesh.n_elements());

}

TEST(GMSHReader, read_mesh_from_file) {
    Profiler::initialize();

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    FilePath mesh_file("mesh/test_input.msh", FilePath::input_file);


    Mesh mesh;
    GmshMeshReader reader(mesh_file);

    reader.read_mesh(&mesh);

    EXPECT_EQ(216, mesh.n_elements());

}

TEST(GMSHReader, find_header) {
    Profiler::initialize();
    
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    FilePath mesh_file("fields/simplest_cube_data.msh", FilePath::input_file);

    GmshMeshReader reader(mesh_file);
    unsigned int n_elements=13;
    unsigned int n_comp=3;
    std::vector<double> data(n_comp*n_elements);
    std::vector<int> element_id_map(n_elements);
    for(unsigned int i=0;i<n_elements;i++)
    	element_id_map[i]=i+1;
    GMSH_DataHeader header;
    header.actual=false;
    header.field_name="vector_fixed";
    header.n_components=n_comp;
    header.n_entities=n_elements;

    header.time=0.0;
    reader.read_element_data(header, &(data[0]), element_id_map);
    EXPECT_EQ(1.0, data[0]);
    EXPECT_EQ(2.0, data[1]);
    EXPECT_EQ(3.0, data[2]);
    EXPECT_EQ(3.0, data[3*8+2]);
    EXPECT_EQ(4.0, data[3*9]);

    header.time=0.1;
    reader.read_element_data(header, &(data[0]), element_id_map);
    EXPECT_EQ(1.0, data[0]);
    EXPECT_EQ(2.0, data[1]);
    EXPECT_EQ(3.0, data[2]);
    EXPECT_EQ(3.0, data[3*8+2]);
    EXPECT_EQ(4.0, data[3*9]);

    header.time=0.9;
    reader.read_element_data(header, &(data[0]), element_id_map);
    EXPECT_EQ(1.0, data[0]);
    EXPECT_EQ(2.0, data[1]);
    EXPECT_EQ(3.0, data[2]);
    EXPECT_EQ(3.0, data[3*8+2]);
    EXPECT_EQ(4.0, data[3*9]);

    header.time=1.0;
    reader.read_element_data(header, &(data[0]), element_id_map);
    EXPECT_EQ(2.0, data[0]);
    EXPECT_EQ(3.0, data[1]);
    EXPECT_EQ(4.0, data[2]);
    EXPECT_EQ(4.0, data[3*8+2]);
    EXPECT_EQ(5.0, data[3*9]);

    header.time=1.1;
    reader.read_element_data(header, &(data[0]), element_id_map);
    EXPECT_EQ(2.0, data[0]);
    EXPECT_EQ(3.0, data[1]);
    EXPECT_EQ(4.0, data[2]);
    EXPECT_EQ(4.0, data[3*8+2]);
    EXPECT_EQ(5.0, data[3*9]);

    header.time=2.1;
    reader.read_element_data(header, &(data[0]), element_id_map);
    EXPECT_EQ(2.0, data[0]);
    EXPECT_EQ(3.0, data[1]);
    EXPECT_EQ(4.0, data[2]);
    EXPECT_EQ(4.0, data[3*8+2]);
    EXPECT_EQ(5.0, data[3*9]);

    header.time=200;
    reader.read_element_data(header, &(data[0]), element_id_map);
    EXPECT_EQ(2.0, data[0]);
    EXPECT_EQ(3.0, data[1]);
    EXPECT_EQ(4.0, data[2]);
    EXPECT_EQ(4.0, data[3*8+2]);
    EXPECT_EQ(5.0, data[3*9]);

}

TEST(ReaderInstances, get_reader) {
	Profiler::initialize();

	// has to introduce some flag for passing absolute path to 'test_units' in source tree
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	{
		FilePath mesh1("mesh/test_input.msh", FilePath::input_file);
		FilePath mesh2("mesh/simplest_cube.msh", FilePath::input_file);
		ReaderInstances::instance()->get_reader(mesh1);
		ReaderInstances::instance()->get_reader(mesh2);
	}

	{
	    Mesh mesh;
		FilePath mesh_file("mesh/test_input.msh", FilePath::input_file);
		std::shared_ptr<GmshMeshReader> mesh_reader = ReaderInstances::instance()->get_reader(mesh_file);

		mesh_reader->read_mesh(&mesh);
		EXPECT_EQ(118, mesh.n_nodes());
	}

	{
	    Mesh mesh;
		FilePath mesh_file("mesh/simplest_cube.msh", FilePath::input_file);
		std::shared_ptr<GmshMeshReader> mesh_reader = ReaderInstances::instance()->get_reader(mesh_file);

		mesh_reader->read_mesh(&mesh);
		EXPECT_EQ(8, mesh.n_nodes());
	}
}
