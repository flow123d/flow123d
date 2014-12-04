/*
 * gmsh_reader_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */


#include <flow_gtest.hh>
#include <sstream>
#include <string>

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
