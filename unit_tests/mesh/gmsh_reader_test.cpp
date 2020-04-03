/*
 * gmsh_reader_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */

#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <sstream>
#include <string>
#include <mesh_constructor.hh>

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"



/*TEST(GMSHReader, read_mesh_from_stream) {
    Profiler::instance();
    
    string fname = string(UNIT_TESTS_SRC_DIR) + "/mesh/simplest_cube.msh";
    ifstream ifs( fname.c_str() );
    if (! ifs) cout << "Can not open file!" << endl;
    stringstream ss;
    ss << ifs.rdbuf();

    GmshMeshReader reader(ss);
    Mesh * mesh = new Mesh();
	reader.read_physical_names(mesh);
	reader.read_raw_mesh(mesh);

    EXPECT_EQ(9, mesh->n_elements());

    delete mesh;
}*/


TEST(GMSHReader, read_mesh_from_file) {
    Profiler::instance();

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	std::string mesh_in_string = "{mesh_file=\"mesh/test_input.msh\"}";
	Mesh * mesh = mesh_constructor(mesh_in_string);
    auto reader = reader_constructor(mesh_in_string);
	reader->read_physical_names(mesh);
	reader->read_raw_mesh(mesh);

    EXPECT_EQ(118, mesh->n_nodes());
    EXPECT_EQ(216, mesh->n_elements());

    delete mesh;
}
