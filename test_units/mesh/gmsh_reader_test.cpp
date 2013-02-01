/*
 * gmsh_reader_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */


#include <gtest/gtest.h>
#include <sstream>
#include <string>

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"



TEST(GMSHReader, read_mesh_from_stream) {
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
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    FilePath mesh_file("mesh/test_input.msh", FilePath::input_file);

    Mesh mesh;
    GmshMeshReader reader(mesh_file);

    reader.read_mesh(&mesh);

    EXPECT_EQ(216, mesh.n_elements());
}

