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

// simplest cube 123d
string gmsh_mesh = R"CODE(
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
8
1 1 1 1
2 -1 1 1
3 -1 -1 1
4 1 -1 1
5 1 -1 -1
6 -1 -1 -1
7 1 1 -1
8 -1 1 -1
$EndNodes
$Elements
9
1 1 2 37 20 7 3
2 2 2 38 34 6 3 7
3 2 2 38 36 3 1 7
4 4 2 39 40 3 7 1 2
5 4 2 39 40 3 7 2 8
6 4 2 39 40 3 7 8 6
7 4 2 39 42 3 7 6 5
8 4 2 39 42 3 7 5 4
9 4 2 39 42 3 7 4 1
$EndElements
)CODE";

string ngh_input = R"CODE(
$NeighbourFormat
1.0 0 8
$EndNeighbourFormat
$Neighbours
10
0   11 2 4 2 5 3 
1   11 2 5 2 6 3 
2   11 2 7 2 8 3 
3   11 2 8 2 9 3 
4   20 1 2 1 1.0
5   20 1 3 2 1.0
6   20 3 4 3 1.0
7   20 2 6 2 1.0
8   20 2 7 3 1.0
9   20 3 9 2 1.0
$EndNeighbours
)CODE";


TEST(GMSHReader, read_mesh_from_stream) {
    string fname = string(UNIT_TESTS_SRC_DIR) + "/mesh/simplest_cube.msh";
    ifstream ifs( fname.c_str() );
    if (! ifs) cout << "Can not open file!" << endl;
    stringstream ss;
    ss << ifs.rdbuf();

    Mesh mesh;
    GmshMeshReader reader(ss);

    reader.read_mesh(&mesh);

    stringstream ngh_ss(ngh_input);
    mesh.setup_topology(&ngh_ss);

    EXPECT_EQ(13, mesh.n_elements());
}


TEST(GMSHReader, read_mesh_from_file) {
    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/test_input.msh", FilePath::input_file);

    Mesh mesh;
    GmshMeshReader reader(mesh_file);

    reader.read_mesh(&mesh);

    EXPECT_EQ(216, mesh.n_elements());
}

