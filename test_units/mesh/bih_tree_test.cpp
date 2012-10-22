/*
 * bih_tree_test.cpp
 *
 *  Created on: Oct 22, 2012
 *      Author: jb
 */




/*
 * gmsh_reader_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */

#define DEBUG

#include <gtest/gtest.h>
#include <sstream>
#include <string>

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "new_mesh/bih_tree.hh"

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



TEST(BIHTree_Test, create_tree) {
    stringstream ss(gmsh_mesh);

    Mesh mesh;
    GmshMeshReader reader;

    reader.read(ss, &mesh);

    // TODO: vytvorit strom a otestovat jeho vysku
    BIHTree bt(&mesh);
}


