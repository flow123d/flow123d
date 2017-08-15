
/*
 *
 *      Author: PE
 */
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

#include "system/global_defs.h"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "mesh_constructor.hh"

#include "intersection/mixed_mesh_intersections.hh"
#include "intersection/intersection_point_aux.hh"
#include "intersection/intersection_aux.hh"
#include "intersection/intersection_local.hh"

#include <dirent.h>

using namespace std;


void compute_intersection(Mesh *mesh)
{
    // compute intersection
//     MixedMeshIntersections ie(mesh);
//     ie.compute_intersections(IntersectionType(IntersectionType::d23
//                                             | IntersectionType::d22));

    mesh->mixed_intersections();
    
    double total_length = mesh->mixed_intersections().measure_22();
    cout << "total_length = " << total_length << endl;
    
    // write computed intersections
//     for(unsigned int i = 0; i < ie.intersection_storage12_.size(); i++)
//     {
//         cout << ie.intersection_storage12_[i];
//     }
}


TEST(intersection_prolongation_23d, all) {
    
//     // directory with testing meshes
    FilePath::set_dirs(UNIT_TESTS_SRC_DIR,"",".");
    string dir_name = "intersection/2d-2d/";
    
    string filename = dir_name + "cube_2f_incomp.msh";

    MessageOut() << "Computing intersection on mesh: " << filename << "\n";
    //FilePath mesh_file(filename, FilePath::input_file);
    string in_mesh_string = "{mesh_file=\"" + filename + "\"}";
    
    Mesh *mesh = mesh_full_constructor(in_mesh_string);
    compute_intersection(mesh);
}
