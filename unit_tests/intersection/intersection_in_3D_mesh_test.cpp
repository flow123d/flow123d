
/*
 *
 *      Author: PE
 */
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

#include "system/global_defs.h"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "mesh_constructor.hh"

#include "intersection/mixed_mesh_intersections.hh"
#include "intersection/intersection_point_aux.hh"
#include "intersection/intersection_aux.hh"
#include "intersection/intersection_local.hh"

#include <dirent.h>

using namespace std;
using namespace computeintersection;


/// Create results for the meshes in directory 'prolong_meshes_23d'.
// void fill_23d_solution(std::vector<std::vector<std::vector<arma::vec3>>> &ils)
// {
//     unsigned int n_files=1;
//     ils.clear();
//     ils.resize(n_files);
//     
//     ils[0].resize(4);
//     ils[0][0] = {   {0.475, 0.325, 0}, 
//                     {0.25, 0.25, 0},
//                     {0, 0.05, 0.05},
//                     {0.2, 0, 0.1},
//                     {0.41, 0, 0.13}};
//     ils[0][1] = {   {0.2, 0, 0.1},
//                     {0.4, -0.05, 0.15},
//                     {0.41, 0, 0.13}};
//     ils[0][2] = {   {0, 0.05, 0.05},
//                     {0.2, 0, 0.1},
//                     {0.04, 0, 0.08}};
//     ils[0][3] = {   {0.2, 0, 0.1},
//                     {0.4, -0.05, 0.15},
//                     {0.2, -0.2 ,0.2},
//                     {0.04, 0, 0.08}};
//     
// 
// }

/// auxiliary function for sorting intersection storage 13d
// bool compare_is23(const computeintersection::IntersectionLocal<2,3>& a,
//                   const computeintersection::IntersectionLocal<2,3>& b)
// {
//     if (a.component_ele_idx() == b.component_ele_idx())
//         return a.bulk_ele_idx() <= b.bulk_ele_idx();
//     else
//         return a.component_ele_idx() < b.component_ele_idx();
// }

void compute_intersection(Mesh *mesh)
{

    // compute intersection
    MixedMeshIntersections ie(mesh);
    ie.compute_intersections(IntersectionType(IntersectionType::d23
                                            | IntersectionType::d22));
    
    double total_length = ie.measure_22();
    cout << "total_length = " << total_length << endl;
    
    // write computed intersections
//     for(unsigned int i = 0; i < ie.intersection_storage12_.size(); i++)
//     {
//         cout << ie.intersection_storage12_[i];
//     }
    
//     //test solution
//     std::vector<computeintersection::IntersectionLocal<2,3>> ilc = ie.intersection_storage23_;
//     
//     // sort the storage, so it is the same for every algorithm (BIH, BB ...)
//     // and we avoid creating the intersection map for exact IPs
//     std::sort(ilc.begin(), ilc.end(),compare_is23);
//     
//     EXPECT_EQ(ilc.size(), il.size());
//     
//     for(unsigned int i=0; i < ilc.size(); i++)
//         for(unsigned int j=0; j < ilc[i].size(); j++)
//     {
//         DBGMSG("---------- check IP[%d] ----------\n",i);
//         arma::vec3 ip = ilc[i][j].coords(mesh->element(ilc[i].component_ele_idx()));
//         EXPECT_NEAR(ip[0], il[i][j][0], 1e-14);
//         EXPECT_NEAR(ip[1], il[i][j][1], 1e-14);
//         EXPECT_NEAR(ip[2], il[i][j][2], 1e-14);
//     }
}


TEST(intersection_prolongation_23d, all) {
    
//     // directory with testing meshes
    FilePath::set_dirs(UNIT_TESTS_SRC_DIR,"",".");
    string dir_name = "intersection/2d-2d/";
    
    string filename = dir_name + "cube_2f_incomp.msh";

    MessageOut() << "Computing intersection on mesh: " << filename << "\n";
    FilePath mesh_file(filename, FilePath::input_file);
    
    Mesh *mesh = mesh_constructor();
    ifstream in(string(mesh_file).c_str());
    mesh->read_gmsh_from_stream(in);
    
    compute_intersection(mesh);
}

//*/
