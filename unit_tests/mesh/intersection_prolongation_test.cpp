/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: VF, PE
 */
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"

#include "mesh/ngh/include/point.h"
#include "mesh/ngh/include/intersection.h"

#include "intersection/inspectelements.h"
#include "intersection/intersectionpoint.h"
#include "intersection/intersectionaux.h"
#include "intersection/intersection_local.h"

#include <dirent.h>

using namespace std;
using namespace computeintersection;

//static const std::string profiler_file = "intersection_profiler.log";
//static const unsigned int profiler_loop = 10000;

//*

// ******************************************************************************************* TEST 1d-3d ****

/// Create results for the meshes in directory 'site_13d'.
void fill_13d_solution(std::vector<computeintersection::IntersectionAux<1,3>> &ils)
{
    ils.clear();
    ils.resize(1);
    
    ils[0].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({1,0}),arma::vec::fixed<4>({1,1,1,1})/4));
    ils[0].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({1,1})/2,arma::vec::fixed<4>({4,3,0,3})*0.1));
}


//Permutes tetrahedron coordinates of IP<1,3> according to given permutation.
// computeintersection::IntersectionLine permute_coords(computeintersection::IntersectionLine il, unsigned int permute[4])
// {
//     computeintersection::IntersectionLine new_il = il;
//     std::vector<computeintersection::IntersectionPoint<1,3>> & points = il.points();
//     for(unsigned int i = 0; i < points.size(); i++)
//     {
//         arma::vec::fixed<4> new_coords;
//         for(unsigned int j = 0; j < 4; j++)
//             new_coords[j] = points[i].local_bcoords_B()[permute[j]];
//         
//         new_il.points()[i].set_coordinates(points[i].local_bcoords_A(), new_coords);
//     }
//     return new_il;
// }

void compute_intersection_13d(Mesh *mesh, const computeintersection::IntersectionAux<1,3> &il)
{
    double length = 0;

    // compute intersection
    DBGMSG("Computing intersection length by NEW algorithm\n");

//     InspectElementsAlgorithm<1> ie(mesh);
//     ie.compute_intersections();
    
    InspectElements ie(mesh);
    ie.compute_intersections();
    
    DBGMSG("N intersections %d\n",ie.intersection_storage13_.size());
    
//     for(unsigned int i = 0; i < ie.intersection_storage13_.size(); i++)
//     {
//         cout << &ie.intersection_storage13_[i] << ie.intersection_storage13_[i];
//     }
    
    // write the first intersection
    FOR_ELEMENTS(mesh, elm){
        
        if( (elm->dim() == 1) && (ie.intersection_map_[elm->index()].size() > 0) )
        {
            computeintersection::IntersectionLocal<1,3>* il13 = 
                static_cast<computeintersection::IntersectionLocal<1,3>*> (ie.intersection_map_[elm->index()][0].second);
            if(il13 != nullptr)
            {
//                 DBGMSG("comp idx %d, bulk idx %d, \n",elm->index(),ie.intersection_map_[elm->index()][0].first);
                cout << *il13;
                break;
            }
        }
    }
//     //test solution
//     std::vector<computeintersection::IntersectionLine> pp = ie.list_intersection_lines(1);
//     computeintersection::IntersectionLine ilc;
//     // component = element index == 1
//     if(pp.size() > 0)
//     {
//         ilc = pp[0];
//         EXPECT_EQ(ilc.size(), il.size());
//     }
//     
//     for(unsigned int i=0; i < ilc.size(); i++)
//     {
//         DBGMSG("---------- check IP[%d] ----------\n",i);
//         EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_A()[0], il[i].local_bcoords_A()[0]);
//         EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_A()[1], il[i].local_bcoords_A()[1]);
//         EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_B()[0], il[i].local_bcoords_B()[0]);
//         EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_B()[1], il[i].local_bcoords_B()[1]);
//         EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_B()[2], il[i].local_bcoords_B()[2]);
//         EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_B()[3], il[i].local_bcoords_B()[3]);
//     }
    
    length = ie.measure_13();
    ie.print_mesh_to_file_13("output_intersection_13");
    DBGMSG("Length of intersection line: (intersections) %.16e\n", length);
//     EXPECT_NEAR(length1, length2, 1e-12);
    EXPECT_DOUBLE_EQ(length,1.5*std::sqrt(0.27)+0.35+0.2
                            + std::sqrt(0.065) + 0.1*std::sqrt(18) + std::sqrt(0.105));
}


TEST(intersection_prolongation_13d, all) {
    Profiler::initialize();
    
//     // directory with testing meshes
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/mesh/meshes_prolongation/";
//     std::vector<string> filenames;
//     
//     // read mesh file names
//     DIR *dir;
//     struct dirent *ent;
//     if ((dir = opendir (dir_name.c_str())) != NULL) {
//         // print all the files and directories within directory 
//         xprintf(Msg,"Testing mesh files: \n");
//         while ((ent = readdir (dir)) != NULL) {
//             string fname = ent->d_name;
//             // test extension ".msh"
//             if(fname.size() >= 4)
//             {
//                 string ext = fname.substr(fname.size()-4);
// //                 xprintf(Msg,"%s\n",ext.c_str());
//                 if(ext == ".msh"){
//                     filenames.push_back(ent->d_name);
//                     xprintf(Msg,"%s\n",ent->d_name);
//                 }
//             }
//         }
//         closedir (dir);
//     } else {
//         ASSERT(0,"Could not open directory with testing meshes.");
//     }
//     
//     std::sort(filenames.begin(), filenames.end(), less<string>());
    
    std::vector<string> filenames = {   "prolongation_13d_01.msh",
                                        "prolongation_13d_02.msh",
                                        "prolongation_13d_03.msh",
                                        "prolongation_13d_04.msh"
    };
    
    std::vector<IntersectionAux<1,3>> solution;
    fill_13d_solution(solution);
    
    // for each mesh, compute intersection area and compare with old NGH
    for(unsigned int s=0; s< filenames.size(); s++)
    {
//         const unsigned int np = 24;
//         unsigned int permutations[np][4] = {{0,1,2,3},
//                                                 {0,1,3,2},  // the tab means permutation with negative jacobian
//                                             {0,3,1,2},
//                                                 {0,3,2,1},
//                                             {0,2,3,1},
//                                                 {0,2,1,3},
//                                                 {1,0,2,3},
//                                             {1,0,3,2},
//                                                 {1,3,0,2},
//                                             {1,3,2,0},
//                                                 {1,2,3,0},
//                                             {1,2,0,3},
//                                                 {2,1,0,3},
//                                             {2,1,3,0},
//                                                 {2,3,1,0},
//                                             {2,3,0,1},
//                                                 {2,0,3,1},
//                                             {2,0,1,3},
//                                                 {3,1,2,0},
//                                             {3,1,0,2},
//                                                 {3,0,1,2},
//                                             {3,0,2,1},
//                                                 {3,2,0,1},
//                                             {3,2,1,0}};
//         for(unsigned int p=0; p<np; p++)
//         {
            xprintf(Msg,"Computing intersection on mesh: %s\n",filenames[s].c_str());
            FilePath::set_io_dirs(".","","",".");
            FilePath mesh_file(dir_name + filenames[s], FilePath::input_file);
            
            Mesh mesh;
            // read mesh with gmshreader
            GmshMeshReader reader(mesh_file);
            reader.read_mesh(&mesh);
        
//             // permute nodes:
//             FOR_ELEMENTS(&mesh,ele)
//             {
//                 if(ele->dim() == 3)
//                 {
//                     Node* tmp[4];
//                     for(unsigned int i=0; i<ele->n_nodes(); i++)
//                     {
//                         tmp[i] = ele->node[permutations[p][i]];
//                     }
//                     for(unsigned int i=0; i<ele->n_nodes(); i++)
//                     {
//                         ele->node[i] = tmp[i];
// //                         ele->node[i]->point().print(cout);
//                     }
// //                     cout << p << ": jac = "  << ele->tetrahedron_jacobian() << endl;
//                 }
//             }
            
            mesh.setup_topology();
            
            xprintf(Msg, "==============\n");
//             for(unsigned int loop = 0; loop < profiler_loop; loop++)
                compute_intersection_13d(&mesh, solution[s]); //permute_coords(solution[s], permutations[p]));
            xprintf(Msg, "==============\n");
//         }
    }
//     std::fstream fs;
//     fs.open(profiler_file.c_str(), std::fstream::out | std::fstream::app);
//     Profiler::instance()->output(PETSC_COMM_WORLD, fs);
    Profiler::uninitialize();
}

//*/




















// ******************************************************************************************* TEST 2d-3d ****

/// Create results for the meshes in directory 'site_13d'.
void fill_23d_solution(std::vector<computeintersection::IntersectionAux<2,3>> &ils)
{
    ils.clear();
    ils.resize(1);
    
    ils[0].points().push_back(computeintersection::IntersectionPoint<2,3>(arma::vec::fixed<3>({1,0,0}),arma::vec::fixed<4>({1,1,1,1})/4));
    ils[0].points().push_back(computeintersection::IntersectionPoint<2,3>(arma::vec::fixed<3>({1,1,1})/3,arma::vec::fixed<4>({4,3,0,3})*0.1));
}


//Permutes tetrahedron coordinates of IP<1,3> according to given permutation.
// computeintersection::IntersectionLine permute_coords(computeintersection::IntersectionLine il, unsigned int permute[4])
// {
//     computeintersection::IntersectionLine new_il = il;
//     std::vector<computeintersection::IntersectionPoint<1,3>> & points = il.points();
//     for(unsigned int i = 0; i < points.size(); i++)
//     {
//         arma::vec::fixed<4> new_coords;
//         for(unsigned int j = 0; j < 4; j++)
//             new_coords[j] = points[i].local_bcoords_B()[permute[j]];
//         
//         new_il.points()[i].set_coordinates(points[i].local_bcoords_A(), new_coords);
//     }
//     return new_il;
// }

void compute_intersection_23d(Mesh *mesh, const computeintersection::IntersectionAux<2,3> &ip)
{
    double area = 0;

    // compute intersection
    DBGMSG("Computing intersection area by NEW algorithm\n");
//     InspectElements ie(mesh);
//     ie.compute_intersections<2,3>();
    InspectElements ie(mesh);
    ie.compute_intersections();
//     //test solution
//     std::vector<computeintersection::IntersectionLine> pp = ie.list_intersection_lines(1);
//     computeintersection::IntersectionLine ilc;
//     // component = element index == 1
//     if(pp.size() > 0)
//     {
//         ilc = pp[0];
//         EXPECT_EQ(ilc.size(), il.size());
//     }
//     
//     for(unsigned int i=0; i < ilc.size(); i++)
//     {
//         DBGMSG("---------- check IP[%d] ----------\n",i);
//         EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_A()[0], il[i].local_bcoords_A()[0]);
//         EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_A()[1], il[i].local_bcoords_A()[1]);
//         EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_B()[0], il[i].local_bcoords_B()[0]);
//         EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_B()[1], il[i].local_bcoords_B()[1]);
//         EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_B()[2], il[i].local_bcoords_B()[2]);
//         EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_B()[3], il[i].local_bcoords_B()[3]);
//     }
    
    area = ie.measure_23();
    ie.print_mesh_to_file_23("output_intersection_23");
    DBGMSG("Area of intersection line: (intersections) %.16e\n", area);
//     EXPECT_NEAR(length1, length2, 1e-12);
//     EXPECT_DOUBLE_EQ(length,1.5*std::sqrt(0.27)+0.35+0.2
//                             + std::sqrt(0.065) + 0.1*std::sqrt(18) + std::sqrt(0.105));
}


TEST(intersection_prolongation_23d, all) {
    Profiler::initialize();
    
//     // directory with testing meshes
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/mesh/meshes_prolongation/";
//     std::vector<string> filenames;
//     
//     // read mesh file names
//     DIR *dir;
//     struct dirent *ent;
//     if ((dir = opendir (dir_name.c_str())) != NULL) {
//         // print all the files and directories within directory 
//         xprintf(Msg,"Testing mesh files: \n");
//         while ((ent = readdir (dir)) != NULL) {
//             string fname = ent->d_name;
//             // test extension ".msh"
//             if(fname.size() >= 4)
//             {
//                 string ext = fname.substr(fname.size()-4);
// //                 xprintf(Msg,"%s\n",ext.c_str());
//                 if(ext == ".msh"){
//                     filenames.push_back(ent->d_name);
//                     xprintf(Msg,"%s\n",ent->d_name);
//                 }
//             }
//         }
//         closedir (dir);
//     } else {
//         ASSERT(0,"Could not open directory with testing meshes.");
//     }
//     
//     std::sort(filenames.begin(), filenames.end(), less<string>());
    
    std::vector<string> filenames = {"prolongation_23d_01.msh"};
    
    std::vector<IntersectionAux<2,3>> solution;
    fill_23d_solution(solution);
    
    // for each mesh, compute intersection area and compare with old NGH
    for(unsigned int s=0; s< filenames.size(); s++)
    {
//         const unsigned int np = 24;
//         unsigned int permutations[np][4] = {{0,1,2,3},
//                                                 {0,1,3,2},  // the tab means permutation with negative jacobian
//                                             {0,3,1,2},
//                                                 {0,3,2,1},
//                                             {0,2,3,1},
//                                                 {0,2,1,3},
//                                                 {1,0,2,3},
//                                             {1,0,3,2},
//                                                 {1,3,0,2},
//                                             {1,3,2,0},
//                                                 {1,2,3,0},
//                                             {1,2,0,3},
//                                                 {2,1,0,3},
//                                             {2,1,3,0},
//                                                 {2,3,1,0},
//                                             {2,3,0,1},
//                                                 {2,0,3,1},
//                                             {2,0,1,3},
//                                                 {3,1,2,0},
//                                             {3,1,0,2},
//                                                 {3,0,1,2},
//                                             {3,0,2,1},
//                                                 {3,2,0,1},
//                                             {3,2,1,0}};
//         for(unsigned int p=0; p<np; p++)
//         {
            xprintf(Msg,"Computing intersection on mesh: %s\n",filenames[s].c_str());
            FilePath::set_io_dirs(".","","",".");
            FilePath mesh_file(dir_name + filenames[s], FilePath::input_file);
            
            Mesh mesh;
            // read mesh with gmshreader
            GmshMeshReader reader(mesh_file);
            reader.read_mesh(&mesh);
        
//             // permute nodes:
//             FOR_ELEMENTS(&mesh,ele)
//             {
//                 if(ele->dim() == 3)
//                 {
//                     Node* tmp[4];
//                     for(unsigned int i=0; i<ele->n_nodes(); i++)
//                     {
//                         tmp[i] = ele->node[permutations[p][i]];
//                     }
//                     for(unsigned int i=0; i<ele->n_nodes(); i++)
//                     {
//                         ele->node[i] = tmp[i];
// //                         ele->node[i]->point().print(cout);
//                     }
// //                     cout << p << ": jac = "  << ele->tetrahedron_jacobian() << endl;
//                 }
//             }
            
            mesh.setup_topology();
            
            xprintf(Msg, "==============\n");
//             for(unsigned int loop = 0; loop < profiler_loop; loop++)
                compute_intersection_23d(&mesh, solution[s]); //permute_coords(solution[s], permutations[p]));
            xprintf(Msg, "==============\n");
//         }
    }
//     std::fstream fs;
//     fs.open(profiler_file.c_str(), std::fstream::out | std::fstream::app);
//     Profiler::instance()->output(PETSC_COMM_WORLD, fs);
    Profiler::uninitialize();
}

//*/