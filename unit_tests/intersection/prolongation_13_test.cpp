/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: VF, PE
 */
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

#include <armadillo>

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
void fill_13d_solution(std::vector<std::vector<std::vector<arma::vec3>>> &ils, std::vector<double> &lengths)
{
    unsigned int n_files=5;
    ils.clear();
    lengths.clear();
    ils.resize(n_files);
    lengths.resize(n_files);
    
    lengths[0] = 1.5*sqrt(0.27)+0.35+0.2 + sqrt(0.065) + 0.1*sqrt(18) + sqrt(0.105);
    lengths[1] = sqrt(0.24*0.24 + 0.0016 + 0.04) + sqrt( 2* 0.0225 + 0.75*0.75) + 
                 sqrt(0.0625 + 2*0.04) + sqrt(0.0225 + 0.09 + 0.01);
    lengths[2] = sqrt(0.24*0.24 + 0.0016 + 0.04) + sqrt( 2* 0.0225 + 0.75*0.75) + 
                 sqrt(0.0625 + 2*0.04) + sqrt(0.0225 + 0.09 + 0.01) +
                 0.1*sqrt(4+9+16) + 0.1*sqrt(4+9+1) + sqrt(1./9 + 0.25 + 1./36);
    lengths[3] = 4*0.25;
    lengths[4] = 0.5 + 0.5 + 2*0.5*sqrt(2);
    
    ils[0].resize(9);
    ils[0][0] = {arma::vec3({1,0,1})*0.3,arma::vec3({35,-25,35})*0.01};
    ils[0][1] = {arma::vec3({1,1,1})/4,arma::vec3({1,0,1})*0.3};
    ils[0][2] = {arma::vec3({1,1,1})/4,arma::vec3({0.2,0.5,0.2})};
    ils[0][3] = {arma::vec3({35,-25,35})*0.01,arma::vec3({0,-25,35})*0.01};
    ils[0][4] = {arma::vec3({0.2,0.5,0.2}),arma::vec3({0,5,2})*0.1};
    
    ils[0][5] = {arma::vec3({1,1,0})/4,arma::vec3({0,1,1})*0.05};
    ils[0][6] = {arma::vec3({0,1,1})*0.05,arma::vec3({2,0,1})*0.1};
    ils[0][7] = {arma::vec3({2,0,1})*0.1,arma::vec3({40,-5,15})*0.01};
    ils[0][8] = {arma::vec3({40,-5,15})*0.01,arma::vec3({1,-1,1})*0.2};

    
    ils[1].resize(7);
    ils[1][0] = {arma::vec3({0.44,0.46,0}),arma::vec3({0.2,0.5,0.2})};
    ils[1][1] = {arma::vec3({0.2,0.5,0.2}),arma::vec3({1,1,1})/4};
    ils[1][2] = {arma::vec3({1,1,1})/4,arma::vec3({0,1,1})*0.05};
    ils[1][3] = {arma::vec3({1,1,1})/4,arma::vec3({0.375,0,1./6})};
    ils[1][4] = {arma::vec3({0.375,0,1./6}),arma::vec3({0.4,-0.05,0.15})};
    ils[1][5] = {arma::vec3({1,1,1})/4,arma::vec3({1,0,1})*0.3};
    ils[1][6] = {arma::vec3({1,0,1})*0.3,arma::vec3({35,-25,35})*0.01};
    
    ils[2].resize(12);
    ils[2][0] = {arma::vec3({0.44,0.46,0}),arma::vec3({0.2,0.5,0.2})};
    ils[2][1] = {arma::vec3({0.2,0.5,0.2}),arma::vec3({1,1,1})/4};
    ils[2][2] = {arma::vec3({1,1,1})/4,arma::vec3({0,1,1})*0.05};
    ils[2][3] = {arma::vec3({1,1,1})/4,arma::vec3({0.375,0,1./6})};
    ils[2][4] = {arma::vec3({0.375,0,1./6}),arma::vec3({0.4,-0.05,0.15})};
    ils[2][5] = {arma::vec3({1,1,1})/4,arma::vec3({1,0,1})*0.3};
    ils[2][6] = {arma::vec3({1,0,1})*0.3,arma::vec3({35,-25,35})*0.01};
    
    ils[2][7] = {arma::vec3({1,-1,1})*0.2,arma::vec3({0,-5,3})*0.1};
    ils[2][8] = {arma::vec3({1,0,3.5})*0.1,arma::vec3({1,-1,1})*0.2};
    ils[2][9] = {arma::vec3({0,2,5})*0.1,arma::vec3({1,0,3.5})*0.1};
    ils[2][10] = {arma::vec3({1./3,0,2./3})};
    ils[2][11] = {arma::vec3({0,-5,5})*0.1,arma::vec3({1./3,0,2./3})};
    
    ils[3].resize(4);
    ils[3][0] = {arma::vec3({0.25,0.25,0.25}),arma::vec3({0.25,0.25,0})};
    ils[3][1] = {arma::vec3({0.25,0.25,-1}),arma::vec3({0.25,0.25,-1.25})};
    ils[3][2] = {arma::vec3({0.25,0.5,0.25}),arma::vec3({0.25,0.5,0})};
    ils[3][3] = {arma::vec3({0.25,0.5,-1}),arma::vec3({0.25,0.5,-1.25})};
    
    ils[4].resize(6);
    ils[4][1] = {arma::vec3({0.4,-0.25,0.25}),arma::vec3({0.4,0,0.25})};
    ils[4][0] = {arma::vec3({0.4,0,0.25}),arma::vec3({0.4,0.25,0.25})};
    
    ils[4][3] = {arma::vec3({0,-0.25,0.25}),arma::vec3({0,0,0.25})};
    ils[4][2] = {arma::vec3({0,0,0.25}),arma::vec3({0,0.25,0.25})};
        
    ils[4][5] = {arma::vec3({0.5,0,0.5}),arma::vec3({0,0,0})};
    ils[4][4] = {arma::vec3({0.5,0,0.5}),arma::vec3({0,0,0})};
}

void compute_intersection_13d(Mesh *mesh, const std::vector<std::vector<arma::vec3>> &il, const double &length)
{
    double computed_length = 0;

    // compute intersection
    DBGMSG("Computing intersection length by NEW algorithm\n");
    
    InspectElements ie(mesh);
    ie.compute_intersections();
    ie.print_mesh_to_file_13("output_intersection_13");
    
    DBGMSG("N intersections %d\n",ie.intersection_storage13_.size());
    
//    // write computed intersections
//     for(unsigned int i = 0; i < ie.intersection_storage13_.size(); i++)
//     {
//         cout << &ie.intersection_storage13_[i] << ie.intersection_storage13_[i];
//     }
    
//     //write the first intersection
//     FOR_ELEMENTS(mesh, elm){
//         
//         if( (elm->dim() == 1) && (ie.intersection_map_[elm->index()].size() > 0) )
//         {
//             computeintersection::IntersectionLocal<1,3>* il13 = 
//                 static_cast<computeintersection::IntersectionLocal<1,3>*> (ie.intersection_map_[elm->index()][0].second);
//             if(il13 != nullptr)
//             {
// //                 DBGMSG("comp idx %d, bulk idx %d, \n",elm->index(),ie.intersection_map_[elm->index()][0].first);
//                 cout << *il13;
//                 break;
//             }
//         }
//     }
    
    //test solution
    std::vector<computeintersection::IntersectionLocal<1,3>> ilc = ie.intersection_storage13_;
    EXPECT_EQ(ilc.size(), il.size());
    
    for(unsigned int i=0; i < ilc.size(); i++)
        for(unsigned int j=0; j < ilc[i].size(); j++)
    {
        DBGMSG("---------- check IP[%d] ----------\n",i);
        arma::vec3 ip = ilc[i][j].coords(mesh->element(ilc[i].component_ele_idx()));
        EXPECT_NEAR(ip[0], il[i][j][0], 1e-14);
        EXPECT_NEAR(ip[1], il[i][j][1], 1e-14);
        EXPECT_NEAR(ip[2], il[i][j][2], 1e-14);
    }
    
    computed_length = ie.measure_13();
    
    DBGMSG("Length of intersection line: (intersections) %.16e\n", computed_length);
    EXPECT_NEAR(computed_length, length, 1e-14);
}


TEST(intersection_prolongation_13d, all) {
    Profiler::initialize();
    
//     // directory with testing meshes
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/prolong_meshes_13d/";
    std::vector<string> filenames;
    
    // read mesh file names
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (dir_name.c_str())) != NULL) {
        // print all the files and directories within directory 
        xprintf(Msg,"Testing mesh files: \n");
        while ((ent = readdir (dir)) != NULL) {
            string fname = ent->d_name;
            // test extension ".msh"
            if(fname.size() >= 4)
            {
                string ext = fname.substr(fname.size()-4);
//                 xprintf(Msg,"%s\n",ext.c_str());
                if(ext == ".msh"){
                    filenames.push_back(ent->d_name);
                    xprintf(Msg,"%s\n",ent->d_name);
                }
            }
        }
        closedir (dir);
    } else {
        ASSERT(0,"Could not open directory with testing meshes.");
    }
    
    std::sort(filenames.begin(), filenames.end(), less<string>());
        
    std::vector<std::vector<std::vector<arma::vec3>>> solution;
    std::vector<double> lengths;
    fill_13d_solution(solution, lengths);
    
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
                compute_intersection_13d(&mesh, solution[s], lengths[s]); //permute_coords(solution[s], permutations[p]));
            xprintf(Msg, "==============\n");
//         }
    }
//     std::fstream fs;
//     fs.open(profiler_file.c_str(), std::fstream::out | std::fstream::app);
//     Profiler::instance()->output(PETSC_COMM_WORLD, fs);
    Profiler::uninitialize();
}

//*/