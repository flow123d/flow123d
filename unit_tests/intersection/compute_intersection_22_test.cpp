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

#include "intersection/inspectelements.h"
#include "intersection/intersectionpoint.h"
#include "intersection/intersection_local.h"

#include "intersection/computeintersection.h"
#include "intersection/intersectionaux.h"

#include <dirent.h>

using namespace std;
using namespace computeintersection;

static const std::string profiler_file = "compute_intersection_22d_profiler.log";
static const unsigned int profiler_loop = 1;


//*

// ******************************************************************************************* TEST 2d-2d ****

/// Create results for the meshes in directory 'simple_meshes_22d'.
void fill_22d_solution(std::vector<computeintersection::IntersectionLocal<2,2>> &ils)
{
    DBGMSG("fill solution\n");
    ils.clear();
    ils.resize(9);
    
    ils[0].points().push_back(computeintersection::IntersectionPoint<2,2>(arma::vec2({0,0.5}),arma::vec2({0.5,0.5})));
    
    ils[1].points().push_back(computeintersection::IntersectionPoint<2,2>(arma::vec2({0,0.5}),arma::vec2({0,1})));
    
    ils[2].points().push_back(computeintersection::IntersectionPoint<2,2>(arma::vec2({0,0}),arma::vec2({0,1})));
    
    ils[3].points().push_back(computeintersection::IntersectionPoint<2,2>(arma::vec2({0.25,0.25}),arma::vec2({0,1})));
    
    ils[4].points().push_back(computeintersection::IntersectionPoint<2,2>(arma::vec2({0,0}),arma::vec2({0,0})));
    ils[4].points().push_back(computeintersection::IntersectionPoint<2,2>(arma::vec2({0,1}),arma::vec2({0,1})));
    
    ils[5].points().push_back(computeintersection::IntersectionPoint<2,2>(arma::vec2({0,0}),arma::vec2({0,0})));
    ils[5].points().push_back(computeintersection::IntersectionPoint<2,2>(arma::vec2({0.25,0.75}),arma::vec2({0,1})));
    
    ils[6].points().push_back(computeintersection::IntersectionPoint<2,2>(arma::vec2({0,0}),arma::vec2({0,0})));
    ils[6].points().push_back(computeintersection::IntersectionPoint<2,2>(arma::vec2({0.25,0.5}),arma::vec2({0,1})));
    
    ils[7].points().push_back(computeintersection::IntersectionPoint<2,2>(arma::vec2({1,1})/4,arma::vec2({0.5,0})));
    ils[7].points().push_back(computeintersection::IntersectionPoint<2,2>(arma::vec2({2,5})/8,arma::vec2({0,0.5})));
    
    ils[8].points().push_back(computeintersection::IntersectionPoint<2,2>(arma::vec2({0,0.5}),arma::vec2({0.4,0.4})));
    ils[8].points().push_back(computeintersection::IntersectionPoint<2,2>(arma::vec2({0.25,0.5}),arma::vec2({0.5,0.5})));
}


//Permutes tetrahedron coordinates of IP<2,2> according to given permutation.
// computeintersection::IntersectionLocal<1,3> permute_coords(computeintersection::IntersectionLocal<1,3> il, unsigned int permute[4])
// {
//     computeintersection::IntersectionLocal<1,3> new_il = il;
//     std::vector<computeintersection::IntersectionPoint<1,3>> & points = il.points();
//     for(unsigned int i = 0; i < points.size(); i++)
//     {
//         arma::vec4 new_coords;
//         arma::vec4 old_coords = {1,1,1,1};
//         for(unsigned int j = 0; j < 3; j++){
//             old_coords[j+1] = points[i].bulk_coords()[j];
//             old_coords[0] = old_coords[0] - points[i].bulk_coords()[j];
//         }
//         if(old_coords[0] < 1e-15) old_coords[0] = 0;
//         
//         for(unsigned int j = 0; j < 4; j++)
//             new_coords[j] = old_coords[permute[j]];
//         
//         new_il.points()[i] = computeintersection::IntersectionPoint<1,3>(points[i].comp_coords(), new_coords.subvec(1,3));
//     }
//     return new_il;
// }

void compute_intersection_22d(Mesh *mesh, const computeintersection::IntersectionLocal<2,2> &il)
{
    Simplex<2> triaA, triaB;
    arma::vec3 *pointsA[3], *pointsB[3];
    for(unsigned int i=0; i < 3; i++){
        pointsA[i]= &(mesh->element(0)->node[i]->point());
        pointsB[i]= &(mesh->element(1)->node[i]->point());
    }
    triaA.set_simplices(pointsA);
    triaB.set_simplices(pointsB);
    
    IntersectionAux<2,2> is;
    std::vector<unsigned int> prolong_table;
    ComputeIntersection< Simplex<2>, Simplex<2>> CI(triaA, triaB);
    CI.init();
    CI.compute(is, prolong_table);
    
    cout << is;
    for(IntersectionPointAux<2,2> &ip: is.points())
    {
        ip.coords(mesh->element(0)).print();
    }

    computeintersection::IntersectionLocal<2,2> ilc(is);
    EXPECT_EQ(ilc.size(), il.size());

    
    for(unsigned int i=0; i < is.size(); i++)
    {
        DBGMSG("---------- check IP[%d] ----------\n",i);
        EXPECT_DOUBLE_EQ(ilc[i].comp_coords()[0], il[i].comp_coords()[0]);
        EXPECT_DOUBLE_EQ(ilc[i].comp_coords()[1], il[i].comp_coords()[1]);
        EXPECT_DOUBLE_EQ(ilc[i].bulk_coords()[0], il[i].bulk_coords()[0]);
        EXPECT_DOUBLE_EQ(ilc[i].bulk_coords()[1], il[i].bulk_coords()[1]);
    }
}


TEST(intersections_22d, all) {
    Profiler::initialize();
    
    // directory with testing meshes
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/simple_meshes_22d/";
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
    
    std::vector<computeintersection::IntersectionLocal<2,2>> solution;
    fill_22d_solution(solution);
    
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
        {
            xprintf(Msg,"Computing intersection on mesh: %s\n",filenames[s].c_str());
            FilePath::set_io_dirs(".","","",".");
            FilePath mesh_file(dir_name + filenames[s], FilePath::input_file);
            
            Mesh mesh;
            // read mesh with gmshreader
            GmshMeshReader reader(mesh_file);
            reader.read_mesh(&mesh);
        
            // permute nodes:
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
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
                compute_intersection_22d(&mesh, solution[s]);//, permute_coords(solution[s], permutations[p]));
            xprintf(Msg, "==============\n");
        }
    }
    std::fstream fs;
    fs.open(profiler_file.c_str(), std::fstream::out | std::fstream::app);
    Profiler::instance()->output(PETSC_COMM_WORLD, fs);
    Profiler::uninitialize();
}

//*/