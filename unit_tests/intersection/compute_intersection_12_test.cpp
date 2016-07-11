/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: VF, PE
 */
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

#include "system/system.hh"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"

#include "intersection/compute_intersection.hh"
#include "intersection/intersection_point_aux.hh"
#include <intersection/intersection_aux.hh>
#include <intersection/intersection_local.hh>

#include <dirent.h>

using namespace std;
using namespace computeintersection;


/// Create results for the meshes in directory 'simple_meshes_12d'.
void fill_12d_solution(std::vector<computeintersection::IntersectionLocal<1,2>> &ils)
{
    ils.clear();
    ils.resize(9);
    // ips[0] is empty
    ils[1].points() = {computeintersection::IntersectionPoint<1,2>({0}, {0, 0})};
    ils[2].points() = {computeintersection::IntersectionPoint<1,2>({1}, {0, 0})};
    ils[3].points() = {computeintersection::IntersectionPoint<1,2>({0.5}, {0, 0})};
    ils[4].points() = {computeintersection::IntersectionPoint<1,2>({0.5}, {0, 0})};
    ils[5].points() = {computeintersection::IntersectionPoint<1,2>({0.5}, {0.5, 0})};
    ils[6].points() = {computeintersection::IntersectionPoint<1,2>({0.5}, {0.25, 0.25})};
    ils[7].points() = {computeintersection::IntersectionPoint<1,2>({1./3}, {0.2, 0}),
                       computeintersection::IntersectionPoint<1,2>({2./3}, {0, 0.4})};
    ils[8].points() = {computeintersection::IntersectionPoint<1,2>({0}, {0, 0}),
                       computeintersection::IntersectionPoint<1,2>({1}, {1, 0})};
}


///Permutes triangle coordinates of IP<1,2> according to given permutation.
std::vector<computeintersection::IntersectionPoint<1,2>> permute_coords(std::vector<computeintersection::IntersectionPoint<1,2>> ips, 
                                                                        unsigned int permute[3])
{
    std::vector<computeintersection::IntersectionPoint<1,2>> new_points(ips.size());
    for(unsigned int i = 0; i < ips.size(); i++)
    {
        arma::vec::fixed<2> new_coords;
        for(unsigned int j = 0; j < 2; j++)
            new_coords[j] = ips[i].bulk_coords()[permute[j]];
        
//         new_points[i].set_coordinates(ips[i].bulk_coords(), new_coords);
        new_points[i] = computeintersection::IntersectionPoint<1,2>(ips[i].comp_coords(), new_coords);
    }
    return new_points;
}

void compute_intersection_12d(Mesh *mesh, const std::vector<computeintersection::IntersectionPoint<1,2>> &ips)
{
    DBGMSG("Computing 1d-2d intersections.\n");
    
    Simplex<1> line;
    Simplex<2> tria;
    arma::vec3 *points_tria[3], *points_line[2];
    for(unsigned int i=0; i < 3; i++)
        points_tria[i]= &(mesh->element(0)->node[i]->point());
    for(unsigned int i=0; i < 2; i++)
        points_line[i]= &(mesh->element(1)->node[i]->point());
    
    tria.set_simplices(points_tria);
    line.set_simplices(points_line);
    
    IntersectionAux<1,2> is(1, 0, 0);
    ComputeIntersection< Simplex<1>, Simplex<2>> CI(line, tria);
    CI.compute_final(is.points());
    
    computeintersection::IntersectionLocal<1,2> ilc(is);
    
//     cout << ilc;
    
    auto ipc = ilc.points();
//     for(computeintersection::IntersectionPoint<1,2> &ip: ipc)
//     {
//         ip.coords(mesh->element(1)).print();
//     }

    
    EXPECT_EQ(ipc.size(), ips.size());
    
    for(unsigned int i=0; i < ipc.size(); i++)
    {
        DBGMSG("---------- check IP[%d] ----------\n",i);
        EXPECT_DOUBLE_EQ(ipc[i].comp_coords()[0], ips[i].comp_coords()[0]);
        EXPECT_DOUBLE_EQ(ipc[i].bulk_coords()[0], ips[i].bulk_coords()[0]);
        EXPECT_DOUBLE_EQ(ipc[i].bulk_coords()[1], ips[i].bulk_coords()[1]);
    }
}


TEST(intersections_12d, all) {
 
    // directory with testing meshes
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/simple_meshes_12d/";
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
        ASSERT(0).error("Could not open directory with testing meshes.");
    }
    
    std::sort(filenames.begin(), filenames.end(), less<string>());
    
    std::vector<computeintersection::IntersectionLocal<1,2>> solution;
    fill_12d_solution(solution);
    
    // for each mesh, compute intersection area and compare with old NGH
    for(unsigned int s=0; s< filenames.size(); s++)
    {
        const unsigned int np = 6;
        unsigned int permutations[np][3] = {{0,1,2},
                                            {1,0,2},
                                            {1,2,0},
                                            {0,2,1},
                                            {2,0,1},
                                            {2,1,0}};
        for(unsigned int p=0; p<1; p++)
        {
            xprintf(Msg,"Computing intersection on mesh: %s\n",filenames[s].c_str());
            FilePath::set_io_dirs(".","","",".");
            FilePath mesh_file(dir_name + filenames[s], FilePath::input_file);
            
            Mesh mesh;
            // read mesh with gmshreader
            GmshMeshReader reader(mesh_file);
            reader.read_mesh(&mesh);
        
            // permute nodes:
            FOR_ELEMENTS(&mesh,ele)
            {
                if(ele->dim() == 2)
                {
                    Node* tmp[3];
                    for(unsigned int i=0; i<ele->n_nodes(); i++)
                    {
                        tmp[i] = ele->node[permutations[p][i]];
                    }
                    for(unsigned int i=0; i<ele->n_nodes(); i++)
                    {
                        ele->node[i] = tmp[i];
//                         ele->node[i]->point().print(cout);
                    }
//                     cout << p << ": jac = "  << ele->tetrahedron_jacobian() << endl;
                }
            }
            
            mesh.setup_topology();
            
            compute_intersection_12d(&mesh, permute_coords(solution[s].points(), permutations[p]));
        }
    }
}