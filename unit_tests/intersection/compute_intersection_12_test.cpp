/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: VF, PE
 */
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>
#include "arma_expect.hh"

#include "system/global_defs.h"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "mesh_constructor.hh"

#include "intersection/compute_intersection.hh"
#include "intersection/intersection_point_aux.hh"
#include <intersection/intersection_aux.hh>
#include <intersection/intersection_local.hh>

#include "compute_intersection_test.hh"

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
                                                                        const std::vector<unsigned int> &permute)
{
    std::vector<computeintersection::IntersectionPoint<1,2>> new_points(ips.size());
    for(unsigned int i = 0; i < ips.size(); i++)
    {
        arma::vec3 new_coords;
        arma::vec3 old_coords = {1,1,1};
        for(unsigned int j = 0; j < 2; j++){
            old_coords[j+1] = ips[i].bulk_coords()[j];
            old_coords[0] = old_coords[0] - ips[i].bulk_coords()[j];
        }
        if(old_coords[0] < 1e-15) old_coords[0] = 0;
        
        for(unsigned int j = 0; j < 3; j++)
            new_coords[j] = old_coords[permute[j]];
        
        new_points[i] = computeintersection::IntersectionPoint<1,2>(ips[i].comp_coords(), new_coords.subvec(1,2));
    }
    return new_points;
}

void compute_intersection_12d(Mesh *mesh, const std::vector<computeintersection::IntersectionPoint<1,2>> &ips)
{
    Simplex<1> line = create_simplex<1>(mesh->element(1));
    Simplex<2> tria = create_simplex<2>(mesh->element(0));
    
    IntersectionAux<1,2> is(1, 0);
    ComputeIntersection< Simplex<1>, Simplex<2>> CI(line, tria);
    CI.compute_final(is.points());
    
    computeintersection::IntersectionLocal<1,2> ilc(is);
    
//     DebugOut() << ilc;
    
    auto ipc = ilc.points();
//     for(computeintersection::IntersectionPoint<1,2> &ip: ipc)
//     {
//         ip.coords(mesh->element(1)).print(DebugOut(),"ip");
//     }

    
    EXPECT_EQ(ipc.size(), ips.size());
    
    for(unsigned int i=0; i < ipc.size(); i++)
    {
        MessageOut().fmt("---------- check IP[{}] ----------\n",i);
        EXPECT_ARMA_EQ(ipc[i].comp_coords(), ips[i].comp_coords());
        EXPECT_ARMA_EQ(ipc[i].bulk_coords(), ips[i].bulk_coords());
    }
}


TEST(intersections_12d, all) {
 
    // directory with testing meshes
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/simple_meshes_12d/";
    std::vector<string> filenames;
    
    read_files_from_dir(dir_name, "msh", filenames);
    
    std::vector<computeintersection::IntersectionLocal<1,2>> solution;
    fill_12d_solution(solution);
    
    // for each mesh, compute intersection area and compare with old NGH
    for(unsigned int s=0; s< filenames.size(); s++)
    {
        const unsigned int np = permutations_triangle.size();
        
        for(unsigned int p=0; p<np; p++)
        {
            MessageOut() << "Computing intersection on mesh: " << filenames[s] << "\n";
            FilePath mesh_file(dir_name + filenames[s], FilePath::input_file);
            
            Mesh *mesh = mesh_constructor();
            // read mesh with gmshreader
            GmshMeshReader reader(mesh_file);
            reader.read_mesh(mesh);
        
            // permute nodes:
            FOR_ELEMENTS(mesh,ele)
            {
                if(ele->dim() == 2)
                    permute_triangle(ele, p);
            }
            
            mesh->setup_topology();
            
            compute_intersection_12d(mesh, permute_coords(solution[s].points(), permutations_triangle[p]));
        }
    }
}
