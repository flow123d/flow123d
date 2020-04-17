/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: VF, PE
 */
#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include "arma_expect.hh"

#include "system/global_defs.h"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "mesh_constructor.hh"

#include "intersection/intersection_point_aux.hh"
#include "intersection/intersection_local.hh"

#include "intersection/compute_intersection.hh"
#include "intersection/intersection_aux.hh"

#include "compute_intersection_test.hh"

using namespace std;

/// auxiliary function for sorting intersection point according to x,y local coordinates of component triangle
bool compare_ip22(const IntersectionPoint<2,2>& a,
                  const IntersectionPoint<2,2>& b)
{
    if (a.comp_coords()[0] == b.comp_coords()[0])
        return a.comp_coords()[1] <= b.comp_coords()[1];
    else
        return a.comp_coords()[0] < b.comp_coords()[0];
}

/// Create results for the meshes in directory 'simple_meshes_22d'.
void fill_22d_solution(std::vector<IntersectionLocal<2,2>> &ils)
{
    ils.clear();
    ils.resize(10);
    
    ils[0].points().push_back(IntersectionPoint<2,2>(arma::vec2({0,0.5}),arma::vec2({0.5,0.5})));
    
    ils[1].points().push_back(IntersectionPoint<2,2>(arma::vec2({0,0.5}),arma::vec2({0,1})));
    
    ils[2].points().push_back(IntersectionPoint<2,2>(arma::vec2({0,0}),arma::vec2({0,1})));
    
    ils[3].points().push_back(IntersectionPoint<2,2>(arma::vec2({0.25,0.25}),arma::vec2({0,1})));
    
    ils[4].points().push_back(IntersectionPoint<2,2>(arma::vec2({0,0}),arma::vec2({0,0})));
    ils[4].points().push_back(IntersectionPoint<2,2>(arma::vec2({0,1}),arma::vec2({0,1})));
    
    ils[5].points().push_back(IntersectionPoint<2,2>(arma::vec2({0,0}),arma::vec2({0,0})));
    ils[5].points().push_back(IntersectionPoint<2,2>(arma::vec2({0.25,0.75}),arma::vec2({0,1})));
    
    ils[6].points().push_back(IntersectionPoint<2,2>(arma::vec2({0,0}),arma::vec2({0,0})));
    ils[6].points().push_back(IntersectionPoint<2,2>(arma::vec2({0.25,0.5}),arma::vec2({0,1})));
    
    ils[7].points().push_back(IntersectionPoint<2,2>(arma::vec2({1,1})/4,arma::vec2({0.5,0})));
    ils[7].points().push_back(IntersectionPoint<2,2>(arma::vec2({2,5})/8,arma::vec2({0,0.5})));
    
    ils[8].points().push_back(IntersectionPoint<2,2>(arma::vec2({0,0.5}),arma::vec2({0.4,0.4})));
    ils[8].points().push_back(IntersectionPoint<2,2>(arma::vec2({0.25,0.5}),arma::vec2({0.5,0.5})));
    
    ils[9].points().push_back(IntersectionPoint<2,2>(arma::vec2({0,0.5}),arma::vec2({0,0})));
    ils[9].points().push_back(IntersectionPoint<2,2>(arma::vec2({2/3.0,1/3.0}),arma::vec2({2.0/30,2.0/30})));
}

// Permutes tetrahedron coordinates of IP<2,2> according to given permutation.
IntersectionLocal<2,2> permute_coords(IntersectionLocal<2,2> il,
                                                           std::vector<unsigned int> permute)
{
    IntersectionLocal<2,2> new_il = il;
    std::vector<IntersectionPoint<2,2>> & points = il.points();
    for(unsigned int i = 0; i < points.size(); i++)
    {
        arma::vec3 new_coords;
        arma::vec3 old_coords = {1,1,1};
        for(unsigned int j = 0; j < 2; j++){
            old_coords[j+1] = points[i].comp_coords()[j];
            old_coords[0] = old_coords[0] - points[i].comp_coords()[j];
        }
        if(old_coords[0] < 1e-15) old_coords[0] = 0;
        
        for(unsigned int j = 0; j < 3; j++)
            new_coords[j] = old_coords[permute[j]];
        
        new_il.points()[i] = IntersectionPoint<2,2>(new_coords.subvec(1,2), points[i].bulk_coords());
    }
    return new_il;
}

void compute_intersection_22d(Mesh *mesh, const IntersectionLocal<2,2> &il)
{
    IntersectionAux<2,2> is;
    ComputeIntersection<2,2> CI(mesh->element_accessor(0), mesh->element_accessor(1), mesh);
    CI.init();
    CI.compute(is);
    
//     DebugOut() << is;
//     for(IntersectionPointAux<2,2> &ip: is.points())
//     {
//         ip.coords(mesh->element(0)).print("ip");
//     }

    IntersectionLocal<2,2> ilc(is);
    EXPECT_EQ(ilc.size(), il.size());
    
    // sort ips in ils according local coordinates of component triangle
    std::sort(ilc.points().begin(), ilc.points().end(),compare_ip22);
    
    for(unsigned int i=0; i < is.size(); i++)
    {
        MessageOut().fmt("---------- check IP[{}] ----------\n",i);
        EXPECT_ARMA_EQ(il[i].comp_coords(), ilc[i].comp_coords());
        EXPECT_ARMA_EQ(il[i].bulk_coords(), ilc[i].bulk_coords());
    }
}


TEST(intersections_22d, all) {   
    // directory with testing meshes
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/simple_meshes_22d/";
    std::vector<string> filenames;
    
    read_files_from_dir(dir_name, "msh", filenames);
    
    std::vector<IntersectionLocal<2,2>> solution;
    fill_22d_solution(solution);
    
    // for each mesh, compute intersection area and compare with old NGH
    for(unsigned int s=0; s< filenames.size(); s++)
    {
        const unsigned int np = permutations_triangle.size();
        for(unsigned int p=0; p<np; p++)
        {
            MessageOut() << "Computing intersection on mesh: " << filenames[s] << "\n";
            string in_mesh_string = "{mesh_file=\"" + dir_name + filenames[s] + "\"}";
            
            Mesh *mesh = mesh_constructor(in_mesh_string);
            // read mesh with gmshreader
            auto reader = reader_constructor(in_mesh_string);
            reader->read_raw_mesh(mesh);
        
            // permute nodes of one triangle:
            mesh->permute_triangle(0, permutations_triangle[p]);
            
            mesh->setup_topology();
            
//             auto il = solution[s];
            auto il = permute_coords(solution[s], permutations_triangle[p]);
            std::sort(il.points().begin(), il.points().end(),compare_ip22);
            compute_intersection_22d(mesh, il);
        }
    }
}
