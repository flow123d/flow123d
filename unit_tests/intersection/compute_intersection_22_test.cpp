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

#include "intersection/intersection_point_aux.hh"
#include "intersection/intersection_local.hh"

#include "intersection/compute_intersection.hh"
#include "intersection/intersection_aux.hh"

#include "compute_intersection_test.hh"

using namespace std;
using namespace computeintersection;

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


// Permutes tetrahedron coordinates of IP<2,2> according to given permutation.
computeintersection::IntersectionLocal<2,2> permute_coords(computeintersection::IntersectionLocal<2,2> il,
                                                           std::vector<unsigned int> permute)
{
    computeintersection::IntersectionLocal<2,2> new_il = il;
    std::vector<computeintersection::IntersectionPoint<2,2>> & points = il.points();
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
        
        new_il.points()[i] = computeintersection::IntersectionPoint<2,2>(new_coords.subvec(1,2), points[i].bulk_coords());
    }
    return new_il;
}

void compute_intersection_22d(Mesh *mesh, const computeintersection::IntersectionLocal<2,2> &il)
{
    Simplex<2> triaA = create_simplex<2>(mesh->element(0)),
               triaB = create_simplex<2>(mesh->element(1));
    
    IntersectionAux<2,2> is;
    std::vector<unsigned int> prolong_table;
    ComputeIntersection< Simplex<2>, Simplex<2>> CI(triaA, triaB);
    CI.init();
    CI.compute(is, prolong_table);
    
//     cout << is;
//     for(IntersectionPointAux<2,2> &ip: is.points())
//     {
//         ip.coords(mesh->element(0)).print();
//     }

    computeintersection::IntersectionLocal<2,2> ilc(is);
    EXPECT_EQ(ilc.size(), il.size());

    
    for(unsigned int i=0; i < is.size(); i++)
    {
        DBGMSG("---------- check IP[%d] ----------\n",i);
        EXPECT_DOUBLE_EQ(il[i].comp_coords()[0], ilc[i].comp_coords()[0]);
        EXPECT_DOUBLE_EQ(il[i].comp_coords()[1], ilc[i].comp_coords()[1]);
        EXPECT_DOUBLE_EQ(il[i].bulk_coords()[0], ilc[i].bulk_coords()[0]);
        EXPECT_DOUBLE_EQ(il[i].bulk_coords()[1], ilc[i].bulk_coords()[1]);
    }
}


TEST(intersections_22d, all) {   
    // directory with testing meshes
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/simple_meshes_22d/";
    std::vector<string> filenames;
    
    read_files_form_dir(dir_name, "msh", filenames);
    
    std::vector<computeintersection::IntersectionLocal<2,2>> solution;
    fill_22d_solution(solution);
    
    // for each mesh, compute intersection area and compare with old NGH
    for(unsigned int s=0; s< filenames.size(); s++)
    {
//         const unsigned int np = permutations_triangle.size();
//         for(unsigned int p=0; p<np; p++)
//         {
            xprintf(Msg,"Computing intersection on mesh: %s\n",filenames[s].c_str());
            FilePath::set_io_dirs(".","","",".");
            FilePath mesh_file(dir_name + filenames[s], FilePath::input_file);
            
            Mesh mesh;
            // read mesh with gmshreader
            GmshMeshReader reader(mesh_file);
            reader.read_mesh(&mesh);
        
            // permute nodes of one triangle:
//             permute_triangle(mesh.element(0),p);
            
            mesh.setup_topology();
            
            compute_intersection_22d(&mesh, solution[s]);//permute_coords(solution[s], permutations_triangle[p]));
//         }
    }
}