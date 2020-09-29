
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
#include "mesh/range_wrapper.hh"
#include "io/msh_gmshreader.h"
#include "mesh_constructor.hh"

#include "intersection/compute_intersection.hh"
#include "intersection/intersection_aux.hh"
#include "intersection/intersection_point_aux.hh"
#include "intersection/intersection_local.hh"

#include "compute_intersection_test.hh"

using namespace std;


/// Create results for the meshes in directory 'simple_meshes_13d'.
void fill_13d_solution(std::vector<IntersectionLocal<1,3>> &ils)
{
    ils.clear();
    ils.resize(12);
    // ils[0] is empty
    ils[1].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({0}),arma::vec3({0,0,0})));
    ils[1].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({1./3}),arma::vec3({1,1,1})/3));
    // only one IP
    ils[2].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({0}),arma::vec3({0,0,0})));
    
    ils[3].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({0}),arma::vec3({0,0,0})));
    ils[3].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({2./3}),arma::vec3({0,1,0})));
    
    ils[4].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({0.25}),arma::vec3({0,0.5,0})));
    ils[4].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({0.5}),arma::vec3({0.5,0.5,0})));
    
    ils[5].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({1./3}),arma::vec3({0,0.5,0})));
    ils[5].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({2./3}),arma::vec3({0.5,0,0})));
    // only one IP
    ils[6].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({0}),arma::vec3({0.25,0,0.25})));
    
    ils[7].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({0}),arma::vec3({1,1,1})/4));
    ils[7].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({0.5}),arma::vec3({3,0,3})*0.1));
    
    ils[8].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({0}),arma::vec3({1,1,1})*0.2));
    ils[8].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({1}),arma::vec3({1,1,1})*0.3));
    
    // ils[9] is empty
    // only one IP
    ils[10].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({0.5}),arma::vec3({1,0,0})/2));
    
    ils[11].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({0.25}),arma::vec3({2,1,0})/4));
    ils[11].points().push_back(IntersectionPoint<1,3>(arma::vec::fixed<1>({0.5}),arma::vec3({2,0,1})/4));
}


//Permutes tetrahedron coordinates of IP<1,3> according to given permutation.
IntersectionLocal<1,3> permute_coords(IntersectionLocal<1,3> il,
                                                           const std::vector<unsigned int> &permute)
{
    IntersectionLocal<1,3> new_il = il;
    std::vector<IntersectionPoint<1,3>> & points = il.points();
    for(unsigned int i = 0; i < points.size(); i++)
    {
        arma::vec4 new_coords;
        arma::vec4 old_coords = {1,1,1,1};
        for(unsigned int j = 0; j < 3; j++){
            old_coords[j+1] = points[i].bulk_coords()[j];
            old_coords[0] = old_coords[0] - points[i].bulk_coords()[j];
        }
        if(old_coords[0] < 1e-15) old_coords[0] = 0;
        
        for(unsigned int j = 0; j < 4; j++)
            new_coords[j] = old_coords[permute[j]];
        
        new_il.points()[i] = IntersectionPoint<1,3>(points[i].comp_coords(), new_coords.subvec(1,3));
    }
    return new_il;
}

void compute_intersection_13d(Mesh *mesh, const IntersectionLocal<1,3> &il)
{
    // compute intersection
    IntersectionAux<1,3> is;
    ComputeIntersection<1,3> CI(mesh->element_accessor(1), mesh->element_accessor(0), mesh);
    CI.init();
    CI.compute(is);
    
//     DebugOut() << is;
//     for(IntersectionPointAux<1,3> &ip: is.points())
//     {
//         ip.coords(mesh->element(0)).print(DebugOut(),"ip");
//     }

    IntersectionLocal<1,3> ilc(is);
    EXPECT_EQ(il.size(), ilc.size());
    
    for(unsigned int i=0; i < ilc.size(); i++)
    {
        MessageOut().fmt("---------- check IP[{}] ----------\n",i);
        EXPECT_ARMA_EQ(il[i].comp_coords(), ilc[i].comp_coords());
        EXPECT_ARMA_EQ(il[i].bulk_coords(), ilc[i].bulk_coords());
    }
}


TEST(intersections_13d, all) {
    
    // directory with testing meshes
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/simple_meshes_13d/";
    std::vector<string> filenames;
    
    read_files_from_dir(dir_name, "msh", filenames);
    
    std::vector<IntersectionLocal<1,3>> solution;
    fill_13d_solution(solution);
    
    // for each mesh, compute intersection area and compare with old NGH
    for(unsigned int s=0; s< filenames.size(); s++)
    {
        const unsigned int np = permutations_tetrahedron.size();
        for(unsigned int p=0; p<np; p++)
        {
            MessageOut().fmt("Computing intersection on mesh #{} permutation #{} :  {}\n", s, p,  filenames[s]);
            string in_mesh_string = "{mesh_file=\"" + dir_name + filenames[s] + "\"}";
            
            Mesh *mesh = mesh_constructor(in_mesh_string);
            // read mesh with gmshreader
            auto reader = reader_constructor(in_mesh_string);
            reader->read_raw_mesh(mesh);
        
            // permute nodes:
            for (auto ele : mesh->elements_range()) {
                if(ele->dim() == 3)
                    mesh->permute_tetrahedron(ele.idx(), permutations_tetrahedron[p]);
            }
            
            mesh->setup_topology();
            
            compute_intersection_13d(mesh, permute_coords(solution[s], permutations_tetrahedron[p]));
        }
    }
}
