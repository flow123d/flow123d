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

#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "mesh/range_wrapper.hh"
#include "io/msh_gmshreader.h"
#include "mesh_constructor.hh"

#include "intersection/compute_intersection.hh"
#include "intersection/intersection_point_aux.hh"
#include <intersection/intersection_aux.hh>
#include <intersection/intersection_local.hh>

#include "compute_intersection_test.hh"

using namespace std;

typedef std::vector<IntersectionPoint<1,2>> TestCaseIPs;
typedef std::pair<std::string, TestCaseIPs> TestCaseResult ;

/// Create results for the meshes in directory 'prolong_meshes_13d'.
void fill_solution(std::vector< TestCaseResult> &c)
{
    // degenerate cases
    // no IP
    c.push_back({ "00_d", {}});

    // 1 IP, A node, T node
    c.push_back({ "01_d", {IntersectionPoint<1,2>({0}, {0, 0})}});
    c.push_back({ "02_d", {IntersectionPoint<1,2>({1}, {0, 0})}});
    c.push_back({ "03_d", {IntersectionPoint<1,2>({1}, {0, 0})}});
    
    // 1 IP, A node, T side
    c.push_back({ "04_d", {IntersectionPoint<1,2>({1}, {0, 0.2})}});
    c.push_back({ "05_d", {IntersectionPoint<1,2>({0}, {0.5, 0.5})}});
    c.push_back({ "06_d", {IntersectionPoint<1,2>({1}, {0, 0.2})}});
    
    // 1 IP, inside A, T node
    c.push_back({ "07_d", {IntersectionPoint<1,2>({2./3}, {0, 1})}});
    c.push_back({ "08_d", {IntersectionPoint<1,2>({1./3}, {0, 0})}});
    
    
    // degenerate cases with 2 IPs
    c.push_back({ "10_d", {IntersectionPoint<1,2>({0}, {0, 0}),
                           IntersectionPoint<1,2>({1}, {1, 0})}});
    c.push_back({ "11_d", {IntersectionPoint<1,2>({1./3}, {0.2, 0}),
                           IntersectionPoint<1,2>({2./3}, {0, 0.4})}});
    c.push_back({ "12_d", {IntersectionPoint<1,2>({0}, {0.5, 0}),
                           IntersectionPoint<1,2>({2./3}, {1, 0})}});
    c.push_back({ "13_d", {IntersectionPoint<1,2>({0}, {0, 0.3}),
                           IntersectionPoint<1,2>({1}, {0.5, 0.3})}});
    c.push_back({ "14_d", {IntersectionPoint<1,2>({0}, {0.2, 0.1}),
                           IntersectionPoint<1,2>({1}, {0.4, 0.2})}});
    c.push_back({ "15_d", {IntersectionPoint<1,2>({0}, {0, 0}),
                           IntersectionPoint<1,2>({2./3}, {0.8, 0.2})}});
    c.push_back({ "16_d", {IntersectionPoint<1,2>({0}, {0, 0}),
                           IntersectionPoint<1,2>({2./3}, {0.8, 0.2})}});
    
    // in 3D ambient space
    // special cases
    c.push_back({ "30_s", {IntersectionPoint<1,2>({0.5}, {0, 0})}});
    c.push_back({ "31_s", {IntersectionPoint<1,2>({0.5}, {0.5, 0})}});
    c.push_back({ "32_s", {}});
    
    // regular cases
    c.push_back({ "50_r", {}});
    c.push_back({ "51_r", {IntersectionPoint<1,2>({0.5}, {0.25, 0.25})}});
}


/// Permutes triangle coordinates of IP<1,2> according to given permutation.
std::vector<IntersectionPoint<1,2>> permute_coords(TestCaseIPs ips, 
                                                   const std::vector<unsigned int> &permute)
{
    TestCaseIPs new_points(ips.size());
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
        
        new_points[i] = IntersectionPoint<1,2>(ips[i].comp_coords(), new_coords.subvec(1,2));
    }
    return new_points;
}


void compute_intersection_12d(Mesh *mesh, const TestCaseIPs &ips, bool degenerate)
{
    IntersectionAux<1,2> is(1, 0);
    ComputeIntersection<1,2> CI(mesh->element_accessor(1), mesh->element_accessor(0), mesh);
    if(degenerate)
        CI.compute_final_in_plane(is.points());
    else
        CI.compute_final(is.points());
    
    IntersectionLocal<1,2> ilc(is);
    
    auto ipc = ilc.points();
//     DebugOut() << ilc;
    
//     for(IntersectionPoint<1,2> &ip: ipc)
//     {
//         ip.coords(mesh->element(1)).print(DebugOut(),"ip");
//     }

    EXPECT_EQ(ips.size(), ipc.size());
    for(unsigned int i=0; i < ipc.size(); i++)
    {
        MessageOut().fmt("---------- check IP[{}] ----------\n",i);
        EXPECT_ARMA_EQ(ips[i].comp_coords(), ipc[i].comp_coords());
        EXPECT_ARMA_EQ(ips[i].bulk_coords(), ipc[i].bulk_coords());
    }
}


TEST(intersections_12d, all) {
 
    // directory with testing meshes
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/simple_meshes_12d/";

    std::vector<TestCaseResult> solution_coords;
    fill_solution(solution_coords);
    
    
    // for each mesh, compute intersection area and compare with old NGH
    unsigned int i_file=0;
    for(auto &test_case : solution_coords)
    {
        string file_name=test_case.first+"_line_triangle.msh";
        TestCaseIPs &case_ips=test_case.second;
        
        bool degenerate = test_case.first.at(3) == 'd';

        FilePath mesh_file(dir_name + file_name, FilePath::input_file);
        ASSERT(mesh_file.exists())(dir_name+file_name);
        
        string in_mesh_string = "{mesh_file=\"" + (string)mesh_file + "\"}";
        
        const unsigned int np = permutations_triangle.size();
//         const unsigned int np = 1;
        for(unsigned int p=0; p<np; p++){
            MessageOut().fmt("## Computing intersection on mesh #{}: {} \n ## triangle permutation: #{}\n",
                                i_file,  file_name, p);
            
            Mesh *mesh = mesh_constructor(in_mesh_string);
            // read mesh with gmshreader
            auto reader = reader_constructor(in_mesh_string);
            reader->read_raw_mesh(mesh);
            
            // permute nodes:
            for (auto ele : mesh->elements_range()) {
                if(ele->dim() == 2)
                	mesh->permute_triangle(ele.idx(), permutations_triangle[p]);
            }
            mesh->setup_topology();
            
//                 compare_with_ngh(mesh);
            // compute both ways
            if(degenerate)
                compute_intersection_12d(mesh, permute_coords(case_ips, permutations_triangle[p]), false);
            
            compute_intersection_12d(mesh, permute_coords(case_ips, permutations_triangle[p]), degenerate);
        }
        i_file++;
    }
}
