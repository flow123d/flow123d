
/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: VF, PE
 */
#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include "system/global_defs.h"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "mesh_constructor.hh"

#include "../../src/intersection/mixed_mesh_intersections.hh"
#include "intersection/intersection_point_aux.hh"
#include "intersection/intersection_aux.hh"
#include "intersection/intersection_local.hh"

#include "arma_expect.hh"
#include "compute_intersection_test.hh"

using namespace std;


/// Create results for the meshes in directory 'prolong_meshes_23d'.
void fill_23d_solution(std::vector<std::vector<std::vector<arma::vec3>>> &ils)
{
    unsigned int n_files=9;
    ils.clear();
    ils.resize(n_files);
    
    ils[0].resize(4);
    ils[0][0] = {   {0.25, 0.25, 0},
                    {0, 0.05, 0.05},
                    {0.2, 0, 0.1},
                    {0.41, 0, 0.13},
                    {0.475, 0.325, 0}};
    ils[0][1] = {   {0.2, 0, 0.1},
                    {0.4, -0.05, 0.15},
                    {0.41, 0, 0.13}};
    ils[0][2] = {   {0, 0.05, 0.05},
                    {0.2, 0, 0.1},
                    {0.04, 0, 0.08}};
    ils[0][3] = {   {0.2, 0, 0.1},
                    {0.4, -0.05, 0.15},
                    {0.2, -0.2 ,0.2},
                    {0.04, 0, 0.08}};
    
    // no prolongation over a gap between two tetrahedra
    ils[1].resize(4);
    ils[1][0] = {   {0.25, 0.25, 0.25},
                    {0.25, 0.25, 0},
                    {0.2, 0.25, 0}};
    ils[1][1] = {   {0.25, 0.25, -1},
                    {0.25, 0.25, -1.25},
                    {0.2, 0.25, -1}};
    ils[1][2] = {   {0.25, 0.5, 0.25},
                    {0.25, 0.5, 0},
                    {0.2, 0.5, 0}};
    ils[1][3] = {   {0.25, 0.5, -1.25},
                    {0.25, 0.5, -1},
                    {0.2, 0.5, -1}};
                    
    // multiple nodes - face
    ils[2].resize(2);
    ils[2][0] = {   {0, 0.25, 0.25},
                    {0, 0.5, 0.25}};
    ils[2][1] = {   {0, 0.25, 0.25},
                    {0, 0.5, 0.25},
                    {0.25, 0.25, 0.25}};

    // node - face
    ils[3].resize(2);
    ils[3][0] = {   {0, 0.25, 0.25}};
    ils[3][1] = {   {0, 0.25, 0.25},
                    {0, 0.375, 0.25},
                    {0.25, 0.25, 0.25}};
                    
    // multiple nodes - edge
    ils[4].resize(2);
    ils[4][0] = {   {0, 0.25, 0},
                    {0, 0.5, 0}};
    ils[4][1] = {   {0, 0.25, 0},
                    {0, 0.5, 0}};
    
    // node - edge
    ils[5].resize(2, {{0, 0.25, 0}});
    
    // node - node
    ils[6].resize(9, {{0, 0, 0}} );
    
    // edge - edge
    ils[7].resize(2);
    ils[7][0] = {   {0, 0.25, 0}};
    ils[7][1] = {   {0, 0.25, 0},
                    {0.25, 0.25, 0},
                    {0.25, 0.25, 0.25},
                    {0, 0.25, 0.25}};
    
    // side - node
    ils[8].resize(3, {{0, 0, 1}} );
}

/// auxiliary function for sorting intersection storage 13d
bool compare_is23(const IntersectionLocal<2,3>& a,
                  const IntersectionLocal<2,3>& b)
{
    if (a.component_ele_idx() == b.component_ele_idx())
        return a.bulk_ele_idx() <= b.bulk_ele_idx();
    else
        return a.component_ele_idx() < b.component_ele_idx();
}

void compute_intersection_23d(Mesh *mesh,
                              const std::vector<std::vector<arma::vec3>> &il)
{
    double area = 0;

    // compute intersection
    MixedMeshIntersections ie(mesh);
    ie.compute_intersections(IntersectionType::d23);
    
    //test solution
    std::vector<IntersectionLocal<2,3>> ilc = ie.intersection_storage23_;
    
    // write computed intersections
    /*
    for(unsigned int i = 0; i < ilc.size(); i++)
    {
        DebugOut() << &ilc[i] << ilc[i];
        for(unsigned int j=0; j < ilc[i].size(); j++)
            DebugOut() << ilc[i][j].coords(mesh->element_accessor(ilc[i].component_ele_idx()));
    }
    */
    
    // sort the storage, so it is the same for every algorithm (BIH, BB ...)
    // and we avoid creating the intersection map for exact IPs
    std::sort(ilc.begin(), ilc.end(),compare_is23);
    
    EXPECT_EQ(ilc.size(), il.size());
    
    for(unsigned int i=0; i < ilc.size(); i++)
        for(unsigned int j=0; j < ilc[i].size(); j++)
    {
        MessageOut().fmt("---------- check Polygon: {} IP: {} ----------\n",i, j);
        arma::vec3 ip = ilc[i][j].coords(mesh->element_accessor(ilc[i].component_ele_idx()));
        EXPECT_ARMA_EQ(il[i][j], ip);
    }
    
    area = ie.measure_23();
    ie.print_mesh_to_file_23("output_intersection_23");
    MessageOut().fmt("Area of intersection line: (intersections) {}\n", area);
//     EXPECT_NEAR(length1, length2, 1e-12);
//     EXPECT_DOUBLE_EQ(length,1.5*std::sqrt(0.27)+0.35+0.2
//                             + std::sqrt(0.065) + 0.1*std::sqrt(18) + std::sqrt(0.105));
}


TEST(intersection_prolongation_23d, all) {
    
//     // directory with testing meshes
    FilePath::set_dirs(UNIT_TESTS_SRC_DIR,"",".");
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/prolong_meshes_23d/";
    std::vector<string> filenames;
    
    read_files_from_dir(dir_name, "msh", filenames);
        
    std::vector<std::vector<std::vector<arma::vec3>>> solution;
    fill_23d_solution(solution);
    
    // for each mesh, compute intersection area and compare with old NGH
    for(unsigned int s=0; s< filenames.size(); s++)
    {
        MessageOut() << "Computing intersection on mesh: " << filenames[s] << "\n";
        string in_mesh_string = "{mesh_file=\"" + dir_name + filenames[s] + "\"}";
        
        Mesh *mesh = mesh_constructor(in_mesh_string);
        // read mesh with gmshreader
        auto reader = reader_constructor(in_mesh_string);
        reader->read_raw_mesh(mesh);
        
        mesh->setup_topology();
        
        compute_intersection_23d(mesh, solution[s]);
    }
}
