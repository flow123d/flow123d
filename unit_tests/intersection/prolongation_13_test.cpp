
/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: VF, PE
 */
#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include <armadillo>

#include "system/global_defs.h"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "mesh_constructor.hh"

#include "../../src/intersection/mixed_mesh_intersections.hh"
#include "intersection/intersection_point_aux.hh"
#include "intersection/intersection_aux.hh"
#include "intersection/intersection_local.hh"

#include "compute_intersection_test.hh"
#include "arma_expect.hh"


using namespace std;

/// Create results for the meshes in directory 'prolong_meshes_13d'.
void fill_13d_solution(std::vector<std::vector<std::vector<arma::vec3>>> &ils,
                       std::vector<double> &lengths
                      )
{
    unsigned int n_files=7;
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
    lengths[5] = 2*0.25*sqrt(3);
    lengths[6] = 2*0.25*sqrt(3) + 2*0.1*sqrt(9+1+1);
    
    ils[0].resize(9);
    ils[0][0] = {arma::vec3({1,1,1})/4,arma::vec3({1,0,1})*0.3};
    ils[0][1] = {arma::vec3({1,0,1})*0.3,arma::vec3({35,-25,35})*0.01};
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
    ils[2][8] = {arma::vec3({0,2,5})*0.1,arma::vec3({1,0,3.5})*0.1};
    ils[2][9] = {arma::vec3({1,0,3.5})*0.1,arma::vec3({1,-1,1})*0.2};
    ils[2][10] = {arma::vec3({1./3,0,2./3})};
    ils[2][11] = {arma::vec3({0,-5,5})*0.1,arma::vec3({1./3,0,2./3})};
    
    ils[3].resize(4);
    ils[3][0] = {arma::vec3({0.25,0.25,0.25}),arma::vec3({0.25,0.25,0})};
    ils[3][1] = {arma::vec3({0.25,0.25,-1}),arma::vec3({0.25,0.25,-1.25})};
    ils[3][2] = {arma::vec3({0.25,0.5,0.25}),arma::vec3({0.25,0.5,0})};
    ils[3][3] = {arma::vec3({0.25,0.5,-1}),arma::vec3({0.25,0.5,-1.25})};
    
    ils[4].resize(6);
    ils[4][0] = {arma::vec3({0.4,0,0.25}),arma::vec3({0.4,0.25,0.25})};
    ils[4][1] = {arma::vec3({0.4,-0.25,0.25}),arma::vec3({0.4,0,0.25})};
    
    ils[4][3] = {arma::vec3({0,-0.25,0.25}),arma::vec3({0,0,0.25})};
    ils[4][2] = {arma::vec3({0,0,0.25}),arma::vec3({0,0.25,0.25})};
    
    ils[4][4] = {arma::vec3({0.5,0,0.5}),arma::vec3({0,0,0})};    
    ils[4][5] = {arma::vec3({0.5,0,0.5}),arma::vec3({0,0,0})};
    
    ils[5].resize(2);
    ils[5][0] = {arma::vec3({0,0,0}),arma::vec3({0.25,0.25,0.25})};
    ils[5][1] = {arma::vec3({-0.25,-0.25,-0.25}),arma::vec3({0,0,0})};
    
    ils[6].resize(9);
    ils[6][0] = {arma::vec3({0,0,0}),arma::vec3({0.25,0.25,0.25})};
    ils[6][1] = {arma::vec3({0,0,0})};
    ils[6][2] = {arma::vec3({-0.25,-0.25,-0.25}),arma::vec3({0,0,0})};
    
    ils[6][3] = {arma::vec3({0,0,0})};
    ils[6][4] = {arma::vec3({0,0,0})};
    ils[6][5] = {arma::vec3({-0.1,-0.3,-0.1}),arma::vec3({0,0,0})};
    ils[6][6] = {arma::vec3({0,0,0}),arma::vec3({0.3,0.1,0.1})};
    ils[6][7] = {arma::vec3({0,0,0})};
    ils[6][8] = {arma::vec3({0,0,0})};
    

}

/// auxiliary function for sorting intersection storage 13d
bool compare_is13(const IntersectionLocal<1,3>& a,
                  const IntersectionLocal<1,3>& b)
{
    if (a.component_ele_idx() == b.component_ele_idx())
        return a.bulk_ele_idx() <= b.bulk_ele_idx();
    else
        return a.component_ele_idx() < b.component_ele_idx();
}

void compute_intersection_13d(Mesh *mesh, const std::vector<std::vector<arma::vec3>> &il,
                              const double &length)
{
    double computed_length = 0;
    
    MixedMeshIntersections ie(mesh);
    ie.compute_intersections(IntersectionType::d13);
    ie.print_mesh_to_file_13("output_intersection_13");
    
    MessageOut().fmt("N intersections {}\n",ie.intersection_storage13_.size());
    
   // write computed intersections
//     for(unsigned int i = 0; i < ie.intersection_storage13_.size(); i++)
//     {
//         DebugOut() << &ie.intersection_storage13_[i] << ie.intersection_storage13_[i];
//     }
    
//     //write the first intersection
//     for (auto elm : mesh->elements_range()) {
//         
//         if( (elm->dim() == 1) && (ie.intersection_map_[elm.index()].size() > 0) )
//         {
//             IntersectionLocal<1,3>* il13 = 
//                 static_cast<IntersectionLocal<1,3>*> (ie.intersection_map_[elm.index()][0].second);
//             if(il13 != nullptr)
//             {
// //                 DebugOut().fmt("comp idx {}, bulk idx {}, \n",elm.index(),ie.intersection_map_[elm.index()][0].first);
//                 DebugOut() << *il13;
//                 break;
//             }
//         }
//     }
    
    //test solution
    std::vector<IntersectionLocal<1,3>> ilc = ie.intersection_storage13_;
    
    // sort the storage, so it is the same for every algorithm (BIH, BB ...)
    // and we avoid creating the intersection map for exact IPs
    std::sort(ilc.begin(), ilc.end(),compare_is13);
    
    EXPECT_EQ(il.size(), ilc.size());
    
    for(unsigned int i=0; i < ilc.size(); i++)
        for(unsigned int j=0; j < ilc[i].size(); j++)
    {
        MessageOut().fmt("---------- check Intersection[{}] el_1d: {} el_3d: {} ----------\n",i,
                ilc[i].component_ele_idx(), ilc[i].bulk_ele_idx());
//         DebugOut()<< "bary: " << ilc[i][j].comp_coords();
        arma::vec3 ip = ilc[i][j].coords(mesh->element_accessor(ilc[i].component_ele_idx()));
        EXPECT_ARMA_EQ(il[i][j], ip);
        //EXPECT_NEAR(il[i][j][0], ip[0], 1e-14);
        //EXPECT_NEAR(il[i][j][1], ip[1], 1e-14);
        //EXPECT_NEAR(il[i][j][2], ip[2], 1e-14);
    }
    
    computed_length = ie.measure_13();
    
    EXPECT_NEAR(computed_length, length, 1e-14);
}


TEST(intersection_prolongation_13d, all) {
//     // directory with testing meshes
    FilePath::set_dirs(UNIT_TESTS_SRC_DIR, "", ".");
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/prolong_meshes_13d/";
    std::vector<string> filenames;
    
    read_files_from_dir(dir_name, "msh", filenames);
        
    std::vector<std::vector<std::vector<arma::vec3>>> solution;
    std::vector<double> lengths;
    fill_13d_solution(solution, lengths);
    
    // for each mesh, compute intersection area and compare with old NGH
    for(unsigned int s=0; s<filenames.size(); s++)
    {
        /*
        // skip failing prolongation
        if (s == 6) continue;
        */

        MessageOut() << "Computing intersection on mesh: " << filenames[s] << "\n";
        string in_mesh_string = "{mesh_file=\"" + dir_name + filenames[s] + "\"}";
        
        Mesh *mesh = mesh_constructor(in_mesh_string);
        // read mesh with gmshreader
        auto reader = reader_constructor(in_mesh_string);
        reader->read_raw_mesh(mesh);
        
        mesh->setup_topology();
        
        compute_intersection_13d(mesh, solution[s], lengths[s]);
    }
}
