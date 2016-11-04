/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: VF, PE
 */
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

#include <armadillo>

#include "system/global_defs.h"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"

#include "intersection/inspect_elements.hh"
#include "intersection/intersection_point_aux.hh"
#include "intersection/intersection_aux.hh"
#include "intersection/intersection_local.hh"

#include "compute_intersection_test.hh"

using namespace std;
using namespace computeintersection;

/// Create results for the meshes in directory 'prolong_meshes_13d'.
void fill_12d_solution(std::vector<std::vector<arma::vec3>> &ils)
{
    unsigned int n_files=1;
    ils.clear();
    ils.resize(n_files);
    
    ils[0].resize(9);
    ils[0][0] = {2.4, 4.4, 10};
    ils[0][1] = {2.4, 4.4, 4.48};
    ils[0][2] = {2.4, 4.4, 0};
    ils[0][3] = {3.0571428571428583, 8.4, 10};
    ils[0][4] = {5.1842105263157876, 8.4, 5.03684210526316};
    ils[0][5] = {7.3428571428571425, 8.4, 0};
    ils[0][6] = {8.4, 4.4, 10};
    ils[0][7] = {8.4, 4.4, 5.68};
    ils[0][8] = {8.4, 4.4, 0};
}

/// auxiliary function for sorting intersection storage 13d
bool compare_is12(const computeintersection::IntersectionLocal<1,2>& a,
                  const computeintersection::IntersectionLocal<1,2>& b)
{
    if (a.component_ele_idx() == b.component_ele_idx())
        return a.bulk_ele_idx() <= b.bulk_ele_idx();
    else
        return a.component_ele_idx() < b.component_ele_idx();
}

void compute_intersection_12d(Mesh *mesh, const std::vector<arma::vec3> &il)
{
    // compute intersection
    InspectElements ie(mesh);
    ie.compute_intersections(computeintersection::IntersectionType::d12_2);
    //ie.print_mesh_to_file_13("output_intersection_13");
    
    MessageOut().fmt("N intersections {}\n",ie.intersection_storage12_.size());
    
//    // write computed intersections
//     for(unsigned int i = 0; i < ie.intersection_storage12_.size(); i++)
//     {
//         DebugOut() << ie.intersection_storage12_[i];
//     }
//     
//     //write the first intersection
//     FOR_ELEMENTS(mesh, elm){
//         
//         if( (elm->dim() == 1) && (ie.intersection_map_[elm->index()].size() > 0) )
//         {
//             computeintersection::IntersectionLocal<1,2>* il12 = 
//                 static_cast<computeintersection::IntersectionLocal<1,2>*> (ie.intersection_map_[elm->index()][0].second);
//             if(il12 != nullptr)
//             {
// //                 DBGMSG("comp idx %d, bulk idx %d, \n",elm->index(),ie.intersection_map_[elm->index()][0].first);
//                 DebugOut() << *il12;
//                 break;
//             }
//         }
//     }
    
    //test solution
    std::vector<computeintersection::IntersectionLocal<1,2>> ilc = ie.intersection_storage12_;
    
    // sort the storage, so it is the same for every algorithm (BIH, BB ...)
    // and we avoid creating the intersection map for exact IPs
    std::sort(ilc.begin(), ilc.end(),compare_is12);
    
    EXPECT_EQ(ilc.size(), il.size());
    
    for(unsigned int i=0; i < ilc.size(); i++)
    {
        MessageOut().fmt("---------- check IP[{}] ----------\n",i);
        arma::vec3 ip = ilc[i][0].coords(mesh->element(ilc[i].component_ele_idx()));
//         ip.print(DebugOut(),"real ip");
        EXPECT_NEAR(ip[0], il[i][0], 1e-14);
        EXPECT_NEAR(ip[1], il[i][1], 1e-14);
        EXPECT_NEAR(ip[2], il[i][2], 1e-14);
    }
}


TEST(intersection_prolongation_12d, all) {
//     // directory with testing meshes
    FilePath::set_dirs(UNIT_TESTS_SRC_DIR,"",".");
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/prolong_meshes_12d/";
    std::vector<string> filenames;
    
    read_files_from_dir(dir_name, "msh", filenames);
        
    std::vector<std::vector<arma::vec3>> solution;
    fill_12d_solution(solution);
    
    // for each mesh, compute intersection area and compare with old NGH
    for(unsigned int s=0; s<filenames.size(); s++)
    {
        MessageOut() << "Computing intersection on mesh: " << filenames[s] << "\n";
        FilePath mesh_file(dir_name + filenames[s], FilePath::input_file);
        
        Mesh mesh;
        // read mesh with gmshreader
        GmshMeshReader reader(mesh_file);
        reader.read_mesh(&mesh);
        
        mesh.setup_topology();
        
        compute_intersection_12d(&mesh, solution[s]);
    }
}
