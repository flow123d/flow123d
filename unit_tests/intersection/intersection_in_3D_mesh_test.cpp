
/*
 *
 *      Author: PE
 */
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

#include "system/global_defs.h"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "mesh_constructor.hh"

#include "intersection/mixed_mesh_intersections.hh"
#include "intersection/intersection_point_aux.hh"
#include "intersection/intersection_aux.hh"
#include "intersection/intersection_local.hh"

#include <dirent.h>

#include "arma_expect.hh"

using namespace std;


void compare_ngh_12d(Mesh *mesh)
{
    // compute intersection
    MixedMeshIntersections ie(mesh);
    ie.compute_intersections(IntersectionType(IntersectionType::d12_1));

    
    MixedMeshIntersections ie_ngh(mesh);
    ie_ngh.compute_intersections(IntersectionType(IntersectionType::d12_ngh));
    
    auto ngh = ie_ngh.intersection_storage12_;
    auto plc = ie.intersection_storage12_;
    
    uint size = plc.size() >= ngh.size() ? ngh.size() : plc.size();
    
    const double expect_tol = 1e-9;
//     unsigned int skip = 0;
    
    for(uint i=0; i<size; i++){
//         IntersectionLocal<1,2> &a = ngh[i+skip];
        IntersectionLocal<1,2> &a = ngh[i];
        IntersectionLocal<1,2> &b = plc[i];
        
        EXPECT_EQ(a.bulk_ele_idx(), b.bulk_ele_idx());
        EXPECT_EQ(a.component_ele_idx(), b.component_ele_idx());
            
        ElementFullIter ele_b = mesh->element(a.bulk_ele_idx());
        ElementFullIter ele_c = mesh->element(a.component_ele_idx());
        ElementFullIter ele_bb = mesh->element(b.bulk_ele_idx());
        ElementFullIter ele_cc = mesh->element(b.component_ele_idx());
        cout << "#######################################   b " << ele_b->id() << "  c " << ele_c->id()
        << "  b " << ele_bb->id() << "  c " << ele_cc->id()<< endl;
        
//         if(a.bulk_ele_idx() != b.bulk_ele_idx()) {skip++;
//             cout << "####     SKIP NGH   ####" << endl;
//         }
        
        EXPECT_EQ(a.size(), b.size());
        if(a.size() != b.size()){
            cout << "different SIZE!!!"<< endl;
            continue;
        }
            
        EXPECT_NEAR(a.compute_measure(), b.compute_measure(), expect_tol);
        for(uint p=0; p<a.size(); p++){
            double bulk = arma::norm(a[p].bulk_coords() - b[p].bulk_coords(),2);
            double comp = arma::norm(a[p].comp_coords() - b[p].comp_coords(),2);
            EXPECT_NEAR(0.0, bulk, expect_tol);
            EXPECT_NEAR(0.0, comp, expect_tol);
//             EXPECT_ARMA_EQ(a[p].bulk_coords(), b[p].bulk_coords());
//             EXPECT_ARMA_EQ(a[p].comp_coords(), b[p].comp_coords());
        }
    }
    
//     cout << "Skip counter = " << skip << endl;
}
    
void compute_intersection(Mesh *mesh)
{
    mesh->mixed_intersections();
//     
//     double total_length = mesh->mixed_intersections().measure_22();
//     cout << "total_length = " << total_length << endl;
    
    // write computed intersections
//     for(unsigned int i = 0; i < ie.intersection_storage12_.size(); i++)
//     {
//         cout << ie.intersection_storage12_[i];
//     }
}


TEST(intersection_test, all) {
    
//     // directory with testing meshes
    FilePath::set_dirs(UNIT_TESTS_SRC_DIR,"",".");
    string dir_name = "intersection/2d-2d/";
    
//     string filename = dir_name + "cube_2f_incomp.msh";
//     string filename = "../tests/00_mesh/square_2x2_frac_nc.msh";
    string filename = "intersection/mesh7.msh";

    MessageOut() << "Computing intersection on mesh: " << filename << "\n";
    //FilePath mesh_file(filename, FilePath::input_file);
    string in_mesh_string = "{mesh_file=\"" + filename + "\"}";
    
    Mesh *mesh = mesh_full_constructor(in_mesh_string);
    compute_intersection(mesh);
}
