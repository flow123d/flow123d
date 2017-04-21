
/*
 *
 *      Author: PE
 */
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

#include "system/global_defs.h"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include <mesh/accessors.hh>
#include "mesh_constructor.hh"

#include "intersection/mixed_mesh_intersections.hh"
#include "intersection/intersection_point_aux.hh"
#include "intersection/intersection_aux.hh"
#include "intersection/intersection_local.hh"

#include <dirent.h>

using namespace std;


void compute_intersection(Mesh *mesh)
{
    // compute intersection
//     MixedMeshIntersections ie(mesh);
//     ie.compute_intersections(IntersectionType(IntersectionType::d23
//                                             | IntersectionType::d22));

    mesh->mixed_intersections();
    
    double total_length = mesh->mixed_intersections().measure_22();
    cout << "total_length = " << total_length << endl;
    
    // write computed intersections
//     for(unsigned int i = 0; i < ie.intersection_storage12_.size(); i++)
//     {
//         cout << ie.intersection_storage12_[i];
//     }
    
//     ElementFullIter t1 = mesh->element.begin();
//     for(uint i=0; i<4; i++){
//         Edge * e = t1->side(i)->edge();
//         for(uint s=0; s < e->n_sides; s++)
//             DBGCOUT(<< "side " << i << "  elements over edge[" << e->side(s)->edge_idx() << "]:  "
//                     << e->side(s)->element()->index() << "\n");
//     }
//     
//     DBGVAR(mesh->n_vb_neighbours());
//     ++t1;
//     ++t1;
//     unsigned int n_neighs = t1->n_neighs_vb;
//     DBGVAR(n_neighs);
//         for (unsigned int i = 0; i < n_neighs; i++) {
//             Neighbour *ngh = t1->neigh_vb[i];
//             DBGCOUT(<< "   ngh: " << ngh->edge_idx() << "\n");
//             ngh->edge_idx();
//         }
//     mesh->element[1];
//     mesh->element[2];
}


TEST(intersection_prolongation_23d, all) {
    
//     // directory with testing meshes
    FilePath::set_dirs(UNIT_TESTS_SRC_DIR,"",".");
    string dir_name = "intersection/";
    
//     string filename = dir_name + "cube_frac_nc23_small.msh";
    string filename = "../tests/00_mesh/square_2x2_frac_nc_small.msh";
    

    MessageOut() << "Computing intersection on mesh: " << filename << "\n";
    FilePath mesh_file(filename, FilePath::input_file);
    
    Mesh *mesh = mesh_constructor();
    ifstream in(string(mesh_file).c_str());
    if (!in.is_open()) {
        string out = "failed to open " + mesh_file.filename() + "\n";
        ASSERT(in.is_open()).error(out);
    }
    else{
        mesh->read_gmsh_from_stream(in);
        compute_intersection(mesh);
    }
}
