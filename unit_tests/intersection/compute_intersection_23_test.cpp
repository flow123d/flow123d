/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: VF, PE
 */
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

#include "system/global_defs.h"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"

#include "mesh/ngh/include/point.h"
#include "mesh/ngh/include/intersection.h"

#include "intersection/inspect_elements.hh"
#include "intersection/intersection_point_aux.hh"
#include "intersection/intersection_local.hh"

#include "compute_intersection_test.hh"

using namespace std;
using namespace computeintersection;


void compute_intersection_area_23d(Mesh *mesh)
{
    double area1, area2 = 0;

    // compute intersection by NGH
    MessageOut() << "Computing polygon area by NGH algorithm\n";
    TTriangle ttr;
    TTetrahedron tte;
    TIntersectionType it = area;

    FOR_ELEMENTS(mesh, elm) {
        if (elm->dim() == 2) {
        ttr.SetPoints(TPoint(elm->node[0]->point()(0), elm->node[0]->point()(1), elm->node[0]->point()(2)),
                     TPoint(elm->node[1]->point()(0), elm->node[1]->point()(1), elm->node[1]->point()(2)),
                     TPoint(elm->node[2]->point()(0), elm->node[2]->point()(1), elm->node[2]->point()(2)) );
        }else if(elm->dim() == 3){
        tte.SetPoints(TPoint(elm->node[0]->point()(0), elm->node[0]->point()(1), elm->node[0]->point()(2)),
                     TPoint(elm->node[1]->point()(0), elm->node[1]->point()(1), elm->node[1]->point()(2)),
                     TPoint(elm->node[2]->point()(0), elm->node[2]->point()(1), elm->node[2]->point()(2)),
                     TPoint(elm->node[3]->point()(0), elm->node[3]->point()(1), elm->node[3]->point()(2)));
//         elm->node[0]->point().print();
//         elm->node[1]->point().print();
//         elm->node[2]->point().print();
//         elm->node[3]->point().print();
//         
//         elm->side(0)->node(0)->point().print();
//         elm->side(0)->node(1)->point().print();
//         elm->side(0)->node(2)->point().print();
        }
    }
    GetIntersection(ttr, tte, it, area2);
    
    
    // compute intersection
    MessageOut() << "Computing polygon area by NEW algorithm\n";
    InspectElements ie(mesh);
    ie.compute_intersections(computeintersection::IntersectionType::d23);
    area1 = ie.measure_23();

//     ie.print_mesh_to_file_23("output_intersection_23");
    
    MessageOut().fmt("Polygon area: (intersections) {},\t(NGH) {}\n", area1, area2);
    EXPECT_NEAR(area1, area2, 1e-14);
//     EXPECT_DOUBLE_EQ(area1,area2);
}


TEST(area_intersections, all) {
    // directory with testing meshes
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/simple_meshes_23d/";
    std::vector<string> filenames;
    
    read_files_form_dir(dir_name, "msh", filenames);

    // for each mesh, compute intersection area and compare with old NGH
    for(auto &fname : filenames)
    {
        const unsigned int np = permutations_triangle.size();
        for(unsigned int p=0; p<np; p++)
        {
            MessageOut() << "Computing intersection on mesh: " << fname << "\n";
            FilePath mesh_file(dir_name + fname, FilePath::input_file);
            
            Mesh mesh;
            // read mesh with gmshreader
            GmshMeshReader reader(mesh_file);
            reader.read_mesh(&mesh);
            
            // permute nodes:
            FOR_ELEMENTS(&mesh,ele)
            {
                if(ele->dim() == 2)
                    permute_triangle(ele,p);
            }
            mesh.setup_topology();
            
            compute_intersection_area_23d(&mesh);
        }
    }
}