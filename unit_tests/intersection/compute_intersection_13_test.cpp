
/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: VF, PE
 */
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

#include "system/system.hh"
#include "system/sys_profiler.hh"
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

static const std::string profiler_file = "compute_intersection_13d_profiler.log";
static const unsigned int profiler_loop = 1;

/// Create results for the meshes in directory 'simple_meshes_13d'.
void fill_13d_solution(std::vector<computeintersection::IntersectionLocal<1,3>> &ils)
{
    DBGMSG("fill solution\n");
    ils.clear();
    ils.resize(12);
    // ils[0] is empty
    ils[1].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({0}),arma::vec3({0,0,0})));
    ils[1].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({1./3}),arma::vec3({1,1,1})/3));
    // only one IP
    ils[2].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({0}),arma::vec3({0,0,0})));
    
    ils[3].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({0}),arma::vec3({0,0,0})));
    ils[3].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({2./3}),arma::vec3({0,1,0})));
    
    ils[4].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({0.25}),arma::vec3({0,0.5,0})));
    ils[4].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({0.5}),arma::vec3({0.5,0.5,0})));
    
    ils[5].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({1./3}),arma::vec3({0,0.5,0})));
    ils[5].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({2./3}),arma::vec3({0.5,0,0})));
    // only one IP
    ils[6].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({0}),arma::vec3({0.25,0,0.25})));
    
    ils[7].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({0}),arma::vec3({1,1,1})/4));
    ils[7].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({0.5}),arma::vec3({3,0,3})*0.1));
    
    ils[8].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({0}),arma::vec3({1,1,1})*0.2));
    ils[8].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({1}),arma::vec3({1,1,1})*0.3));
    
    // ils[9] is empty
    // only one IP
    ils[10].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({0.5}),arma::vec3({1,0,0})/2));
    
    ils[11].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({0.25}),arma::vec3({2,1,0})/4));
    ils[11].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<1>({0.5}),arma::vec3({2,0,1})/4));

    DBGMSG("fill solution\n");
}


//Permutes tetrahedron coordinates of IP<1,3> according to given permutation.
computeintersection::IntersectionLocal<1,3> permute_coords(computeintersection::IntersectionLocal<1,3> il,
                                                           const std::vector<unsigned int> &permute)
{
    computeintersection::IntersectionLocal<1,3> new_il = il;
    std::vector<computeintersection::IntersectionPoint<1,3>> & points = il.points();
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
        
        new_il.points()[i] = computeintersection::IntersectionPoint<1,3>(points[i].comp_coords(), new_coords.subvec(1,3));
    }
    return new_il;
}

void compute_intersection_13d(Mesh *mesh, const computeintersection::IntersectionLocal<1,3> &il)
{
    double length1, length2 = 0;

    // compute intersection
    DBGMSG("Computing intersection length by NEW algorithm\n");
    InspectElements ie(mesh);
    ie.compute_intersections(computeintersection::IntersectionType::d13);
    
    //test solution
    std::vector<computeintersection::IntersectionLocal<1,3>> pp = ie.intersection_storage13_;
    computeintersection::IntersectionLocal<1,3> ilc;
    // component = element index == 1
    if(pp.size() > 0)
    {
        ilc = pp[0];
        EXPECT_EQ(ilc.size(), il.size());
    }
    
    for(unsigned int i=0; i < ilc.size(); i++)
    {
        DBGMSG("---------- check IP[%d] ----------\n",i);
        EXPECT_DOUBLE_EQ(ilc[i].comp_coords()[0], il[i].comp_coords()[0]);
        EXPECT_DOUBLE_EQ(ilc[i].bulk_coords()[0], il[i].bulk_coords()[0]);
        EXPECT_DOUBLE_EQ(ilc[i].bulk_coords()[1], il[i].bulk_coords()[1]);
        EXPECT_DOUBLE_EQ(ilc[i].bulk_coords()[2], il[i].bulk_coords()[2]);
    }
    
    length1 = ie.measure_13();
        ie.print_mesh_to_file_13("output_intersection_13");
    
    //TODO: delete comparison with NGH
    // compute intersection by NGH
    DBGMSG("Computing intersection length by NGH algorithm\n");
    START_TIMER("OLD intersections 1D-3D");
    TAbscissa tabs;
    TTetrahedron tte;
    TIntersectionType it = line;

    FOR_ELEMENTS(mesh, elm) {
        if (elm->dim() == 1) {
        tabs.SetPoints(TPoint(elm->node[0]->point()(0), elm->node[0]->point()(1), elm->node[0]->point()(2)),
                       TPoint(elm->node[1]->point()(0), elm->node[1]->point()(1), elm->node[1]->point()(2)));
        }
        else if(elm->dim() == 3){
        tte.SetPoints(TPoint(elm->node[0]->point()(0), elm->node[0]->point()(1), elm->node[0]->point()(2)),
                     TPoint(elm->node[1]->point()(0), elm->node[1]->point()(1), elm->node[1]->point()(2)),
                     TPoint(elm->node[2]->point()(0), elm->node[2]->point()(1), elm->node[2]->point()(2)),
                     TPoint(elm->node[3]->point()(0), elm->node[3]->point()(1), elm->node[3]->point()(2)));
        }
    }
    GetIntersection(tabs, tte, it, length2); // get only relative length of the intersection to the abscissa
    length2 *= tabs.Length(); 
    END_TIMER("OLD intersections 1D-3D");
    
    DBGMSG("Length of intersection line: (intersections) %.16e,\t(NGH) %.16e\n", length1, length2);
//     EXPECT_NEAR(length1, length2, 1e-12);
    EXPECT_DOUBLE_EQ(length1,length2);
}


TEST(intersections_13d, all) {
    Profiler::initialize();
    
    // directory with testing meshes
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/simple_meshes_13d/";
    std::vector<string> filenames;
    
    read_files_form_dir(dir_name, "msh", filenames);
    
    std::vector<computeintersection::IntersectionLocal<1,3>> solution;
    fill_13d_solution(solution);
    
    // for each mesh, compute intersection area and compare with old NGH
    for(unsigned int s=0; s< filenames.size(); s++)
    {
        const unsigned int np = permutations_tetrahedron.size();
        for(unsigned int p=0; p<np; p++)
        {
            xprintf(Msg,"Computing intersection on mesh: %s\n",filenames[s].c_str());
            FilePath::set_io_dirs(".","","",".");
            FilePath mesh_file(dir_name + filenames[s], FilePath::input_file);
            
            Mesh mesh;
            // read mesh with gmshreader
            GmshMeshReader reader(mesh_file);
            reader.read_mesh(&mesh);
        
            // permute nodes:
            FOR_ELEMENTS(&mesh,ele)
            {
                if(ele->dim() == 3)
                    permute_tetrahedron(ele,p);
            }
            
            mesh.setup_topology();
            
            xprintf(Msg, "==============\n");
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
                compute_intersection_13d(&mesh, permute_coords(solution[s], permutations_tetrahedron[p]));
            xprintf(Msg, "==============\n");
        }
    }
    std::fstream fs;
    fs.open(profiler_file.c_str(), std::fstream::out);
    Profiler::instance()->output(PETSC_COMM_WORLD, fs);
    Profiler::uninitialize();
}