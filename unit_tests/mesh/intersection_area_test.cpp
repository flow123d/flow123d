/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: VF
 */
#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/file_path.hh"
#include "mesh/mesh.h"

#include "mesh/ngh/include/point.h"
#include "mesh/ngh/include/intersection.h"

#include "intersection/inspectelements.h"

#include <dirent.h>

using namespace std;
using namespace computeintersection;


void compute_intersection_area(Mesh *mesh)
{
    double area1, area2 = 0;

    // compute intersection
    xprintf(Msg, "Computing polygon area by NEW algorithm\n");
    InspectElements ie(mesh);
    ie.compute_intersections<2,3>();
    area1 = ie.polygonArea();

    // compute intersection by NGH
    xprintf(Msg, "Computing polygon area by NGH algorithm\n");
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
        }
    }
    GetIntersection(ttr, tte, it, area2);
    
    xprintf(Msg,"Polygon area: (intersections) %.16e,\t(NGH) %.16e\n", area1, area2);
//     EXPECT_NEAR(area1, area2, 1e-12);
    EXPECT_DOUBLE_EQ(area1,area2);
}

TEST(area_intersections, all) {
    Profiler::initialize();
    
    // directory with testing meshes
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/mesh/site/";
    std::vector<string> filenames;
    
    // read mesh file names
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (dir_name.c_str())) != NULL) {
        /* print all the files and directories within directory */
        xprintf(Msg,"Testing mesh files: \n");
        while ((ent = readdir (dir)) != NULL) {
            string fname = ent->d_name;
            // test extension ".msh"
            if(fname.size() >= 4)
            {
                string ext = fname.substr(fname.size()-4);
//                 xprintf(Msg,"%s\n",ext.c_str());
                if(ext == ".msh"){
                    filenames.push_back(ent->d_name);
                    xprintf(Msg,"%s\n",ent->d_name);
                }
            }
        }
        closedir (dir);
    } else {
        ASSERT(0,"Could not open directory with testing meshes.");
    }
    
    // for each mesh, compute intersection area and compare with old NGH
    for(auto &fname : filenames)
    {
        xprintf(Msg,"Computing intersection on mesh: %s\n",fname.c_str());
        string mesh_file = dir_name + fname;
        
        Mesh mesh;
        ifstream ifs(mesh_file.c_str());
        mesh.read_gmsh_from_stream(ifs);
        
        xprintf(Msg, "==============\n");
        compute_intersection_area(&mesh);
        xprintf(Msg, "==============\n");
    }

    
    Profiler::uninitialize();
}




