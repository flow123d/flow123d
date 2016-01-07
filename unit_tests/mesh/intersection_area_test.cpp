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

#include "intersection/inspectelements.h"
#include "intersection/intersectionpoint.h"
#include "intersection/intersectionline.h"

#include <dirent.h>

using namespace std;
using namespace computeintersection;


void compute_intersection_area_23d(Mesh *mesh)
{
    double area1, area2 = 0;

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
    xprintf(Msg, "Computing polygon area by NEW algorithm\n");
    InspectElements ie(mesh);
    ie.compute_intersections<2,3>();
    area1 = ie.polygon_area();

//     ie.print_mesh_to_file("output_intersection");
    
    xprintf(Msg,"Polygon area: (intersections) %.16e,\t(NGH) %.16e\n", area1, area2);
    EXPECT_NEAR(area1, area2, 1e-14);
//     EXPECT_DOUBLE_EQ(area1,area2);
}


TEST(area_intersections, all) {
    Profiler::initialize();
    
    // directory with testing meshes
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/mesh/site/";//notfunctional/
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
    
    std::sort(filenames.begin(), filenames.end(), less<string>());

    // for each mesh, compute intersection area and compare with old NGH
    for(auto &fname : filenames)
    {
        unsigned int permutations[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
        for(unsigned int p=0; p<6; p++)
        {
            xprintf(Msg,"Computing intersection on mesh: %s\n",fname.c_str());
            FilePath::set_io_dirs(".","","",".");
            FilePath mesh_file(dir_name + fname, FilePath::input_file);
            
            Mesh mesh;
            // read mesh with gmshreader
            GmshMeshReader reader(mesh_file);
            reader.read_mesh(&mesh);
            
            // permute nodes:
            FOR_ELEMENTS(&mesh,ele)
            {
                if(ele->dim() == 2)
                {
                    Node* tmp[3];
                    for(unsigned int i=0; i<ele->n_nodes(); i++)
                    {
                        tmp[i] = ele->node[permutations[p][i]];
                    }
                    for(unsigned int i=0; i<ele->n_nodes(); i++)
                    {
                        ele->node[i] = tmp[i];
//                         ele->node[i]->point().print(cout);
                    }
                }
            }
            mesh.setup_topology();
            
            xprintf(Msg, "==============\n");
            compute_intersection_area_23d(&mesh);
            xprintf(Msg, "==============\n");
        }
    }
    Profiler::uninitialize();
}






/******************************************************************************************** TEST 1d-3d ****/

/// Create results for the meshes in directory 'site_13d'.
void fill_13d_solution(std::vector<computeintersection::IntersectionLine> &ils)
{
    ils.clear();
    ils.resize(12);
    // ils[0] is empty
    ils[1].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({1,0}),arma::vec::fixed<4>({1,0,0,0})));
    ils[1].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({1-1.0/3,1.0/3}),arma::vec::fixed<4>({0,1,1,1})/3));
    // only one IP
    ils[2].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({1,0}),arma::vec::fixed<4>({1,0,0,0})));
    
    ils[3].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({1,0}),arma::vec::fixed<4>({1,0,0,0})));
    ils[3].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({1,2})/3,arma::vec::fixed<4>({0,0,1,0})));
    
    ils[4].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({3,1})/4,arma::vec::fixed<4>({0.5,0,0.5,0})));
    ils[4].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({1,1})/2,arma::vec::fixed<4>({0,0.5,0.5,0})));
    
    ils[5].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({2,1})/3,arma::vec::fixed<4>({0.5,0,0.5,0})));
    ils[5].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({1,2})/3,arma::vec::fixed<4>({0.5,0.5,0,0})));
    // only one IP
    ils[6].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({1,0}),arma::vec::fixed<4>({0.5,0.25,0,0.25})));
    
    ils[7].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({1,0}),arma::vec::fixed<4>({1,1,1,1})/4));
    ils[7].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({1,1})/2,arma::vec::fixed<4>({4,3,0,3})*0.1));
    
    ils[8].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({1,0}),arma::vec::fixed<4>({2,1,1,1})*0.2));
    ils[8].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({0,1}),arma::vec::fixed<4>({1./3,1,1,1})*0.3));
    
    // ils[9] is empty
    // only one IP
    ils[10].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({1,1})/2,arma::vec::fixed<4>({1,1,0,0})/2));
    
    ils[11].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({3,1})/4,arma::vec::fixed<4>({1,2,1,0})/4));
    ils[11].points().push_back(computeintersection::IntersectionPoint<1,3>(arma::vec::fixed<2>({1,1})/2,arma::vec::fixed<4>({1,2,0,1})/4));
}

/**
 * Permutes tetrahedron coordinates of IP<1,3> according to given permutation.
 */
computeintersection::IntersectionLine permute_coords(computeintersection::IntersectionLine il, unsigned int permute[4])
{
    computeintersection::IntersectionLine new_il = il;
    std::vector<computeintersection::IntersectionPoint<1,3>> & points = il.points();
    for(unsigned int i = 0; i < points.size(); i++)
    {
        arma::vec::fixed<4> new_coords;
        for(unsigned int j = 0; j < 4; j++)
            new_coords[j] = points[i].local_bcoords_B()[permute[j]];
        
        new_il.points()[i].set_coordinates(points[i].local_bcoords_A(), new_coords);
    }
    return new_il;
}

void compute_intersection_area_13d(Mesh *mesh, const computeintersection::IntersectionLine &il)
{
    double length1, length2 = 0;

    // compute intersection
    xprintf(Msg, "Computing intersection length by NEW algorithm\n");
    InspectElements ie(mesh);
    ie.compute_intersections<1,3>();
    
    //test solution
    std::vector<computeintersection::IntersectionLine> pp = ie.list_intersection_lines(1);
    computeintersection::IntersectionLine ilc;
    // component = element index == 1
    if(pp.size() > 0)
    {
        ilc = pp[0];
        EXPECT_EQ(ilc.size(), il.size());
    }
    
    for(unsigned int i=0; i < ilc.size(); i++)
    {
        cout << "---------- check IP[" << i << "] ----------" << endl;
        EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_A()[0], il[i].local_bcoords_A()[0]);
        EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_A()[1], il[i].local_bcoords_A()[1]);
        EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_B()[0], il[i].local_bcoords_B()[0]);
        EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_B()[1], il[i].local_bcoords_B()[1]);
        EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_B()[2], il[i].local_bcoords_B()[2]);
        EXPECT_DOUBLE_EQ(ilc[i].local_bcoords_B()[3], il[i].local_bcoords_B()[3]);
    }
    
    length1 = ie.line_length();

    //TODO: delete comparison with NGH
    // compute intersection by NGH
    xprintf(Msg, "Computing intersection length by NGH algorithm\n");
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
    
    xprintf(Msg,"Length of intersection line: (intersections) %.16e,\t(NGH) %.16e\n", length1, length2);
//     EXPECT_NEAR(length1, length2, 1e-12);
    EXPECT_DOUBLE_EQ(length1,length2);
}


TEST(area_intersections_13d, all) {
    Profiler::initialize();
    
    // directory with testing meshes
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/mesh/site_13d/";
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
    
    std::sort(filenames.begin(), filenames.end(), less<string>());
    
    std::vector<IntersectionLine> solution;
    fill_13d_solution(solution);
    
    // for each mesh, compute intersection area and compare with old NGH
    for(unsigned int s=0; s< filenames.size(); s++)
    {
        const unsigned int np = 24;
        unsigned int permutations[np][4] = {{0,1,2,3},
                                                {0,1,3,2},  // the tab means permutation with negative jacobian
                                            {0,3,1,2},
                                                {0,3,2,1},
                                            {0,2,3,1},
                                                {0,2,1,3},
                                                {1,0,2,3},
                                            {1,0,3,2},
                                                {1,3,0,2},
                                            {1,3,2,0},
                                                {1,2,3,0},
                                            {1,2,0,3},
                                                {2,1,0,3},
                                            {2,1,3,0},
                                                {2,3,1,0},
                                            {2,3,0,1},
                                                {2,0,3,1},
                                            {2,0,1,3},
                                                {3,1,2,0},
                                            {3,1,0,2},
                                                {3,0,1,2},
                                            {3,0,2,1},
                                                {3,2,0,1},
                                            {3,2,1,0}};
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
                {
                    Node* tmp[4];
                    for(unsigned int i=0; i<ele->n_nodes(); i++)
                    {
                        tmp[i] = ele->node[permutations[p][i]];
                    }
                    for(unsigned int i=0; i<ele->n_nodes(); i++)
                    {
                        ele->node[i] = tmp[i];
//                         ele->node[i]->point().print(cout);
                    }
//                     cout << p << ": jac = "  << ele->tetrahedron_jacobian() << endl;
                }
            }
            
            mesh.setup_topology();
            
            xprintf(Msg, "==============\n");
            compute_intersection_area_13d(&mesh, permute_coords(solution[s], permutations[p]));
            xprintf(Msg, "==============\n");
        }
    }
    Profiler::uninitialize();
}



