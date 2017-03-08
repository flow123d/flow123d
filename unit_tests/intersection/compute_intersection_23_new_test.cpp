
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
#include "mesh_constructor.hh"

#include "mesh/ngh/include/point.h"
#include "mesh/ngh/include/intersection.h"

#include "../../src/intersection/mixed_mesh_intersections.hh"
#include "intersection/compute_intersection.hh"
#include "intersection/intersection_point_aux.hh"
#include "intersection/intersection_aux.hh"
#include "intersection/intersection_local.hh"

#include "arma_expect.hh"

#include "compute_intersection_test.hh"

using namespace std;

/**
 * Ideas for testing:
 * - load files, sort alphabetically => analytical results must follow that order
 * - mesh prerequisite: triangle must have index 1 (0), tetrahedron 2 (1)
 * - do all permutations of triangle
 * - do all permutations of tetrahedron (with positive Jacobian)
 * - analytical real coordinates of IPs: due to permutations, we must have a fixed order
 *      therefore we sort real points, see compare_coords()
 * - due to sorting, we have no control of correct tracing
 *      (at the moment the comparison of final area with old NGH sort of checks that)
 * - in future, we could have analytical objects in std::vector<IntersectionAux<2,3>> solution
 *      and check also for topologic position, correct tracing order, orientation of IP, barycentric coords...
 *      - problem here is, we have to always pay attention to given permutations...
 * 
 * - TODO: one test case fails due to some rounding error; look into, if it is harmless or not...
 *         (possibly add tolerance into armadillo expect)
 * - TODO: when finished, remove NGH comparison
 */

typedef std::vector< arma::vec3>    TestCaseIPs;
typedef std::pair<std::string, TestCaseIPs>   TestCaseResult ;

/// Create results for the meshes in directory 'prolong_meshes_13d'.
void fill_solution(std::vector< TestCaseResult> &c)
{
    c.clear();

    // regular cases
    c.push_back({ "31_r", {
                {0.2, 0.5, 0.2},
                {0.2, 0, 0.2},
                {0.4, 0, 0.3}}});

    c.push_back({ "32_r", {
                {0, 0.3, 0.1},
                {0.2, 0.1, 0},
                {0, 0.4, 0}}});
    c.push_back({ "33_r", {
                {0, 0.2, 0},
                {0, 0.5, 0.5},
                {0.2, 0.8, 0}}});
    c.push_back({ "34_r", {
                {0.1, 0.15, 0.3},
                {0.4, 0.05, 0.1},
                {0.2, 0.25, 0.15}}});

    c.push_back({ "40_r", {
                {0.6, 0, 0},
                {0, 0, 0.6},
                {0, 0.2, 0.8},
                {0.8, 0.2, 0}}});
    c.push_back({ "41_r", {
                {0.2, 0.5, 0.2},
                {0, 0.375, 0.125},
                {0, 0, 0.14},     //FIXME: rounding error, should be zero
                {0.4, 0, 0.3}}});
    c.push_back({ "42_r", {
                {0.05, 0.3, 0.2},
                {0.1, 0.2, 0},
                {0.3, 0.1, 0},
                {0.35, 0.15, 0.2}}});
    
    c.push_back({ "50_r", {
                {0.035, 0.57, 0.15},
                {0, 0.4, 0.25},
                {0, 0, 0.475},
                {0.6, 0, 0.4},
                {0.5, 0.2, 0.3}
                }});

    // no intersection
    c.push_back({ "90_n", {
                }});
    c.push_back({ "91_n", {
                }});


    // special cases:
    // s3 touch cases (intersection on boundary of s3)
    c.push_back({ "00_s", {{0, 0, 1}} }); // s2 corner - s3 vertex
    c.push_back({ "01_s", {{0.5, 0, 0.5}} }); // s2 corner - s3 edge

    c.push_back({ "05_s", {{0.2, 0, 0.8}, {0.8,0, 0.2}} }); // s2 plane - s3 edge
    c.push_back({ "06_s", {{0, 0, 0}, {1, 0, 0}}}); // s2 side identical with s3 edge
    c.push_back({ "07_s", {{0, 0, 0}, {0, 0, 1}}}); // s3 edge in s2 plane
    c.push_back({ "08_s", {{0, 0, 0}, {0, 0, 0.5}}}); // s3 edge in s2 plane
    c.push_back({ "09_s", {{0.5, 0.5, 0}, {1, 0, 0}}}); // s3 edge in s2 plane
    c.push_back({ "10_s", {{1, 0, 0}, {0.5, 0.5, 0} }}); // s3 edge in s2 plane

    c.push_back({ "11_s", { // s2 in face of s3
                {0, 0.5, 0},
                {0.5, 0, 0},
                {0.8, 0, 0},
                {0.8, 0.2, 0},
                {0.3, 0.7, 0},
                {0, 0.7, 0}
                }});

    c.push_back({ "12_s", { // s2 identical with face of s3
                {0, 0, 0},
                {1, 0, 0},
                {0, 1, 0}}});
    c.push_back({ "13_s", { // s2 side on s3 edge; s2 corner on s3 edge
                {0.2, 0, 0},
                {0.2, 0.8, 0},
                {0.8, 0.2, 0}}});
    c.push_back({ "14_s", { // s2 side on s3 edge; s2 corner on s3 edge
                {0, 0, 0},
                {0.5, 0.5, 0},
                {0.4, 0.6, 0},
                {0.2, 0.6, 0},
                {0, 0.2, 0}}});
    
    // special polygon, on S3 boundary

    c.push_back({ "22_s", {
                {0.2, 0.2, 0.4},
                {0, 0, 0.4},
                {0.2, 0, 0.4}
                }});
    c.push_back({ "24_s", { // s2 corners on : s3 vertex, edge , face
                {0.5, 0, 0.125},
                {0.5, 0.25, 0},
                {1.0, 0, 0}}});
    c.push_back({ "25_s", { // s2 corners on : s3 vertex, edge , face
                {0.5, 0.5, 0},
                {0, 0, 0.5},
                {0, 0, 1}}});
    c.push_back({ "26_s", { // s2 corners on : s3 vertex, edge , face; two s2 sides in s3 faces
                {0.5, 0.5, 0},
                {0, 0.25, 0.25},
                {0, 0, 1}}});

}

/// auxiliary function for sorting intersection point according to x,y local coordinates of component triangle
bool compare_coords(const arma::vec3& a,
                    const arma::vec3& b)
{
    double eps = std::numeric_limits<double>::epsilon();
    if(std::fabs(a[0] - b[0]) <= eps){
        if(std::fabs(a[1] - b[1]) < eps)
            return a[2] <= b[2];
        else
            return a[1] < b[1];
    }
    else
        return a[0] < b[0];
}

void compute_intersection_23d(Mesh *mesh, const std::vector<arma::vec3> &il){
    
    // fixed element indices in the tests
    unsigned int triangle_ele_idx = 0,
                 tetra_ele_idx = 1;
    Simplex<2> triangle = create_simplex<2>(mesh->element(triangle_ele_idx));
    Simplex<3> tetra = create_simplex<3>(mesh->element(tetra_ele_idx));
    
    IntersectionAux<2,3> is;
    ComputeIntersection< Simplex<2>, Simplex<3>> CI(triangle, tetra);
    CI.init();
    CI.compute(is);
    
//     cout << is;
//    for(IntersectionPointAux<2,3> &ip: is.points())
//    {
//        ip.coords(mesh->element(0)).print(std::cout,"ip");
//    }

    EXPECT_EQ(il.size(), is.size());
    
    std::vector<arma::vec3> coords(is.size());
    IntersectionLocal<2,3> temp_ilc(is);
    //     std::cout << temp_ilc;
    for(unsigned int i=0; i < is.size(); i++)
    {
        coords[i] = temp_ilc[i].coords(mesh->element(triangle_ele_idx));
    }
    // sort computed coords according to real coordinates
    //std::sort(coords.begin(), coords.end(),compare_coords);
    
    
    for(unsigned int i=0; i < coords.size(); i++)
    {
        std::cout << "---------- check IP[" << i << "] ----------\n";
        
        EXPECT_ARMA_EQ(il[i], coords[i]);
//         EXPECT_ARMA_EQ(il[i].comp_coords(), ilc[i].comp_coords());
//         EXPECT_ARMA_EQ(il[i].bulk_coords(), ilc[i].bulk_coords());
    }
}

void compare_with_ngh(Mesh *mesh)
{
    double area1, area2 = 0;

    // compute intersection by NGH
    MessageOut() << "Computing polygon area by NGH algorithm\n";
    ngh::TTriangle ttr;
    ngh::TTetrahedron tte;
    ngh::TIntersectionType it = ngh::area;

    FOR_ELEMENTS(mesh, elm) {
        if (elm->dim() == 2) {
        ttr.SetPoints(
                ngh::TPoint(elm->node[0]->point()(0), elm->node[0]->point()(1), elm->node[0]->point()(2)),
                ngh::TPoint(elm->node[1]->point()(0), elm->node[1]->point()(1), elm->node[1]->point()(2)),
                ngh::TPoint(elm->node[2]->point()(0), elm->node[2]->point()(1), elm->node[2]->point()(2)) );
        }else if(elm->dim() == 3){
        tte.SetPoints(
                ngh::TPoint(elm->node[0]->point()(0), elm->node[0]->point()(1), elm->node[0]->point()(2)),
                ngh::TPoint(elm->node[1]->point()(0), elm->node[1]->point()(1), elm->node[1]->point()(2)),
                ngh::TPoint(elm->node[2]->point()(0), elm->node[2]->point()(1), elm->node[2]->point()(2)),
                ngh::TPoint(elm->node[3]->point()(0), elm->node[3]->point()(1), elm->node[3]->point()(2)));
        }
    }
    ngh::GetIntersection(ttr, tte, it, area2);
    
    
    // compute intersection
    MessageOut() << "Computing polygon area by NEW algorithm\n";
    MixedMeshIntersections ie(mesh);
    ie.compute_intersections(IntersectionType::d23);
    area1 = ie.measure_23();

//     ie.print_mesh_to_file_23("output_intersection_23");
    
    MessageOut().fmt("Polygon area: (intersections) {},\t(NGH) {}\n", area1, area2);
    EXPECT_NEAR(area1, area2, 1e-14);
//     EXPECT_DOUBLE_EQ(area1,area2);
}


TEST(area_intersections, all) {
    // directory with testing meshes
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/simple_meshes_23d_new/";
    std::vector<string> filenames;
    
    read_files_from_dir(dir_name, "msh", filenames);
    
    //std::vector<IntersectionAux<2,3>> solution;
    std::vector<TestCaseResult> solution_coords;
    fill_solution(solution_coords);
    //EXPECT_EQ(solution_coords.size(), filenames.size());
    
    // for each mesh, compute intersection area and compare with old NGH
    unsigned int i_file=0;
    for(auto &test_case : solution_coords)
    {
        string file_name=test_case.first+"_triangle_tetrahedron.msh";
        TestCaseIPs &case_ips=test_case.second;

        FilePath mesh_file(dir_name + file_name, FilePath::input_file);
        
        const unsigned int np = 1;//permutations_triangle.size();
        for(unsigned int p=0; p<np; p++){
//             if (p>0) break; //FIXME dont forget to remove later
            
            const unsigned int npt = 1;//permutations_tetrahedron.size();
            for(unsigned int pt=0; pt<npt; pt++){
//                 if (pt>0) break; //FIXME dont forget to remove later
                MessageOut().fmt("## Computing intersection on mesh #{}: {} \n ## permutation:  triangle #{}, tetrahedron #{}\n", i_file,  file_name, p, pt);
                
                Mesh *mesh = mesh_constructor();
                // read mesh with gmshreader
                GmshMeshReader reader(mesh_file);
                reader.read_mesh(mesh);
                
                // permute nodes:
                FOR_ELEMENTS(mesh,ele)
                {
                    if(ele->dim() == 2)
                        permute_triangle(ele,p);
                    if(ele->dim() == 3)
                        permute_tetrahedron(ele,pt);
                }
                mesh->setup_topology();
                
                compare_with_ngh(mesh);
                compute_intersection_23d(mesh, case_ips);
            }
        }
        i_file++;
    }
}
