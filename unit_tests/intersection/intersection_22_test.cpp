
/*
 *
 *      Author: PE
 */
#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include "system/global_defs.h"
#include "system/file_path.hh"
#include "system/sys_profiler.hh"
#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "mesh_constructor.hh"

#include "intersection/mixed_mesh_intersections.hh"
#include "intersection/intersection_point_aux.hh"

#include <dirent.h>

using namespace std;

typedef std::pair<unsigned int, double>    TestCaseResult;  // number of components, total length
typedef std::pair<std::string, TestCaseResult>   TestCase;

/// Create results for the meshes in directory 'prolong_meshes_13d'.
void fill_solution(std::vector< TestCase> &c)
{
    c.clear();
    
    // TODO: seems that intersection coumputing is terribly slow
    // This test for all meshes (up to about 1000 elements) runs nearly 20s.

    c.push_back({"cube_2f_comp_coarse", {1,0}});
    c.push_back({"cube_2f_comp_fine", {1,0}});
    //c.push_back({"cube_2f_incomp", {2,1}});
    c.push_back({"cube_2f_incomp_SurfaceComp", {2,1}});
    //c.push_back({"cube_mult_compXincomp", {3,2}});
    //c.push_back({"cube_mult_compXincomp_triangle", {4,6.840032688952172}});
    c.push_back({"cube_mult_compXincomp_2triangles", {5,10.988973338817276}});
}

void compute_intersection(Mesh *mesh, TestCaseResult result)
{

    // compute intersection
    MixedMeshIntersections ie(mesh);
    ie.compute_intersections(IntersectionType(IntersectionType::d23
                                            | IntersectionType::d22));
    
    double total_length = ie.measure_22();
    cout << "total_length = " << setprecision(17) << total_length << endl;
    EXPECT_EQ(result.first, ie.number_of_components(2));
//     EXPECT_DOUBLE_EQ(result.second, total_length);
    EXPECT_NEAR(result.second, total_length, geometry_epsilon*result.second);
}


TEST(intersection_prolongation_23d, all) {
    DebugOut() << "start";    
  
    Profiler::instance();

    // directory with testing meshes
    FilePath::set_dirs(UNIT_TESTS_SRC_DIR,"",".");
    string dir_name = "intersection/2d-2d/";
    
    std::vector<TestCase> solution_coords;
    fill_solution(solution_coords);
    
    // for each mesh, compute intersection area and compare with old NGH
    for(auto &test_case : solution_coords)
    {
        string file_name = test_case.first+".msh";
        TestCaseResult &case_result = test_case.second;

        MessageOut() << "Computing intersection on mesh: " << file_name << "\n";
        
        string in_mesh_string = "{ mesh_file=\"" + dir_name + file_name + "\", optimize_mesh=false }";
        Mesh *mesh = mesh_full_constructor(in_mesh_string);
        
        compute_intersection(mesh, case_result);
    }
    Profiler::uninitialize();
}
