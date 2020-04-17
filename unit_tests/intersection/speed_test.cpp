/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: VF, PE
 */
#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include "system/global_defs.h"

#ifdef FLOW123D_RUN_UNIT_BENCHMARKS

#include <armadillo>

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"

#include "mesh_constructor.hh"

// #include "intersection/intersection_point_aux.hh"
// #include "intersection/intersection_aux.hh"
#include "intersection/intersection_local.hh"
#include "intersection/mixed_mesh_intersections.hh"

#include "compute_intersection_test.hh"

#include <dirent.h>

using namespace std;

static const std::string profiler_file = "speed_test_profiler.log";
static const unsigned int profiler_loop = 10;

// bool compare_il_idx(const ILpair &ilA,
//                     const ILpair &ilB)
// {
//     return ilA.first < ilB.first;
// }

void compute_intersection(Mesh *mesh)
{
    // compute intersection
    
    MixedMeshIntersections ie(mesh);
    ie.compute_intersections(IntersectionType::all);
//     ie.compute_intersections((IntersectionType)(IntersectionType::d23 | IntersectionType::d22));
//     ie.print_mesh_to_file_13("output_intersection_speed_13");
//     ie.print_mesh_to_file_23("output_intersection_speed_23");
    
    MessageOut().fmt("N intersections 13({}), 23({}), 22({})\n",ie.intersection_storage13_.size(),
                                                                ie.intersection_storage23_.size(),
                                                                ie.intersection_storage22_.size());
    
    // write computed intersections
//     for(unsigned int i = 0; i < ie.intersection_storage13_.size(); i++)
//     {
//         cout << ie.intersection_storage13_[i];
//     }
//     cout << "----------------------------------------------------------------" << endl;
    
//     for (auto ele : mesh->elements_range()) {
//         if(ele->dim() == 2){
//             int ele_idx = ele.index();
//             auto vec = ie.intersection_map_[ele_idx];
//             std::sort(vec.begin(), vec.end(), compare_il_idx);
//             for(unsigned int i = 0; i < vec.size(); i++)
//                 cout << "ele2d x ele3d:  " << ele_idx << " " << vec[i].first << "\n";
//         }
//     }
    
//     for(unsigned int i = 0; i < ie.intersection_storage23_.size(); i++)
//     {
//         cout << ie.intersection_storage23_[i];
//     }
//     cout << "----------------------------------------------------------------" << endl;
//     for(unsigned int i = 0; i < ie.intersection_storage22_.size(); i++)
//     {
//         cout << ie.intersection_storage22_[i];
//     }
}


/// This test is currently turned off.
/// Its purpose is to compute intersections on large meshes (bedrichov, melechov, cube etc. on bacula).

// TEST(benchmark_meshes, all) {
//     Profiler::instance();
//     
//     // directory with testing meshes
//     string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/benchmarks/";
//     std::vector<string> filenames;
//     read_files_from_dir(dir_name, "msh", filenames, false);
//     
//     
//     if(filenames.size() == 0)
//     {
//         WarningOut().fmt("No benchmark meshes were found in directory: '{}'\n", dir_name.c_str());
//         EXPECT_EQ(1,1);
//         return;
//     }
//     
//     Profiler::instance()->set_task_info("Speed test Inspect Elements Algorithm. "+filenames[0],2);
// 
//     // for each mesh, compute intersection area and compare with old NGH
//     for(unsigned int s=0; s< filenames.size(); s++)
//     {
//             MessageOut().fmt("Computing intersection on mesh: {}\n",filenames[s].c_str());
//             FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
//             string in_mesh_string = "{mesh_file=\"" + dir_name + filenames[s] + "\"}";
//             
//             Mesh *mesh = mesh_constructor(in_mesh_string);
//             // read mesh with gmshreader
//             auto reader = reader_constructor(in_mesh_string);
//             reader->read_raw_mesh(mesh);
//             mesh->setup_topology();
//             
//             MessageOut() << "==============\n";
//             for(unsigned int loop = 0; loop < profiler_loop; loop++)
//                 compute_intersection(mesh);
//             MessageOut() <<  "==============\n";
//     }
//     std::fstream fs;
//     fs.open(profiler_file.c_str(), std::fstream::out /*| std::fstream::app*/);
//     Profiler::instance()->output(PETSC_COMM_WORLD, fs);
//     Profiler::uninitialize();
// }

#endif // FLOW123D_RUN_UNIT_BENCHMARKS
