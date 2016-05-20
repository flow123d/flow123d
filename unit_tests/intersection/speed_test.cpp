/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: VF, PE
 */
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

#include <armadillo>

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"

#include "mesh/ngh/include/point.h"
#include "mesh/ngh/include/intersection.h"

#include "intersection/inspectelements.h"
#include "intersection/intersectionpoint.h"
#include "intersection/intersectionaux.h"
#include "intersection/intersection_local.h"

#include <dirent.h>

using namespace std;
using namespace computeintersection;

static const std::string profiler_file = "speed_test_profiler.log";
static const unsigned int profiler_loop = 10;

void compute_intersection(Mesh *mesh)
{
    // compute intersection
    
    InspectElements ie(mesh);
    ie.compute_intersections(computeintersection::IntersectionType::d23);
//     ie.compute_intersections((computeintersection::IntersectionType)(computeintersection::IntersectionType::d23 | computeintersection::IntersectionType::d22));
//     ie.print_mesh_to_file_13("output_intersection_speed_13");
//     ie.print_mesh_to_file_23("output_intersection_speed_23");
    
    xprintf(Msg,"N intersections 13(%d), 23(%d), 22(%d)\n",ie.intersection_storage13_.size(),
                                                      ie.intersection_storage23_.size(),
                                                      ie.intersection_storage22_.size());
    
    // write computed intersections
//     for(unsigned int i = 0; i < ie.intersection_storage13_.size(); i++)
//     {
//         cout << ie.intersection_storage13_[i];
//     }
//     cout << "----------------------------------------------------------------" << endl;
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


std::vector<string> read_filenames(string dir_name)
{
    std::vector<string> filenames;
    
    // read mesh file names
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (dir_name.c_str())) != NULL) {
        // print all the files and directories within directory 
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
    return filenames;
}

TEST(benchmark_meshes, all) {
    Profiler::initialize();
    
    DBGMSG("tolerances: %e\t%e\n", rounding_epsilon, geometry_epsilon);
    
    // directory with testing meshes
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/benchmarks/";
    std::vector<string> filenames = read_filenames(dir_name);
    
    Profiler::instance()->set_task_info("Speed test Inspect Elements Algorithm. "+filenames[0],2);

    // for each mesh, compute intersection area and compare with old NGH
    for(unsigned int s=0; s< filenames.size(); s++)
    {
            xprintf(Msg,"Computing intersection on mesh: %s\n",filenames[s].c_str());
            FilePath::set_io_dirs(".","","",".");
            FilePath mesh_file(dir_name + filenames[s], FilePath::input_file);
            
            Mesh mesh;
            // read mesh with gmshreader
            GmshMeshReader reader(mesh_file);
            reader.read_mesh(&mesh);
            mesh.setup_topology();
            
            xprintf(Msg, "==============\n");
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
                compute_intersection(&mesh);
            xprintf(Msg, "==============\n");
    }
    std::fstream fs;
    fs.open(profiler_file.c_str(), std::fstream::out /*| std::fstream::app*/);
    Profiler::instance()->output(PETSC_COMM_WORLD, fs);
    Profiler::uninitialize();
}

//*/