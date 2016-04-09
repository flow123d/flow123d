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

static const std::string profiler_file = "speed_test_13d_profiler.log";
static const unsigned int profiler_loop = 1;

//*
// ******************************************************************************************* TEST 1d-3d ****

void compute_intersection_13d(Mesh *mesh)
{
    // compute intersection
    
    InspectElements ie(mesh);
    ie.compute_intersections(computeintersection::IntersectionType::d13);
    ie.print_mesh_to_file_13("output_intersection_13");
    
    DBGMSG("N intersections %d\n",ie.intersection_storage13_.size());
    
    // write computed intersections
    for(unsigned int i = 0; i < ie.intersection_storage13_.size(); i++)
    {
        cout << ie.intersection_storage13_[i];
    }
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

TEST(intersection_prolongation_13d, all) {
    Profiler::initialize();
    
    // directory with testing meshes
    string dir_name = string(UNIT_TESTS_SRC_DIR) + "/intersection/benchmarks/";
    std::vector<string> filenames = read_filenames(dir_name);
        
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
                compute_intersection_13d(&mesh);
            xprintf(Msg, "==============\n");
    }
    std::fstream fs;
    fs.open(profiler_file.c_str(), std::fstream::out | std::fstream::app);
    Profiler::instance()->output(PETSC_COMM_WORLD, fs);
    Profiler::uninitialize();
}

//*/