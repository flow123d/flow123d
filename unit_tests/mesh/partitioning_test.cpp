/*
 * partitioning_test.cpp
 *
 *  Created on: Jun 24, 2013
 *      Author: jb
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>

#include "mesh/partitioning.hh"
#include "la/distribution.hh"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"
#include "system/index_types.hh"
#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"


// Test input for mesh
const string mesh_input = R"JSON(
{ 
  mesh_file="mesh/simplest_cube.msh",
  partitioning={
    tool="METIS",
    graph_type="any_neighboring"
  }
}
)JSON";


TEST(Partitioning, all) {
    Profiler::instance();

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Mesh * mesh = mesh_full_constructor(mesh_input);

    const Distribution * init_ds = mesh->get_part()->get_init_distr();

    cout << *init_ds;
    const LongIdx * part = mesh->get_part()->get_loc_part();
    // print partitioning
    for(unsigned int i=0; i < init_ds->lsize(); i++) cout << "proc: " << init_ds->myp() << " i: " << i << " part: " << part[i] << endl;

    vector<LongIdx> old_ids(mesh->n_elements());
    for(unsigned int i=0; i < mesh->n_elements(); i++) old_ids[i] = 2*i;

    Distribution * new_ds;
    int * id_4_loc, * new_4_id;
    mesh->get_part()->id_maps( 2*old_ids.size(), &(old_ids[0]), new_ds, id_4_loc, new_4_id);

    cout << *new_ds;
    for(unsigned int i=0; i < new_ds->lsize(); i++) {
        cout << "proc: " << new_ds->myp() << " loc i: " << i << " id: " << id_4_loc[i] << endl;
        EXPECT_EQ(i+new_ds->begin(), new_4_id[ id_4_loc[i] ] );
    }

    /*
    for(unsigned int i=0; i < old_ids.size(); i++) cout << "proc: " << new_ds->myp() << " id: " << old_ids[i] << " new: " << new_4_id[old_ids[i]] << endl;
    vector<int> &global_part = mesh->get_part()->subdomain_id_field_data().get();
    if (global_part.size() > 1) {
        for(unsigned int i=0; i< old_ids.size(); i++) {
            EXPECT_EQ( new_ds->get_proc( new_4_id[ old_ids[i] ]),  global_part[i] );
        }
    }*/

    delete mesh;
}
