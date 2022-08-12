/*
 * raw_mesh_test.cpp
 *
 *  Created on: Jan 26, 2013
 *      Author: jb
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>   // we need get_record_accessor method

#include "mesh/raw_mesh.hh"
#include "mesh/mesh.h"
#include "input/accessors.hh"
#include "system/file_path.hh"
#include "system/sys_profiler.hh"


TEST(RawMesh, construct) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Profiler::instance();

	Input::Record in_rec = get_record_accessor("{mesh_file=\"mesh/simplest_cube.msh\"}", Input::FileFormat::format_JSON);
	Mesh* source_mesh = BaseMeshReader::mesh_factory( in_rec, true );
	RawMesh *mesh = new RawMesh(source_mesh);

    delete mesh;
	delete source_mesh;
}
