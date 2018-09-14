/*
 * range_test.cpp
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest.hh>
#include <mesh_constructor.hh>

#include "mesh/side_impl.hh"
#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "system/sys_profiler.hh"
#include "mesh/range_wrapper.hh"
#include "mesh/accessors.hh"

TEST(RangeWrapper, range) {
    Profiler::initialize();

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	std::string mesh_in_string = "{mesh_file=\"fields/simplest_cube_data.msh\"}";
	Mesh * mesh = mesh_constructor(mesh_in_string);
    auto reader = reader_constructor(mesh_in_string);
    reader->read_physical_names(mesh);
    reader->read_raw_mesh(mesh);
    mesh->check_and_finish();

    unsigned int counter=0;
    std::vector<unsigned int> expected_dims = {1,2,2,3,3,3,3,3,3};
	Range<ElementAccessor<3>, Mesh> range(const_cast<const Mesh *>(mesh), 0, 9);
	for (auto v : range) {
		EXPECT_EQ(v.idx(), counter);
		EXPECT_EQ(v.dim(), expected_dims[counter]);
		counter++;
	}

    delete mesh;
}

