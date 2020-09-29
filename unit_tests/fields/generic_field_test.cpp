
#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>


#include "fields/generic_field.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "io/msh_gmshreader.h"
#include "system/sys_profiler.hh"
#include "fields/field_flag.hh"


TEST(GenericField, all) {
    Profiler::instance();
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Mesh * mesh = mesh_full_constructor("{mesh_file=\"mesh/simplest_cube.msh\"}");

    GenericField<3>::IndexField subdomain;
    subdomain.flags(FieldFlag::input_copy);
    GenericField<3>::IndexField region_id;
    region_id.flags(FieldFlag::input_copy);

    subdomain = GenericField<3>::subdomain(*mesh);
    subdomain.set_time(TimeGovernor().step(), LimitSide::right);
    // TODO: After we have support to read partitioning form the MSH file.

    region_id = GenericField<3>::region_id(*mesh);
    region_id.set_time(TimeGovernor().step(), LimitSide::right);
    for (auto ele : mesh->elements_range())
    	EXPECT_EQ( ele.region().id(),
    			   region_id.value(ele.centre(), ele)
    			   );

    //delete mesh;
}
