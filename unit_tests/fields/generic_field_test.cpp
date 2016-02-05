
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>


#include "fields/generic_field.hh"
#include "mesh/mesh.h"
#include "system/sys_profiler.hh"
#include "fields/field_flag.hh"


TEST(GenericField, all) {
    Profiler::initialize();

    FilePath mesh_file( string(UNIT_TESTS_SRC_DIR) + "/mesh/simplest_cube.msh", FilePath::input_file);
    Mesh mesh;
    ifstream in(string(mesh_file).c_str());
    mesh.read_gmsh_from_stream(in);

    GenericField<3>::IndexField subdomain;
    subdomain.flags(FieldFlag::input_copy);
    GenericField<3>::IndexField region_id;
    region_id.flags(FieldFlag::input_copy);

    subdomain = GenericField<3>::subdomain(mesh);
    subdomain.set_time(TimeGovernor().step(), LimitSide::right);
    // TODO: After we have support to read partitioning form the MSH file.

    region_id = GenericField<3>::region_id(mesh);
    region_id.set_time(TimeGovernor().step(), LimitSide::right);
    FOR_ELEMENTS(&mesh, ele)
    	EXPECT_EQ( ele->region().id(),
    			   region_id.value(ele->centre(), ele->element_accessor())
    			   );
}
