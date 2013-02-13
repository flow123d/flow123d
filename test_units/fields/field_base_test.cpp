/*
 * field_base_test.cpp
 *
 *  Created on: Feb 3, 2013
 *      Author: jb
 */



#include <gtest/gtest.h>


#include "fields/field_base.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/json_to_storage.hh"
#include "fields/field_constant.hh"

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"

string input = R"INPUT(
{   
   init_conc={ # formula on 2d 
       TYPE="FieldFormula",
       value=["x", "x*y", "y+t"]
   },
   conductivity_3d={ #3x3 tensor
       TYPE="FieldFormula",
       value=["sin(x)+cos(y)","exp(x)+y^2", "base:=(x+y); base+base^2"]
   }
}
)INPUT";

/* Regions in the test mesh:
 * $PhysicalNames
    6
    1       37      "1D diagonal"
    2       38      "2D XY diagonal"
    2       101     ".top side"
    2       102     ".bottom side"
    3       39      "3D back"
    3       40      "3D front"
    $EndPhysicalNames
 */
TEST(Field, init_from_default) {
    DBGMSG("here\n");
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Mesh mesh;
    ifstream in(string( FilePath("mesh/simplest_cube.msh", FilePath::input_file) ).c_str());
    mesh.read_gmsh_from_stream(in);

    Point<3> p("1 2 3");

    {
        Field<3, FieldValue<3>::Scalar > scalar_field;

        // test default initialization of scalar field
        scalar_field.set_default( Input::Type::Default("45") );
        scalar_field.set_mesh(&mesh);
        scalar_field.set_time(0.0);

        EXPECT_EQ( 45.0, scalar_field.value(p, mesh.element_accessor(0)) );
        EXPECT_EQ( 45.0, scalar_field.value(p, mesh.element_accessor(6)) );
        EXPECT_DEATH( { scalar_field.value(p, mesh.element_accessor(0,true)); }, "Null field ptr " );
    }

    {
        BCField<3, FieldValue<3>::Scalar > scalar_field;

        // test death of set_time without default value
        scalar_field.set_mesh(&mesh);
        EXPECT_DEATH( {scalar_field.set_time(0.0);} , "Missing value of the field");
    }
    //
    {
        BCField<3, FieldValue<3>::Enum > enum_field;
        Input::Type::Selection sel("TestType");
        sel.add_value(0, "none")
           .add_value(1,"dirichlet")
           .close();

        enum_field.set_selection(&sel);
        enum_field.set_default( Input::Type::Default("none") );
        enum_field.set_mesh(&mesh);
        enum_field.set_time(0.0);

        EXPECT_EQ( 0 , enum_field.value(p, mesh.element_accessor(0, true)) );

    }
    Field<3, FieldValue<3>::Vector > vector_field;

}


TEST(Field, no_check) {
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";

    enum {
        dirichlet,
        neumann,
        robin
    };
    // test optional checking in the set_time method
    BCField<3, FieldValue<3>::Enum > bc_type;
    bc_type.set_name("bc_type");

    std::vector<FieldEnum> list;
    BCField<3, FieldValue<3>::Scalar > bc_value;
    bc_value.set_name("bc_value");
    list.clear(); list.push_back(neumann);
    bc_value.disable_where( &bc_type, list );

    BCField<3, FieldValue<3>::Scalar > bc_flux;
    bc_flux.set_name("bc_flux");
    list.clear(); list.push_back(dirichlet); list.push_back(robin);
    bc_flux.disable_where( &bc_type, list );

    BCField<3, FieldValue<3>::Scalar > bc_sigma;
    bc_sigma.set_name("bc_sigma");
    list.clear(); list.push_back(dirichlet); list.push_back(neumann);
    bc_sigma.disable_where( &bc_type, list );

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Mesh mesh;
    ifstream in(string( FilePath("mesh/simplest_cube.msh", FilePath::input_file) ).c_str());
    mesh.read_gmsh_from_stream(in);

    bc_type.set_mesh(&mesh);
    bc_flux.set_mesh(&mesh);
    bc_value.set_mesh(&mesh);
    bc_sigma.set_mesh(&mesh);

    /*
    1       37      "1D diagonal"
    2       38      "2D XY diagonal"
    2       101     ".top side"
    2       102     ".bottom side"
    3       39      "3D back"
    3       40      "3D front"
     */

    typedef FieldConstant<3, FieldValue<3>::Scalar > SConst;
    typedef FieldConstant<3, FieldValue<3>::Enum > EConst;
    auto neumann_type = EConst().set_value(neumann);
    auto robin_type = EConst().set_value(robin);
    auto one = SConst().set_value(1.0);

    bc_type.set_field(mesh.region_db().find_id(101), & neumann_type );
    bc_flux.set_field(mesh.region_db().find_id(101), & one );

    bc_type.set_field(mesh.region_db().find_id(102), & robin_type );
    bc_value.set_field(mesh.region_db().find_id(102), & one );
    bc_sigma.set_field(mesh.region_db().find_id(102), & one );

    bc_type.set_field(mesh.region_db().find_id(-3), & neumann_type );
    bc_flux.set_field(mesh.region_db().find_id(-3), & one );

    bc_type.set_time(0.0);
    bc_flux.set_time(0.0);
    bc_value.set_time(0.0);
    bc_sigma.set_time(0.0);
}
