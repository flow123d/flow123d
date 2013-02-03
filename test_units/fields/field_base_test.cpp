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


TEST(Field, init_from_default) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Mesh mesh;
    GmshMeshReader( FilePath("mesh/simplest_cube.msh", FilePath::input_file) ).read_mesh(&mesh);
    Point<3> p("1 2 3");

    {
        Field<3, FieldValue<3>::Scalar > scalar_field;

        // test default initialization of scalar field
        scalar_field.set_default( Input::Type::Default("45") );
        scalar_field.set_mesh(&mesh);
        scalar_field.set_time(0.0);

        EXPECT_EQ( 45.0, scalar_field.value(p, mesh.element_accessor(0)) );
        EXPECT_EQ( 45.0, scalar_field.value(p, mesh.element_accessor(6)) );
        EXPECT_DEATH( { scalar_field.value(p, mesh.element_accessor(0,true)); }, "Null field ptr on region id: 101, field:" );
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

