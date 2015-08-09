/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id: interpolation_main.cc 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>

#include "system/system.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/json_to_storage.hh"
#include "system/sys_profiler.hh"
#include "mesh/region.hh"
#include "input/type_output.hh"

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"

#include "fields/field_interpolated_p0.hh"
#include "fields/field_interpolated_p0.impl.hh"



// tests are started from 'build/test_units'
string input = R"CODE(
{   
   scalar={
       TYPE="FieldInterpolatedP0",
       gmsh_file="fields/simplest_cube_3d.msh",
       field_name="scalar"
   },
   scalar_large={
       TYPE="FieldInterpolatedP0",
       gmsh_file="fields/bigger_3d_cube_0.5.msh",
       field_name="scalar"
   },
   vector_fixed={
       TYPE="FieldInterpolatedP0",
       gmsh_file="fields/simplest_cube_3d.msh",
       field_name="vector_fixed"
   },
   vector={
       TYPE="FieldInterpolatedP0",
       gmsh_file="fields/simplest_cube_3d.msh",
       field_name="vector_fixed"
   },
   tensor_fixed={
       TYPE="FieldInterpolatedP0",
       gmsh_file="fields/simplest_cube_3d.msh",
       field_name="tensor_fixed"
   }
}
)CODE";


// simplest cube 123d
string gmsh_mesh = R"CODE(
$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
1
0       101     ".point_boundary"
$EndPhysicalNames
$Nodes
8
1 1 1 1
2 -1 1 1
3 -1 -1 1
4 1 -1 1
5 1 -1 -1
6 -1 -1 -1
7 1 1 -1
8 -1 1 -1
$EndNodes
$Elements
12
1 1 2 37 20 7 3
2 2 2 38 34 6 3 7
3 2 2 38 36 3 1 7
4 1 2 37 20 8 4
5 2 2 38 34 6 4 8
6 2 2 38 36 4 1 8
7 15 2 101 101 3
# ------ interpolation to following 3D elements not implemented
8 4 2 39 40 3 7 2 8
9 4 2 39 40 3 7 8 6
10 4 2 39 42 3 7 6 5
11 4 2 39 42 3 7 5 4
12 4 2 39 42 3 7 4 1
$EndElements
)CODE";


class FieldInterpolatedP0Test : public testing::Test {
public:
    typedef FieldInterpolatedP0<3, FieldValue<3>::Scalar > ScalarField;
    typedef FieldInterpolatedP0<3, FieldValue<3>::Enum > EnumField;
    typedef FieldInterpolatedP0<3, FieldValue<3>::VectorFixed > VecFixField;
    typedef FieldInterpolatedP0<3, FieldValue<3>::Vector > VecField;
    typedef FieldInterpolatedP0<3, FieldValue<2>::TensorFixed > TensorField;
    typedef FieldInterpolatedP0<3, FieldValue<3>::EnumVector > EnumVector;

    virtual void SetUp() {
        // setup FilePath directories
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

        Profiler::initialize();

        //FilePath mesh_file( "mesh/simplest_cube.msh", FilePath::input_file);
        mesh = new Mesh;
        stringstream in(gmsh_mesh.c_str());
        mesh->read_gmsh_from_stream(in);

        Input::Type::Record rec_type = Input::Type::Record("Test","")
            .declare_key("scalar", ScalarField::get_input_type(FieldAlgorithmBase<3, FieldValue<3>::Scalar>::get_input_type(), NULL), Input::Type::Default::obligatory(),"" )
            .declare_key("scalar_large", ScalarField::get_input_type(FieldAlgorithmBase<3, FieldValue<3>::Scalar>::get_input_type(), NULL), Input::Type::Default::obligatory(),"" )
            .declare_key("vector_fixed", VecFixField::get_input_type(FieldAlgorithmBase<3, FieldValue<3>::VectorFixed>::get_input_type(), NULL), Input::Type::Default::obligatory(),"" )
            .declare_key("vector", VecField::get_input_type(FieldAlgorithmBase<3, FieldValue<3>::Vector>::get_input_type(), NULL), Input::Type::Default::obligatory(),"" )
            .declare_key("tensor_fixed", TensorField::get_input_type(FieldAlgorithmBase<3, FieldValue<2>::TensorFixed>::get_input_type(), NULL), Input::Type::Default::obligatory(),"" )
            .close();

        Input::JSONToStorage reader( input, rec_type, Input::FileFormat::format_JSON );
        rec=reader.get_root_interface<Input::Record>();

        test_time[0] = 0.0;
        test_time[1] = 1.0;
    }
    virtual void TearDown() {

    }

    Mesh *mesh;
    Input::Record rec;
    Space<3>::Point point;
    double test_time[2];

};


TEST_F(FieldInterpolatedP0Test, 1d_2d_elements_small) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("scalar"));

    for (unsigned int j=1; j<3; j++) {
    	field.set_time(test_time[j-1]);

    	EXPECT_DOUBLE_EQ( j*0.650, field.value(point, mesh->element_accessor(0)) );
        EXPECT_DOUBLE_EQ( j*0.650, field.value(point, mesh->element_accessor(1)) );
        EXPECT_DOUBLE_EQ( j*0.650, field.value(point, mesh->element_accessor(2)) );
        EXPECT_DOUBLE_EQ( j*0.700, field.value(point, mesh->element_accessor(3)) );
        EXPECT_DOUBLE_EQ( j*0.675, field.value(point, mesh->element_accessor(4)) );
        EXPECT_DOUBLE_EQ( j*0.675, field.value(point, mesh->element_accessor(5)) );
        EXPECT_DOUBLE_EQ( j*0.650, field.value(point, mesh->element_accessor(0, true)) );
    }

}

TEST_F(FieldInterpolatedP0Test, 1d_2d_elements_large) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("scalar_large"));
    field.set_time(0.0);

    //EXPECT_DOUBLE_EQ( 0.650, field.value(point, mesh->element_accessor(0)) );
}

