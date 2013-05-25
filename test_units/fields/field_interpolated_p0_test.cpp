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

#include <gtest/gtest.h>

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
#include "fields/field_interpolated_p0_impl.hh"



// tests are started from 'build/test_units'
string input = R"CODE(
{   
   scalar={
       TYPE="FieldInterpolatedP0",
       gmsh_file="fields/simplest_cube_data.msh",
       field_name="scalar"
   },
   vector_fixed={
       TYPE="FieldInterpolatedP0",
       gmsh_file="fields/simplest_cube_data.msh",
       field_name="vector_fixed"
   },
   vector={
       TYPE="FieldInterpolatedP0",
       gmsh_file="fields/simplest_cube_data.msh",
       field_name="vector_fixed"
   },
   tensor_fixed={
       TYPE="FieldInterpolatedP0",
       gmsh_file="fields/simplest_cube_data.msh",
       field_name="tensor_fixed"
   }
}
)CODE";


// simplest cube 123d
string gmsh_mesh = R"CODE(
$MeshFormat
2.2 0 8
$EndMeshFormat
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
9
1 1 2 37 20 7 3
2 2 2 38 34 6 3 7
3 2 2 38 36 3 1 7
4 4 2 39 40 3 7 1 2
5 4 2 39 40 3 7 2 8
6 4 2 39 40 3 7 8 6
7 4 2 39 42 3 7 6 5
8 4 2 39 42 3 7 5 4
9 4 2 39 42 3 7 4 1
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

        FilePath mesh_file( "mesh/simplest_cube.msh", FilePath::input_file);
        mesh= new Mesh;
        ifstream in(string( mesh_file ).c_str());
        mesh->read_gmsh_from_stream(in);

        Input::Type::Record  rec_type("Test","");
        rec_type.declare_key("scalar", ScalarField::input_type, Input::Type::Default::obligatory(),"" );
        rec_type.declare_key("vector_fixed", VecFixField::input_type, Input::Type::Default::obligatory(),"" );
        rec_type.declare_key("vector", VecField::input_type, Input::Type::Default::obligatory(),"" );
        rec_type.declare_key("tensor_fixed", TensorField::input_type, Input::Type::Default::obligatory(),"" );
        rec_type.finish();

        std::stringstream ss(input);
        Input::JSONToStorage reader;
        reader.read_stream( ss, rec_type );
        rec=reader.get_root_interface<Input::Record>();

    }
    virtual void TearDown() {

    }

    Mesh *mesh;
    Input::Record rec;
    Point<3> point;

};


TEST_F(FieldInterpolatedP0Test, 2d_elements) {
    ScalarField field;
    field.init_from_input(rec.val<Input::Record>("scalar"));
    field.set_mesh(mesh);
    field.set_time(0.0);

    const Point<3> p;
    const ElementAccessor<3> ele(mesh, 2, 0);
    /*Element ele( 2, mesh, RegionIdx() );
    ele.node[0]= new Node(0.01, 0.01, 0.00);
    ele.node[1]= new Node(0.16, 0.16, 0.00);
    ele.node[2]= new Node(0.02, 0.02, 0.05); // */

    field.value(p, ele);
    field.value(p, ele);

    /*
    // setup FilePath directories
    FilePath::set_io_dirs(".","/",UNIT_TESTS_SRC_DIR,".");
    Profiler::initialize();

    // initialize Input:Types
    //FieldBase< 3, FieldValue<3>::Scalar >::get_input_type();

    // read input string
    std::stringstream ss(input.c_str());
    Input::JSONToStorage reader;
    reader.read_stream( ss, FieldInterpolatedP0< 3, FieldValue<3>::Scalar >::input_type );
    Input::Record in_rec=reader.get_root_interface<Input::Record>();

    // load mesh
    Mesh mesh;
    std::stringstream msh_ss(gmsh_mesh);
    mesh.read_gmsh_from_stream(msh_ss);

    FieldInterpolatedP0< 3, FieldValue<3>::Scalar > func;
    func.init_from_input(in_rec);
    Point<3> p;
    Element ele( 2, &mesh, RegionIdx() );
    ele.node[0]= new Node(0.01, 0.01, 0.00);
    ele.node[1]= new Node(0.16, 0.16, 0.00);
    ele.node[2]= new Node(0.02, 0.02, 0.05);
    //func.set_element(&ele);
    //EXPECT_EQ(3.5 , func.value(p));*/
}

/*
int main(int argc, char **argv) {

    TPoint pointA(-0.10, -0.10, 0.00);
    TPoint pointB(1.60, 1.60, 0.00);
    TPoint pointC(0.10, 0.10, 0.50);
    TTriangle triangle(pointA, pointB, pointC);

    /*TPoint* point0 = new TPoint(0.00, 0.00, 0.00);
    TPoint* point1 = new TPoint(3.00, 0.00, 0.00);
    TPoint* point2 = new TPoint(0.00, 3.00, 0.00);
    TPoint* point3 = new TPoint(0.00, 0.00, 3.00);
    TTetrahedron tetrahedron(point0, point1, point2, point3);*/

    /*TPoint* point1 = new TPoint( 1.0, 1.0, 1.0);
    TPoint* point2 = new TPoint( 2.0, 6.0, 3.0);
    TPoint* point3 = new TPoint( 2.0, 2.0, 1.0);
    TPoint* point4 = new TPoint(-1.0, 3.0, 3.0);
    TPoint* point5 = new TPoint(-1.0, 5.0, 4.0);
    TPoint* point6 = new TPoint( 3.0, 4.0, 1.5);
    TPoint* point7 = new TPoint( 0.0, 7.0, 4.5);
    TPolygon* pol = new TPolygon();
    pol->Add(point1);
    pol->Add(point2);
    pol->Add(point3);
    pol->Add(point4);
    pol->Add(point5);
    pol->Add(point6);
    pol->Add(point7);
    pol->Write();
    xprintf(Msg, "Polygon: %f\n", pol->GetArea());*/



