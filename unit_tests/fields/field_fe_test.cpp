/*
 * field_elementwise_test.cpp
 *
 *  Created on: Jan 25, 2013
 *      Author: jb
 */



#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>


#include "fields/field_fe.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"

#include "fem/dofhandler.hh"
#include "fem/fe_p.hh"
#include "quadrature/quadrature.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_values.hh"
#include "fem/fe_rt.hh"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"




class FieldFETest : public testing::Test {
public:
    typedef FieldFE<3, FieldValue<3>::Scalar > ScalarField;
    typedef FieldFE<3, FieldValue<3>::VectorFixed > VecField;

    virtual void SetUp() {
        // setup FilePath directories
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

        Profiler::initialize();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);
        
        FilePath mesh_file( "fields/one_element_2d.msh", FilePath::input_file);
        mesh= new Mesh;
        ifstream in(string( mesh_file ).c_str());
        mesh->read_gmsh_from_stream(in);
        dh = new DOFHandlerMultiDim(*mesh);
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, 3, dof_values, &v);
    }
    virtual void TearDown() {
    	delete dh;
    	delete mesh;
    }

    Mesh *mesh;
    DOFHandlerMultiDim *dh;
    double dof_values[3];
    Vec v;

	MappingP1<1,3> map1;
	MappingP1<2,3> map2;
	MappingP1<3,3> map3;

};


TEST_F(FieldFETest, scalar) {
	FE_P_disc<1,1,3> fe1;
	FE_P_disc<1,2,3> fe2;
	FE_P_disc<1,3,3> fe3;
    ScalarField field;

    dh->distribute_dofs(fe1, fe2, fe3);
    field.set_fe_data(dh, &map1, &map2, &map3, &v);
    field.set_time(0.0);

    dof_values[0] = 1;
    dof_values[1] = 2;
    dof_values[2] = 3;

    vector<double> values(3);

    // test values at vertices of the triangle
    field.value_list( { { 1, 1, 5 }, { 4, 0, 5 }, { 2, 3, 5 } }, mesh->element_accessor(0), values );
    EXPECT_DOUBLE_EQ( dof_values[0], values[0] );
    EXPECT_DOUBLE_EQ( dof_values[1], values[1] );
    EXPECT_DOUBLE_EQ( dof_values[2], values[2] );

    // test value at barycenter
    EXPECT_DOUBLE_EQ( (dof_values[0]+dof_values[1]+dof_values[2])/3, field.value({ 7./3, 4./3, 5 }, mesh->element_accessor(0)) );
}


TEST_F(FieldFETest, vector) {
	FE_RT0<1,3> fe1;
	FE_RT0<2,3> fe2;
	FE_RT0<3,3> fe3;
    VecField field;

    dh->distribute_dofs(fe1, fe2, fe3);
    field.set_fe_data(dh, &map1, &map2, &map3, &v);
    field.set_time(0.0);

    // The Raviart-Thomas function given by the following dofs
    // is 3/7*(x-7/3, y-4/3, 0).
    dof_values[0] = 0;
    dof_values[1] = 1;
    dof_values[2] = 0;

    arma::vec3 result = { 2./7, 1./14, 0 };

    EXPECT_NEAR( 0, arma::norm(result - field.value({ 3, 1.5, 5 }, mesh->element_accessor(0)), 2), 1e-15 );
}




