/*
 * field_evaluate_test.cpp
 *
 *  Created on: Jan 25, 2013
 *      Author: jb
 */



#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>

#include "la/vector_mpi.hh"
#include "fields/fe_value_handler.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_p.hh"
#include "mesh/mesh.h"
#include "quadrature/quadrature_lib.hh"

class FieldEvaluateTest : public testing::Test {
public:

    virtual void SetUp() {
        // setup FilePath directories
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    	this->mesh = mesh_full_constructor("{mesh_file=\"fields/one_element_2d.msh\"}");
    }

    virtual void TearDown() {
    	dh.reset();
    	if (mesh != nullptr) delete mesh;
    }

    void create_dof_handler(double val1, double val2, double val3) {
        dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
        v.resize(3);
        v[0] = val1;
        v[1] = val2;
        v[2] = val3;
    }

    FEValueInitData fill_init_data(unsigned int n_comp) {
    	FEValueInitData init_data;
    	init_data.dh = dh;
    	init_data.data_vec = v;
    	init_data.ndofs = dh->max_elem_dofs();
    	init_data.n_comp = 3;
    	init_data.comp_index = 0;
        return init_data;
    }

    Mesh *mesh;
    std::shared_ptr<DOFHandlerMultiDim> dh;
    VectorMPI v;

};


TEST_F(FieldEvaluateTest, value_handler_scalar) {
	create_dof_handler(1, 2, 3);

	FE_P_disc<0> fe0(1);
	FE_P_disc<1> fe1(1);
	FE_P_disc<2> fe2(1);
	FE_P_disc<3> fe3(1);
	std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, &fe0, &fe1, &fe2, &fe3);
    dh->distribute_dofs(ds);

    FEValueHandler<2, 3, FieldValue<3>::Scalar> val_handler; // Mesh with only one 2D element
	val_handler.initialize( fill_init_data(1) );

	auto acc = mesh->element_accessor(0);
	auto val = val_handler.value(acc.centre(), acc);
	std::cout << "val: " << val << std::endl;
}


TEST_F(FieldEvaluateTest, field_evaluate_scalar) {
	create_dof_handler(1, 2, 3);

	FE_P_disc<0> fe0(1);
	FE_P_disc<1> fe1(1);
	FE_P_disc<2> fe2(1);
	FE_P_disc<3> fe3(1);
	std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, &fe0, &fe1, &fe2, &fe3);
    dh->distribute_dofs(ds);

    Quadrature<2> *quad = new QGauss<2>(1);
    FieldEvaluate<2, 3, FieldValue<3>::Scalar> field_eval(quad); // Mesh with only one 2D element
	field_eval.initialize( fill_init_data(1) );

	auto val = field_eval.value(mesh->element_accessor(0));
	std::cout << "size: " << val.size() << std::endl;
	std::cout << "val: " << val[0] << std::endl;
}


TEST_F(FieldEvaluateTest, value_handler_vector) {
	create_dof_handler(1, 2, 3);

	FE_RT0<0> fe0;
    FE_RT0<1> fe1;
    FE_RT0<2> fe2;
    FE_RT0<3> fe3;
	std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, &fe0, &fe1, &fe2, &fe3);
    dh->distribute_dofs(ds);

    FEValueHandler<2, 3, FieldValue<3>::VectorFixed> val_handler; // Mesh with only one 2D element
	val_handler.initialize( fill_init_data(3) );

	auto acc = mesh->element_accessor(0);
	auto val = val_handler.value(acc.centre(), acc);
	std::cout << "size: " << val.size() << std::endl;
	std::cout << "val: " << val[0] << ", " << val[1] << ", " << val[2] << std::endl;
}


TEST_F(FieldEvaluateTest, field_evaluate_vector) {
	create_dof_handler(1, 2, 3);

    FE_RT0<0> fe0;
    FE_RT0<1> fe1;
    FE_RT0<2> fe2;
    FE_RT0<3> fe3;
	std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, &fe0, &fe1, &fe2, &fe3);
    dh->distribute_dofs(ds);

    Quadrature<2> *quad = new QGauss<2>(1);
    FieldEvaluate<2, 3, FieldValue<3>::VectorFixed> field_eval(quad); // Mesh with only one 2D element
	field_eval.initialize( fill_init_data(3) );

	auto val = field_eval.value(mesh->element_accessor(0));
	std::cout << "size: " << val.size() << std::endl;
	std::cout << "val: " << val[0] << std::endl;
}


