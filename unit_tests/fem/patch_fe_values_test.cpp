/*
 * fe_values_test.cpp
 *
 *  Created on: Sep 9, 2012
 *      Author: jb
 */

#define FEAL_OVERRIDE_ASSERTS

#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>
#include <cmath>
#include "arma_expect.hh"
#include "armadillo"
#include <mesh_constructor.hh>

#include "system/armadillo_tools.hh"
#include "system/sys_profiler.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/fe_p.hh"
#include "fem/fe_system.hh"
#include "fem/fe_values.hh"
#include "fem/mapping_p1.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "mesh/region.hh"
#include "mesh/accessors.hh"


class PatchFEValuesTest : public testing::Test {
protected:
    static constexpr unsigned int DIM = 2;
    static constexpr unsigned int ORDER = 2;

    Mesh *mesh_;
    Quadrature *quad_;
    UpdateFlags u;
    shared_ptr<FiniteElement<DIM>> fe_;
    shared_ptr<FiniteElement<DIM>> fe_vector_;
    shared_ptr<FiniteElement<DIM>> fe_tensor_;

public:
    PatchFEValuesTest() {
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        mesh_ = mesh_full_constructor("{ mesh_file=\"fem/small_mesh.msh\", optimize_mesh=false }");
        quad_ = new QGauss(DIM, ORDER);
        u = update_values | update_gradients | update_JxW_values | update_jacobians | update_quadrature_points;

        fe_ = std::make_shared< FE_P_disc<DIM> >(ORDER);
        fe_vector_ = std::make_shared<FESystem<DIM>>(fe_, FEVector, 3);
        fe_tensor_ = std::make_shared<FESystem<DIM>>(fe_, FETensor, 9);
    }

    void compute_and_check_values( std::vector< unsigned int > elm_idx_vec, shared_ptr<FiniteElement<DIM>> fe) {
        PatchFEValues<3> patch_fe_values(mesh_->n_elements());
        FEValues<3> fe_values;
        patch_fe_values.initialize(*quad_, *fe, u);
        fe_values.initialize(*quad_, *fe, u);

        // Reinitializes FE data of all 2D elements in one step
        patch_fe_values.reinit(mesh_, elm_idx_vec);
        EXPECT_EQ(patch_fe_values.used_size(), 5);
        EXPECT_EQ(patch_fe_values.max_size(), 5);

        // Iterates over elements, reinitializes FEValues, compares values with FEPatchValues
        for (uint i=0; i<elm_idx_vec.size(); ++i)
        {
            ElementAccessor<3> elm(mesh_, elm_idx_vec[i]);
            if (elm.dim() != 2) continue;
            fe_values.reinit(elm);
            patch_fe_values.get_cell(i);

            // Loops are swapped, JxW and determinant don't depend on dofs
            for (unsigned int k=0; k<patch_fe_values.n_points(); k++)
            {
                for (unsigned int j=0; j<patch_fe_values.n_dofs(); j++)
                {
                    EXPECT_EQ(fe_values.shape_value(j,k), patch_fe_values.shape_value(j,k));
                    EXPECT_ARMA_EQ(fe_values.shape_grad(j,k), patch_fe_values.shape_grad(j,k));
                }
                EXPECT_EQ(fe_values.JxW(k), patch_fe_values.JxW(k));
                //EXPECT_EQ(fe_values.determinant(k), patch_fe_values.determinant(k));
            }
        }
    }
};


TEST_F(PatchFEValuesTest, patch_values_test) {
    // Fills vector of element ids
    std::vector< unsigned int > elm_idx_vec;
    for (auto elm : mesh_->elements_range()) {
        elm_idx_vec.push_back(elm.idx());
    }

    compute_and_check_values(elm_idx_vec, this->fe_);
    compute_and_check_values(elm_idx_vec, this->fe_vector_);
    compute_and_check_values(elm_idx_vec, this->fe_tensor_);
}

