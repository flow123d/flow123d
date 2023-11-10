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
#include "system/armadillo_tools.hh"
#include "system/sys_profiler.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/mapping_p1.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "mesh/region.hh"
#include "mesh/accessors.hh"


TEST(FeValues, patch_test_basic) {
    PatchFEValues<3> patch_fe_values(15);
    FEValues<3> fe_values;

    patch_fe_values.resize(10);
    EXPECT_EQ(patch_fe_values.used_size(), 10);
    EXPECT_EQ(patch_fe_values.max_size(), 15);
    patch_fe_values.resize(7);
    EXPECT_EQ(patch_fe_values.used_size(), 7);
    EXPECT_EQ(patch_fe_values.max_size(), 15);
    patch_fe_values.resize(15);
    EXPECT_EQ(patch_fe_values.used_size(), 15);
    EXPECT_EQ(patch_fe_values.max_size(), 15);
}

