/*
 * elasticity_asm_bench.cpp
 *
 *  Created on: Oct 11, 2024
 *      Author: David Flanderka
 *
 *  Speed tests of elasticity assembly
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include "elasticity_mockup.hh"
#include "elasticity_mockup_assembly.hh"
//#include "DG_mockup_meshes.hh"


TEST_F(ElasticityMockupTest, simple_asm) {
    this->profiler_output("elasticity_asm");
    Profiler::uninitialize();
}
