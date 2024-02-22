/*
 * assembly_speed_test.cpp
 *
 *  Created on: Aug 16, 2023
 *      Author: David Flanderka
 *
 *  Speed tests of assembly
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include "DG_mockup.impl.hh"
#include "DG_mockup_meshes.hh"


/****************************************************************************************
 *                 Speed test of Assembly
 *
 * Results:
 * Mesh 27936 elements, 100 assemblation loops - TODO not actual
 * Checked GenericAssembly with active bulk integral only vs. with all active integrals
 *
 *                           bulk      all
 * add_integrals_to_patch   19.12    43.95
 * create_patch              3.00    18.59
 * cache_update              8.43    26.38
 *
 ****************************************************************************************/


TEST_F(DGMockupTest, simple_asm) {
    string eq_data_input = R"YAML(
    solver: !Petsc
      a_tol: 1.0e-12
      r_tol: 1.0e-12
    input_fields:
      - region: .BOUNDARY
        bc_type: diffusive_flux
      - region: BULK
        porosity: !FieldFormula
          value: 0.1*X[0]*X[1]
        init_conc:
          - !FieldConstant
            value: 0.31
          - !FieldFormula
            value: 1-X[0]*X[1]
        diff_m:
          - !FieldConstant
            value: [ [ 0.04, 0.02, 0 ], [ 0.02, 0.01, 0 ], [ 0, 0, 0 ] ]
          - !FieldFormula
            value: "[ [ 0.01*X[0], 0.2*X[1], 1 ], [ 0.2*X[1], 0.01*X[0], 2 ], [ 1, 2, 3 ] ]"
    )YAML";

    std::vector< std::shared_ptr<CodePoint> > cp_vec;
    for (uint i=0; i<meshes_table.size(); ++i)
    {
        // replace START_TIMER tag, we can't set constexpr string converted from meshes_table[i]
        cp_vec.emplace_back( new CODE_POINT(meshes_table[i].c_str()) );
        TimerFrame timer = TimerFrame( *cp_vec[i] );

        DGMockup<Mass_FullAssembly, Stiffness_FullAssembly, Sources_FullAssembly> test;
        test.create_and_set_mesh( "../../benchmark_meshes/" + meshes_table[i] + ".msh");
        test.initialize( eq_data_input, {"A", "B"} );
        test.eq_fields_->init_field_constants(1, 0.5, 0.75, 1, 0.25, 0.5, arma::vec3("1 2 3"), arma::mat33("0.5 0 0, 0 0.75 0, 0 0 1"));
        test.run_simulation();

        // replace END_TIMER equivalent as START_TIMER
        Profiler::instance()->stop_timer( *cp_vec[i] );
    }
    this->profiler_output("const_simple");
    Profiler::uninitialize();
}
