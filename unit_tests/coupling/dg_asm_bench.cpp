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

        std::string mesh_file = "../../benchmark_meshes/" + meshes_table[i] + ".msh";
        uint n_repeats = 1;
        if (meshes_table[i].find("small") != std::string::npos) {
            // 10 repeats of simulation for small meshes
            n_repeats = 10;
        }

        for (uint j=0; j<n_repeats; ++j) {
            // FullAssembly + field_const
            START_TIMER("FullAssembly_const");
            DGMockup<Mass_FullAssembly, Stiffness_FullAssembly, Sources_FullAssembly> test_full_asm_const;
            test_full_asm_const.create_and_set_mesh(mesh_file);
            test_full_asm_const.initialize( eq_data_input, {"A", "B"} );
            test_full_asm_const.eq_fields_->init_field_constants(1, 0.5, 0.75, 1, 0.25, 0.5, arma::vec3("1 2 3"), arma::mat33("0.5 0 0, 0 0.75 0, 0 0 1"));
            test_full_asm_const.run_simulation();
            END_TIMER("FullAssembly_const");

            // FullAssembly + field_model
            START_TIMER("FullAssembly_model");
            DGMockup<Mass_FullAssembly, Stiffness_FullAssembly, Sources_FullAssembly> test_full_asm_model;
            test_full_asm_model.create_and_set_mesh(mesh_file);
            test_full_asm_model.initialize( eq_data_input, {"A", "B"} );
            test_full_asm_model.eq_fields_->init_field_models();
            test_full_asm_model.run_simulation();
            END_TIMER("FullAssembly_model");

//            // ComputeLocal + field_const
//            START_TIMER("ComputeLocal_const");
//            DGMockup<Mass_ComputeLocal, Stiffness_ComputeLocal, Sources_ComputeLocal> test_comp_local_const;
//            test_comp_local_const.create_and_set_mesh(mesh_file);
//            test_comp_local_const.initialize( eq_data_input, {"A", "B"} );
//            test_comp_local_const.eq_fields_->init_field_constants(1, 0.5, 0.75, 1, 0.25, 0.5, arma::vec3("1 2 3"), arma::mat33("0.5 0 0, 0 0.75 0, 0 0 1"));
//            test_comp_local_const.run_simulation();
//            END_TIMER("ComputeLocal_const");
//
//            // EvalFields + field_const
//            START_TIMER("EvalFields_const");
//            DGMockup<Mass_EvalFields, Stiffness_EvalFields, Sources_EvalFields> test_eval_fields_const;
//            test_eval_fields_const.create_and_set_mesh(mesh_file);
//            test_eval_fields_const.initialize( eq_data_input, {"A", "B"} );
//            test_eval_fields_const.eq_fields_->init_field_constants(1, 0.5, 0.75, 1, 0.25, 0.5, arma::vec3("1 2 3"), arma::mat33("0.5 0 0, 0 0.75 0, 0 0 1"));
//            test_eval_fields_const.run_simulation();
//            END_TIMER("EvalFields_const");
//
//            // EvalFields + field_model
//            START_TIMER("EvalFields_model");
//            DGMockup<Mass_EvalFields, Stiffness_EvalFields, Sources_EvalFields> test_eval_fields_model;
//            test_eval_fields_model.create_and_set_mesh(mesh_file);
//            test_eval_fields_model.initialize( eq_data_input, {"A", "B"} );
//            test_eval_fields_model.eq_fields_->init_field_models();
//            test_eval_fields_model.run_simulation();
//            END_TIMER("EvalFields_model");

        }

        // replace END_TIMER equivalent as START_TIMER
        Profiler::instance()->stop_timer( *cp_vec[i] );
    }
    this->profiler_output("dg_asm");
    Profiler::uninitialize();
}
