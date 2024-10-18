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
#include "bench_meshes_handler.hh"


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

    BenchMeshesHandler mesh_handler;
    std::vector<std::string> meshes_table = mesh_handler.get_mesh_names("dg_asm");
    std::vector<std::string> meshes_sizes = mesh_handler.get_mesh_sizes();

    std::vector< std::shared_ptr<CodePoint> > cp_vec;
    std::vector< std::shared_ptr<CodePoint> > cp_vec_in;
    for (uint i=0; i<meshes_table.size(); ++i)
    {
        // replace START_TIMER tag, we can't set constexpr string converted from meshes_table[i]
    	cp_vec.emplace_back( new CODE_POINT(meshes_table[i].c_str()) );
        TimerFrame timer = TimerFrame( *cp_vec[i] );

        for (uint j=0; j<meshes_sizes.size(); ++j)
        {
        	cp_vec_in.emplace_back( new CODE_POINT(meshes_sizes[j].c_str()) );
            TimerFrame timer = TimerFrame( *cp_vec_in[i*meshes_sizes.size()+j] );

            std::string mesh_file = "../../benchmark_meshes/" + meshes_table[i] + "_" + meshes_sizes[j] + ".msh";
            uint n_repeats = 1;
            if (meshes_sizes[j] == "small") {
                // 10 repeats of simulation for small meshes
                n_repeats = 10;
            }

            for (uint k=0; k<n_repeats; ++k) {
                // measure needs call of assembly algorithm in separate methods
                this->run_fullassembly_model(eq_data_input, mesh_file);
                this->run_computelocal_model(eq_data_input, mesh_file);
                this->run_evalfields_const(eq_data_input, mesh_file);
                this->run_evalfields_model(eq_data_input, mesh_file);
            }
            Profiler::instance()->stop_timer( *cp_vec_in[i*meshes_sizes.size()+j] );

        }
        // replace END_TIMER equivalent as START_TIMER
        Profiler::instance()->stop_timer( *cp_vec[i] );
    }
    this->profiler_output("dg_asm");
    Profiler::uninitialize();
}
