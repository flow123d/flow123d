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

#include "elasticity_mockup.impl.hh"
#include "DG_mockup_meshes.hh"


TEST_F(ElasticityMockupTest, simple_asm) {
    string eq_data_input = R"YAML(
    solver: !Petsc
      a_tol: 1.0e-12
      r_tol: 1.0e-12
    input_fields:
      - region: .BOUNDARY
        bc_type: displacement
        bc_displacement: 0
        # bc_type: traction
        # bc_traction: [0,-20,0]      
      - region: BULK
        young_modulus: !FieldFormula
          value: 1e8*X[1]
        poisson_ratio: !FieldFormula
          value: 0.1*X[0]
    )YAML";

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
    this->profiler_output("elasticity_asm");
    Profiler::uninitialize();
}
