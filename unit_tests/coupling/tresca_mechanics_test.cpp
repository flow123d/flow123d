#define TEST_USE_PETSC
#define TEST_HAS_MAIN
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "mechanics/elasticity.hh"
#include "la/linsys_PERMON.hh"
#include "fields/field_constant.hh"
#include "tools/time_governor.hh"
#include "permonsys.h"
#include "system/sys_profiler.hh"

TEST(TrescaMechanics, reference_is_frozen_during_coupling_iterations)
{
    Profiler::instance();
    FilePath::set_io_dirs(".", UNIT_TESTS_SRC_DIR, "", UNIT_TESTS_BIN_DIR "/coupling");
    {
        std::unique_ptr<Mesh> mesh(mesh_full_constructor(
            "{mesh_file=\"../tests/00_mesh/square_1x1_frac_chopped.msh\", optimize_mesh=false}"));
        const std::string input = R"YAML(
contact: true
tresca_friction: true
solver: !Petsc
  r_tol: 1e-10
  a_tol: 1e-12
  max_it: 10000
output_stream:
  file: tresca_history.pvd
  format: !vtk
    variant: ascii
output:
  fields: [displacement]
input_fields:
  - region: bulk
    young_modulus: 1
    poisson_ratio: 0.25
  - region: fracture
    young_modulus: 1e-6
    poisson_ratio: 0.25
    cross_section_min: 0.01
    friction_threshold: 0.001
  - region: .left_x
    bc_type: displacement
    bc_displacement: 0
  - region: .right_x
    bc_type: traction
    bc_traction: [-0.2, 0.2, 0]
  - region: .right_x
    time: 1
    bc_traction: [-0.2, 0.1, 0]
  - region: .right_x
    time: 2
    bc_traction: [-0.2, -0.2, 0]
)YAML";
        auto &input_type = const_cast<Input::Type::Record &>(Elasticity::get_input_type());
        input_type.finish();
        Input::ReaderToStorage reader(input, input_type, Input::FileFormat::format_YAML);
        TimeGovernor time(0.0, 1.0);
        Elasticity mechanics(*mesh, reader.get_root_interface<Input::Record>(), &time);
        auto scalar = [](double value) {
            auto field = std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
            field->set_value(value);
            return field;
        };
        mechanics.eq_fields().cross_section.set(scalar(1.0), 0.0, {"bulk", ".BOUNDARY"});
        mechanics.eq_fields().cross_section.set(scalar(0.01), 0.0, {"fracture"});
        Field<3, FieldValue<3>::Scalar> potential, reference_potential;
        potential.name("pressure_potential"); reference_potential.name("reference_pressure_potential");
        potential.set_mesh(*mesh); reference_potential.set_mesh(*mesh);
        potential.set(scalar(0.0), 0.0); reference_potential.set(scalar(0.0), 0.0);
        potential.set_time(time.step(), LimitSide::right);
        reference_potential.set_time(time.step(), LimitSide::right);
        mechanics.set_potential_load(potential, reference_potential);
        mechanics.initialize();
        mechanics.zero_time_step();
        auto &data = mechanics.eq_data();
        auto *solver = static_cast<LinSys_PERMON *>(data.ls);
        Vec previous, difference;
        chkerr(VecDuplicate(data.tresca_reference_jump, &previous));
        chkerr(VecDuplicate(previous, &difference));
        chkerr(MatMult(data.tangential_jump_matrix, mechanics.get_solution(), previous));
        PetscReal slip;
        chkerr(VecNorm(previous, NORM_2, &slip));
        EXPECT_GT(slip, 1e-8);

        // A second solve at t=0 must not turn the current iterate into slip history.
        mechanics.solve_linear_system();
        PetscReal reference_norm;
        chkerr(VecNorm(data.tresca_reference_jump, NORM_2, &reference_norm));
        EXPECT_EQ(0.0, reference_norm);

        for (unsigned int step = 0; step < 2; ++step) {
            chkerr(MatMult(data.tangential_jump_matrix, mechanics.get_solution(), previous));
            time.next_time(); // The shared governor advances without Elasticity::next_time().
            for (unsigned int iteration = 0; iteration < 2; ++iteration) {
                mechanics.solve_linear_system();
                chkerr(VecWAXPY(difference, -1.0, previous, data.tresca_reference_jump));
                chkerr(VecNorm(difference, NORM_2, &reference_norm));
                EXPECT_LT(reference_norm, 1e-12);
                EXPECT_LT(solver->compute_residual(), 1e-6);
            }
        }
        chkerr(VecDestroy(&difference)); chkerr(VecDestroy(&previous));
    }
    Profiler::uninitialize();
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (PermonInitialize(&argc, &argv, nullptr, nullptr)) return 1;
    auto &listeners = ::testing::UnitTest::GetInstance()->listeners();
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new ::testing::internal::MPI_PrettyUnitTestResultPrinter);
    const int result = RUN_ALL_TESTS();
    if (PermonFinalize()) return 1;
    return result;
}
