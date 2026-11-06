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

namespace {

/// Gather a small test vector so expected field values are independent of PETSc row ownership.
std::vector<PetscScalar> all_values(Vec source)
{
    Vec local;
    VecScatter scatter;
    chkerr(VecScatterCreateToAll(source, &scatter, &local));
    chkerr(VecScatterBegin(scatter, source, local, INSERT_VALUES, SCATTER_FORWARD));
    chkerr(VecScatterEnd(scatter, source, local, INSERT_VALUES, SCATTER_FORWARD));
    PetscInt size;
    chkerr(VecGetSize(local, &size));
    const PetscScalar *data;
    chkerr(VecGetArrayRead(local, &data));
    std::vector<PetscScalar> result(size);
    if (size) std::copy_n(data, size, result.begin());
    chkerr(VecRestoreArrayRead(local, &data));
    chkerr(VecScatterDestroy(&scatter));
    chkerr(VecDestroy(&local));
    return result;
}

/// Compare real mechanics output with retained physical multipliers after every coupling solve.
void expect_reaction_output(Elasticity &mechanics)
{
    mechanics.update_output_fields();
    auto &data = mechanics.eq_data();
    auto *solver = static_cast<LinSys_PERMON *>(data.ls);
    auto normal = all_values(solver->dual_ineq_multipliers());
    auto tangential = all_values(solver->dual_box_multipliers());
    auto normal_field = mechanics.eq_fields().output_contact_normal_stress_ptr;
    auto tangent_field = mechanics.eq_fields().output_contact_tangential_stress_ptr;
    for (auto cell : data.dh_->own_range()) {
        double expected_normal = 0.0;
        auto normal_row = data.constraint_idx.find(cell.elm_idx());
        if (normal_row != data.constraint_idx.end())
            expected_normal = -normal[normal_row->second] / cell.elm().measure();
        auto scalar_cell = cell.cell_with_other_dh(normal_field->get_dofhandler().get());
        EXPECT_NEAR(expected_normal, normal_field->vec().get(scalar_cell.get_loc_dof_indices()[0]), 1e-10);
        arma::vec3 expected_tangent;
        expected_tangent.zeros();
        auto tangent_row = data.tangential_jump_idx.find(cell.elm_idx());
        if (tangent_row != data.tangential_jump_idx.end()) {
            arma::vec3 direction = *cell.elm().node(1) - *cell.elm().node(0);
            direction /= arma::norm(direction, 2);
            expected_tangent = -direction * tangential[tangent_row->second] / cell.elm().measure();
            EXPECT_LE(arma::norm(expected_tangent, 2), 0.001 + 1e-8);
        }
        auto vector_cell = cell.cell_with_other_dh(tangent_field->get_dofhandler().get());
        auto dofs = vector_cell.get_loc_dof_indices();
        for (unsigned int i = 0; i < dofs.size(); ++i)
            EXPECT_NEAR(arma::dot(vector_cell.cell_dof(i).coefs, expected_tangent),
                        tangent_field->vec().get(dofs[i]), 1e-10);
    }
}

} // namespace

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
  fields: [displacement, contact_normal_stress, contact_tangential_stress]
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
        expect_reaction_output(mechanics);
        Vec previous, difference;
        chkerr(VecDuplicate(data.tresca_reference_jump, &previous));
        chkerr(VecDuplicate(previous, &difference));
        chkerr(MatMult(data.tangential_jump_matrix, mechanics.get_solution(), previous));
        PetscReal slip;
        chkerr(VecNorm(previous, NORM_2, &slip));
        EXPECT_GT(slip, 1e-8);

        // A second solve at t=0 must not turn the current iterate into slip history.
        mechanics.solve_linear_system();
        expect_reaction_output(mechanics);
        PetscReal reference_norm;
        chkerr(VecNorm(data.tresca_reference_jump, NORM_2, &reference_norm));
        EXPECT_EQ(0.0, reference_norm);

        for (unsigned int step = 0; step < 2; ++step) {
            chkerr(MatMult(data.tangential_jump_matrix, mechanics.get_solution(), previous));
            time.next_time(); // The shared governor advances without Elasticity::next_time().
            for (unsigned int iteration = 0; iteration < 2; ++iteration) {
                mechanics.solve_linear_system();
                expect_reaction_output(mechanics);
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
