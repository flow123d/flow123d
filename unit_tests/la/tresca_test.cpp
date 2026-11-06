#define TEST_USE_PETSC
#define TEST_HAS_MAIN
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include "la/distribution.hh"
#include "la/linsys_PERMON.hh"
#include "permonsys.h"
#include "system/sys_profiler.hh"
#include <algorithm>

/// Independent supported normal/tangential spring pair on every MPI rank.
class TrescaSolverTest : public testing::Test {
protected:
    /// Initialize profiling for the algebraic solver's timers.
    void SetUp() override { Profiler::instance(); }
    /// Release profiling after local solver objects have been destroyed.
    void TearDown() override { Profiler::uninitialize(); }

    /// Check exact solutions of a coupled 2x2 primal energy under prescribed loading histories.
    void check_loading(bool open_contact, bool general_bounds = false, bool reverse_jump = false) {
        Distribution distribution(2, PETSC_COMM_WORLD);
        double solution[2] = {0.0, 0.0};
        LinSys_PERMON ls(&distribution);
        ls.set_solution(solution);
        ls.set_tolerances(1e-12, 1e-12, 1e6, 1000);
        int dofs[2] = {static_cast<int>(distribution.begin()), static_cast<int>(distribution.begin()+1)};
        double K[4] = {2.0, 0.4, 0.4, 3.0};
        ls.start_allocation();
        ls.mat_set_values(2, dofs, 2, dofs, K);
        ls.start_add_assembly();
        ls.mat_set_values(2, dofs, 2, dofs, K);
        ls.finish_assembly();

        int rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        Mat N, T;
        Vec gap, bound, lower, reference;
        chkerr(MatCreateAIJ(PETSC_COMM_WORLD, 1, 2, PETSC_DECIDE, PETSC_DECIDE, 1, nullptr, 0, nullptr, &N));
        chkerr(MatCreateAIJ(PETSC_COMM_WORLD, 1, 2, PETSC_DECIDE, PETSC_DECIDE, 1, nullptr, 0, nullptr, &T));
        chkerr(MatSetValue(N, rank, dofs[0], 1.0, INSERT_VALUES));
        const double jump_sign = reverse_jump ? -1.0 : 1.0;
        chkerr(MatSetValue(T, rank, dofs[1], jump_sign, INSERT_VALUES));
        chkerr(MatAssemblyBegin(N, MAT_FINAL_ASSEMBLY)); chkerr(MatAssemblyEnd(N, MAT_FINAL_ASSEMBLY));
        chkerr(MatAssemblyBegin(T, MAT_FINAL_ASSEMBLY)); chkerr(MatAssemblyEnd(T, MAT_FINAL_ASSEMBLY));
        chkerr(MatCreateVecs(N, nullptr, &gap));
        chkerr(MatCreateVecs(T, nullptr, &bound));
        chkerr(VecDuplicate(bound, &reference));
        chkerr(VecDuplicate(bound, &lower));
        chkerr(VecSet(gap, 0.0));
        chkerr(VecSet(reference, 0.0));
        ls.set_inequality(N, gap);
        ls.set_dual_box(T, reference, lower, bound);

        // Rows: {normal load, tangential load, lower, upper, reference, expected u_t, expected q}.
        const double closed[][7] = {
            {2., 0.3, -0.5, 0.5, 0., 0., 0.3},        // stick
            {2., 2., -0.5, 0.5, 0., 0.5, 0.5},        // slide forward
            {2., 1.7, -0.5, 0.5, 0., 0.4, 0.5},       // repeated load-step solve, same reference
            {2., 1.7, -0.5, 0.5, 0.5, 0.5, 0.2},      // accepted previous step, unload and stick
            {2., -1., -0.5, 0.5, 0.5, -1./6., -0.5},  // reversed sliding
            {2., 2., 0., 0., 0., 2./3., 0.},           // zero threshold
            {2., 2., -5., 5., 0., 0., 2.}              // large threshold
        };
        const double open[][7] = {{-2., 0.2, -1., 1., 0., 0., 0.6}};
        const double boxes[][7] = {
            {2., 0.3, 0.4, 0.8, 0., -1./30., 0.4}, // positive lower bound
            {2., 2., -0.8, -0.4, 0., 0.8, -0.4},   // negative upper bound
            {2., 2., 0.25, 0.25, 0., 7./12., 0.25}, // fixed nonzero multiplier
            {2., 0.3, PETSC_NINFINITY, PETSC_INFINITY, 0., 0., 0.3}
        };
        unsigned int n_steps = general_bounds ? sizeof(boxes)/sizeof(boxes[0])
                                             : open_contact ? 1 : sizeof(closed)/sizeof(closed[0]);
        for (unsigned int step = 0; step < n_steps; ++step) {
            double load[7];
            std::copy_n(general_bounds ? boxes[step] : open_contact ? open[step] : closed[step], 7, load);
            // An empty active-friction block on odd ranks exercises distributed
            // row extraction and reconstruction with unequal multiplier layouts.
            if (!open_contact && !general_bounds && step == 1 && rank % 2) {
                load[2] = load[3] = 0.0; load[5] = 2.0/3.0; load[6] = 0.0;
            }
            if (general_bounds && step == 0 && rank % 2) {
                load[2] = load[3] = 0.25; load[5] = 1./60.; load[6] = 0.25;
            }
            ls.start_add_assembly();
            ls.rhs_zero_entries();
            double f[2] = {load[0], load[1]};
            ls.rhs_set_values(2, dofs, f);
            ls.finish_assembly();
            chkerr(VecSetValue(bound, rank, load[3], INSERT_VALUES));
            chkerr(VecAssemblyBegin(bound)); chkerr(VecAssemblyEnd(bound));
            chkerr(VecSetValue(lower, rank, load[2], INSERT_VALUES));
            chkerr(VecAssemblyBegin(lower)); chkerr(VecAssemblyEnd(lower));
            chkerr(VecSet(reference, jump_sign * load[4]));
            auto result = ls.solve();
            EXPECT_GT(result.converged_reason, 0);
            EXPECT_NEAR(open_contact ? -1.0 : 0.0, solution[0], 1e-9);
            EXPECT_NEAR(load[5], solution[1], 1e-9);
            PetscScalar q, normal;
            chkerr(VecGetValues(ls.dual_box_multipliers(), 1, &rank, &q));
            chkerr(VecGetValues(ls.dual_ineq_multipliers(), 1, &rank, &normal));
            EXPECT_NEAR(jump_sign * load[6], q, 1e-9);
            EXPECT_NEAR(open_contact ? 0.0 : load[0] - 0.4*load[5], normal, 1e-9);
            EXPECT_LT(ls.compute_residual(), 1e-9);
        }

        chkerr(VecDestroy(&gap)); chkerr(VecDestroy(&bound)); chkerr(VecDestroy(&lower));
        chkerr(VecDestroy(&reference));
        chkerr(MatDestroy(&N)); chkerr(MatDestroy(&T));
    }
};

TEST_F(TrescaSolverTest, stick_slide_unload_reverse) { check_loading(false); }
TEST_F(TrescaSolverTest, prescribed_tresca_strength_on_open_contact) { check_loading(true); }
TEST_F(TrescaSolverTest, asymmetric_dual_box_and_fixed_multipliers) { check_loading(false, true); }
TEST_F(TrescaSolverTest, reversing_jump_preserves_displacement_and_physical_reactions) {
    check_loading(false, false, true);
    check_loading(true, false, true);
}

TEST_F(TrescaSolverTest, contact_only_retains_original_inequality_multipliers)
{
    Distribution distribution(1, PETSC_COMM_WORLD);
    double solution = 0.0;
    LinSys_PERMON ls(&distribution);
    ls.set_solution(&solution);
    ls.set_tolerances(1e-12, 1e-12, 1e6, 1000);
    int dof = distribution.begin();
    double stiffness = 2.0;
    ls.start_allocation();
    ls.mat_set_values(1, &dof, 1, &dof, &stiffness);
    ls.start_add_assembly();
    ls.mat_set_values(1, &dof, 1, &dof, &stiffness);
    ls.finish_assembly();

    Mat normal;
    Vec gap;
    chkerr(MatCreateAIJ(PETSC_COMM_WORLD, 1, 1, PETSC_DECIDE, PETSC_DECIDE,
                       1, nullptr, 0, nullptr, &normal));
    chkerr(MatSetValue(normal, dof, dof, 1.0, INSERT_VALUES));
    chkerr(MatAssemblyBegin(normal, MAT_FINAL_ASSEMBLY));
    chkerr(MatAssemblyEnd(normal, MAT_FINAL_ASSEMBLY));
    chkerr(MatCreateVecs(normal, nullptr, &gap));
    chkerr(VecSet(gap, 0.0));
    ls.set_inequality(normal, gap);

    for (double load : {2.0, -2.0, 2.0}) {
        ls.start_add_assembly();
        ls.rhs_zero_entries();
        ls.rhs_set_values(1, &dof, &load);
        ls.finish_assembly();
        EXPECT_GT(ls.solve().converged_reason, 0);
        ASSERT_NE(nullptr, ls.dual_ineq_multipliers());
        PetscInt row = dof;
        PetscScalar multiplier;
        chkerr(VecGetValues(ls.dual_ineq_multipliers(), 1, &row, &multiplier));
        EXPECT_NEAR(load > 0.0 ? 2.0 : 0.0, multiplier, 1e-9);
        EXPECT_NEAR(load > 0.0 ? 0.0 : -1.0, solution, 1e-9);
    }

    ls.set_inequality(nullptr, nullptr);
    EXPECT_GT(ls.solve().converged_reason, 0);
    EXPECT_EQ(nullptr, ls.dual_ineq_multipliers());
    EXPECT_EQ(nullptr, ls.dual_box_multipliers());
    EXPECT_NEAR(1.0, solution, 1e-9);
    chkerr(VecDestroy(&gap));
    chkerr(MatDestroy(&normal));
}

TEST_F(TrescaSolverTest, original_dual_solution_survives_objective_scaling)
{
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    Mat normal, stiffness;
    Vec rhs, gap, displacement;
    chkerr(MatCreateAIJ(PETSC_COMM_WORLD, 1, 1, PETSC_DECIDE, PETSC_DECIDE,
                       1, nullptr, 0, nullptr, &normal));
    chkerr(MatSetValue(normal, rank, rank, 1.0, INSERT_VALUES));
    chkerr(MatAssemblyBegin(normal, MAT_FINAL_ASSEMBLY));
    chkerr(MatAssemblyEnd(normal, MAT_FINAL_ASSEMBLY));
    chkerr(MatDuplicate(normal, MAT_COPY_VALUES, &stiffness));
    chkerr(MatScale(stiffness, 2.0));
    chkerr(MatCreateVecs(normal, &displacement, &gap));
    chkerr(VecDuplicate(displacement, &rhs));
    chkerr(VecSet(displacement, 0.0));
    chkerr(VecSet(gap, 0.0));
    chkerr(VecSet(rhs, 2.0));
    QP primal, unscaled_dual;
    QPS solver;
    chkerr(QPCreate(PETSC_COMM_WORLD, &primal));
    chkerr(QPSetOperator(primal, stiffness));
    chkerr(QPSetRhs(primal, rhs));
    chkerr(QPSetInitialVector(primal, displacement));
    chkerr(QPSetIneq(primal, normal, gap));
    chkerr(QPTDualize(primal, MAT_INV_MONOLITHIC, MAT_REG_NONE));
    chkerr(QPGetChild(primal, &unscaled_dual));
    chkerr(QPTScaleObjectiveByScalar(unscaled_dual, 2.0, 3.0));
    chkerr(QPSCreate(PETSC_COMM_WORLD, &solver));
    chkerr(QPSSetQP(solver, primal));
    chkerr(QPSSetTolerances(solver, 1e-12, 1e-12, 1e6, 1000));
    chkerr(QPSSetFromOptions(solver));
    chkerr(QPSSolve(solver));
    Vec multiplier;
    chkerr(QPGetSolutionVector(unscaled_dual, &multiplier));
    PetscInt row = rank;
    PetscScalar value;
    chkerr(VecGetValues(multiplier, 1, &row, &value));
    EXPECT_NEAR(2.0, value, 1e-9);
    chkerr(VecGetValues(displacement, 1, &row, &value));
    EXPECT_NEAR(0.0, value, 1e-9);
    chkerr(QPSDestroy(&solver));
    chkerr(QPDestroy(&primal));
    chkerr(VecDestroy(&rhs));
    chkerr(VecDestroy(&gap));
    chkerr(VecDestroy(&displacement));
    chkerr(MatDestroy(&normal));
    chkerr(MatDestroy(&stiffness));
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
