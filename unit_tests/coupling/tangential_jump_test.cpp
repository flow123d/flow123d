#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>

#include "mechanics/assembly_elasticity.hh"
#include "fields/field_constant.hh"
#include "tools/time_governor.hh"

/// Exercise the assembled trace operator on exact affine bulk displacement fields.
class TangentialJumpTest : public testing::Test {
protected:
    /// Initialize profiling and fixture paths before constructing a mesh.
    void SetUp() override {
        FilePath::set_io_dirs(".", UNIT_TESTS_SRC_DIR, "", ".");
        Profiler::instance();
    }

    /// Release profiling after all local mesh and assembly objects have been destroyed.
    void TearDown() override {
        Profiler::uninitialize();
    }

    /// Verify force units/measure scaling and collective rejection of unsupported 3D friction.
    void check_threshold(const std::string &mesh_file, bool unsupported_surface) {
        std::unique_ptr<Mesh> mesh(mesh_full_constructor(
            "{mesh_file=\"coupling/" + mesh_file + "\", optimize_mesh=false}"));
        auto fields = std::make_shared<Elasticity::EqFields>();
        fields->set_mesh(*mesh);
        auto strength = std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
        strength->set_value(2.0);
        fields->friction_threshold.set(strength, 0.0);
        TimeGovernor time(0.0, 1.0);
        fields->friction_threshold.set_time(time.step(), LimitSide::right);
        Elasticity::EqData data(fields);
        data.create_dh(mesh.get());
        data.allocate_tangential_jump_matrix();
        chkerr(MatCreateVecs(data.tangential_jump_matrix, nullptr, &data.tresca_force_bound));
        GenericAssembly<Elasticity::TrescaThresholdAssemblyDim> assembly(&data, data.dh_.get());
        if (unsupported_surface) {
            EXPECT_THROW(assembly.assemble(data.dh_), std::exception);
        } else {
            assembly.assemble(data.dh_);
            for (const auto &entry : data.tangential_jump_idx) {
                PetscInt row = entry.second;
                PetscScalar value;
                chkerr(VecGetValues(data.tresca_force_bound, 1, &row, &value));
                EXPECT_NEAR(6.0, value, 1e-12); // Tilted line length is 3; traction threshold is 2.
            }
        }
    }

    /// Check geometry, signs, averaging, fracture-DOF exclusion, and distributed ownership.
    void check_jump(const std::string &mesh_file, unsigned int fracture_dim, bool two_sided = true) {
        std::unique_ptr<Mesh> mesh(mesh_full_constructor(
            "{mesh_file=\"coupling/" + mesh_file + "\", optimize_mesh=false}"));
        auto fields = std::make_shared<Elasticity::EqFields>();
        fields->set_mesh(*mesh);
        fields->set_components({"displacement"});
        Elasticity::EqData data(fields);
        data.create_dh(mesh.get());
        data.allocate_tangential_jump_matrix();
        GenericAssembly<Elasticity::TangentialJumpAssemblyDim> assembly(&data, data.dh_.get());
        chkerr(MatZeroEntries(data.tangential_jump_matrix));
        assembly.assemble(data.dh_);
        chkerr(MatAssemblyBegin(data.tangential_jump_matrix, MAT_FINAL_ASSEMBLY));
        chkerr(MatAssemblyEnd(data.tangential_jump_matrix, MAT_FINAL_ASSEMBLY));

        PetscInt n_rows, n_cols, first_row, last_row;
        chkerr(MatGetSize(data.tangential_jump_matrix, &n_rows, &n_cols));
        EXPECT_EQ(two_sided ? fracture_dim : 0, static_cast<unsigned int>(n_rows));
        EXPECT_EQ(data.dh_->n_global_dofs(), static_cast<unsigned int>(n_cols));
        chkerr(MatGetOwnershipRange(data.tangential_jump_matrix, &first_row, &last_row));
        for (const auto &entry : data.tangential_jump_idx) {
            EXPECT_EQ(10, mesh->element_accessor(entry.first).input_id());
            EXPECT_GE(entry.second, first_row);
            EXPECT_LE(entry.second + fracture_dim, last_row);
        }

        arma::vec3 tangents[2], centre, normal;
        for (auto fracture : mesh->elements_range()) {
            if (fracture.input_id() != 10) continue;
            centre = fracture.centre();
            tangents[0] = *fracture.node(1) - *fracture.node(0);
            tangents[0] /= arma::norm(tangents[0], 2);
            if (fracture_dim == 2) {
                normal = arma::cross(tangents[0], *fracture.node(2) - *fracture.node(0));
                normal /= arma::norm(normal, 2);
                tangents[1] = arma::cross(normal, tangents[0]);
            } else {
                normal = arma::cross(tangents[0], arma::vec3("0 0 1"));
                normal /= arma::norm(normal, 2);
            }
        }

        Vec displacement, jump;
        chkerr(MatCreateVecs(data.tangential_jump_matrix, &displacement, &jump));
        arma::mat33 gradient("1 2 -1; -2 0 3; 1 -1 2");
        for (unsigned int motion = 0; motion < 4; ++motion) {
            arma::vec3 plus("0.2 -0.4 0.5"), minus("-0.1 0.3 -0.2");
            if (motion == 0) minus = plus;
            if (motion == 1) { plus = normal; minus.zeros(); }
            chkerr(VecSet(displacement, 0.0));
            for (auto cell : data.dh_->own_range()) {
                std::vector<LongIdx> indices(cell.n_dofs());
                cell.get_dof_indices(indices);
                std::vector<PetscScalar> values(indices.size());
                for (unsigned int i = 0; i < indices.size(); ++i) {
                    const auto &dof = cell.cell_dof(i);
                    arma::vec3 value;
                    if (cell.dim() == fracture_dim + 1) {
                        bool plus_side = cell.elm().input_id() == 20;
                        value = plus_side ? plus : minus;
                        if (motion == 3)
                            value += (plus_side ? 1.0 : -0.5) * gradient * (*cell.elm().node(dof.n_face_idx));
                    } else {
                        value.fill(99.0); // Fracture DOFs must not contribute to either bulk trace.
                    }
                    values[i] = arma::dot(dof.coefs, value);
                }
                chkerr(VecSetValues(displacement, indices.size(), indices.data(), values.data(), INSERT_VALUES));
            }
            chkerr(VecAssemblyBegin(displacement));
            chkerr(VecAssemblyEnd(displacement));
            chkerr(MatMult(data.tangential_jump_matrix, displacement, jump));
            arma::vec3 expected = plus - minus;
            if (motion == 3) expected += 1.5 * gradient * centre;
            for (const auto &entry : data.tangential_jump_idx) {
                for (unsigned int k = 0; k < fracture_dim; ++k) {
                    PetscInt row = entry.second + k;
                    PetscScalar value;
                    chkerr(VecGetValues(jump, 1, &row, &value));
                    EXPECT_NEAR(arma::dot(tangents[k], expected), value, 1e-12);
                }
            }
        }
        chkerr(VecDestroy(&jump));
        chkerr(VecDestroy(&displacement));
    }
};

TEST_F(TangentialJumpTest, tilted_line) {
    check_jump("tangential_jump_2d.msh", 1);
}

TEST_F(TangentialJumpTest, tilted_triangle) {
    check_jump("tangential_jump_3d.msh", 2);
}

TEST_F(TangentialJumpTest, one_sided_fracture_has_no_jump_rows) {
    check_jump("tangential_jump_one_sided.msh", 1, false);
}

TEST_F(TangentialJumpTest, integrated_tresca_threshold) {
    check_threshold("tangential_jump_2d.msh", false);
}

TEST_F(TangentialJumpTest, rejects_isotropic_tresca_on_triangle) {
    check_threshold("tangential_jump_3d.msh", true);
}
