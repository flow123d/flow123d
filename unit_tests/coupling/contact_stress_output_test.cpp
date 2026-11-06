#define TEST_USE_PETSC
#define TEST_HAS_MAIN
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "mechanics/elasticity.hh"
#include "mechanics/assembly_elasticity.hh"
#include "la/linsys_PERMON.hh"
#include "fields/field_constant.hh"
#include "tools/time_governor.hh"
#include "permonsys.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

/// Supply physical-row results to the real mechanics output path without solving a synthetic load.
class StoredMultiplierSolver : public LinSys_PERMON {
public:
    using LinSys_PERMON::LinSys_PERMON;

    /// Own copies in the original A and T row layouts, just as after a PERMON solve.
    void store(Vec normal, Vec tangential) {
        chkerr(VecDuplicate(normal, &ineq_multipliers_));
        chkerr(VecCopy(normal, ineq_multipliers_));
        if (tangential) {
            chkerr(VecDuplicate(tangential, &dual_box_multipliers_));
            chkerr(VecCopy(tangential, dual_box_multipliers_));
        }
    }
};

/// Borrow output fields and a primal force vector for independent interface-work integration.
struct PhysicalReactionData {
    using EqFields = Elasticity::EqFields;
    /// Share mechanics fields for the lifetime of the independent quadrature assembly.
    std::shared_ptr<EqFields> eq_fields_;
    /// Borrow the distributed force vector in original global primal DOF numbering.
    Vec force;
};

/// Integrate stress times bulk trace shape functions, without undoing multiplier normalization.
template<unsigned int dim>
class PhysicalReactionAssembly : public AssemblyBasePatch<dim> {
public:
    using EqData = PhysicalReactionData;
    using EqFields = Elasticity::EqFields;
    /// Identify the independent quadrature calculation in profiler output.
    static constexpr const char *name() { return "Test_physical_interface_work"; }

    /// Register physical trace weights, side normals, and adjacent bulk cross-sections.
    PhysicalReactionAssembly(EqData *data, PatchInternals *internals)
    : AssemblyBasePatch<dim>(1, internals), eq_fields_(data->eq_fields_.get()), data_(data),
      integral_(this->create_coupling_integral(this->quad_)),
      weights_(integral_->JxW()), normals_(integral_->normal_vector()), shapes_(integral_->vector_shape()) {
        used_fields_ += eq_fields_->cross_section;
    }

    /// Allocate insertion buffers only when a dimensional join is visited.
    void initialize() {}

    /// Integrate normal and signed tangential output tractions against bulk trace basis functions.
    void dimjoin_intergral(DHCellAccessor cell, DHCellSide neighbour) override {
        if (dim == 3 || !cell.is_own()) return;
        auto normal_field = eq_fields_->output_contact_normal_stress_ptr;
        auto tangent_field = eq_fields_->output_contact_tangential_stress_ptr;
        auto scalar = cell.cell_with_other_dh(normal_field->get_dofhandler().get());
        double normal_stress = normal_field->vec().get(scalar.get_loc_dof_indices()[0]);
        auto vector = cell.cell_with_other_dh(tangent_field->get_dofhandler().get());
        arma::vec3 tangent_stress;
        tangent_stress.zeros();
        auto vector_dofs = vector.get_loc_dof_indices();
        for (unsigned int i = 0; i < vector_dofs.size(); ++i)
            tangent_stress += vector.cell_dof(i).coefs * tangent_field->vec().get(vector_dofs[i]);
        double sign = 1.0;
        if (cell.elm()->n_neighs_vb() == 2) {
            auto first = cell.elm()->neigh_vb[0]->side()->element();
            auto second = cell.elm()->neigh_vb[1]->side()->element();
            sign = neighbour.element().input_id() == std::min(first.input_id(), second.input_id()) ? 1.0 : -1.0;
        }
        std::vector<LongIdx> dofs;
        neighbour.cell().get_dof_indices(dofs);
        std::vector<PetscScalar> forces(dofs.size(), 0.0);
        for (auto p : integral_->points(neighbour)) {
            arma::vec3 traction = normal_stress * eq_fields_->cross_section(p) * normals_(p)
                                  + sign * tangent_stress;
            for (unsigned int i = 0; i < dofs.size(); ++i)
                forces[i] += weights_(p) * arma::dot(shapes_.shape(i)(p), traction);
        }
        chkerr(VecSetValues(data_->force, dofs.size(), dofs.data(), forces.data(), ADD_VALUES));
    }

private:
    EqFields *eq_fields_; ///< Borrow mechanics fields, including output FieldFE vectors.
    EqData *data_; ///< Borrow the force-vector handle and fields for this assembly.
    FieldSet used_fields_; ///< Register the bulk cross-section cache dependency.
    std::shared_ptr<CouplingIntegralAcc<dim>> integral_; ///< Physical interface quadrature.
    FeQ<Scalar> weights_; ///< Physical interface weights, without average-jump normalization.
    ElQ<Vector> normals_; ///< Outward normals of the adjacent bulk side.
    FeQArray<Vector> shapes_; ///< Original displacement basis traces on the bulk side.
    template <template<IntDim...> class DimAssembly>
    friend class GenericAssembly;
};

template<unsigned int dim> using PhysicalReactionDim = PhysicalReactionAssembly<dim>;

/// Verify original-row imports and P0 conversion on a tilted line of length three.
class ContactStressOutputTest : public testing::Test {
protected:
    /// Initialize profiling and isolate output under the unit-test binary directory.
    void SetUp() override {
        Profiler::instance();
        FilePath::set_io_dirs(".", UNIT_TESTS_SRC_DIR, "", UNIT_TESTS_BIN_DIR "/coupling");
    }

    /// Release mechanics and PETSc objects before ending profiling.
    void TearDown() override {
        mechanics_.reset();
        time_.reset();
        mesh_.reset();
        Profiler::uninitialize();
    }

    /// Initialize the production fields and output assembly without a mechanics solve.
    void initialize(bool contact, bool friction) {
        mesh_.reset(mesh_full_constructor(
            "{mesh_file=\"coupling/tangential_jump_2d.msh\", optimize_mesh=false}"));
        std::string input = std::string("contact: ") + (contact ? "true\n" : "false\n")
            + "tresca_friction: " + (friction ? "true\n" : "false\n") + R"YAML(
solver: !Petsc
output_stream:
  file: contact_stress_output.pvd
  precision: 17
  format: !vtk
    variant: ascii
output:
  fields: [contact_normal_stress, contact_tangential_stress]
input_fields:
  - region: bulk
    young_modulus: 1
    poisson_ratio: 0.25
  - region: fracture
    young_modulus: 1e-6
    poisson_ratio: 0.25
    friction_threshold: 2
  - region: .BOUNDARY
    bc_type: displacement
    bc_displacement: 0
)YAML";
        auto &type = const_cast<Input::Type::Record &>(Elasticity::get_input_type());
        type.finish();
        Input::ReaderToStorage reader(input, type, Input::FileFormat::format_YAML);
        time_.reset(new TimeGovernor(0.0, 1.0));
        mechanics_.reset(new Elasticity(*mesh_, reader.get_root_interface<Input::Record>(), time_.get()));
        auto constant = [](double value) {
            auto result = std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
            result->set_value(value);
            return result;
        };
        mechanics_->eq_fields().cross_section.set(constant(1.0), 0.0, {"bulk", ".BOUNDARY"});
        mechanics_->eq_fields().cross_section.set(constant(0.01), 0.0, {"fracture"});
        Field<3, FieldValue<3>::Scalar> potential, reference;
        potential.name("test_potential");
        reference.name("test_reference_potential");
        potential.set_mesh(*mesh_); reference.set_mesh(*mesh_);
        potential.set(constant(0.0), 0.0); reference.set(constant(0.0), 0.0);
        potential.set_time(time_->step(), LimitSide::right);
        reference.set_time(time_->step(), LimitSide::right);
        mechanics_->set_potential_load(potential, reference);
        mechanics_->initialize();
        EXPECT_TRUE(mechanics_->eq_fields().output_contact_normal_stress.units() == UnitSI::Pa());
        EXPECT_TRUE(mechanics_->eq_fields().output_contact_tangential_stress.units() == UnitSI::Pa());
        EXPECT_EQ(OutputTime::ELEM_DATA, mechanics_->eq_fields().output_contact_normal_stress.get_output_type());
        EXPECT_EQ(OutputTime::ELEM_DATA, mechanics_->eq_fields().output_contact_tangential_stress.get_output_type());
    }

    /// Place all normal rows on a rank different from an owned fracture to force remote imports.
    void seed(bool friction, bool remote_normal_rows) {
        auto &data = mechanics_->eq_data();
        auto *stored = new StoredMultiplierSolver(data.dh_->distr().get());
        delete data.ls;
        data.ls = stored;
        if (remote_normal_rows) {
            int rank, size, first_owner, local_owner;
            MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
            MPI_Comm_size(PETSC_COMM_WORLD, &size);
            local_owner = size;
            for (auto cell : data.dh_->own_range())
                if (data.constraint_idx.count(cell.elm_idx())) local_owner = rank;
            MPI_Allreduce(&local_owner, &first_owner, 1, MPI_INT, MPI_MIN, PETSC_COMM_WORLD);
            PetscInt rows, columns;
            chkerr(MatGetSize(data.constraint_matrix, &rows, &columns));
            chkerr(MatDestroy(&data.constraint_matrix));
            chkerr(MatCreateAIJ(PETSC_COMM_WORLD, rank == (first_owner + 1) % size ? rows : 0,
                               data.dh_->lsize(), rows, columns, 0, nullptr, 0, nullptr,
                               &data.constraint_matrix));
        }
        Vec normal, tangent = nullptr;
        chkerr(MatCreateVecs(data.constraint_matrix, nullptr, &normal));
        PetscInt begin, end;
        chkerr(VecGetOwnershipRange(normal, &begin, &end));
        PetscScalar *values;
        chkerr(VecGetArray(normal, &values));
        for (PetscInt row = begin; row < end; ++row) values[row - begin] = 6.0 + row;
        chkerr(VecRestoreArray(normal, &values));
        if (friction) {
            chkerr(MatCreateVecs(data.tangential_jump_matrix, nullptr, &tangent));
            chkerr(VecSet(tangent, 6.0));
        }
        stored->store(normal, tangent);
        chkerr(VecDestroy(&normal));
        chkerr(VecDestroy(&tangent));
    }

    /// Check every owned output cell, including cells without the corresponding constraints.
    void expect_stresses(bool contact, bool friction, bool zero = false) {
        auto &data = mechanics_->eq_data();
        auto &fields = mechanics_->eq_fields();
        auto normal = fields.output_contact_normal_stress_ptr;
        auto tangent = fields.output_contact_tangential_stress_ptr;
        for (auto cell : data.dh_->own_range()) {
            auto scalar_cell = cell.cell_with_other_dh(normal->get_dofhandler().get());
            auto vector_cell = cell.cell_with_other_dh(tangent->get_dofhandler().get());
            ASSERT_EQ(1, scalar_cell.n_dofs());
            ASSERT_EQ(3, vector_cell.n_dofs());
            double expected_normal = 0.0;
            auto row = data.constraint_idx.find(cell.elm_idx());
            if (contact && !zero && row != data.constraint_idx.end())
                expected_normal = -(6.0 + row->second) / cell.elm().measure();
            EXPECT_NEAR(expected_normal, normal->vec().get(scalar_cell.get_loc_dof_indices()[0]), 1e-12);
            arma::vec3 expected_tangent;
            expected_tangent.zeros();
            if (friction && !zero && data.tangential_jump_idx.count(cell.elm_idx())) {
                // Use an independent direction for this fixture: (2,1,2), length three.
                expected_tangent = arma::vec3("2 1 2") * (-2.0 / 3.0);
            }
            auto dofs = vector_cell.get_loc_dof_indices();
            for (unsigned int i = 0; i < dofs.size(); ++i)
                EXPECT_NEAR(arma::dot(vector_cell.cell_dof(i).coefs, expected_tangent),
                            tangent->vec().get(dofs[i]), 1e-12);
        }
    }

    std::unique_ptr<Mesh> mesh_; ///< Own the fixture geometry until mechanics is destroyed.
    std::unique_ptr<TimeGovernor> time_; ///< Own the time governor borrowed by mechanics.
    std::unique_ptr<Elasticity> mechanics_; ///< Own initialized production output machinery.
};

TEST_F(ContactStressOutputTest, initializes_zero_stresses_before_first_solve)
{
    initialize(true, true);
    mechanics_->update_output_fields();
    expect_stresses(false, false);
    seed(true, false);
    mechanics_->update_output_fields();
    expect_stresses(true, true);
}

TEST_F(ContactStressOutputTest, converts_original_rows_and_replaces_previous_values)
{
    initialize(true, true);
    seed(true, false);
    mechanics_->update_output_fields();
    expect_stresses(true, true);
    auto *solver = static_cast<LinSys_PERMON *>(mechanics_->eq_data().ls);
    chkerr(VecSet(solver->dual_ineq_multipliers(), 0.0));
    chkerr(VecSet(solver->dual_box_multipliers(), 0.0));
    mechanics_->update_output_fields();
    expect_stresses(true, true, true);
}

TEST_F(ContactStressOutputTest, imports_remote_normal_rows_with_empty_rank_selections)
{
    initialize(true, true);
    seed(true, true);
    mechanics_->update_output_fields();
    expect_stresses(true, true);
}

TEST_F(ContactStressOutputTest, contact_without_friction_has_zero_tangential_output)
{
    initialize(true, false);
    seed(false, false);
    mechanics_->update_output_fields();
    expect_stresses(true, false);
}

TEST_F(ContactStressOutputTest, disabled_contact_has_zero_output)
{
    initialize(false, false);
    mechanics_->update_output_fields();
    expect_stresses(false, false);
}

TEST_F(ContactStressOutputTest, stresses_reproduce_virtual_work_with_bulk_cross_section_weights)
{
    initialize(true, true);
    seed(true, false);
    auto &data = mechanics_->eq_data();
    auto &fields = mechanics_->eq_fields();
    auto section = std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
    section->set_value(2.0);
    fields.cross_section.set(section, 0.0, {"bulk"});
    fields.set_time(time_->step(), LimitSide::right);

    // Assemble the production trace weights independently of the output conversion.
    chkerr(VecSet(data.constraint_vec, 0.0));
    GenericAssembly<Elasticity::ConstraintAssemblyDim> normal_assembly(&data, data.dh_.get());
    GenericAssembly<Elasticity::TangentialJumpAssemblyDim> tangent_assembly(&data, data.dh_.get());
    normal_assembly.assemble(data.dh_);
    tangent_assembly.assemble(data.dh_);
    chkerr(MatAssemblyBegin(data.constraint_matrix, MAT_FINAL_ASSEMBLY));
    chkerr(MatAssemblyEnd(data.constraint_matrix, MAT_FINAL_ASSEMBLY));
    chkerr(MatAssemblyBegin(data.tangential_jump_matrix, MAT_FINAL_ASSEMBLY));
    chkerr(MatAssemblyEnd(data.tangential_jump_matrix, MAT_FINAL_ASSEMBLY));
    chkerr(VecAssemblyBegin(data.constraint_vec));
    chkerr(VecAssemblyEnd(data.constraint_vec));
    mechanics_->update_output_fields();

    auto *solver = static_cast<LinSys_PERMON *>(data.ls);
    Vec expected, actual, work;
    chkerr(MatCreateVecs(data.constraint_matrix, &expected, nullptr));
    chkerr(VecDuplicate(expected, &actual)); chkerr(VecDuplicate(expected, &work));
    chkerr(MatMultTranspose(data.constraint_matrix, solver->dual_ineq_multipliers(), expected));
    chkerr(MatMultTranspose(data.tangential_jump_matrix, solver->dual_box_multipliers(), work));
    chkerr(VecAXPY(expected, 1.0, work));
    chkerr(VecScale(expected, -1.0));
    chkerr(VecSet(actual, 0.0));
    PhysicalReactionData physical_data{data.eq_fields_, actual};
    GenericAssembly<PhysicalReactionDim> physical_assembly(&physical_data, data.dh_.get());
    physical_assembly.assemble(data.dh_);
    chkerr(VecAssemblyBegin(actual)); chkerr(VecAssemblyEnd(actual));
    PetscReal norm;
    chkerr(VecNorm(expected, NORM_2, &norm));
    EXPECT_GT(norm, 1.0);
    chkerr(VecAXPY(actual, -1.0, expected));
    chkerr(VecNorm(actual, NORM_2, &norm));
    EXPECT_LT(norm, 1e-10);
    chkerr(VecDestroy(&expected)); chkerr(VecDestroy(&actual)); chkerr(VecDestroy(&work));
}

TEST_F(ContactStressOutputTest, writes_scalar_and_vector_cell_arrays_with_field_values)
{
    initialize(true, true);
    seed(true, false);
    mechanics_->output_data();
    auto &data = mechanics_->eq_data();
    auto &fields = mechanics_->eq_fields();
    std::vector<double> local_records;
    for (auto cell : data.dh_->own_range()) {
        auto scalar = cell.cell_with_other_dh(fields.output_contact_normal_stress_ptr->get_dofhandler().get());
        auto vector = cell.cell_with_other_dh(fields.output_contact_tangential_stress_ptr->get_dofhandler().get());
        auto centre = cell.elm().centre();
        for (unsigned int i = 0; i < 3; ++i) local_records.push_back(centre[i]);
        local_records.push_back(fields.output_contact_normal_stress_ptr->vec().get(scalar.get_loc_dof_indices()[0]));
        arma::vec3 traction;
        traction.zeros();
        auto dofs = vector.get_loc_dof_indices();
        for (unsigned int i = 0; i < dofs.size(); ++i)
            traction += vector.cell_dof(i).coefs * fields.output_contact_tangential_stress_ptr->vec().get(dofs[i]);
        for (unsigned int i = 0; i < 3; ++i) local_records.push_back(traction[i]);
    }
    auto gather = [](const std::vector<double> &local) {
        int size;
        MPI_Comm_size(PETSC_COMM_WORLD, &size);
        int count = local.size();
        std::vector<int> counts(size), offsets(size, 0);
        MPI_Allgather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, PETSC_COMM_WORLD);
        for (int i = 1; i < size; ++i) offsets[i] = offsets[i - 1] + counts[i - 1];
        std::vector<double> result(offsets.back() + counts.back());
        MPI_Allgatherv(local.data(), count, MPI_DOUBLE, result.data(), counts.data(), offsets.data(),
                       MPI_DOUBLE, PETSC_COMM_WORLD);
        return result;
    };
    auto expected = gather(local_records);
    mechanics_.reset(); // Flush the final time frame and close the serial VTK file.
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank != 0) return;
    FilePath path("contact_stress_output/contact_stress_output-000000.vtu", FilePath::output_file);
    std::ifstream input(static_cast<std::string>(path));
    ASSERT_TRUE(input.good());
    boost::property_tree::ptree xml;
    ASSERT_NO_THROW(boost::property_tree::read_xml(input, xml));
    std::map<std::string, std::vector<double>> arrays;
    for (const auto &array : xml.get_child("VTKFile.UnstructuredGrid.Piece.CellData")) {
        if (array.first != "DataArray") continue;
        std::string name = array.second.get<std::string>("<xmlattr>.Name");
        if (name != "contact_normal_stress" && name != "contact_tangential_stress") continue;
        EXPECT_EQ(name == "contact_normal_stress" ? 1 : 3,
                  array.second.get<int>("<xmlattr>.NumberOfComponents", 1));
        std::istringstream values(array.second.data());
        double value;
        while (values >> value) arrays[name].push_back(value);
    }
    const std::size_t n_cells = expected.size() / 7;
    ASSERT_EQ(n_cells, arrays["contact_normal_stress"].size());
    ASSERT_EQ(3 * n_cells, arrays["contact_tangential_stress"].size());

    // Identify each written cell by geometry rather than assuming MPI/output element order.
    auto numbers = [](const boost::property_tree::ptree &node) {
        std::vector<double> result;
        std::istringstream values(node.data());
        double value;
        while (values >> value) result.push_back(value);
        return result;
    };
    auto points = numbers(xml.get_child("VTKFile.UnstructuredGrid.Piece.Points.DataArray"));
    std::vector<double> connectivity, offsets;
    for (const auto &array : xml.get_child("VTKFile.UnstructuredGrid.Piece.Cells")) {
        if (array.first != "DataArray") continue;
        auto name = array.second.get<std::string>("<xmlattr>.Name");
        if (name == "connectivity") connectivity = numbers(array.second);
        if (name == "offsets") offsets = numbers(array.second);
    }
    ASSERT_EQ(n_cells, offsets.size());
    std::size_t first = 0;
    std::vector<bool> seen(n_cells, false);
    for (std::size_t cell = 0; cell < n_cells; ++cell) {
        std::size_t end = offsets[cell];
        ASSERT_GT(end, first);
        ASSERT_LE(end, connectivity.size());
        arma::vec3 centre;
        centre.zeros();
        for (std::size_t i = first; i < end; ++i) {
            std::size_t node = connectivity[i];
            ASSERT_LT(3 * node + 2, points.size());
            for (unsigned int j = 0; j < 3; ++j) centre[j] += points[3 * node + j];
        }
        centre /= end - first;
        first = end;
        std::size_t match = n_cells;
        for (std::size_t i = 0; i < n_cells; ++i) {
            arma::vec3 reference;
            for (unsigned int j = 0; j < 3; ++j) reference[j] = expected[7 * i + j];
            if (arma::norm(centre - reference, 2) < 1e-12) match = i;
        }
        ASSERT_LT(match, n_cells);
        EXPECT_FALSE(seen[match]);
        seen[match] = true;
        EXPECT_NEAR(expected[7 * match + 3], arrays["contact_normal_stress"][cell], 1e-12);
        for (unsigned int j = 0; j < 3; ++j)
            EXPECT_NEAR(expected[7 * match + 4 + j], arrays["contact_tangential_stress"][3 * cell + j], 1e-12);
    }
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
