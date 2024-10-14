#ifndef ELASTICITY_MOCKUP_HH_
#define ELASTICITY_MOCKUP_HH_


#include <mesh_constructor.hh>
#include "arma_expect.hh"
#include <rev_num.h>

#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fields/field_value_cache.hh"
#include "fields/field_values.hh"
#include "fields/field_set.hh"
#include "fields/field_fe.hh"
#include "fields/generic_field.hh"
#include "fields/bc_field.hh"
#include "fields/equation_output.hh"
#include "fields/field_model.hh"
#include "fields/field_constant.hh"
#include "coupling/equation.hh"
#include "tools/unit_si.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "fem/fe_p.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "la/linsys.hh"
#include "la/linsys_PETSC.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"


class ElasticityMockupTest : public testing::Test {
public:
	ElasticityMockupTest()
    {
		string root_dir=string(UNIT_TESTS_BIN_DIR) + "/coupling";
		string build = string(__DATE__) + ", " + string(__TIME__)
	            + " flags: (unknown compiler flags)";

        FilePath::set_io_dirs(".",root_dir,"",".");
        Profiler::instance();
        Profiler::instance()->set_program_info("Flow123d",
                string(FLOW123D_VERSION_NAME_), string(FLOW123D_GIT_BRANCH_), string(FLOW123D_GIT_REVISION_), build);
        Profiler::set_memory_monitoring(false, false);
    }

    ~ElasticityMockupTest() {}

    /// Perform profiler output.
    void profiler_output(std::string file_name) {
		FilePath fp(file_name + "_profiler.json", FilePath::output_file);
		Profiler::instance()->output(MPI_COMM_WORLD, fp.filename());
	}
};

namespace equation_data {

class EqFields : public FieldSet {
public:
    enum Bc_types {
      bc_type_displacement,
      bc_type_displacement_normal,
      bc_type_traction,
	  bc_type_stress,
    };

    static const Input::Type::Selection & get_bc_type_selection() {
        return Input::Type::Selection("Elasticity_BC_Type", "Types of boundary conditions for mechanics.")
                .add_value(bc_type_displacement, "displacement",
                      "Prescribed displacement.")
                .add_value(bc_type_displacement_normal, "displacement_n",
                      "Prescribed displacement in the normal direction to the boundary.")
                .add_value(bc_type_traction, "traction",
                      "Prescribed traction.")
                .add_value(bc_type_stress, "stress",
                      "Prescribed stress tensor.")
                .close();
    }

    /// Constructor
    EqFields()
    {
        *this+=bc_type
            .name("bc_type")
            .description(
            "Type of boundary condition.")
            .units( UnitSI::dimensionless() )
            .input_default("\"traction\"")
            .input_selection( get_bc_type_selection() )
            .flags_add(FieldFlag::in_rhs & FieldFlag::in_main_matrix);

        *this+=bc_displacement
            .name("bc_displacement")
            .description("Prescribed displacement on boundary.")
            .units( UnitSI().m() )
            .input_default("0.0")
            .flags_add(in_rhs);

        *this+=bc_traction
            .name("bc_traction")
            .description("Prescribed traction on boundary.")
            .units( UnitSI().Pa() )
            .input_default("0.0")
            .flags_add(in_rhs);

        *this+=bc_stress
            .name("bc_stress")
            .description("Prescribed stress on boundary.")
            .units( UnitSI().Pa() )
            .input_default("0.0")
            .flags_add(in_rhs);

        *this+=load
            .name("load")
            .description("Prescribed bulk load.")
            .units( UnitSI().N().m(-3) )
            .input_default("0.0")
            .flags_add(in_rhs);

        *this+=young_modulus
            .name("young_modulus")
            .description("Young's modulus.")
            .units( UnitSI().Pa() )
             .input_default("0.0")
            .flags_add(in_main_matrix & in_rhs);

        *this+=poisson_ratio
            .name("poisson_ratio")
            .description("Poisson's ratio.")
            .units( UnitSI().dimensionless() )
             .input_default("0.0")
            .flags_add(in_main_matrix & in_rhs);

        *this+=fracture_sigma
                .name("fracture_sigma")
                .description(
                "Coefficient of transfer of forces through fractures.")
                .units( UnitSI::dimensionless() )
                .input_default("1.0")
                .flags_add(in_main_matrix & in_rhs);

        *this+=initial_stress
            .name("initial_stress")
            .description("Initial stress tensor.")
            .units( UnitSI().Pa() )
            .input_default("0.0")
            .flags_add(in_rhs);

        *this += region_id.name("region_id")
        	        .units( UnitSI::dimensionless())
        	        .flags(FieldFlag::equation_external_output);

        *this += subdomain.name("subdomain")
          .units( UnitSI::dimensionless() )
          .flags(FieldFlag::equation_external_output);

        *this+=cross_section
          .name("cross_section")
          .units( UnitSI().m(3).md() )
          .flags(input_copy & in_time_term & in_main_matrix & in_rhs);

        *this+=cross_section_min
          .name("cross_section_min")
          .description("Minimal cross-section of fractures.")
          .units( UnitSI().m(3).md() )
          .input_default("0.0");

        *this+=potential_load
          .name("potential_load")
          .units( UnitSI().m() )
          .flags(input_copy & in_rhs);

        *this+=output_field
                .name("displacement")
                .description("Displacement vector field output.")
                .units( UnitSI().m() )
                .flags(equation_result);

        *this += output_stress
                .name("stress")
                .description("Stress tensor output.")
                .units( UnitSI().Pa() )
                .flags(equation_result);

        *this += output_von_mises_stress
                .name("von_mises_stress")
                .description("von Mises stress output.")
                .units( UnitSI().Pa() )
                .flags(equation_result);

        *this += output_mean_stress
                .name("mean_stress")
                .description("mean stress output.")
                .units( UnitSI().Pa() )
                .flags(equation_result);

        *this += output_cross_section
                .name("cross_section_updated")
                .description("Cross-section after deformation - output.")
                .units( UnitSI().m() )
                .flags(equation_result);

        *this += output_divergence
                .name("displacement_divergence")
                .description("Displacement divergence output.")
                .units( UnitSI().dimensionless() )
                .flags(equation_result);

        *this += lame_mu.name("lame_mu")
                .description("Field lame_mu.")
                .input_default("0.0")
                .units( UnitSI().Pa() );

        *this += lame_lambda.name("lame_lambda")
                .description("Field lame_lambda.")
                .input_default("0.0")
                .units( UnitSI().Pa() );

        *this += dirichlet_penalty.name("dirichlet_penalty")
                .description("Field dirichlet_penalty.")
                .input_default("0.0")
                .units( UnitSI().Pa() );

        // add all input fields to the output list
        output_fields += *this;

        this->add_coords_field();
        this->set_default_fieldset();
    }


    // Definition of Fields
    BCField<3, FieldValue<3>::Enum > bc_type;
    BCField<3, FieldValue<3>::VectorFixed> bc_displacement;
    BCField<3, FieldValue<3>::VectorFixed> bc_traction;
    BCField<3, FieldValue<3>::TensorFixed> bc_stress;
    Field<3, FieldValue<3>::VectorFixed> load;
    Field<3, FieldValue<3>::Scalar> young_modulus;
    Field<3, FieldValue<3>::Scalar> poisson_ratio;
    Field<3, FieldValue<3>::Scalar> fracture_sigma;    ///< Transition parameter for diffusive transfer on fractures.
    Field<3, FieldValue<3>::TensorFixed> initial_stress;

    /// Pointer to DarcyFlow field cross_section
    Field<3, FieldValue<3>::Scalar > cross_section;
    Field<3, FieldValue<3>::Scalar > cross_section_min;
    Field<3, FieldValue<3>::Scalar > potential_load;   ///< Potential of an additional (external) load.
    Field<3, FieldValue<3>::Scalar > ref_potential_load; ///< Potential of reference external load on boundary. TODO: Switch to BCField when possible.
    Field<3, FieldValue<3>::Scalar> region_id;
    Field<3, FieldValue<3>::Scalar> subdomain;

    Field<3, FieldValue<3>::VectorFixed> output_field;
    Field<3, FieldValue<3>::TensorFixed> output_stress;
    Field<3, FieldValue<3>::Scalar> output_von_mises_stress;
    Field<3, FieldValue<3>::Scalar> output_mean_stress;
    Field<3, FieldValue<3>::Scalar> output_cross_section;
    Field<3, FieldValue<3>::Scalar> output_divergence;

    /// @name Instances of FieldModel used in assembly methods
    // @{

    Field<3, FieldValue<3>::Scalar > lame_mu;
    Field<3, FieldValue<3>::Scalar > lame_lambda;
    Field<3, FieldValue<3>::Scalar > dirichlet_penalty;

    // @}

    std::shared_ptr<FieldFE<3, FieldValue<3>::VectorFixed> > output_field_ptr;
    std::shared_ptr<FieldFE<3, FieldValue<3>::TensorFixed> > output_stress_ptr;
    std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > output_von_mises_stress_ptr;
    std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > output_mean_stress_ptr;
    std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > output_cross_section_ptr;
    std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > output_div_ptr;

    EquationOutput output_fields;

};

class EqData {
public:
    EqData()
    : ls(nullptr) {}

    ~EqData() {
        if (ls!=nullptr) delete ls;
    }

    /// Create DOF handler objects
    void create_dh(Mesh * mesh, unsigned int fe_order) {
        ASSERT_EQ(fe_order, 1)(fe_order).error("Unsupported polynomial order for finite elements in Elasticity");
        MixedPtr<FE_P> fe_p(1);
        MixedPtr<FiniteElement> fe = mixed_fe_system(fe_p, FEVector, 3);

        std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe);
        dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh);
        dh_->distribute_dofs(ds);
    }

    /// Object for distribution of dofs.
    std::shared_ptr<DOFHandlerMultiDim> dh_;

    /// Linear algebraic system.
    LinSys *ls;

};

} // end of namespace equation_data

#endif /* ELASTICITY_MOCKUP_HH_ */
