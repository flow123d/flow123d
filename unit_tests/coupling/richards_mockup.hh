#ifndef RICHARDS_MOCKUP_HH_
#define RICHARDS_MOCKUP_HH_


#include <mesh_constructor.hh>
#include "arma_expect.hh"
#include <rev_num.h>

#include "fem/eval_points.hh"
#include "fem/eval_subset.hh"
#include "fem/element_cache_map.hh"
#include "fields/field_values.hh"
#include "fields/field_set.hh"
#include "fields/field_fe.hh"
#include "fields/generic_field.hh"
#include "fields/multi_field.hh"
#include "fields/bc_multi_field.hh"
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
#include "flow/assembly_richards.hh"
#include "flow/assembly_models.hh"
#include "richards_mockup_assembly.hh"
//#include "darcy_mockup.hh"
#include "darcy_mockup.impl.hh"
#include "balance_null.hh"

class GenericAssemblyBase;
template< template<IntDim...> class DimAssembly> class GenericAssembly;



namespace it = Input::Type;

namespace equation_data_richards {

class EqFields : public equation_data::EqFields {
public:
    /**
     * For compatibility with old BCD file we have to assign integer codes starting from 1.
     */
    enum BC_Type {
        none=0,
        dirichlet=1,
        total_flux=4,
        seepage=5,
        river=6
    };

    /// Return a Selection corresponding to enum BC_Type.
    static const Input::Type::Selection & get_bc_type_selection() {
    	return it::Selection("Flow_Darcy_BC_Type")
            .add_value(none, "none",
                "Homogeneous Neumann boundary condition\n(zero normal flux over the boundary).")
            .add_value(dirichlet, "dirichlet",
                "Dirichlet boundary condition. "
                "Specify the pressure head through the ``bc_pressure`` field "
                "or the piezometric head through the ``bc_piezo_head`` field.")
            .add_value(total_flux, "total_flux", "Flux boundary condition (combines Neumann and Robin type). "
                "Water inflow equal to (($ \\delta_d(q_d^N + \\sigma_d (h_d^R - h_d) )$)). "
                "Specify the water inflow by the ``bc_flux`` field, the transition coefficient by ``bc_robin_sigma`` "
                "and the reference pressure head or piezometric head through ``bc_pressure`` or ``bc_piezo_head`` respectively.")
            .add_value(seepage, "seepage",
                "Seepage face boundary condition. Pressure and inflow bounded from above. Boundary with potential seepage flow "
                "is described by the pair of inequalities: "
                "(($h_d \\le h_d^D$)) and (($ -\\boldsymbol q_d\\cdot\\boldsymbol n \\le \\delta q_d^N$)), where the equality holds in at least one of them. "
                "Caution: setting (($q_d^N$)) strictly negative "
                "may lead to an ill posed problem since a positive outflow is enforced. "
                "Parameters (($h_d^D$)) and (($q_d^N$)) are given by the fields ``bc_switch_pressure`` (or ``bc_switch_piezo_head``) and ``bc_flux`` respectively."
                )
            .add_value(river, "river",
                "River boundary condition. For the water level above the bedrock, (($H_d > H_d^S$)), the Robin boundary condition is used with the inflow given by: "
                "(( $ \\delta_d(q_d^N + \\sigma_d(H_d^D - H_d) )$)). For the water level under the bedrock, constant infiltration is used: "
                "(( $ \\delta_d(q_d^N + \\sigma_d(H_d^D - H_d^S) )$)). Parameters: ``bc_pressure``, ``bc_switch_pressure``, "
                " ``bc_sigma``, ``bc_flux``."
                )
            .close();
    }

    /// Creation of all fields.
    EqFields() {
        *this += water_content.name("water_content")
                .units(UnitSI::dimensionless())
                .flags(FieldFlag::equation_result)
                .description(R"(Water content.
                    It is a fraction of water volume to the whole volume.)");
        *this += conductivity_richards.name("conductivity_richards")
                .units( UnitSI().m().s(-1) )
                .flags(FieldFlag::equation_result)
                .description("Computed isotropic scalar conductivity by the soil model.");

        *this += water_content_saturated.name("water_content_saturated")
                .description(R"(Saturated water content (($ \theta_s $)).
                    Relative volume of water in a reference volume of a saturated porous media.)")
                .input_default("0.0")
                .units( UnitSI::dimensionless() );

        *this += water_content_residual.name("water_content_residual")
                .description(R"(Residual water content (($ \theta_r $)).
                    Relative volume of water in a reference volume of an ideally dry porous media.)")
                .input_default("0.0")
                .units( UnitSI::dimensionless() );

        *this += genuchten_p_head_scale.name("genuchten_p_head_scale")
                .description(R"(The van Genuchten pressure head scaling parameter (($ \alpha $)).
                    It is related to the inverse of the air entry pressure, i.e. the pressure
                    where the relative water content starts to decrease below 1.)")
                .input_default("0.0")
                .units( UnitSI().m(-1) );

        *this += genuchten_n_exponent.name("genuchten_n_exponent")
                .description("The van Genuchten exponent parameter (($ n $)).")
                .input_default("2.0")
                .units( UnitSI::dimensionless() );

        this->set_default_fieldset();
    }

    // input fields
    Field<3, FieldValue<3>::Scalar > water_content_saturated;   // corresponds to the porosity (theta_s = Vw/V = porosity)
    Field<3, FieldValue<3>::Scalar > water_content_residual;
    Field<3, FieldValue<3>::Scalar > genuchten_p_head_scale;
    Field<3, FieldValue<3>::Scalar > genuchten_n_exponent;

    //output fields
    Field<3, FieldValue<3>::Scalar > water_content;
    std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> water_content_ptr;

    Field<3, FieldValue<3>::Scalar > conductivity_richards;
//         FieldFE<3, FieldValue<3>::Scalar > conductivity_richards;
    std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> conductivity_ptr;
};

class EqData : public equation_data::EqData {
public:
    typedef equation_data_richards::EqFields EqFields;

    EqData(std::shared_ptr<EqFields> eq_fields)
    : equation_data::EqData::EqData(eq_fields), eq_fields_(eq_fields) {}

	/// Shared pointer of EqFields
	std::shared_ptr<EqFields> eq_fields_;

    // Auxiliary assembly fields.
    VectorMPI water_content_previous_time;
    VectorMPI capacity;

    std::shared_ptr<SoilModelBase> soil_model_;
};

} // end of namespace equation_data_richards


/// Test class
class RichardsMockupTest : public testing::Test {
public:
    template<unsigned int dim> using MHMatrixAssemblyRichardsDim = MHMatrixAssemblyRichards<dim, equation_data_richards::EqData>;
    template<unsigned int dim> using MHMatrixEvalFieldsRichardsDim = MHMatrixEvalFieldsRichards<dim, equation_data_richards::EqData>;

    RichardsMockupTest()
    {
		string root_dir=string(UNIT_TESTS_BIN_DIR) + "/coupling";
		string build = string(__DATE__) + ", " + string(__TIME__)
	            + " flags: (unknown compiler flags)";

        FilePath::set_io_dirs(".",root_dir,"",".");
        Profiler::instance();
        Profiler::instance()->set_program_info("Flow123d",
                string(FLOW123D_VERSION_NAME_), string(FLOW123D_GIT_BRANCH_), string(FLOW123D_GIT_REVISION_), build);
        Profiler::set_memory_monitoring(false);
    }

    ~RichardsMockupTest() {}

    /// Run assembly algorithms with different type of assembly and type of field
    void run_fullassembly_const(const string &eq_data_input, const std::string &mesh_file);
    void run_fullassembly_model(const string &eq_data_input, const std::string &mesh_file);
    void run_computelocal_const(const string &eq_data_input, const std::string &mesh_file);
    void run_computelocal_model(const string &eq_data_input, const std::string &mesh_file);
    void run_evalfields_const(const string &eq_data_input, const std::string &mesh_file);
    void run_evalfields_model(const string &eq_data_input, const std::string &mesh_file);

	/// Perform profiler output.
    void profiler_output(std::string file_name) {
		FilePath fp(file_name + "_profiler.json", FilePath::output_file);
		Profiler::instance()->output(MPI_COMM_WORLD, fp.filename());
	}
};


/*******************************************************************************
 * Equivalent to TransportDG class
 */
template<template<IntDim...> class MhMatrix>
class RichardsMockup : public DarcyMockup<MhMatrix> {
public:
    it::Record & get_input_type() {
        namespace it=Input::Type;
        std::string equation_name = "TestEquation";

        it::Record field_descriptor = it::Record(equation_name + "_Data",
            FieldCommon::field_descriptor_record_description(equation_name + "_Data"))
        .copy_keys( DarcyMockup<MhMatrix>::type_field_descriptor() )
        .copy_keys( equation_data_richards::EqFields().make_field_descriptor_type(equation_name + "_Data_aux") )
        .close();

        auto model_selection = it::Selection("Soil_Model_Type", "")
                .add_value(SoilModelBase::van_genuchten, "van_genuchten", "Van Genuchten soil model with cutting near zero.")
                .add_value(SoilModelBase::irmay, "irmay", "Irmay model for conductivity, Van Genuchten model for the water content. Suitable for bentonite.")
                .close();

        auto soil_rec = it::Record("SoilModel", "Soil model settings.")
            .allow_auto_conversion("model_type")
            .declare_key("model_type", model_selection, it::Default("\"van_genuchten\""),
                "Selection of the globally applied soil model. In future we replace this key by a field for selection of the model."
                "That will allow usage of different soil model in a single simulation.")
            .declare_key("cut_fraction", it::Double(0.0,1.0), it::Default("0.999"),
                    "Fraction of the water content where we cut  and rescale the curve.")
            .close();

        equation_data_richards::EqFields eq_fields;

        return it::Record(equation_name, "Lumped Mixed-Hybrid solver for unsteady unsaturated Darcy flow.")
            .derive_from(DarcyFlowInterface::get_input_type())
            .copy_keys(DarcyMockup<MhMatrix>::get_input_type())
            .declare_key("input_fields", it::Array( field_descriptor ), it::Default::obligatory(),
                    "Input data for Darcy flow model.")
//            .declare_key("output", DarcyFlowMHOutput::get_input_type(eq_fields, equation_name()),
//                    IT::Default("{ \"fields\": [ \"pressure_p0\", \"velocity_p0\" ] }"),
//                    "Specification of output fields and output times.")
            .declare_key("soil_model", soil_rec, it::Default("\"van_genuchten\""),
                    "Soil model settings.")
            .close();
    }

    RichardsMockup(bool use_linsys)
    : DarcyMockup<MhMatrix>(use_linsys)
    {
        eq_fields_ = make_shared<equation_data_richards::EqFields>();
        eq_data_ = make_shared<equation_data_richards::EqData>(eq_fields_);
        DarcyMockup<MhMatrix>::eq_data_ = eq_data_;
        DarcyMockup<MhMatrix>::eq_fields_ = eq_fields_;
        this->eq_fieldset_ = eq_fields_;
    }

    ~RichardsMockup() {}

    bool zero_time_term(bool time_global=false) {
        if (time_global) {
            return (eq_fields_->storativity.input_list_size() == 0)
                    && (eq_fields_->water_content_saturated.input_list_size() == 0);

        } else {
            return (eq_fields_->storativity.field_result(this->mesh_->region_db().get_region_set("BULK"))
                    == result_zeros)
                    && (eq_fields_->water_content_saturated.field_result(this->mesh_->region_db().get_region_set("BULK"))
                    == result_zeros);
        }
    }

    void accept_time_step() override {
    	eq_data_->p_edge_solution_previous_time.copy_from(eq_data_->p_edge_solution);
        VectorMPI water_content_vec = eq_fields_->water_content_ptr->vec();
        eq_data_->water_content_previous_time.copy_from(water_content_vec);

        eq_data_->p_edge_solution_previous_time.local_to_ghost_begin();
        eq_data_->p_edge_solution_previous_time.local_to_ghost_end();
    }

    /// Assembly or update whole linear system.
    void assembly_linear_system() override;

    /// Create and initialize assembly objects
    void initialize_asm() override;



	std::shared_ptr<equation_data_richards::EqFields> eq_fields_;
	std::shared_ptr<equation_data_richards::EqData> eq_data_;
};


#endif /* RICHARDS_MOCKUP_HH_ */
