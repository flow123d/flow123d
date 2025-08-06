/*
 * richards_lmh.cc
 *
 *  Created on: Sep 16, 2015
 *      Author: jb
 */


#include "system/global_defs.h"
#include "system/sys_profiler.hh"
#include "system/asserts.hh"

#include "coupling/balance.hh"

#include "input/input_type.hh"
#include "input/factory.hh"
#include "flow/richards_lmh.hh"
#include "flow/soil_models.hh"
#include "flow/darcy_flow_mh_output.hh"
#include "flow/assembly_richards.hh"
#include "tools/time_governor.hh"

#include "petscmat.h"
#include "petscviewer.h"
#include "petscerror.h"
#include <armadillo>

#include "la/schur.hh"
#include "la/vector_mpi.hh"


#include "tools/include_fadbad.hh" // for "fadbad.h", "badiff.h", "fadiff.h"


FLOW123D_FORCE_LINK_IN_CHILD(richards_lmh)



RichardsLMH::EqFields::EqFields()
{
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


RichardsLMH::EqData::EqData()
: DarcyLMH::EqData::EqData() {}


const Input::Type::Record & RichardsLMH::get_input_type() {
	namespace it=Input::Type;

	it::Record field_descriptor = it::Record(equation_name() + "_Data",
        FieldCommon::field_descriptor_record_description(equation_name() + "_Data"))
    .copy_keys( DarcyLMH::type_field_descriptor() )
    .copy_keys( RichardsLMH::EqFields().make_field_descriptor_type(equation_name() + "_Data_aux") )
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

    RichardsLMH::EqFields eq_fields;
    
    return it::Record(equation_name(), "Lumped Mixed-Hybrid solver for unsteady unsaturated Darcy flow.")
        .derive_from(DarcyFlowInterface::get_input_type())
        .copy_keys(DarcyLMH::get_input_type())
        .declare_key("input_fields", it::Array( field_descriptor ), it::Default::obligatory(),
                "Input data for Darcy flow model.")
        .declare_key("output", DarcyFlowMHOutput::get_input_type(eq_fields, equation_name()),
                IT::Default("{ \"fields\": [ \"pressure_p0\", \"velocity_p0\" ] }"),
                "Specification of output fields and output times.")
        .declare_key("soil_model", soil_rec, it::Default("\"van_genuchten\""),
                "Soil model settings.")
        .close();
}


const int RichardsLMH::registrar =
        Input::register_class< RichardsLMH, Mesh &, const Input::Record >(equation_name()) +
        RichardsLMH::get_input_type().size();



RichardsLMH::RichardsLMH(Mesh &mesh_in, const  Input::Record in_rec, TimeGovernor *tm)
    : DarcyLMH(mesh_in, in_rec, tm),
    init_cond_postprocess_assembly_(nullptr)
{
    eq_fields_ = make_shared<EqFields>();
    eq_data_ = make_shared<EqData>();
    DarcyLMH::eq_data_ = eq_data_;
    DarcyLMH::eq_fields_ = eq_fields_;
    this->eq_fieldset_ = eq_fields_;
    //eq_data_->edge_new_local_4_mesh_idx_ = &(this->edge_new_local_4_mesh_idx_);

    eq_fields_->set_mesh(*mesh_);
}


void RichardsLMH::initialize_specific() {

    auto model_rec = input_record_.val<Input::Record>("soil_model");
    auto model_type = model_rec.val<SoilModelBase::SoilModelType>("model_type");
    double fraction= model_rec.val<double>("cut_fraction");
    if (model_type == SoilModelBase::van_genuchten)
        eq_data_->soil_model_ = std::make_shared<SoilModel_VanGenuchten>(fraction);
    else if (model_type == SoilModelBase::irmay)
        eq_data_->soil_model_ = std::make_shared<SoilModel_Irmay>(fraction);
    else
        ASSERT_PERMANENT(false);

    // create edge vectors
    eq_data_->water_content_previous_time = eq_data_->dh_cr_disc_->create_vector();
    eq_data_->capacity = eq_data_->dh_cr_disc_->create_vector();

    ASSERT_PTR(mesh_);
    eq_data_->mesh = mesh_;

    eq_fields_->water_content_ptr = create_field_fe< 3, FieldValue<3>::Scalar >(eq_data_->dh_cr_disc_);
    eq_fields_->water_content.set(eq_fields_->water_content_ptr, 0.0);
    
    eq_fields_->conductivity_ptr = create_field_fe< 3, FieldValue<3>::Scalar >(eq_data_->dh_p_);
    eq_fields_->conductivity_richards.set(eq_fields_->conductivity_ptr, 0.0);


    //eq_data_->multidim_assembler = AssemblyFlowBase::create< AssemblyRichards >(eq_fields_, eq_data_);
}


//void RichardsLMH::initial_condition_postprocess()
//{
//    // modify side fluxes in parallel
//    // for every local edge take time term on diagonal and add it to the corresponding flux
//
//    for ( DHCellAccessor dh_cell : eq_data_->dh_->own_range() ) {
//    	eq_data_->multidim_assembler[dh_cell.elm().dim()-1]->update_water_content(dh_cell);
//    }
//}


void RichardsLMH::accept_time_step()
{
	eq_data_->p_edge_solution_previous_time.copy_from(eq_data_->p_edge_solution);
    VectorMPI water_content_vec = eq_fields_->water_content_ptr->vec();
    eq_data_->water_content_previous_time.copy_from(water_content_vec);

    eq_data_->p_edge_solution_previous_time.local_to_ghost_begin();
    eq_data_->p_edge_solution_previous_time.local_to_ghost_end();
}

bool RichardsLMH::zero_time_term(bool time_global) {
    if (time_global) {
        return (eq_fields_->storativity.input_list_size() == 0)
                && (eq_fields_->water_content_saturated.input_list_size() == 0);

    } else {
        return (eq_fields_->storativity.field_result(mesh_->region_db().get_region_set("BULK"))
                == result_zeros)
                && (eq_fields_->water_content_saturated.field_result(mesh_->region_db().get_region_set("BULK"))
                == result_zeros);
    }
}


void RichardsLMH::assembly_linear_system()
{
    START_TIMER("RicharsLMH::assembly_linear_system");

    eq_data_->p_edge_solution.local_to_ghost_begin();
    eq_data_->p_edge_solution.local_to_ghost_end();

    eq_data_->is_linear = eq_fields_->genuchten_p_head_scale.field_result(mesh_->region_db().get_region_set("BULK")) == result_zeros;

    //DebugOut() << "Assembly linear system\n";
        START_TIMER("full assembly");
//         if (typeid(*schur0) != typeid(LinSys_BDDC)) {
//             schur0->start_add_assembly(); // finish allocation and create matrix
//             schur_compl->start_add_assembly();
//         }
        
        lin_sys_schur().start_add_assembly();
            
        eq_data_->time_step_ = time_->dt();
        
        lin_sys_schur().mat_zero_entries();
        lin_sys_schur().rhs_zero_entries();

//        assembly_mh_matrix( eq_data_->multidim_assembler ); // fill matrix
        START_TIMER("RichardsLMH::assembly_steady_mh_matrix");
        this->mh_matrix_assembly_->assemble(eq_data_->dh_); // fill matrix
        END_TIMER("RichardsLMH::assembly_steady_mh_matrix");

        lin_sys_schur().finish_assembly();
        lin_sys_schur().set_matrix_changed();
}


void RichardsLMH::initialize_asm() {
    this->read_init_cond_assembly_ = new GenericAssembly< ReadInitCondAssemblyLMH >(eq_fields_.get(), eq_data_.get());
    this->init_cond_postprocess_assembly_ = new GenericAssembly< InitCondPostprocessAssembly >(this->eq_fields_.get(), this->eq_data_.get());
    this->mh_matrix_assembly_ = new GenericAssembly< MHMatrixAssemblyRichards >(this->eq_fields_.get(), this->eq_data_.get());
    this->reconstruct_schur_assembly_ = new GenericAssembly< ReconstructSchurAssemblyRichards >(this->eq_fields_.get(), this->eq_data_.get());
}


void RichardsLMH::read_init_cond_asm() {
    this->read_init_cond_assembly_->assemble(eq_data_->dh_cr_);
    this->init_cond_postprocess_assembly_->assemble(eq_data_->dh_cr_);
}


RichardsLMH::~RichardsLMH() {
    if (init_cond_postprocess_assembly_!=nullptr) {
        delete init_cond_postprocess_assembly_;
        init_cond_postprocess_assembly_ = nullptr;
    }
}
