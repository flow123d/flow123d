/*
 * richards_lmh.cc
 *
 *  Created on: Sep 16, 2015
 *      Author: jb
 */


#include "system/global_defs.h"
#include "system/sys_profiler.hh"

#include "coupling/balance.hh"

#include "input/input_type.hh"
#include "input/factory.hh"
#include "flow/richards_lmh.hh"
#include "flow/soil_models.hh"
#include "flow/assembly_richards.hh"
#include "flow/darcy_flow_mh_output.hh"
#include "tools/time_governor.hh"

#include "petscmat.h"
#include "petscviewer.h"
#include "petscerror.h"
#include <armadillo>

#include "la/schur.hh"
#include "la/vector_mpi.hh"

// in the third_party/FADBAD++ dir, namespace "fadbad"
#include "fadbad.h"
#include "badiff.h"
#include "fadiff.h"



FLOW123D_FORCE_LINK_IN_CHILD(richards_lmh);


namespace it=Input::Type;


RichardsLMH::EqData::EqData()
{
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
}


const it::Record & RichardsLMH::get_input_type() {
    it::Record field_descriptor = it::Record("RichardsLMH_Data",FieldCommon::field_descriptor_record_description("RichardsLMH_Data"))
    .copy_keys( DarcyLMH::type_field_descriptor() )
    .copy_keys( RichardsLMH::EqData().make_field_descriptor_type("RichardsLMH_Data_aux") )
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

    RichardsLMH::EqData eq_data;
    
    return it::Record("Flow_Richards_LMH", "Lumped Mixed-Hybrid solver for unsteady unsaturated Darcy flow.")
        .derive_from(DarcyFlowInterface::get_input_type())
        .copy_keys(DarcyLMH::get_input_type())
        .declare_key("input_fields", it::Array( field_descriptor ), it::Default::obligatory(),
                "Input data for Darcy flow model.")
        .declare_key("output", DarcyFlowMHOutput::get_input_type(eq_data, "Flow_Richards_LMH"),
                IT::Default("{ \"fields\": [ \"pressure_p0\", \"velocity_p0\" ] }"),
                "Specification of output fields and output times.")
        .declare_key("soil_model", soil_rec, it::Default("\"van_genuchten\""),
                "Soil model settings.")
        .close();
}


const int RichardsLMH::registrar =
        Input::register_class< RichardsLMH, Mesh &, const Input::Record >("Flow_Richards_LMH") +
        RichardsLMH::get_input_type().size();



RichardsLMH::RichardsLMH(Mesh &mesh_in, const  Input::Record in_rec)
    : DarcyLMH(mesh_in, in_rec)
{
    data_ = make_shared<EqData>();
    DarcyLMH::data_ = data_;
    EquationBase::eq_data_ = data_.get();
    //data_->edge_new_local_4_mesh_idx_ = &(this->edge_new_local_4_mesh_idx_);
}

void RichardsLMH::initialize_specific() {

    auto model_rec = input_record_.val<Input::Record>("soil_model");
    auto model_type = model_rec.val<SoilModelBase::SoilModelType>("model_type");
    double fraction= model_rec.val<double>("cut_fraction");
    if (model_type == SoilModelBase::van_genuchten)
        data_->soil_model_ = std::make_shared<SoilModel_VanGenuchten>(fraction);
    else if (model_type == SoilModelBase::irmay)
        data_->soil_model_ = std::make_shared<SoilModel_Irmay>(fraction);
    else
        ASSERT(false);

    // create edge vectors
    data_->phead_edge_ = data_->dh_cr_->create_vector();
    data_->water_content_previous_it = data_->dh_cr_disc_->create_vector();
    data_->water_content_previous_time = data_->dh_cr_disc_->create_vector();
    data_->capacity = data_->dh_cr_disc_->create_vector();
}


void RichardsLMH::initial_condition_postprocess()
{
    // set water_content
    // pretty ugly since postprocess change fluxes, which cause bad balance, so we must set them back
    data_->previous_solution.copy_from(data_->data_vec_); // store solution vector
    postprocess();
//     data_->previous_solution.swap(data_->data_vec_); // currently does not mimic VecSwap
    VecSwap(schur0->get_solution(), data_->previous_solution.petsc_vec()); // restore solution vector
}


void RichardsLMH::prepare_new_time_step()
{
    data_->previous_solution.copy_from(data_->data_vec_);
    data_->water_content_previous_time.copy_from(data_->water_content_previous_it);
}

bool RichardsLMH::zero_time_term(bool time_global) {
    if (time_global) {
        return (data_->storativity.input_list_size() == 0)
                && (data_->water_content_saturated.input_list_size() == 0);

    } else {
        return (data_->storativity.field_result(mesh_->region_db().get_region_set("BULK"))
                == result_zeros)
                && (data_->water_content_saturated.field_result(mesh_->region_db().get_region_set("BULK"))
                == result_zeros);
    }
}


void RichardsLMH::assembly_linear_system()
{

    START_TIMER("RicharsLMH::assembly_linear_system");
    
    // update the subvector with edge pressure for solution vector
    data_->dh_cr_->update_subvector(data_->data_vec_, data_->phead_edge_);
    data_->phead_edge_.local_to_ghost_begin();
    data_->phead_edge_.local_to_ghost_end();

    data_->is_linear = data_->genuchten_p_head_scale.field_result(mesh_->region_db().get_region_set("BULK")) == result_zeros;

    bool is_steady = zero_time_term();
    //DebugOut() << "Assembly linear system\n";
        START_TIMER("full assembly");
        if (typeid(*schur0) != typeid(LinSys_BDDC)) {
            schur0->start_add_assembly(); // finish allocation and create matrix
            schur_compl->start_add_assembly();
        }
        data_->time_step_ = time_->dt();
        auto multidim_assembler = AssemblyBase::create< AssemblyRichards >(data_);


        schur0->mat_zero_entries();
        schur0->rhs_zero_entries();
        
        schur_compl->mat_zero_entries();
        schur_compl->rhs_zero_entries();

        assembly_mh_matrix( multidim_assembler ); // fill matrix

            //MatView( *const_cast<Mat*>(schur0->get_matrix()), PETSC_VIEWER_STDOUT_WORLD  );
            //VecView( *const_cast<Vec*>(schur0->get_rhs()),   PETSC_VIEWER_STDOUT_WORLD);

        schur0->finish_assembly();
        schur_compl->finish_assembly();
        schur0->set_matrix_changed();
        schur_compl->set_matrix_changed();


        if (! is_steady) {
            START_TIMER("fix time term");
            //DebugOut() << "setup time term\n";
            // assembly time term and rhs
            solution_changed_for_scatter=true;
        }
}



void RichardsLMH::postprocess() {

    // update structures for balance of water volume
    assembly_linear_system();
    
    // update the subvector with edge pressure for solution vector
    data_->dh_cr_->update_subvector(data_->data_vec_, data_->phead_edge_);
    data_->phead_edge_.local_to_ghost_begin();
    data_->phead_edge_.local_to_ghost_end();

    // modify side fluxes in parallel
    // for every local edge take time term on diagonal and add it to the corresponding flux
    auto multidim_assembler = AssemblyBase::create< AssemblyRichards >(data_);
    
    for ( DHCellAccessor dh_cell : data_->dh_->own_range() ) {
        multidim_assembler[dh_cell.elm().dim()-1]->postprocess_velocity(dh_cell);
    }
}
