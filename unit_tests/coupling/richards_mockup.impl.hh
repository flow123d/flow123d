#ifndef RICHARDS_MOCKUP_IMPL_HH_
#define RICHARDS_MOCKUP_IMPL_HH_

#include "richards_mockup.hh"
#include "linsys_null.hh"
#include "flow/assembly_richards.hh"
#include "fields/field_add_potential.hh"


void RichardsMockupTest::run_fullassembly_const(const string &eq_data_input, const std::string &mesh_file) {
    // FullAssembly + field_const
    START_TIMER("FullAssembly_const");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    RichardsMockup<MHMatrixAssemblyRichardsDim, equation_data_richards::EqData> test_full_asm_const(true);
    test_full_asm_const.create_and_set_mesh(mesh_file);
    test_full_asm_const.initialize( eq_data_input );
    test_full_asm_const.eq_fields_->init_field_constants(0.5, arma::vec3("0.8 0.6 0"));
    test_full_asm_const.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("FullAssembly_const");
}

void RichardsMockupTest::run_fullassembly_model(const string &eq_data_input, const std::string &mesh_file) {
    // FullAssembly + field_model
    START_TIMER("FullAssembly_model");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    RichardsMockup<MHMatrixAssemblyRichardsDim, equation_data_richards::EqData> test_full_asm_model(true);
    test_full_asm_model.create_and_set_mesh(mesh_file);
    test_full_asm_model.initialize( eq_data_input );
    test_full_asm_model.eq_fields_->init_field_models();
    test_full_asm_model.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("FullAssembly_model");
}

void RichardsMockupTest::run_computelocal_const(const string &eq_data_input, const std::string &mesh_file) {
    // ComputeLocal + field_const
    START_TIMER("ComputeLocal_const");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    RichardsMockup<MHMatrixAssemblyRichardsDim, equation_data_richards::EqData> test_comp_local_const(false);
    test_comp_local_const.create_and_set_mesh(mesh_file);
    test_comp_local_const.initialize( eq_data_input );
    test_comp_local_const.eq_fields_->init_field_constants(0.5, arma::vec3("0.8 0.6 0"));
    test_comp_local_const.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("ComputeLocal_const");
}

void RichardsMockupTest::run_computelocal_model(const string &eq_data_input, const std::string &mesh_file) {
    // ComputeLocal + field_model
    START_TIMER("ComputeLocal_model");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    RichardsMockup<MHMatrixAssemblyRichardsDim, equation_data_richards::EqData> test_comp_local_model(false);
    test_comp_local_model.create_and_set_mesh(mesh_file);
    test_comp_local_model.initialize( eq_data_input );
    test_comp_local_model.eq_fields_->init_field_models();
    test_comp_local_model.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("ComputeLocal_model");
}

void RichardsMockupTest::run_evalfields_const(const string &eq_data_input, const std::string &mesh_file) {
    // EvalFields + field_const
    START_TIMER("EvalFields_const");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    RichardsMockup<MHMatrixAssemblyRichardsDim, equation_data_richards::EqData> test_eval_fields_const(false);
    test_eval_fields_const.create_and_set_mesh(mesh_file);
    test_eval_fields_const.initialize( eq_data_input );
    test_eval_fields_const.eq_fields_->init_field_constants(0.5, arma::vec3("0.8 0.6 0"));
    test_eval_fields_const.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("EvalFields_const");
}

void RichardsMockupTest::run_evalfields_model(const string &eq_data_input, const std::string &mesh_file) {
    // EvalFields + field_model
    START_TIMER("EvalFields_model");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    RichardsMockup<MHMatrixAssemblyRichardsDim, equation_data_richards::EqData> test_eval_fields_model(false);
    test_eval_fields_model.create_and_set_mesh(mesh_file);
    test_eval_fields_model.initialize( eq_data_input );
    test_eval_fields_model.eq_fields_->init_field_models();
    test_eval_fields_model.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("EvalFields_model");
}


template<template<IntDim...> class MhMatrix, class TEqData>
void RichardsMockup<MhMatrix, TEqData>::initialize_asm() {
    this->mh_matrix_assembly_ = new GenericAssembly< MhMatrix >(this->eq_data_.get());
}

template<template<IntDim...> class MhMatrix, class TEqData>
void RichardsMockup<MhMatrix, TEqData>::initialize_specific() {
    // create edge vectors
    this->eq_data_->water_content_previous_time = this->eq_data_->dh_cr_disc_->create_vector();
    this->eq_data_->capacity = this->eq_data_->dh_cr_disc_->create_vector();

    this->eq_fields_->water_content_ptr = create_field_fe< 3, FieldValue<3>::Scalar >(this->eq_data_->dh_cr_disc_);
    this->eq_fields_->water_content.set(this->eq_fields_->water_content_ptr, 0.0);

    this->eq_fields_->conductivity_ptr = create_field_fe< 3, FieldValue<3>::Scalar >(this->eq_data_->dh_p_);
    this->eq_fields_->conductivity_richards.set(this->eq_fields_->conductivity_ptr, 0.0);
}

template<template<IntDim...> class MhMatrix, class TEqData>
void RichardsMockup<MhMatrix, TEqData>::assembly_linear_system()
{
    this->eq_data_->p_edge_solution.local_to_ghost_begin();
    this->eq_data_->p_edge_solution.local_to_ghost_end();

    this->eq_data_->is_linear = this->eq_fields_->genuchten_p_head_scale.field_result(this->mesh_->region_db().get_region_set("BULK")) == result_zeros;

    //DebugOut() << "Assembly linear system\n";
    START_TIMER("full assembly");
    this->eq_data_->lin_sys_schur->start_add_assembly();

    this->eq_data_->time_step_ = this->time_->dt();

    this->eq_data_->lin_sys_schur->mat_zero_entries();
    this->eq_data_->lin_sys_schur->rhs_zero_entries();

    this->mh_matrix_assembly_->assemble(this->eq_data_->dh_); // fill matrix

    this->eq_data_->lin_sys_schur->finish_assembly();
    this->eq_data_->lin_sys_schur->set_matrix_changed();
}

#endif /*RICHARDS_MOCKUP_IMPL_HH_ */
