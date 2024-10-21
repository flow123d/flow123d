#ifndef ELASTICITY_MOCKUP_IMPL_HH_
#define ELASTICITY_MOCKUP_IMPL_HH_

#include "elasticity_mockup.hh"
#include "elasticity_mockup_assembly.hh"


void ElasticityMockupTest::run_fullassembly_const(const string &eq_data_input, const std::string &mesh_file) {
    // FullAssembly + field_const
    START_TIMER("FullAssembly_const");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    ElasticityMockup<Stiffness_FullAssembly, Rhs_FullAssembly> test_full_asm_const;
    test_full_asm_const.create_and_set_mesh(mesh_file);
    test_full_asm_const.initialize( eq_data_input );
    test_full_asm_const.eq_fields_->init_field_constants(0.5, 0.75, 1);
    test_full_asm_const.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("FullAssembly_const");
}

void ElasticityMockupTest::run_fullassembly_model(const string &eq_data_input, const std::string &mesh_file) {
    // FullAssembly + field_model
    START_TIMER("FullAssembly_model");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    ElasticityMockup<Stiffness_FullAssembly, Rhs_FullAssembly> test_full_asm_model;
    test_full_asm_model.create_and_set_mesh(mesh_file);
    test_full_asm_model.initialize( eq_data_input );
    test_full_asm_model.eq_fields_->init_field_models();
    test_full_asm_model.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("FullAssembly_model");
}

void ElasticityMockupTest::run_computelocal_const(const string &eq_data_input, const std::string &mesh_file) {
    // FullAssembly + field_const
    START_TIMER("ComputeLocal_const");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    ElasticityMockup<Stiffness_ComputeLocal, Rhs_ComputeLocal> test_comp_local_const;
    test_comp_local_const.create_and_set_mesh(mesh_file);
    test_comp_local_const.initialize( eq_data_input );
    test_comp_local_const.eq_fields_->init_field_constants(0.5, 0.75, 1);
    test_comp_local_const.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("ComputeLocal_const");
}

void ElasticityMockupTest::run_computelocal_model(const string &eq_data_input, const std::string &mesh_file) {
    // FullAssembly + field_model
    START_TIMER("ComputeLocal_model");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    ElasticityMockup<Stiffness_ComputeLocal, Rhs_ComputeLocal> test_comp_local_model;
    test_comp_local_model.create_and_set_mesh(mesh_file);
    test_comp_local_model.initialize( eq_data_input );
    test_comp_local_model.eq_fields_->init_field_models();
    test_comp_local_model.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("ComputeLocal_model");
}

void ElasticityMockupTest::run_evalfields_const(const string &eq_data_input, const std::string &mesh_file) {
    // FullAssembly + field_const
    START_TIMER("EvalFields_const");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    ElasticityMockup<Stiffness_EvalFields, Rhs_EvalFields> test_eval_fields_const;
    test_eval_fields_const.create_and_set_mesh(mesh_file);
    test_eval_fields_const.initialize( eq_data_input );
    test_eval_fields_const.eq_fields_->init_field_constants(0.5, 0.75, 1);
    test_eval_fields_const.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("EvalFields_const");
}

void ElasticityMockupTest::run_evalfields_model(const string &eq_data_input, const std::string &mesh_file) {
    // FullAssembly + field_model
    START_TIMER("EvalFields_model");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    ElasticityMockup<Stiffness_EvalFields, Rhs_EvalFields> test_eval_fields_model;
    test_eval_fields_model.create_and_set_mesh(mesh_file);
    test_eval_fields_model.initialize( eq_data_input );
    test_eval_fields_model.eq_fields_->init_field_models();
    test_eval_fields_model.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("FullAssembly_model");
}



template<template<IntDim...> class Stiffness, template<IntDim...> class Rhs>
void ElasticityMockup<Stiffness, Rhs>::initialize(const string &input) {
    Input::ReaderToStorage reader( input, get_input_type(), Input::FileFormat::format_YAML );
    input_rec = reader.get_root_interface<Input::Record>();

    this->time_ = new TimeGovernor(0.0, 0.5);

    eq_fields_->set_components({"displacement"});
    eq_fields_->set_input_list( input_rec.val<Input::Array>("input_fields"), time() );

    // create shared pointer to a FieldFE, pass FE data and push this FieldFE to output_field on all regions
    eq_fields_->output_field_ptr = create_field_fe<3, FieldValue<3>::VectorFixed>(eq_data_->dh_);
    eq_fields_->output_field.set(eq_fields_->output_field_ptr, 0.);

    eq_fields_->output_fields.set_mesh(*mesh_);
    eq_fields_->output_field.output_type(OutputTime::CORNER_DATA);

    // set time marks for writing the output
    //eq_fields_->output_fields.initialize(output_stream_, mesh_, input_rec.val<Input::Record>("output"), this->time());

    // equation default PETSc solver options
    std::string petsc_default_opts;
    petsc_default_opts = "-ksp_type cg -pc_type hypre -pc_hypre_type boomeramg";

    // allocate matrix and vector structures
    LinSys *ls = new LinSys_PETSC(eq_data_->dh_->distr().get(), petsc_default_opts);
    ((LinSys_PETSC*)ls)->set_initial_guess_nonzero();
    ls->set_from_input( input_rec.val<Input::Record>("solver") );
    ls->set_solution(eq_fields_->output_field_ptr->vec().petsc_vec());
    eq_data_->ls = ls;

    stiffness_assembly_ = new GenericAssembly< Stiffness >(eq_fields_.get(), eq_data_.get(), eq_data_->dh_.get());
    rhs_assembly_ = new GenericAssembly< Rhs >(eq_fields_.get(), eq_data_.get(), eq_data_->dh_.get());
}

template<template<IntDim...> class Stiffness, template<IntDim...> class Rhs>
void ElasticityMockup<Stiffness, Rhs>::zero_time_step() {
    START_TIMER("ZERO-TIME STEP");

    DebugOut().fmt("zero_time_step time {}\n", this->time_->t());

    START_TIMER("data initialize");
	eq_fields_->mark_input_times( *time_ );
	eq_fields_->set_time(time_->step(), LimitSide::right);
    END_TIMER("data initialize");

    START_TIMER("assembly");
	eq_data_->ls->start_allocation();
    stiffness_assembly_->assemble(eq_data_->dh_);
    rhs_assembly_->assemble(eq_data_->dh_);

    eq_data_->ls->start_add_assembly();
    MatSetOption(*eq_data_->ls->get_matrix(), MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    eq_data_->ls->mat_zero_entries();
    eq_data_->ls->rhs_zero_entries();
    stiffness_assembly_->assemble(eq_data_->dh_);
    rhs_assembly_->assemble(eq_data_->dh_);
    END_TIMER("assembly");
    eq_data_->ls->finish_assembly();
    //eq_data_->ls->solve();

    END_TIMER("ZERO-TIME STEP");
}

template<template<IntDim...> class Stiffness, template<IntDim...> class Rhs>
void ElasticityMockup<Stiffness, Rhs>::update_solution() {
    START_TIMER("MECH-ONE STEP");

    time_->next_time();

    START_TIMER("data reinit");
    eq_fields_->set_time(time_->step(), LimitSide::right);
    END_TIMER("data reinit");

    // assemble stiffness matrix
    if (eq_data_->ls->get_matrix() == NULL
        || eq_fields_->subset(FieldFlag::in_main_matrix).changed())
    {
        DebugOut() << "Mechanics: Assembling matrix.\n";
        eq_data_->ls->start_add_assembly();
        eq_data_->ls->mat_zero_entries();
        stiffness_assembly_->assemble(eq_data_->dh_);
        eq_data_->ls->finish_assembly();
    }

    // assemble right hand side (due to sources and boundary conditions)
    if (eq_data_->ls->get_rhs() == NULL
        || eq_fields_->subset(FieldFlag::in_rhs).changed())
    {
        DebugOut() << "Mechanics: Assembling right hand side.\n";
        eq_data_->ls->start_add_assembly();
        eq_data_->ls->rhs_zero_entries();
        rhs_assembly_->assemble(eq_data_->dh_);
        eq_data_->ls->finish_assembly();
    }

//    START_TIMER("solve");
//    eq_data_->ls->solve();
//    END_TIMER("solve");

    END_TIMER("MECH-ONE STEP");
}


#endif /* ELASTICITY_MOCKUP_IMPL_HH_ */
