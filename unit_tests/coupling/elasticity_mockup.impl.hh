#ifndef ELASTICITY_MOCKUP_IMPL_HH_
#define ELASTICITY_MOCKUP_IMPL_HH_

#include "elasticity_mockup.hh"
#include "linsys_null.hh"
#include "bench_fieldalgo_init.hh"
#include "mechanics/assembly_elasticity.hh"


namespace equation_data {

/// Initialize selected fields as FieldModels
template <>
void EqFields::init_fields<InitFieldModel>()
{
    lame_mu.set(Model<3, FieldValue<3>::Scalar>::create(fn_lame_mu(), young_modulus, poisson_ratio), 0.0);
    lame_lambda.set(Model<3, FieldValue<3>::Scalar>::create(fn_lame_lambda(), young_modulus, poisson_ratio), 0.0);
    dirichlet_penalty.set(Model<3, FieldValue<3>::Scalar>::create(fn_dirichlet_penalty(), lame_mu, lame_lambda), 0.0);
}

/// Initialize selected fields as FieldConstant
template <>
void EqFields::init_fields<InitFieldConst>()
{
    /* Definition of const values */
    double val_lame_mu = 0.5;
    double val_lame_lambda = 0.75;
    double val_dirichl_penalty = 1;

    {
	    FieldValue<3>::Scalar f_value(val_lame_mu);
        auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
        field_algo->set_value(f_value);
        lame_mu.set(field_algo, 0.0);
    }
    {
	    FieldValue<3>::Scalar f_value(val_lame_lambda);
        auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
        field_algo->set_value(f_value);
        lame_lambda.set(field_algo, 0.0);
    }
    {
	    FieldValue<3>::Scalar f_value(val_dirichl_penalty);
        auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
        field_algo->set_value(f_value);
        dirichlet_penalty.set(field_algo, 0.0);
    }
}

/// Initialize selected fields as FieldFE
template <>
void EqFields::init_fields<InitFieldFE>()
{
	uint seed = 0;

    MixedPtr<FE_P_disc> fe_p(0);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh_ptr_, fe_p);
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh_ptr_);
    dh->distribute_dofs(ds);

    {
        auto field_algo=std::make_shared<FieldFE<3, FieldValue<3>::Scalar>>();
        VectorMPI data_vec = InitFieldFE::create_data_vec(dh, seed++);
        field_algo->set_fe_data(dh, data_vec);
        lame_mu.set(field_algo, 0.0);
    }
    {
        auto field_algo=std::make_shared<FieldFE<3, FieldValue<3>::Scalar>>();
        VectorMPI data_vec = InitFieldFE::create_data_vec(dh, seed++);
        field_algo->set_fe_data(dh, data_vec);
        lame_lambda.set(field_algo, 0.0);
    }
    {
        auto field_algo=std::make_shared<FieldFE<3, FieldValue<3>::Scalar>>();
        VectorMPI data_vec = InitFieldFE::create_data_vec(dh, seed++);
        field_algo->set_fe_data(dh, data_vec);
        dirichlet_penalty.set(field_algo, 0.0);
    }
}

} // end of namespace equation_data


template <template<IntDim...> class StiffnessAsm, template<IntDim...> class RhsAsm, class FT>
void ElasticityMockupTest::run_mockup(const string &eq_data_input, const std::string &mesh_file, const std::string tag_name, bool use_linsys) {
    START_TIMER(tag_name);
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    ElasticityMockup<StiffnessAsm, RhsAsm> test_mockup(use_linsys);
    test_mockup.create_and_set_mesh(mesh_file);
    test_mockup.initialize( eq_data_input );
    test_mockup.eq_fields_->template init_fields<FT>();
    test_mockup.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER(tag_name);
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

    if (use_linsys_) {
        // allocate matrix and vector structures
        std::string petsc_default_opts = "-ksp_type cg -pc_type hypre -pc_hypre_type boomeramg";
        LinSys *ls = new LinSys_PETSC(eq_data_->dh_->distr().get(), petsc_default_opts);
        ((LinSys_PETSC*)ls)->set_initial_guess_nonzero();
        ls->set_from_input( input_rec.val<Input::Record>("solver") );
        ls->set_solution(eq_fields_->output_field_ptr->vec().petsc_vec());
        eq_data_->ls = ls;
    } else { // use linSysEmpty
        eq_data_->ls = new LinSysNull(eq_data_->dh_->distr().get());
    }

    stiffness_assembly_ = new GenericAssembly< Stiffness >(eq_data_.get(), eq_data_->dh_.get());
    rhs_assembly_ = new GenericAssembly< Rhs >(eq_data_.get(), eq_data_->dh_.get());
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
    if (use_linsys_) MatSetOption(*eq_data_->ls->get_matrix(), MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
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
