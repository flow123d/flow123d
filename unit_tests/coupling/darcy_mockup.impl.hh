#ifndef DARCY_MOCKUP_IMPL_HH_
#define DARCY_MOCKUP_IMPL_HH_

#include "darcy_mockup.hh"
#include "linsys_null.hh"
#include "flow/assembly_lmh.hh"
#include "fields/field_add_potential.hh"


void DarcyMockupTest::run_fullassembly_const(const string &eq_data_input, const std::string &mesh_file) {
    // FullAssembly + field_const
    START_TIMER("FullAssembly_const");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    DarcyMockup<MHMatrixAssemblyLMHDim> test_full_asm_const(true);
    test_full_asm_const.create_and_set_mesh(mesh_file);
    test_full_asm_const.initialize( eq_data_input );
    test_full_asm_const.eq_fields_->init_field_constants(0.5, arma::vec3("0.8 0.6 0"));
    test_full_asm_const.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("FullAssembly_const");
}

void DarcyMockupTest::run_fullassembly_model(const string &eq_data_input, const std::string &mesh_file) {
    // FullAssembly + field_model
    START_TIMER("FullAssembly_model");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    DarcyMockup<MHMatrixAssemblyLMHDim> test_full_asm_model(true);
    test_full_asm_model.create_and_set_mesh(mesh_file);
    test_full_asm_model.initialize( eq_data_input );
    test_full_asm_model.eq_fields_->init_field_models();
    test_full_asm_model.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("FullAssembly_model");
}

void DarcyMockupTest::run_computelocal_const(const string &eq_data_input, const std::string &mesh_file) {
    // ComputeLocal + field_const
    START_TIMER("ComputeLocal_const");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    DarcyMockup<MHMatrixAssemblyLMHDim> test_comp_local_const(false);
    test_comp_local_const.create_and_set_mesh(mesh_file);
    test_comp_local_const.initialize( eq_data_input );
    test_comp_local_const.eq_fields_->init_field_constants(0.5, arma::vec3("0.8 0.6 0"));
    test_comp_local_const.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("ComputeLocal_const");
}

void DarcyMockupTest::run_computelocal_model(const string &eq_data_input, const std::string &mesh_file) {
    // ComputeLocal + field_model
    START_TIMER("ComputeLocal_model");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    DarcyMockup<MHMatrixAssemblyLMHDim> test_comp_local_model(false);
    test_comp_local_model.create_and_set_mesh(mesh_file);
    test_comp_local_model.initialize( eq_data_input );
    test_comp_local_model.eq_fields_->init_field_models();
    test_comp_local_model.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("ComputeLocal_model");
}

void DarcyMockupTest::run_evalfields_const(const string &eq_data_input, const std::string &mesh_file) {
    // EvalFields + field_const
    START_TIMER("EvalFields_const");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    DarcyMockup<MHMatrixEvalFieldsDim> test_eval_fields_const(false);
    test_eval_fields_const.create_and_set_mesh(mesh_file);
    test_eval_fields_const.initialize( eq_data_input );
    test_eval_fields_const.eq_fields_->init_field_constants(0.5, arma::vec3("0.8 0.6 0"));
    test_eval_fields_const.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("EvalFields_const");
}

void DarcyMockupTest::run_evalfields_model(const string &eq_data_input, const std::string &mesh_file) {
    // EvalFields + field_model
    START_TIMER("EvalFields_model");
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    DarcyMockup<MHMatrixEvalFieldsDim> test_eval_fields_model(false);
    test_eval_fields_model.create_and_set_mesh(mesh_file);
    test_eval_fields_model.initialize( eq_data_input );
    test_eval_fields_model.eq_fields_->init_field_models();
    test_eval_fields_model.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER("EvalFields_model");
}


template<template<IntDim...> class MhMatrix>
void DarcyMockup<MhMatrix>::initialize(const string &input) {
    Input::ReaderToStorage reader( input, get_input_type(), Input::FileFormat::format_YAML );
    in_rec_ = reader.get_root_interface<Input::Record>();

    // init eq_data
    this->time_ = new TimeGovernor(0.0, 0.5);
	eq_data_->mesh = mesh_;

    auto gravity_array = in_rec_.val<Input::Array>("gravity");
    std::vector<double> gvec;
    gravity_array.copy_to(gvec);
    gvec.push_back(0.0); // zero pressure shift
    eq_data_->gravity_ =  arma::vec(gvec);
    eq_data_->gravity_vec_ = eq_data_->gravity_.subvec(0,2);

    FieldValue<3>::VectorFixed gvalue(eq_data_->gravity_vec_);
    auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::VectorFixed>>();
    field_algo->set_value(gvalue);
    eq_fields_->gravity_field.set(field_algo, 0.0);
    eq_fields_->bc_gravity.set(field_algo, 0.0);

    eq_fields_->bc_pressure.add_factory(
            std::make_shared<AddPotentialFactory<3, FieldValue<3>::Scalar> >
            (eq_fields_->bc_gravity, eq_fields_->X(), eq_fields_->bc_piezo_head) );
    eq_fields_->bc_switch_pressure.add_factory(
            std::make_shared<AddPotentialFactory<3, FieldValue<3>::Scalar> >
            (eq_fields_->bc_gravity, eq_fields_->X(), eq_fields_->bc_switch_piezo_head) );
    eq_fields_->init_pressure.add_factory(
            std::make_shared<AddPotentialFactory<3, FieldValue<3>::Scalar> >
            (eq_fields_->gravity_field, eq_fields_->X(), eq_fields_->init_piezo_head) );


    eq_fields_->set_input_list( this->in_rec_.val<Input::Array>("input_fields"), *time_ );

    // Check that the time step was set for the transient simulation.
    if (! zero_time_term(true) && time_->is_default() ) {
        //THROW(ExcAssertMsg());
        //THROW(ExcMissingTimeGovernor() << input_record_.ei_address());
        MessageOut() << "Missing the key 'time', obligatory for the transient problems." << endl;
        ASSERT_PERMANENT(false);
    }

    eq_fields_->mark_input_times(*time_);

    // init FieldsCoords
    eq_fields_->add_coords_field();

    { // construct flux and pressure fields
		uint rt_component = 0;
		eq_data_->full_solution = eq_data_->dh_->create_vector();
        auto ele_flux_ptr = create_field_fe<3, FieldValue<3>::VectorFixed>(eq_data_->dh_, &eq_data_->full_solution, rt_component);
        eq_fields_->flux.set(ele_flux_ptr, 0.0);

		uint p_ele_component = 1;
        auto ele_pressure_ptr = create_field_fe<3, FieldValue<3>::Scalar>(eq_data_->dh_, &eq_data_->full_solution, p_ele_component);
        eq_fields_->field_ele_pressure.set(ele_pressure_ptr, 0.0);

        uint p_edge_component = 2;
        auto edge_pressure_ptr = create_field_fe<3, FieldValue<3>::Scalar>(eq_data_->dh_, &eq_data_->full_solution, p_edge_component);
        eq_fields_->field_edge_pressure.set(edge_pressure_ptr, 0.0);
    }


    { // init DOF handlers represents element pressure DOFs
        uint p_element_component = 1;
        eq_data_->dh_p_ = std::make_shared<SubDOFHandlerMultiDim>(eq_data_->dh_,p_element_component);
    }

    { // init DOF handlers represents edge DOFs
        uint p_edge_component = 2;
        eq_data_->dh_cr_ = std::make_shared<SubDOFHandlerMultiDim>(eq_data_->dh_,p_edge_component);
    }

    { // init DOF handlers represents side DOFs
		MixedPtr<FE_CR_disc> fe_cr_disc;
		std::shared_ptr<DiscreteSpace> ds_cr_disc = std::make_shared<EqualOrderDiscreteSpace>( mesh_, fe_cr_disc);
		eq_data_->dh_cr_disc_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
		eq_data_->dh_cr_disc_->distribute_dofs(ds_cr_disc);
    }

    eq_data_->init();

    // create solution vector for 2. Schur complement linear system
    eq_data_->p_edge_solution = eq_data_->dh_cr_->create_vector();
    eq_data_->p_edge_solution_previous = eq_data_->dh_cr_->create_vector();
    eq_data_->p_edge_solution_previous_time = eq_data_->dh_cr_->create_vector();

    // Initialize bc_switch_dirichlet to size of global boundary.
    eq_data_->bc_switch_dirichlet.resize(mesh_->n_elements()+mesh_->bc_mesh()->n_elements(), 1);

    eq_data_->nonlinear_iteration_=0;
    Input::AbstractRecord rec = this->in_rec_
            .val<Input::Record>("nonlinear_solver")
            .val<Input::AbstractRecord>("linear_solver");

    // auxiliary set_time call  since allocation assembly evaluates fields as well
    data_changed_ = eq_fields_->set_time(time_->step(), LimitSide::right) || data_changed_;
    if (use_linsys_) {
        eq_data_->lin_sys_schur = std::make_shared<LinSys_PETSC>( &(*eq_data_->dh_cr_->distr()) );
        eq_data_->lin_sys_schur->set_from_input(rec);
        eq_data_->lin_sys_schur->set_positive_definite();
        eq_data_->lin_sys_schur->set_solution( eq_data_->p_edge_solution.petsc_vec() );
        eq_data_->lin_sys_schur->set_symmetric();
        ((LinSys_PETSC *)eq_data_->lin_sys_schur.get())->set_initial_guess_nonzero(true);

        eq_data_->lin_sys_schur->start_allocation();

        allocate_mh_matrix();

        eq_data_->full_solution.zero_entries();
        eq_data_->p_edge_solution.zero_entries();
    } else {
        eq_data_->lin_sys_schur = std::make_shared<LinSysNull>( &(*eq_data_->dh_cr_->distr()) );
    }

    // initialization of balance object
    eq_data_->balance_ = std::make_shared<BalanceNull>("water", mesh_);

    this->initialize_asm();
}

template<template<IntDim...> class MhMatrix>
void DarcyMockup<MhMatrix>::initialize_asm() {
    this->mh_matrix_assembly_ = new GenericAssembly< MhMatrix >(eq_data_.get());
}

template<template<IntDim...> class MhMatrix>
void DarcyMockup<MhMatrix>::zero_time_step() {
    START_TIMER("ZERO-TIME STEP");
    data_changed_ = eq_fields_->set_time(time_->step(), LimitSide::right) || data_changed_;

    eq_data_->p_edge_solution.zero_entries();

    MessageOut() << "Flow zero time step - unsteady case\n";
    eq_data_->time_step_ = time_->dt();

    assembly_linear_system();

    accept_time_step(); // accept zero time step, i.e. initial condition

    END_TIMER("ZERO-TIME STEP");
}


template<template<IntDim...> class MhMatrix>
void DarcyMockup<MhMatrix>::update_solution()
{
    START_TIMER("SIMULATION-ONE STEP");

    time_->next_time();

    time_->view("DARCY"); //time governor information output

    data_changed_ = eq_fields_->set_time(time_->step(), LimitSide::left) || data_changed_;
    bool zero_time_term_from_left=zero_time_term();

    bool jump_time = eq_fields_->storativity.is_jump_time();
    if (! zero_time_term_from_left) {
        MessageOut() << "Flow time step - unsteady case\n";
        // time term not treated as zero
        // Unsteady solution up to the T.

        // this flag is necesssary for switching BC to avoid setting zero neumann on the whole boundary in the steady case
    	eq_data_->use_steady_assembly_ = false;

        solve_nonlinear(); // with left limit data
    }

    data_changed_ = eq_fields_->set_time(time_->step(), LimitSide::right) || data_changed_;
    bool zero_time_term_from_right=zero_time_term();
    if (zero_time_term_from_right) {
        MessageOut() << "Flow time step - steady case\n";
        // this flag is necesssary for switching BC to avoid setting zero neumann on the whole boundary in the steady case
    	eq_data_->use_steady_assembly_ = true;
        solve_nonlinear(); // with right limit data
    } else if (! zero_time_term_from_left && jump_time) {
    	WarningOut() << "Discontinuous time term not supported yet.\n";
    }

    eq_data_->full_solution.local_to_ghost_begin();
    eq_data_->full_solution.local_to_ghost_end();
    END_TIMER("SIMULATION-ONE STEP");
}

template<template<IntDim...> class MhMatrix>
void DarcyMockup<MhMatrix>::assembly_linear_system()
{
    eq_data_->p_edge_solution.local_to_ghost_begin();
    eq_data_->p_edge_solution.local_to_ghost_end();

    eq_data_->is_linear=true;
    //DebugOut() << "Assembly linear system\n";

    eq_data_->lin_sys_schur->start_add_assembly();
    eq_data_->lin_sys_schur->mat_zero_entries();
    eq_data_->lin_sys_schur->rhs_zero_entries();

    eq_data_->time_step_ = time_->dt();

    this->mh_matrix_assembly_->assemble(eq_data_->dh_);; // fill matrix
    eq_data_->lin_sys_schur->finish_assembly();
    eq_data_->lin_sys_schur->set_matrix_changed();
}

template<template<IntDim...> class MhMatrix>
void DarcyMockup<MhMatrix>::solve_nonlinear()
{
    assembly_linear_system();
    double residual_norm = eq_data_->lin_sys_schur->compute_residual();
    eq_data_->nonlinear_iteration_ = 0;
    MessageOut().fmt("[nonlinear solver] norm of initial residual: {}\n", residual_norm);

    // Reduce is_linear flag.
    int is_linear_common;
    MPI_Allreduce(&(eq_data_->is_linear), &is_linear_common,1, MPI_INT ,MPI_MIN,PETSC_COMM_WORLD);

    Input::Record nl_solver_rec = in_rec_.val<Input::Record>("nonlinear_solver");
    this->tolerance_ = nl_solver_rec.val<double>("tolerance");
    this->max_n_it_  = nl_solver_rec.val<unsigned int>("max_it");
    this->min_n_it_  = nl_solver_rec.val<unsigned int>("min_it");
    if (this->min_n_it_ > this->max_n_it_) this->min_n_it_ = this->max_n_it_;

    if (! is_linear_common) {
        // set tolerances of the linear solver unless they are set by user.
        eq_data_->lin_sys_schur->set_tolerances(0.1*this->tolerance_, 0.01*this->tolerance_, 10000, 100);
    }
    vector<double> convergence_history;

    while (eq_data_->nonlinear_iteration_ < this->min_n_it_ ||
           (residual_norm > this->tolerance_ &&  eq_data_->nonlinear_iteration_ < this->max_n_it_ )) {
        ASSERT_EQ( convergence_history.size(), eq_data_->nonlinear_iteration_ );
        convergence_history.push_back(residual_norm);

        // print_matlab_matrix("matrix_" + std::to_string(time_->step().index()) + "_it_" + std::to_string(nonlinear_iteration_));
        // stagnation test
        if (convergence_history.size() >= 5 &&
            convergence_history[ convergence_history.size() - 1]/convergence_history[ convergence_history.size() - 2] > 0.9 &&
            convergence_history[ convergence_history.size() - 1]/convergence_history[ convergence_history.size() - 5] > 0.8) {
            // stagnation
            if (input_record_.val<Input::Record>("nonlinear_solver").val<bool>("converge_on_stagnation")) {
            	WarningOut().fmt("Accept solution on stagnation. Its: {} Residual: {}\n", eq_data_->nonlinear_iteration_, residual_norm);
                break;
            } else {
                ASSERT(false).error("Nonlinear solver did not converge. Reason: Stagnation.");
            }
        }

        if (! is_linear_common){
        	eq_data_->p_edge_solution_previous.copy_from(eq_data_->p_edge_solution);
        	eq_data_->p_edge_solution_previous.local_to_ghost_begin();
        	eq_data_->p_edge_solution_previous.local_to_ghost_end();
        }

        LinSys::SolveInfo si = eq_data_->lin_sys_schur->solve();
        MessageOut().fmt("[schur solver] lin. it: {}, reason: {}, residual: {}\n",
        		si.n_iterations, si.converged_reason, eq_data_->lin_sys_schur->compute_residual());

        eq_data_->nonlinear_iteration_++;

        // hack to make BDDC work with empty compute_residual
        if (is_linear_common){
            // we want to print this info in linear (and steady) case
            residual_norm = eq_data_->lin_sys_schur->compute_residual();
            MessageOut().fmt("[nonlinear solver] lin. it: {}, reason: {}, residual: {}\n",
        		si.n_iterations, si.converged_reason, residual_norm);
            break;
        }
        data_changed_=true; // force reassembly for non-linear case

        double alpha = 1; // how much of new solution
        VecAXPBY(eq_data_->p_edge_solution.petsc_vec(), (1-alpha), alpha, eq_data_->p_edge_solution_previous.petsc_vec());

        assembly_linear_system();

        residual_norm = eq_data_->lin_sys_schur->compute_residual();
        MessageOut().fmt("[nonlinear solver] it: {} lin. it: {}, reason: {}, residual: {}\n",
                eq_data_->nonlinear_iteration_, si.n_iterations, si.converged_reason, residual_norm);
    }

    // adapt timestep
    if (! this->zero_time_term()) {
        double mult = 1.0;
        if (eq_data_->nonlinear_iteration_ < 3) mult = 1.6;
        if (eq_data_->nonlinear_iteration_ > 7) mult = 0.7;
        time_->set_upper_constraint(time_->dt() * mult, "Darcy adaptivity.");
    }
}

#endif /*DARCY_MOCKUP_IMPL_HH_ */
