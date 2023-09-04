#ifndef ASSEMBLY_BENCHMARK_IMPL_HH_
#define ASSEMBLY_BENCHMARK_IMPL_HH_

#include "assembly_benchmark.hh"
#include "assembly_speed_test.hh"

void AssemblyBenchmarkTest::initialize(const string &input, std::vector<std::string> substances) {
    Input::ReaderToStorage reader( input, get_input_type(), Input::FileFormat::format_YAML );
    in_rec_ = reader.get_root_interface<Input::Record>();

    this->time_ = new TimeGovernor(0.0, 0.5);
    eq_data_->substances_ = substances;

    // Initialize output stream.
    //std::shared_ptr<OutputTime> stream = OutputTime::create_output_stream("solute", in_rec.val<Input::Record>("output_stream"), time().get_unit_conversion());
    //this->set_output_stream(stream);

    eq_fields_->set_components(eq_data_->subst_names());
    eq_fields_->set_input_list( in_rec_.val<Input::Array>("input_fields"), *this->time_ );

    // Resize coefficient array
    eq_data_->max_edg_sides = max(this->mesh_->max_edge_sides(1), max(this->mesh_->max_edge_sides(2), this->mesh_->max_edge_sides(3)));

    eq_data_->output_vec.resize(eq_data_->n_substances());
    eq_fields_->output_field.set_components(eq_data_->subst_names());
    eq_fields_->output_field.set_mesh(*this->mesh_);
    eq_fields_->output_fields.set_mesh(*this->mesh_);
    eq_fields_->output_type(OutputTime::CORNER_DATA);

    eq_fields_->output_field.setup_components();
    for (unsigned int sbi=0; sbi<eq_data_->n_substances(); sbi++)
    {
        // create shared pointer to a FieldFE, pass FE data and push this FieldFE to output_field on all regions
        auto output_field_ptr = create_field_fe< 3, FieldValue<3>::Scalar >(eq_data_->dh_);
        eq_fields_->output_field[sbi].set(output_field_ptr, 0);
        eq_data_->output_vec[sbi] = output_field_ptr->vec();
    }

    // set time marks for writing the output
    //eq_fields_->output_fields.initialize(Model::output_stream_, this->mesh_, in_rec_.val<Input::Record>("output"), this->time_);

    // allocate matrix and vector structures
    eq_data_->ls    = new LinSys*[eq_data_->n_substances()];
    eq_data_->ls_dt = new LinSys*[eq_data_->n_substances()];
    eq_data_->conc_fe.resize(eq_data_->n_substances());
    std::string petsc_default_opts = "-ksp_type bcgs -pc_type ilu -pc_factor_levels 2 -ksp_diagonal_scale_fix -pc_factor_fill 6.0";

    MixedPtr<FE_P_disc> fe(0);
    shared_ptr<DiscreteSpace> ds = make_shared<EqualOrderDiscreteSpace>(this->mesh_, fe);
    eq_data_->dh_p0 = make_shared<DOFHandlerMultiDim>(*this->mesh_);
    eq_data_->dh_p0->distribute_dofs(ds);

    stiffness_matrix.resize(eq_data_->n_substances(), nullptr);
    mass_matrix.resize(eq_data_->n_substances(), nullptr);
    rhs.resize(eq_data_->n_substances(), nullptr);
    mass_vec.resize(eq_data_->n_substances(), nullptr);
    eq_data_->ret_vec.resize(eq_data_->n_substances(), nullptr);

    for (unsigned int sbi = 0; sbi < eq_data_->n_substances(); sbi++) {
        eq_data_->ls[sbi] = new LinSys_PETSC(eq_data_->dh_->distr().get(), petsc_default_opts);
        ( (LinSys_PETSC *)eq_data_->ls[sbi] )->set_from_input( in_rec_.val<Input::Record>("solver") );
        eq_data_->ls[sbi]->set_solution(eq_data_->output_vec[sbi].petsc_vec());

        eq_data_->ls_dt[sbi] = new LinSys_PETSC(eq_data_->dh_->distr().get(), petsc_default_opts);
        ( (LinSys_PETSC *)eq_data_->ls_dt[sbi] )->set_from_input( in_rec_.val<Input::Record>("solver") );
        eq_data_->conc_fe[sbi] = create_field_fe< 3, FieldValue<3>::Scalar >(eq_data_->dh_p0);

        VecDuplicate(eq_data_->ls[sbi]->get_solution(), &eq_data_->ret_vec[sbi]);
    }

    // create assemblation object, finite element structures and distribute DOFs
    mass_assembly_ = new GenericAssembly< MassAssembly >(eq_fields_.get(), eq_data_.get());
    stiffness_assembly_ = new GenericAssembly< StiffnessAssembly >(eq_fields_.get(), eq_data_.get());
    sources_assembly_ = new GenericAssembly< SourcesAssembly >(eq_fields_.get(), eq_data_.get());

    int qsize = mass_assembly_->eval_points()->max_size();
    eq_data_->dif_coef.resize(eq_data_->n_substances());
    for (unsigned int sbi=0; sbi<eq_data_->n_substances(); sbi++)
    {
        eq_data_->dif_coef[sbi].resize(qsize);
    }
}

void AssemblyBenchmarkTest::zero_time_step() {
    START_TIMER("ZERO-TIME STEP");

    DebugOut().fmt("zero_time_step time {}\n", this->time_->t());

    START_TIMER("data initialize");
    eq_fields_->mark_input_times( *this->time_ );
    eq_fields_->set_time(this->time_->step(), LimitSide::left);

    for (unsigned int sbi = 0; sbi < eq_data_->n_substances(); sbi++)
        ( (LinSys_PETSC *)eq_data_->ls[sbi] )->set_initial_guess_nonzero();

    // preallocate system matrix
    for (unsigned int i=0; i<eq_data_->n_substances(); i++)
    {
        // preallocate system matrix
    	eq_data_->ls[i]->start_allocation();
        stiffness_matrix[i] = NULL;
        rhs[i] = NULL;

        // preallocate mass matrix
        eq_data_->ls_dt[i]->start_allocation();
        mass_matrix[i] = NULL;
        VecZeroEntries(eq_data_->ret_vec[i]);
    }
    END_TIMER("data initialize");

    START_TIMER("assembly");
    stiffness_assembly_->assemble(eq_data_->dh_);
    mass_assembly_->assemble(eq_data_->dh_);
    sources_assembly_->assemble(eq_data_->dh_);
    END_TIMER("assembly");

    for (unsigned int i=0; i<eq_data_->n_substances(); i++)
    {
      VecAssemblyBegin(eq_data_->ret_vec[i]);
      VecAssemblyEnd(eq_data_->ret_vec[i]);
    }

    //output_data();
    END_TIMER("ZERO-TIME STEP");
}


void AssemblyBenchmarkTest::update_solution()
{
    START_TIMER("SIMULATION-ONE STEP");

    this->time_->next_time();
    DebugOut().fmt("update_solution time {}\n", this->time_->t());

    START_TIMER("data reinit");
    eq_fields_->set_time(this->time_->step(), LimitSide::left);
    END_TIMER("data reinit");

    START_TIMER("assembly");

    // assemble mass matrix
    for (unsigned int i=0; i<eq_data_->n_substances(); i++)
    {
        eq_data_->ls_dt[i]->start_add_assembly();
        eq_data_->ls_dt[i]->mat_zero_entries();
        VecZeroEntries(eq_data_->ret_vec[i]);
    }
    DebugOut() << " - mass_assembly \n";
    mass_assembly_->assemble(eq_data_->dh_);
    for (unsigned int i=0; i<eq_data_->n_substances(); i++)
    {
        eq_data_->ls_dt[i]->finish_assembly();
        VecAssemblyBegin(eq_data_->ret_vec[i]);
        VecAssemblyEnd(eq_data_->ret_vec[i]);
    }

    // assemble stiffness matrix
    // new fluxes can change the location of Neumann boundary,
    // thus stiffness matrix must be reassembled
    for (unsigned int i=0; i<eq_data_->n_substances(); i++)
    {
        eq_data_->ls[i]->start_add_assembly();
        eq_data_->ls[i]->mat_zero_entries();
    }
    DebugOut() << " - stiffness_assembly \n";
    stiffness_assembly_->assemble(eq_data_->dh_);
    for (unsigned int i=0; i<eq_data_->n_substances(); i++)
    {
        eq_data_->ls[i]->finish_assembly();

        if (stiffness_matrix[i] == NULL)
            MatConvert(*( eq_data_->ls[i]->get_matrix() ), MATSAME, MAT_INITIAL_MATRIX, &stiffness_matrix[i]);
        else
            MatCopy(*( eq_data_->ls[i]->get_matrix() ), stiffness_matrix[i], DIFFERENT_NONZERO_PATTERN);
    }

    // assemble right hand side (due to sources and boundary conditions)
    for (unsigned int i=0; i<eq_data_->n_substances(); i++)
    {
        eq_data_->ls[i]->start_add_assembly();
        eq_data_->ls[i]->rhs_zero_entries();
    }
    DebugOut() << " - sources_assembly \n";
    sources_assembly_->assemble(eq_data_->dh_);
    for (unsigned int i=0; i<eq_data_->n_substances(); i++)
    {
        eq_data_->ls[i]->finish_assembly();

        if (rhs[i] == nullptr) VecDuplicate(*( eq_data_->ls[i]->get_rhs() ), &rhs[i]);
        VecCopy(*( eq_data_->ls[i]->get_rhs() ), rhs[i]);
    }

    END_TIMER("assembly");

    END_TIMER("SIMULATION-ONE STEP");
}

#endif /* ASSEMBLY_BENCHMARK_IMPL_HH_ */
