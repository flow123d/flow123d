#ifndef DG_MOCKUP_IMPL_HH_
#define DG_MOCKUP_IMPL_HH_

#include "DG_mockup.hh"
#include "linsys_null.hh"
#include "bench_fieldalgo_init.hh"
#include "transport/assembly_dg.hh"


namespace equation_data {

/// Initialize selected fields as FieldModels
template <>
void EqFields::init_fields<InitFieldModel>()
{
    setup_mf_components();

    v_norm.set(Model<3, FieldValue<3>::Scalar>::create(fn_conc_v_norm(), flow_flux), 0.0);
    mass_matrix_coef.set(Model<3, FieldValue<3>::Scalar>::create(fn_conc_mass_matrix(), cross_section, water_content), 0.0);
    retardation_coef.set(
        Model<3, FieldValue<3>::Scalar>::create_multi(fn_conc_retardation(), cross_section, porosity, rock_density, sorption_coefficient),
        0.0
    );
    sources_density_out.set(Model<3, FieldValue<3>::Scalar>::create_multi(fn_conc_sources_dens(), cross_section, sources_density), 0.0);
    sources_sigma_out.set(Model<3, FieldValue<3>::Scalar>::create_multi(fn_conc_sources_sigma(), cross_section, sources_sigma), 0.0);
    sources_conc_out.set(Model<3, FieldValue<3>::Scalar>::create_multi(fn_conc_sources_conc(), sources_conc), 0.0);
    std::vector<typename Field<3, FieldValue<3>::VectorFixed>::FieldBasePtr> ad_coef_ptr_vec;
    for (unsigned int sbi=0; sbi<sorption_coefficient.size(); sbi++)
        ad_coef_ptr_vec.push_back( Model<3, FieldValue<3>::VectorFixed>::create(fn_conc_ad_coef(), flow_flux) );
    advection_coef.set(ad_coef_ptr_vec, 0.0);
    diffusion_coef.set(
        Model<3, FieldValue<3>::TensorFixed>::create_multi(
            fn_conc_diff_coef(), diff_m, flow_flux, v_norm, disp_l, disp_t, water_content, porosity, cross_section
        ),
        0.0
    );
}

/// Initialize selected fields as FieldConstant
template <>
void EqFields::init_fields<InitFieldConst>()
{
    /* Definition of const values */
	double vnorm = 1;
    double mass_matrix = 0.5;
    double ret_coef = 0.75;
    double sources_dens = 1;
    double sourc_sigma = 0.25;
    double sourc_conc = 0.5;
    arma::vec3 advec_coef = arma::vec3("1 2 3");
    arma::mat33 diff_coef = arma::mat33("0.5 0 0; 0 0.75 0; 0 0 1");

	setup_mf_components();

    {
	    FieldValue<3>::Scalar f_value(vnorm);
        auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
        field_algo->set_value(f_value);
        v_norm.set(field_algo, 0.0);
    }
    {
	    FieldValue<3>::Scalar f_value(mass_matrix);
        auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
        field_algo->set_value(f_value);
        mass_matrix_coef.set(field_algo, 0.0);
    }
    {
        FieldValue<3>::Scalar f_value(ret_coef);
        std::vector<typename Field<3, FieldValue<3>::Scalar>::FieldBasePtr> field_vec;
        for (unsigned int sbi=0; sbi<sorption_coefficient.size(); sbi++) {
            auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
            field_algo->set_value(f_value);
            field_vec.push_back(field_algo);
        }
        retardation_coef.set(field_vec, 0.0);
    }
    {
        FieldValue<3>::Scalar f_value(sources_dens);
        std::vector<typename Field<3, FieldValue<3>::Scalar>::FieldBasePtr> field_vec;
        for (unsigned int sbi=0; sbi<sources_density.size(); sbi++) {
            auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
            field_algo->set_value(f_value);
            field_vec.push_back(field_algo);
        }
        sources_density_out.set(field_vec, 0.0);
    }
    {
        FieldValue<3>::Scalar f_value(sourc_sigma);
        std::vector<typename Field<3, FieldValue<3>::Scalar>::FieldBasePtr> field_vec;
        for (unsigned int sbi=0; sbi<sources_sigma.size(); sbi++) {
            auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
            field_algo->set_value(f_value);
            field_vec.push_back(field_algo);
        }
        sources_sigma_out.set(field_vec, 0.0);
    }
    {
        FieldValue<3>::Scalar f_value(sourc_conc);
        std::vector<typename Field<3, FieldValue<3>::Scalar>::FieldBasePtr> field_vec;
        for (unsigned int sbi=0; sbi<sources_conc.size(); sbi++) {
            auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
            field_algo->set_value(f_value);
            field_vec.push_back(field_algo);
        }
        sources_conc_out.set(field_vec, 0.0);
    }
    {
        FieldValue<3>::TensorFixed f_value(diff_coef);
        std::vector<typename Field<3, FieldValue<3>::TensorFixed>::FieldBasePtr> field_vec;
        for (unsigned int sbi=0; sbi<diff_m.size(); sbi++) {
            auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::TensorFixed>>();
            field_algo->set_value(f_value);
            field_vec.push_back(field_algo);
        }
        diffusion_coef.set(field_vec, 0.0);
    }
}

/// Initialize selected fields as FieldFE
template <>
void EqFields::init_fields<InitFieldFE>()
{
	setup_mf_components();
	uint seed = 0;

    MixedPtr<FE_P_disc> fe_base(0);
    MixedPtr<FiniteElement> fe_tens = mixed_fe_system(fe_base, FEType::FETensor, 9);
    MixedPtr<FiniteElement> fe_rt = MixedPtr<FE_RT0_disc>();
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh_ptr_, fe_base);
    std::shared_ptr<DiscreteSpace> ds_tens = std::make_shared<EqualOrderDiscreteSpace>(mesh_ptr_, fe_tens);
    std::shared_ptr<DiscreteSpace> ds_rt = std::make_shared<EqualOrderDiscreteSpace>(mesh_ptr_, fe_rt);
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh_ptr_);
    std::shared_ptr<DOFHandlerMultiDim> dh_tens = std::make_shared<DOFHandlerMultiDim>(*mesh_ptr_);
    std::shared_ptr<DOFHandlerMultiDim> dh_rt = std::make_shared<DOFHandlerMultiDim>(*mesh_ptr_);
    dh->distribute_dofs(ds);
    dh_tens->distribute_dofs(ds_tens);
    dh_rt->distribute_dofs(ds_rt);

    {
        auto field_algo=std::make_shared<FieldFE<3, FieldValue<3>::VectorFixed>>();
        VectorMPI data_vec = InitFieldFE::create_data_vec(dh_rt, seed++);
        field_algo->set_fe_data(dh_rt, data_vec);
        flow_flux.set(field_algo, 0.0);
    }
    v_norm.set(Model<3, FieldValue<3>::Scalar>::create(fn_conc_v_norm(), flow_flux), 0.0);
    {
        auto field_algo=std::make_shared<FieldFE<3, FieldValue<3>::Scalar>>();
        VectorMPI data_vec = InitFieldFE::create_data_vec(dh, seed++);
        field_algo->set_fe_data(dh, data_vec);
        mass_matrix_coef.set(field_algo, 0.0);
    }
    {
        std::vector<typename Field<3, FieldValue<3>::Scalar>::FieldBasePtr> field_vec;
        for (unsigned int sbi=0; sbi<sorption_coefficient.size(); sbi++) {
            auto field_algo=std::make_shared<FieldFE<3, FieldValue<3>::Scalar>>();
            VectorMPI data_vec = InitFieldFE::create_data_vec(dh, seed++);
            field_algo->set_fe_data(dh, data_vec);
            field_vec.push_back(field_algo);
        }
        retardation_coef.set(field_vec, 0.0);
    }
    {
        std::vector<typename Field<3, FieldValue<3>::Scalar>::FieldBasePtr> field_vec;
        for (unsigned int sbi=0; sbi<sources_density.size(); sbi++) {
            auto field_algo=std::make_shared<FieldFE<3, FieldValue<3>::Scalar>>();
            VectorMPI data_vec = InitFieldFE::create_data_vec(dh, seed++);
            field_algo->set_fe_data(dh, data_vec);
            field_vec.push_back(field_algo);
        }
        sources_density_out.set(field_vec, 0.0);
    }
    {
        std::vector<typename Field<3, FieldValue<3>::Scalar>::FieldBasePtr> field_vec;
        for (unsigned int sbi=0; sbi<sources_sigma.size(); sbi++) {
            auto field_algo=std::make_shared<FieldFE<3, FieldValue<3>::Scalar>>();
            VectorMPI data_vec = InitFieldFE::create_data_vec(dh, seed++);
            field_algo->set_fe_data(dh, data_vec);
            field_vec.push_back(field_algo);
        }
        sources_sigma_out.set(field_vec, 0.0);
    }
    {
        std::vector<typename Field<3, FieldValue<3>::Scalar>::FieldBasePtr> field_vec;
        for (unsigned int sbi=0; sbi<sources_conc.size(); sbi++) {
            auto field_algo=std::make_shared<FieldFE<3, FieldValue<3>::Scalar>>();
            VectorMPI data_vec = InitFieldFE::create_data_vec(dh, seed++);
            field_algo->set_fe_data(dh, data_vec);
            field_vec.push_back(field_algo);
        }
        sources_conc_out.set(field_vec, 0.0);
    }
    {
        std::vector<typename Field<3, FieldValue<3>::TensorFixed>::FieldBasePtr> field_vec;
        for (unsigned int sbi=0; sbi<diff_m.size(); sbi++) {
            auto field_algo=std::make_shared<FieldFE<3, FieldValue<3>::TensorFixed>>();
            VectorMPI data_vec = InitFieldFE::create_data_vec(dh_tens, seed++);
            field_algo->set_fe_data(dh_tens, data_vec);
            field_vec.push_back(field_algo);
        }
        diffusion_coef.set(field_vec, 0.0);
    }
}

} // end of namespace equation_data


template <template<IntDim...> class MassAsm, template<IntDim...> class StiffnessAsm, template<IntDim...> class SourcesAsm, class FT>
void DGMockupTest::run_mockup(const string &eq_data_input, const std::string &mesh_file, const std::string tag_name, bool use_linsys) {
    START_TIMER(tag_name);
    START_TIMER("full_mesh"); // necessary for correct process of profiler output
    DGMockup<MassAsm, StiffnessAsm, SourcesAsm> test_mockup(use_linsys);
    test_mockup.create_and_set_mesh(mesh_file);
    test_mockup.initialize( eq_data_input, {"A", "B"} );
    test_mockup.eq_fields_->template init_fields<FT>();
    test_mockup.run_simulation();
    END_TIMER("full_mesh");
    END_TIMER(tag_name);
}


template<template<IntDim...> class Mass, template<IntDim...> class Stiffness, template<IntDim...> class Sources>
void DGMockup<Mass, Stiffness, Sources>::initialize(const string &input, std::vector<std::string> substances) {
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
    eq_data_->subst_idx_.resize(eq_data_->n_substances(), 0);
    eq_data_->balance_ = std::make_shared<BalanceNull>("mass", this->mesh_);

    MixedPtr<FE_P_disc> fe(0);
    shared_ptr<DiscreteSpace> ds = make_shared<EqualOrderDiscreteSpace>(this->mesh_, fe);
    eq_data_->dh_p0 = make_shared<DOFHandlerMultiDim>(*this->mesh_);
    eq_data_->dh_p0->distribute_dofs(ds);

    stiffness_matrix.resize(eq_data_->n_substances(), nullptr);
    mass_matrix.resize(eq_data_->n_substances(), nullptr);
    rhs.resize(eq_data_->n_substances(), nullptr);
    mass_vec.resize(eq_data_->n_substances(), nullptr);
    eq_data_->ret_vec.resize(eq_data_->n_substances(), nullptr);

    if (use_linsys_) {
        std::string petsc_default_opts = "-ksp_type bcgs -pc_type ilu -pc_factor_levels 2 -ksp_diagonal_scale_fix -pc_factor_fill 6.0";
        for (unsigned int sbi = 0; sbi < eq_data_->n_substances(); sbi++) {
            eq_data_->ls[sbi] = new LinSys_PETSC(eq_data_->dh_->distr().get(), petsc_default_opts);
            ( (LinSys_PETSC *)eq_data_->ls[sbi] )->set_from_input( in_rec_.val<Input::Record>("solver") );
            eq_data_->ls[sbi]->set_solution(eq_data_->output_vec[sbi].petsc_vec());

            eq_data_->ls_dt[sbi] = new LinSys_PETSC(eq_data_->dh_->distr().get(), petsc_default_opts);
            ( (LinSys_PETSC *)eq_data_->ls_dt[sbi] )->set_from_input( in_rec_.val<Input::Record>("solver") );
            eq_data_->conc_fe[sbi] = create_field_fe< 3, FieldValue<3>::Scalar >(eq_data_->dh_p0);

            VecDuplicate(eq_data_->ls[sbi]->get_solution(), &eq_data_->ret_vec[sbi]);
        }
    } else {
        for (unsigned int sbi = 0; sbi < eq_data_->n_substances(); sbi++) {
            eq_data_->ls[sbi] = new LinSysNull(eq_data_->dh_->distr().get());
            eq_data_->ls_dt[sbi] = new LinSysNull(eq_data_->dh_->distr().get());
            VecDuplicate(eq_data_->output_vec[sbi].petsc_vec(), &eq_data_->ret_vec[sbi]);
        }
    }

    // create assemblation object, finite element structures and distribute DOFs
    mass_assembly_ = new GenericAssembly< Mass >(eq_data_.get(), eq_data_->dh_.get());
    stiffness_assembly_ = new GenericAssembly< Stiffness >(eq_data_.get(), eq_data_->dh_.get());
    sources_assembly_ = new GenericAssembly< Sources >(eq_data_.get(), eq_data_->dh_.get());

    int qsize = mass_assembly_->eval_points()->max_size();
    eq_data_->dif_coef.resize(eq_data_->n_substances());
    for (unsigned int sbi=0; sbi<eq_data_->n_substances(); sbi++)
    {
        eq_data_->dif_coef[sbi].resize(qsize);
    }
}

template<template<IntDim...> class Mass, template<IntDim...> class Stiffness, template<IntDim...> class Sources>
void DGMockup<Mass, Stiffness, Sources>::zero_time_step() {
    START_TIMER("ZERO-TIME STEP");

    DebugOut().fmt("zero_time_step time {}\n", this->time_->t());

    START_TIMER("data initialize");
    eq_fields_->mark_input_times( *this->time_ );
    eq_fields_->set_time(this->time_->step(), LimitSide::left);

    if (use_linsys_) {
        for (unsigned int sbi = 0; sbi < eq_data_->n_substances(); sbi++)
            ( (LinSys_PETSC *)eq_data_->ls[sbi] )->set_initial_guess_nonzero();
    }

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


template<template<IntDim...> class Mass, template<IntDim...> class Stiffness, template<IntDim...> class Sources>
void DGMockup<Mass, Stiffness, Sources>::update_solution()
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
    if (use_linsys_)
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
    if (use_linsys_)
        for (unsigned int i=0; i<eq_data_->n_substances(); i++)
        {
            eq_data_->ls[i]->finish_assembly();

            if (rhs[i] == nullptr) VecDuplicate(*( eq_data_->ls[i]->get_rhs() ), &rhs[i]);
            VecCopy(*( eq_data_->ls[i]->get_rhs() ), rhs[i]);
        }

    END_TIMER("assembly");

    END_TIMER("SIMULATION-ONE STEP");
}

#endif /*DG_MOCKUP_IMPL_HH_ */
