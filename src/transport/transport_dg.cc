/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    transport_dg.cc
 * @brief   Discontinuous Galerkin method for equation of transport with dispersion.
 * @author  Jan Stebel
 */

#include "system/sys_profiler.hh"
#include "transport/transport_dg.hh"

#include "io/output_time.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_values.hh"
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "fields/field_fe.hh"
#include "la/linsys_PETSC.hh"
#include "transport/advection_diffusion_model.hh"
#include "transport/concentration_model.hh"
#include "transport/heat_model.hh"
#include "coupling/balance.hh"

#include "fields/multi_field.hh"
#include "fields/generic_field.hh"
#include "input/factory.hh"
#include "io/equation_output.hh"

FLOW123D_FORCE_LINK_IN_CHILD(concentrationTransportModel);
FLOW123D_FORCE_LINK_IN_CHILD(heatModel);



using namespace Input::Type;

template<class Model>
const Selection & TransportDG<Model>::get_dg_variant_selection_input_type() {
	return Selection("DG_variant", "Type of penalty term.")
		.add_value(non_symmetric, "non-symmetric", "non-symmetric weighted interior penalty DG method")
		.add_value(incomplete,    "incomplete",    "incomplete weighted interior penalty DG method")
		.add_value(symmetric,     "symmetric",     "symmetric weighted interior penalty DG method")
		.close();
}

/*
 *  Should be removed
template<class Model>
const Selection & TransportDG<Model>::EqData::get_output_selection() {
	return Model::ModelEqData::get_output_selection_input_type(
	        "DG",
	        "Implicit in time Discontinuous Galerkin solver")
		.copy_values(EqData().make_output_field_selection("").close())
		ConvectionTransport::EqData().output_fields
		                            .make_output_field_selection(
		                                "ConvectionTransport_output_fields",
		                                "Selection of output fields for Convection Solute Transport model.")
		                            .close()),
		.close();
}
*/

template<class Model>
const Record & TransportDG<Model>::get_input_type() {
    std::string equation_name = std::string(Model::ModelEqData::name()) + "_DG";
	return Model::get_input_type("DG", "DG solver")
		.declare_key("solver", LinSys_PETSC::get_input_type(), Default::obligatory(),
				"Linear solver for MH problem.")
		.declare_key("input_fields", Array(
		        TransportDG<Model>::EqData()
		            .make_field_descriptor_type(equation_name)),
		        IT::Default::obligatory(),
		        "Input fields of the equation.")
		.declare_key("dg_variant", TransportDG<Model>::get_dg_variant_selection_input_type(), Default("\"non-symmetric\""),
				"Variant of interior penalty discontinuous Galerkin method.")
		.declare_key("dg_order", Integer(0,3), Default("1"),
				"Polynomial order for finite element in DG method (order 0 is suitable if there is no diffusion/dispersion).")
        .declare_key("output",
                EqData().output_fields.make_output_type(equation_name, ""),
                IT::Default("{ \"fields\": [ " + Model::ModelEqData::default_output_field() + "] }"),
                "Setting of the field output.")
		.close();
}

template<class Model>
const int TransportDG<Model>::registrar =
		Input::register_class< TransportDG<Model>, Mesh &, const Input::Record>(std::string(Model::ModelEqData::name()) + "_DG") +
		TransportDG<Model>::get_input_type().size();




FEObjects::FEObjects(Mesh *mesh_, unsigned int fe_order)
{
    unsigned int q_order;

	switch (fe_order)
	{
	case 0:
		q_order = 0;
		fe1_ = new FE_P_disc<0,1,3>;
		fe2_ = new FE_P_disc<0,2,3>;
		fe3_ = new FE_P_disc<0,3,3>;
		break;

	case 1:
		q_order = 2;
		fe1_ = new FE_P_disc<1,1,3>;
		fe2_ = new FE_P_disc<1,2,3>;
		fe3_ = new FE_P_disc<1,3,3>;
		break;

	case 2:
		q_order = 4;
		fe1_ = new FE_P_disc<2,1,3>;
		fe2_ = new FE_P_disc<2,2,3>;
		fe3_ = new FE_P_disc<2,3,3>;
		break;

	case 3:
		q_order = 6;
		fe1_ = new FE_P_disc<3,1,3>;
		fe2_ = new FE_P_disc<3,2,3>;
		fe3_ = new FE_P_disc<3,3,3>;
		break;

	default:
	    q_order=0;
		xprintf(PrgErr, "Unsupported polynomial order %d for finite elements in TransportDG ", fe_order);
		break;
	}

	fe_rt1_ = new FE_RT0<1,3>;
	fe_rt2_ = new FE_RT0<2,3>;
	fe_rt3_ = new FE_RT0<3,3>;

	q0_ = new QGauss<0>(q_order);
	q1_ = new QGauss<1>(q_order);
	q2_ = new QGauss<2>(q_order);
	q3_ = new QGauss<3>(q_order);

	map0_ = new MappingP1<0,3>;
	map1_ = new MappingP1<1,3>;
	map2_ = new MappingP1<2,3>;
	map3_ = new MappingP1<3,3>;

	dh_ = new DOFHandlerMultiDim(*mesh_);

	dh_->distribute_dofs(*fe1_, *fe2_, *fe3_);
}


FEObjects::~FEObjects()
{
	delete fe1_;
	delete fe2_;
	delete fe3_;
	delete fe_rt1_;
	delete fe_rt2_;
	delete fe_rt3_;
	delete q0_;
	delete q1_;
	delete q2_;
	delete q3_;
	delete map0_;
	delete map1_;
	delete map2_;
	delete map3_;
	delete dh_;
}

template<> FiniteElement<0,3> *FEObjects::fe<0>() { return 0; }
template<> FiniteElement<1,3> *FEObjects::fe<1>() { return fe1_; }
template<> FiniteElement<2,3> *FEObjects::fe<2>() { return fe2_; }
template<> FiniteElement<3,3> *FEObjects::fe<3>() { return fe3_; }

template<> FiniteElement<0,3> *FEObjects::fe_rt<0>() { return 0; }
template<> FiniteElement<1,3> *FEObjects::fe_rt<1>() { return fe_rt1_; }
template<> FiniteElement<2,3> *FEObjects::fe_rt<2>() { return fe_rt2_; }
template<> FiniteElement<3,3> *FEObjects::fe_rt<3>() { return fe_rt3_; }

template<> Quadrature<0> *FEObjects::q<0>() { return q0_; }
template<> Quadrature<1> *FEObjects::q<1>() { return q1_; }
template<> Quadrature<2> *FEObjects::q<2>() { return q2_; }
template<> Quadrature<3> *FEObjects::q<3>() { return q3_; }

template<> Mapping<0,3> *FEObjects::mapping<0>() { return map0_; }
template<> Mapping<1,3> *FEObjects::mapping<1>() { return map1_; }
template<> Mapping<2,3> *FEObjects::mapping<2>() { return map2_; }
template<> Mapping<3,3> *FEObjects::mapping<3>() { return map3_; }

DOFHandlerMultiDim *FEObjects::dh() { return dh_; }


template<class Model>
TransportDG<Model>::EqData::EqData() : Model::ModelEqData()
{
    *this+=fracture_sigma
            .name("fracture_sigma")
            .description(
            "Coefficient of diffusive transfer through fractures (for each substance).")
            .units( UnitSI::dimensionless() )
            .input_default("1.0")
            .flags_add(FieldFlag::in_main_matrix);

    *this+=dg_penalty
            .name("dg_penalty")
            .description(
            "Penalty parameter influencing the discontinuity of the solution (for each substance). "
            "Its default value 1 is sufficient in most cases. Higher value diminishes the inter-element jumps.")
            .units( UnitSI::dimensionless() )
            .input_default("1.0")
            .flags_add(FieldFlag::in_rhs & FieldFlag::in_main_matrix);

    *this += region_id.name("region_id")
    	        .units( UnitSI::dimensionless())
    	        .flags(FieldFlag::equation_external_output);


    // add all input fields to the output list
    output_fields += *this;

}

template<class Model>
TransportDG<Model>::TransportDG(Mesh & init_mesh, const Input::Record in_rec)
        : Model(init_mesh, in_rec),
		  input_rec(in_rec),
		  allocation_done(false)
{
	// Can not use name() + "constructor" here, since START_TIMER only accepts const char *
	// due to constexpr optimization.
	START_TIMER(Model::ModelEqData::name());
	// Check that Model is derived from AdvectionDiffusionModel.
	static_assert(std::is_base_of<AdvectionDiffusionModel, Model>::value, "");

	this->eq_data_ = &data_;


    // Set up physical parameters.
    data_.set_mesh(init_mesh);
    data_.region_id = GenericField<3>::region_id(*Model::mesh_);


    // DG variant and order
    dg_variant = in_rec.val<DGVariant>("dg_variant");
    dg_order = in_rec.val<unsigned int>("dg_order");
    
    Model::init_from_input(in_rec);

    // create finite element structures and distribute DOFs
    feo = new FEObjects(Model::mesh_, dg_order);
    //DebugOut().fmt("TDG: solution size {}\n", feo->dh()->n_global_dofs());

}


template<class Model>
void TransportDG<Model>::initialize()
{
    data_.set_components(Model::substances_.names());
    data_.set_input_list( input_rec.val<Input::Array>("input_fields") );

    // DG stabilization parameters on boundary edges
    gamma.resize(Model::n_substances());
    for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
    	gamma[sbi].resize(Model::mesh_->boundary_.size());

    // Resize coefficient arrays
    int qsize = max(feo->q<0>()->size(), max(feo->q<1>()->size(), max(feo->q<2>()->size(), feo->q<3>()->size())));
    int max_edg_sides = max(Model::mesh_->max_edge_sides(1), max(Model::mesh_->max_edge_sides(2), Model::mesh_->max_edge_sides(3)));
    mm_coef.resize(qsize);
    ret_coef.resize(Model::n_substances());
    sorption_sources.resize(Model::n_substances());
    ad_coef.resize(Model::n_substances());
    dif_coef.resize(Model::n_substances());
    for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
    {
    	ret_coef[sbi].resize(qsize);
    	ad_coef[sbi].resize(qsize);
    	dif_coef[sbi].resize(qsize);
    }
    ad_coef_edg.resize(max_edg_sides);
    dif_coef_edg.resize(max_edg_sides);
    for (int sd=0; sd<max_edg_sides; sd++)
    {
    	ad_coef_edg[sd].resize(Model::n_substances());
    	dif_coef_edg[sd].resize(Model::n_substances());
    	for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
    	{
    		ad_coef_edg[sd][sbi].resize(qsize);
    		dif_coef_edg[sd][sbi].resize(qsize);
    	}
    }

	output_vec.resize(Model::n_substances());
	//output_solution.resize(Model::n_substances());
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	unsigned int output_vector_size= (rank==0)?feo->dh()->n_global_dofs():0;
	for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
	{
		// for each substance we allocate output array and vector
		//output_solution[sbi] = new double[feo->dh()->n_global_dofs()];
		VecCreateSeq(PETSC_COMM_SELF, output_vector_size, &output_vec[sbi]);
	}
	data_.output_field.set_components(Model::substances_.names());
	data_.output_field.set_mesh(*Model::mesh_);
    data_.output_type(OutputTime::CORNER_DATA);

    data_.output_field.setup_components();
	for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
	{
		// create shared pointer to a FieldFE, pass FE data and push this FieldFE to output_field on all regions
		std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > output_field_ptr(new FieldFE<3, FieldValue<3>::Scalar>);
		output_field_ptr->set_fe_data(feo->dh(), feo->mapping<1>(), feo->mapping<2>(), feo->mapping<3>(), &output_vec[sbi]);
		data_.output_field[sbi].set_field(Model::mesh_->region_db().get_region_set("ALL"), output_field_ptr, 0);
	}

    // set time marks for writing the output
    data_.output_fields.initialize(Model::output_stream_, input_rec.val<Input::Record>("output"), this->time());

    // allocate matrix and vector structures
    ls    = new LinSys*[Model::n_substances()];
    ls_dt = new LinSys*[Model::n_substances()];
    solution_elem_ = new double*[Model::n_substances()];
    for (unsigned int sbi = 0; sbi < Model::n_substances(); sbi++) {
    	ls[sbi] = new LinSys_PETSC(feo->dh()->distr());
    	( (LinSys_PETSC *)ls[sbi] )->set_from_input( input_rec.val<Input::Record>("solver") );
    	ls[sbi]->set_solution(NULL);

    	ls_dt[sbi] = new LinSys_PETSC(feo->dh()->distr());
    	( (LinSys_PETSC *)ls_dt[sbi] )->set_from_input( input_rec.val<Input::Record>("solver") );
    	solution_elem_[sbi] = new double[Model::mesh_->get_el_ds()->lsize()];
    }
    stiffness_matrix = new Mat[Model::n_substances()];
    mass_matrix = new Mat[Model::n_substances()];
    rhs = new Vec[Model::n_substances()];
    mass_vec = new Vec[Model::n_substances()];


    // initialization of balance object
    if (Model::balance_ != nullptr)
    {
	    Model::balance_->allocate(feo->dh()->distr()->lsize(),
	    		max(feo->fe<1>()->n_dofs(), max(feo->fe<2>()->n_dofs(), feo->fe<3>()->n_dofs())));
    }

}


template<class Model>
TransportDG<Model>::~TransportDG()
{
    delete Model::time_;

    if (gamma.size() > 0) {
        // initialize called

        for (auto &vec : output_vec) VecDestroy(&vec);

        for (unsigned int i=0; i<Model::n_substances(); i++)
        {
            delete ls[i];
            delete[] solution_elem_[i];
            delete ls_dt[i];
            MatDestroy(&stiffness_matrix[i]);
            MatDestroy(&mass_matrix[i]);
            VecDestroy(&rhs[i]);
            VecDestroy(&mass_vec[i]);
        }
        delete[] ls;
        delete[] solution_elem_;
        delete[] ls_dt;
        delete[] stiffness_matrix;
        delete[] mass_matrix;
        delete[] rhs;
        delete[] mass_vec;
        delete feo;

    }

}


template<class Model>
void TransportDG<Model>::output_vector_gather()
{
    VecScatter output_scatter;
    VecScatterCreateToZero(ls[0]->get_solution(), &output_scatter, PETSC_NULL);
	for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
	{
		// gather solution to output_vec[sbi]
	    VecScatterBegin(output_scatter, ls[sbi]->get_solution(), output_vec[sbi], INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(output_scatter, ls[sbi]->get_solution(), output_vec[sbi], INSERT_VALUES, SCATTER_FORWARD);
	}
    VecScatterDestroy(&(output_scatter));

}



template<class Model>
void TransportDG<Model>::zero_time_step()
{
	START_TIMER(Model::ModelEqData::name());
	data_.mark_input_times( *(Model::time_) );
	data_.set_time(Model::time_->step(), LimitSide::left);


    // set initial conditions
    set_initial_condition();
    for (unsigned int sbi = 0; sbi < Model::n_substances(); sbi++)
    	( (LinSys_PETSC *)ls[sbi] )->set_initial_guess_nonzero();

    // check first time assembly - needs preallocation
    if (!allocation_done) preallocate();

    // after preallocation we assemble the matrices and vectors required for mass balance
    if (Model::balance_ != nullptr)
    {
		for (unsigned int sbi=0; sbi<Model::n_substances(); ++sbi)
		{
			Model::balance_->calculate_mass(Model::subst_idx[sbi], ls[sbi]->get_solution());
			Model::balance_->calculate_source(Model::subst_idx[sbi], ls[sbi]->get_solution());
			Model::balance_->calculate_flux(Model::subst_idx[sbi], ls[sbi]->get_solution());

			// add sources due to sorption
			vector<double> masses(Model::mesh_->region_db().bulk_size());
			double mass = 0;
			Model::balance_->calculate_mass(Model::subst_idx[sbi], ls[sbi]->get_solution(), masses);
			for (auto reg_mass : masses)
				mass += reg_mass;

			sorption_sources[sbi] = mass;
		}
    }

    output_data();
}


template<class Model>
void TransportDG<Model>::preallocate()
{
	// preallocate system matrix
	for (unsigned int i=0; i<Model::n_substances(); i++)
	{
		// preallocate system matrix
		ls[i]->start_allocation();
		stiffness_matrix[i] = NULL;
		rhs[i] = NULL;

		// preallocate mass matrix
		ls_dt[i]->start_allocation();
		mass_matrix[i] = NULL;
	}
	assemble_stiffness_matrix();
	assemble_mass_matrix();
	set_sources();
	set_boundary_conditions();

	allocation_done = true;
}



template<class Model>
void TransportDG<Model>::update_solution()
{
	START_TIMER("DG-ONE STEP");

	Model::time_->next_time();
	Model::time_->view("TDG");
    
    START_TIMER("data reinit");
    data_.set_time(Model::time_->step(), LimitSide::left);
    END_TIMER("data reinit");
    
	// assemble mass matrix
    if (mass_matrix[0] == NULL || data_.subset(FieldFlag::in_time_term).changed() )
	{
    	for (unsigned int i=0; i<Model::n_substances(); i++)
        {
            ls_dt[i]->start_add_assembly();
            ls_dt[i]->mat_zero_entries();
        }
        assemble_mass_matrix();
        for (unsigned int i=0; i<Model::n_substances(); i++)
        {
            ls_dt[i]->finish_assembly();
            // construct mass_vec for initial time
            if (mass_matrix[i] == NULL)
            {
                VecDuplicate(ls[i]->get_solution(), &mass_vec[i]);
                MatMult(*(ls_dt[i]->get_matrix()), ls[i]->get_solution(), mass_vec[i]);
                if (Model::balance_ != nullptr && Model::balance_->cumulative())
                {
                    double total_mass = 0;
                    VecSum(mass_vec[i], &total_mass);
                    sorption_sources[i] -= total_mass;
                }
                MatConvert(*( ls_dt[i]->get_matrix() ), MATSAME, MAT_INITIAL_MATRIX, &mass_matrix[i]);
            }
            else
                MatCopy(*( ls_dt[i]->get_matrix() ), mass_matrix[i], DIFFERENT_NONZERO_PATTERN);
        }
    }

	// assemble stiffness matrix
    if (stiffness_matrix[0] == NULL
    		|| data_.subset(FieldFlag::in_main_matrix).changed()
    		|| Model::flux_changed)
    {
        // new fluxes can change the location of Neumann boundary,
        // thus stiffness matrix must be reassembled
    	for (unsigned int i=0; i<Model::n_substances(); i++)
    	{
    		ls[i]->start_add_assembly();
    		ls[i]->mat_zero_entries();
    	}
        assemble_stiffness_matrix();
        for (unsigned int i=0; i<Model::n_substances(); i++)
        {
        	ls[i]->finish_assembly();

        	if (stiffness_matrix[i] == NULL)
        		MatConvert(*( ls[i]->get_matrix() ), MATSAME, MAT_INITIAL_MATRIX, &stiffness_matrix[i]);
        	else
        		MatCopy(*( ls[i]->get_matrix() ), stiffness_matrix[i], DIFFERENT_NONZERO_PATTERN);
        }
    }

    // assemble right hand side (due to sources and boundary conditions)
    if (rhs[0] == NULL
    		|| data_.subset(FieldFlag::in_rhs).changed()
    		|| Model::flux_changed)
    {
    	for (unsigned int i=0; i<Model::n_substances(); i++)
    	{
    		ls[i]->start_add_assembly();
    		ls[i]->rhs_zero_entries();
    	}
    	set_sources();
    	set_boundary_conditions();
    	for (unsigned int i=0; i<Model::n_substances(); i++)
    	{
    		ls[i]->finish_assembly();

    		if (rhs[i] == nullptr) VecDuplicate(*( ls[i]->get_rhs() ), &rhs[i]);
    		VecCopy(*( ls[i]->get_rhs() ), rhs[i]);
    	}
    }

    Model::flux_changed = false;


    /* Apply backward Euler time integration.
     *
     * Denoting A the stiffness matrix and M the mass matrix, the algebraic system at the k-th time level reads
     *
     *   (1/dt M + A)u^k = f + 1/dt M.u^{k-1}
     *
     * Hence we modify at each time level the right hand side:
     *
     *   f^k = f + 1/dt M u^{k-1},
     *
     * where f stands for the term stemming from the force and boundary conditions.
     * Accordingly, we set
     *
     *   A^k = A + 1/dt M.
     *
     */
    Mat m;
    START_TIMER("solve");
    for (unsigned int i=0; i<Model::n_substances(); i++)
    {
    	MatConvert(stiffness_matrix[i], MATSAME, MAT_INITIAL_MATRIX, &m);
		MatAXPY(m, 1./Model::time_->dt(), mass_matrix[i], SUBSET_NONZERO_PATTERN);
		ls[i]->set_matrix(m, DIFFERENT_NONZERO_PATTERN);
		Vec w;
		VecDuplicate(rhs[i], &w);
		VecWAXPY(w, 1./Model::time_->dt(), mass_vec[i], rhs[i]);
		ls[i]->set_rhs(w);

		VecDestroy(&w);
		MatDestroy(&m);

		ls[i]->solve();

		// update mass_vec due to possible changes in mass matrix
		MatMult(*(ls_dt[i]->get_matrix()), ls[i]->get_solution(), mass_vec[i]);
    }
    END_TIMER("solve");

    calculate_cumulative_balance();

    END_TIMER("DG-ONE STEP");
}


template<class Model>
void TransportDG<Model>::calculate_concentration_matrix()
{
    // calculate element averages of solution
    for (unsigned int i_cell=0; i_cell<Model::mesh_->get_el_ds()->lsize(); i_cell++)
    {
    	typename DOFHandlerBase::CellIterator elem = Model::mesh_->element(feo->dh()->el_index(i_cell));

    	unsigned int n_dofs;
    	switch (elem->dim())
    	{
    	case 1:
    		n_dofs = feo->fe<1>()->n_dofs();
    		break;
    	case 2:
			n_dofs = feo->fe<2>()->n_dofs();
			break;
    	case 3:
			n_dofs = feo->fe<3>()->n_dofs();
			break;
    	}

    	unsigned int dof_indices[n_dofs];
    	feo->dh()->get_dof_indices(elem, dof_indices);

    	for (unsigned int sbi=0; sbi<Model::n_substances(); ++sbi)
    	{
    		solution_elem_[sbi][i_cell] = 0;

    		for (unsigned int j=0; j<n_dofs; ++j)
    			solution_elem_[sbi][i_cell] += ls[sbi]->get_solution_array()[dof_indices[j]-feo->dh()->distr()->begin()];

			solution_elem_[sbi][i_cell] = max(solution_elem_[sbi][i_cell]/n_dofs, 0.);
    	}
    }
}




template<class Model>
void TransportDG<Model>::output_data()
{
    //if (!Model::time_->is_current( Model::time_->marks().type_output() )) return;


    START_TIMER("DG-OUTPUT");

    // gather the solution from all processors
    data_.output_fields.set_time( this->time().step(), LimitSide::left);
    //if (data_.output_fields.is_field_output_time(data_.output_field, this->time().step()) )
    output_vector_gather();
    data_.output_fields.output(this->time().step());

	Model::output_data();

    END_TIMER("DG-OUTPUT");
}


template<class Model>
void TransportDG<Model>::calculate_cumulative_balance()
{
    if (Model::balance_ != nullptr && Model::balance_->cumulative())
    {
    	for (unsigned int sbi=0; sbi<Model::n_substances(); ++sbi)
    	{
    		Model::balance_->calculate_cumulative_sources(Model::subst_idx[sbi], ls[sbi]->get_solution(), Model::time_->dt());
    		Model::balance_->calculate_cumulative_fluxes(Model::subst_idx[sbi], ls[sbi]->get_solution(), Model::time_->dt());

    		// add sources due to sorption
    		vector<double> masses(Model::mesh_->region_db().bulk_size());
    		double mass = 0;
    		Model::balance_->calculate_mass(Model::subst_idx[sbi], ls[sbi]->get_solution(), masses);
    		for (auto reg_mass : masses)
    			mass += reg_mass;
    		double total_mass = 0;
    		VecSum(mass_vec[sbi], &total_mass);

    		Model::balance_->add_cumulative_source(Model::subst_idx[sbi], mass-total_mass-sorption_sources[sbi]);
    		sorption_sources[sbi] = mass - total_mass;
    	}
    }
}


template<class Model>
void TransportDG<Model>::calculate_instant_balance()
{
    if (Model::balance_ != nullptr)
    {
    	for (unsigned int sbi=0; sbi<Model::n_substances(); ++sbi)
    	{
    		Model::balance_->calculate_mass(Model::subst_idx[sbi], ls[sbi]->get_solution());
			Model::balance_->calculate_source(Model::subst_idx[sbi], ls[sbi]->get_solution());
			Model::balance_->calculate_flux(Model::subst_idx[sbi], ls[sbi]->get_solution());
    	}
    }
}


template<class Model>
void TransportDG<Model>::assemble_mass_matrix()
{
  START_TIMER("assemble_mass");
  	if (Model::balance_ != nullptr)
  		Model::balance_->start_mass_assembly(Model::subst_idx);
	assemble_mass_matrix<1>();
	assemble_mass_matrix<2>();
	assemble_mass_matrix<3>();
	if (Model::balance_ != nullptr)
		Model::balance_->finish_mass_assembly(Model::subst_idx);
  END_TIMER("assemble_mass");
}


template<class Model> template<unsigned int dim>
void TransportDG<Model>::assemble_mass_matrix()
{
    FEValues<dim,3> fe_values(*feo->mapping<dim>(), *feo->q<dim>(), *feo->fe<dim>(), update_values | update_JxW_values | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
    vector<int> dof_indices(ndofs);
    PetscScalar local_mass_matrix[ndofs*ndofs];
    vector<PetscScalar> local_mass_balance_vector(ndofs);

    // assemble integral over elements
    for (unsigned int i_cell=0; i_cell<Model::mesh_->get_el_ds()->lsize(); i_cell++)
    {
    	typename DOFHandlerBase::CellIterator cell = Model::mesh_->element(feo->dh()->el_index(i_cell));
        if (cell->dim() != dim) continue;

        fe_values.reinit(cell);
        feo->dh()->get_dof_indices(cell, (unsigned int*)&(dof_indices[0]));
        ElementAccessor<3> ele_acc = cell->element_accessor();

        Model::compute_mass_matrix_coefficient(fe_values.point_list(), ele_acc, mm_coef);
        Model::compute_retardation_coefficient(fe_values.point_list(), ele_acc, ret_coef);

        for (unsigned int sbi=0; sbi<Model::n_substances(); ++sbi)
        {
            // assemble the local mass matrix
            for (unsigned int i=0; i<ndofs; i++)
            {
                for (unsigned int j=0; j<ndofs; j++)
                {
                    local_mass_matrix[i*ndofs+j] = 0;
                    for (unsigned int k=0; k<qsize; k++)
                        local_mass_matrix[i*ndofs+j] += (mm_coef[k]+ret_coef[sbi][k])*fe_values.shape_value(j,k)*fe_values.shape_value(i,k)*fe_values.JxW(k);
                }
            }

            if (Model::balance_ != nullptr)
            {
                for (unsigned int i=0; i<ndofs; i++)
                {
                        local_mass_balance_vector[i] = 0;
                        for (unsigned int k=0; k<qsize; k++)
                            local_mass_balance_vector[i] += mm_coef[k]*fe_values.shape_value(i,k)*fe_values.JxW(k);
                }
            }
            
            Model::balance_->add_mass_matrix_values(Model::subst_idx[sbi], ele_acc.region().bulk_idx(), dof_indices, local_mass_balance_vector);
			ls_dt[sbi]->mat_set_values(ndofs, &(dof_indices[0]), ndofs, &(dof_indices[0]), local_mass_matrix);
        }
    }
}




template<class Model>
void TransportDG<Model>::assemble_stiffness_matrix()
{
  START_TIMER("assemble_stiffness");
   START_TIMER("assemble_volume_integrals");
    assemble_volume_integrals<1>();
    assemble_volume_integrals<2>();
    assemble_volume_integrals<3>();
   END_TIMER("assemble_volume_integrals");

   START_TIMER("assemble_fluxes_boundary");
    assemble_fluxes_boundary<1>();
    assemble_fluxes_boundary<2>();
    assemble_fluxes_boundary<3>();
   END_TIMER("assemble_fluxes_boundary");

   START_TIMER("assemble_fluxes_elem_elem");
    assemble_fluxes_element_element<1>();
    assemble_fluxes_element_element<2>();
    assemble_fluxes_element_element<3>();
   END_TIMER("assemble_fluxes_elem_elem");

   START_TIMER("assemble_fluxes_elem_side");
    assemble_fluxes_element_side<1>();
    assemble_fluxes_element_side<2>();
    assemble_fluxes_element_side<3>();
   END_TIMER("assemble_fluxes_elem_side");
  END_TIMER("assemble_stiffness");
}



template<class Model>
template<unsigned int dim>
void TransportDG<Model>::assemble_volume_integrals()
{
    FEValues<dim,3> fv_rt(*feo->mapping<dim>(), *feo->q<dim>(), *feo->fe_rt<dim>(),
    		update_values | update_gradients);
    FEValues<dim,3> fe_values(*feo->mapping<dim>(), *feo->q<dim>(), *feo->fe<dim>(),
    		update_values | update_gradients | update_JxW_values | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
    unsigned int dof_indices[ndofs];
    vector<arma::vec3> velocity(qsize);
    vector<arma::vec> sources_sigma(qsize, arma::vec(Model::n_substances()));
    PetscScalar local_matrix[ndofs*ndofs];

	// assemble integral over elements
    for (unsigned int i_cell=0; i_cell<Model::mesh_->get_el_ds()->lsize(); i_cell++)
    {
    	typename DOFHandlerBase::CellIterator cell = Model::mesh_->element(feo->dh()->el_index(i_cell));
        if (cell->dim() != dim) continue;

        fe_values.reinit(cell);
        fv_rt.reinit(cell);
        ElementAccessor<3> ele_acc = cell->element_accessor();
        feo->dh()->get_dof_indices(cell, dof_indices);

        calculate_velocity(cell, velocity, fv_rt);
        Model::compute_advection_diffusion_coefficients(fe_values.point_list(), velocity, ele_acc, ad_coef, dif_coef);
        Model::compute_sources_sigma(fe_values.point_list(), ele_acc, sources_sigma);

        // assemble the local stiffness matrix
        for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
        {
        	for (unsigned int i=0; i<ndofs; i++)
        		for (unsigned int j=0; j<ndofs; j++)
        			local_matrix[i*ndofs+j] = 0;

        	for (unsigned int k=0; k<qsize; k++)
        	{
        		for (unsigned int i=0; i<ndofs; i++)
        		{
        			arma::vec3 Kt_grad_i = dif_coef[sbi][k].t()*fe_values.shape_grad(i,k);
        			double ad_dot_grad_i = arma::dot(ad_coef[sbi][k], fe_values.shape_grad(i,k));

        			for (unsigned int j=0; j<ndofs; j++)
						local_matrix[i*ndofs+j] += (arma::dot(Kt_grad_i, fe_values.shape_grad(j,k))
												   -fe_values.shape_value(j,k)*ad_dot_grad_i
												   +sources_sigma[k][sbi]*fe_values.shape_value(j,k)*fe_values.shape_value(i,k))*fe_values.JxW(k);
				}
			}
			ls[sbi]->mat_set_values(ndofs, (int *)dof_indices, ndofs, (int *)dof_indices, local_matrix);
        }
    }
}


template<class Model>
void TransportDG<Model>::set_sources()
{
  START_TIMER("assemble_sources");
    if (Model::balance_ != nullptr)
    	Model::balance_->start_source_assembly(Model::subst_idx);
	set_sources<1>();
	set_sources<2>();
	set_sources<3>();
	if (Model::balance_ != nullptr)
		Model::balance_->finish_source_assembly(Model::subst_idx);
  END_TIMER("assemble_sources");
}

template<class Model>
template<unsigned int dim>
void TransportDG<Model>::set_sources()
{
    FEValues<dim,3> fe_values(*feo->mapping<dim>(), *feo->q<dim>(), *feo->fe<dim>(),
    		update_values | update_JxW_values | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
    vector<arma::vec> sources_conc(qsize, arma::vec(Model::n_substances())),
    		sources_density(qsize, arma::vec(Model::n_substances())),
			sources_sigma(qsize, arma::vec(Model::n_substances()));
    vector<int> dof_indices(ndofs);
    PetscScalar local_rhs[ndofs];
    vector<PetscScalar> local_source_balance_vector(ndofs), local_source_balance_rhs(ndofs);
    double source;

	// assemble integral over elements
    for (unsigned int i_cell=0; i_cell<Model::mesh_->get_el_ds()->lsize(); i_cell++)
    {
    	typename DOFHandlerBase::CellIterator cell = Model::mesh_->element(feo->dh()->el_index(i_cell));
        if (cell->dim() != dim) continue;

        fe_values.reinit(cell);
        feo->dh()->get_dof_indices(cell, (unsigned int *)&(dof_indices[0]));

        Model::compute_source_coefficients(fe_values.point_list(), cell->element_accessor(), sources_conc, sources_density, sources_sigma);

        // assemble the local stiffness matrix
        for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
        {
        	fill_n(local_rhs, ndofs, 0);
        	local_source_balance_vector.assign(ndofs, 0);
        	local_source_balance_rhs.assign(ndofs, 0);

        	// compute sources
        	for (unsigned int k=0; k<qsize; k++)
        	{
        		source = (sources_density[k][sbi] + sources_conc[k][sbi]*sources_sigma[k][sbi])*fe_values.JxW(k);

        		for (unsigned int i=0; i<ndofs; i++)
        			local_rhs[i] += source*fe_values.shape_value(i,k);
        	}
        	ls[sbi]->rhs_set_values(ndofs, &(dof_indices[0]), local_rhs);

        	if (Model::balance_ != nullptr)
        	{
        		for (unsigned int i=0; i<ndofs; i++)
        		{
        			for (unsigned int k=0; k<qsize; k++)
        				local_source_balance_vector[i] -= sources_sigma[k][sbi]*fe_values.shape_value(i,k)*fe_values.JxW(k);

        			local_source_balance_rhs[i] += local_rhs[i];
        		}
        		Model::balance_->add_source_matrix_values(Model::subst_idx[sbi], cell->region().bulk_idx(), dof_indices, local_source_balance_vector);
        		Model::balance_->add_source_vec_values(Model::subst_idx[sbi], cell->region().bulk_idx(), dof_indices, local_source_balance_rhs);
        	}
        }
    }
}



template<class Model>
template<unsigned int dim>
void TransportDG<Model>::assemble_fluxes_element_element()
{
    vector<FESideValues<dim,3>*> fe_values;
    FESideValues<dim,3> fsv_rt(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe_rt<dim>(),
    		update_values);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim-1>()->size(),
    		n_max_sides = ad_coef_edg.size();
    vector<unsigned int*> side_dof_indices;
    PetscScalar local_matrix[ndofs*ndofs];
    vector<vector<arma::vec3> > side_velocity(n_max_sides);
    vector<arma::vec> dg_penalty(n_max_sides);
    double gamma_l, omega[2], transport_flux;

    for (unsigned int sid=0; sid<n_max_sides; sid++)
    {
    	side_dof_indices.push_back(new unsigned int[ndofs]);
    	fe_values.push_back(new FESideValues<dim,3>(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe<dim>(),
    			update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points));
    }

    // assemble integral over sides
    for (unsigned int iedg=0; iedg<feo->dh()->n_loc_edges(); iedg++)
    {
    	Edge *edg = &Model::mesh_->edges[feo->dh()->edge_index(iedg)];
        if (edg->n_sides < 2 || edg->side(0)->element()->dim() != dim) continue;

		for (int sid=0; sid<edg->n_sides; sid++)
		{
			typename DOFHandlerBase::CellIterator cell = edg->side(sid)->element();
			ElementAccessor<3> ele_acc = cell->element_accessor();
			feo->dh()->get_dof_indices(cell, side_dof_indices[sid]);
			fe_values[sid]->reinit(cell, edg->side(sid)->el_idx());
			fsv_rt.reinit(cell, edg->side(sid)->el_idx());
			calculate_velocity(cell, side_velocity[sid], fsv_rt);
			Model::compute_advection_diffusion_coefficients(fe_values[sid]->point_list(), side_velocity[sid], ele_acc, ad_coef_edg[sid], dif_coef_edg[sid]);
			dg_penalty[sid] = data_.dg_penalty.value(cell->centre(), ele_acc);
		}

        // fluxes and penalty
		for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
		{
			vector<double> fluxes(edg->n_sides);
			for (int sid=0; sid<edg->n_sides; sid++)
			{
				fluxes[sid] = 0;
				for (unsigned int k=0; k<qsize; k++)
					fluxes[sid] += arma::dot(ad_coef_edg[sid][sbi][k], fe_values[sid]->normal_vector(k))*fe_values[sid]->JxW(k);
				fluxes[sid] /= edg->side(sid)->measure();
			}

			for (int s1=0; s1<edg->n_sides; s1++)
			{
				for (int s2=s1+1; s2<edg->n_sides; s2++)
				{
					OLD_ASSERT(edg->side(s1)->valid(), "Invalid side of edge.");
					if (!feo->dh()->el_is_local(edg->side(s1)->element().index())
							&& !feo->dh()->el_is_local(edg->side(s2)->element().index())) continue;

					arma::vec3 nv = fe_values[s1]->normal_vector(0);

					// set up the parameters for DG method
					set_DG_parameters_edge(*edg, s1, s2, qsize, dif_coef_edg[s1][sbi], dif_coef_edg[s2][sbi], fluxes, fe_values[0]->normal_vector(0), dg_penalty[s1][sbi], dg_penalty[s2][sbi], gamma_l, omega, transport_flux);

					int sd[2];
					sd[0] = s1;
					sd[1] = s2;

#define AVERAGE(i,k,side_id)  (fe_values[sd[side_id]]->shape_value(i,k)*0.5)
#define WAVERAGE(i,k,side_id) (arma::dot(dif_coef_edg[sd[side_id]][sbi][k]*fe_values[sd[side_id]]->shape_grad(i,k),nv)*omega[side_id])
#define JUMP(i,k,side_id)     ((side_id==0?1:-1)*fe_values[sd[side_id]]->shape_value(i,k))

					// For selected pair of elements:
					for (int n=0; n<2; n++)
					{
						if (!feo->dh()->el_is_local(edg->side(sd[n])->element().index())) continue;

						for (int m=0; m<2; m++)
						{
							for (unsigned int i=0; i<fe_values[sd[n]]->n_dofs(); i++)
								for (unsigned int j=0; j<fe_values[sd[m]]->n_dofs(); j++)
									local_matrix[i*fe_values[sd[m]]->n_dofs()+j] = 0;

							for (unsigned int k=0; k<qsize; k++)
							{
								double flux_times_JxW = transport_flux*fe_values[0]->JxW(k);
								double gamma_times_JxW = gamma_l*fe_values[0]->JxW(k);

								for (unsigned int i=0; i<fe_values[sd[n]]->n_dofs(); i++)
								{
									double flux_JxW_jump_i = flux_times_JxW*JUMP(i,k,n);
									double gamma_JxW_jump_i = gamma_times_JxW*JUMP(i,k,n);
									double JxW_jump_i = fe_values[0]->JxW(k)*JUMP(i,k,n);
									double JxW_var_wavg_i = fe_values[0]->JxW(k)*WAVERAGE(i,k,n)*dg_variant;

									for (unsigned int j=0; j<fe_values[sd[m]]->n_dofs(); j++)
									{
										int index = i*fe_values[sd[m]]->n_dofs()+j;

										// flux due to transport (applied on interior edges) (average times jump)
										local_matrix[index] += flux_JxW_jump_i*AVERAGE(j,k,m);

										// penalty enforcing continuity across edges (applied on interior and Dirichlet edges) (jump times jump)
										local_matrix[index] += gamma_JxW_jump_i*JUMP(j,k,m);

										// terms due to diffusion
										local_matrix[index] -= WAVERAGE(j,k,m)*JxW_jump_i;
										local_matrix[index] -= JUMP(j,k,m)*JxW_var_wavg_i;
									}
								}
							}
							ls[sbi]->mat_set_values(fe_values[sd[n]]->n_dofs(), (int *)side_dof_indices[sd[n]], fe_values[sd[m]]->n_dofs(), (int *)side_dof_indices[sd[m]], local_matrix);
						}
					}
#undef AVERAGE
#undef WAVERAGE
#undef JUMP
				}
			}
		}
    }

    for (unsigned int i=0; i<n_max_sides; i++)
    {
        delete fe_values[i];
        delete[] side_dof_indices[i];
    }
}


template<class Model>
template<unsigned int dim>
void TransportDG<Model>::assemble_fluxes_boundary()
{
    FESideValues<dim,3> fe_values_side(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe<dim>(),
    		update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
    FESideValues<dim,3> fsv_rt(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe_rt<dim>(),
    		update_values);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim-1>()->size();
    unsigned int side_dof_indices[ndofs];
    PetscScalar local_matrix[ndofs*ndofs];
    vector<arma::vec3> side_velocity;
    vector<double> robin_sigma(qsize);
    vector<double> csection(qsize);
    arma::vec dg_penalty;
    double gamma_l;

    // assemble boundary integral
    for (unsigned int iedg=0; iedg<feo->dh()->n_loc_edges(); iedg++)
    {
    	Edge *edg = &Model::mesh_->edges[feo->dh()->edge_index(iedg)];
    	if (edg->n_sides > 1) continue;
    	// check spatial dimension
    	if (edg->side(0)->dim() != dim-1) continue;
    	// skip edges lying not on the boundary
    	if (edg->side(0)->cond() == NULL) continue;

    	SideIter side = edg->side(0);
        typename DOFHandlerBase::CellIterator cell = side->element();
        ElementAccessor<3> ele_acc = cell->element_accessor();
        feo->dh()->get_dof_indices(cell, side_dof_indices);
        fe_values_side.reinit(cell, side->el_idx());
        fsv_rt.reinit(cell, side->el_idx());

        calculate_velocity(cell, side_velocity, fsv_rt);
        Model::compute_advection_diffusion_coefficients(fe_values_side.point_list(), side_velocity, ele_acc, ad_coef, dif_coef);
        dg_penalty = data_.dg_penalty.value(cell->centre(), ele_acc);
        arma::uvec bc_type;
        Model::get_bc_type(side->cond()->element_accessor(), bc_type);
        data_.cross_section.value_list(fe_values_side.point_list(), ele_acc, csection);

        for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
        {
        	for (unsigned int i=0; i<ndofs; i++)
        		for (unsigned int j=0; j<ndofs; j++)
        			local_matrix[i*ndofs+j] = 0;

			// On Neumann boundaries we have only term from integrating by parts the advective term,
			// on Dirichlet boundaries we additionally apply the penalty which enforces the prescribed value.
			double side_flux = 0;
			for (unsigned int k=0; k<qsize; k++)
				side_flux += arma::dot(ad_coef[sbi][k], fe_values_side.normal_vector(k))*fe_values_side.JxW(k);
			double transport_flux = side_flux/side->measure();

			if (bc_type[sbi] == AdvectionDiffusionModel::abc_dirichlet)
			{
				// set up the parameters for DG method
				set_DG_parameters_boundary(side, qsize, dif_coef[sbi], transport_flux, fe_values_side.normal_vector(0), dg_penalty[sbi], gamma_l);
				gamma[sbi][side->cond_idx()] = gamma_l;
				transport_flux += gamma_l;
			}

			// fluxes and penalty
			for (unsigned int k=0; k<qsize; k++)
			{
				double flux_times_JxW;
				if (bc_type[sbi] == AdvectionDiffusionModel::abc_total_flux)
                {
                    Model::get_flux_bc_sigma(sbi, fe_values_side.point_list(), side->cond()->element_accessor(), robin_sigma);
					flux_times_JxW = csection[k]*robin_sigma[k]*fe_values_side.JxW(k);
                }
				else if (bc_type[sbi] == AdvectionDiffusionModel::abc_diffusive_flux)
                {
                    Model::get_flux_bc_sigma(sbi, fe_values_side.point_list(), side->cond()->element_accessor(), robin_sigma);
					flux_times_JxW = (transport_flux + csection[k]*robin_sigma[k])*fe_values_side.JxW(k);
                }
				else if (bc_type[sbi] == AdvectionDiffusionModel::abc_inflow && side_flux < 0)
					flux_times_JxW = 0;
				else
					flux_times_JxW = transport_flux*fe_values_side.JxW(k);

				for (unsigned int i=0; i<ndofs; i++)
				{
					for (unsigned int j=0; j<ndofs; j++)
					{
					    // flux due to advection and penalty
						local_matrix[i*ndofs+j] += flux_times_JxW*fe_values_side.shape_value(i,k)*fe_values_side.shape_value(j,k);

						// flux due to diffusion (only on dirichlet and inflow boundary)
						if (bc_type[sbi] == AdvectionDiffusionModel::abc_dirichlet)
							local_matrix[i*ndofs+j] -= (arma::dot(dif_coef[sbi][k]*fe_values_side.shape_grad(j,k),fe_values_side.normal_vector(k))*fe_values_side.shape_value(i,k)
									+ arma::dot(dif_coef[sbi][k]*fe_values_side.shape_grad(i,k),fe_values_side.normal_vector(k))*fe_values_side.shape_value(j,k)*dg_variant
									)*fe_values_side.JxW(k);
					}
				}
			}

			ls[sbi]->mat_set_values(ndofs, (int *)side_dof_indices, ndofs, (int *)side_dof_indices, local_matrix);
        }
    }
}


template<class Model>
template<unsigned int dim>
void TransportDG<Model>::assemble_fluxes_element_side()
{

	if (dim == 1) return;
    FEValues<dim-1,3> fe_values_vb(*feo->mapping<dim-1>(), *feo->q<dim-1>(), *feo->fe<dim-1>(),
    		update_values | update_gradients | update_JxW_values | update_quadrature_points);
    FESideValues<dim,3> fe_values_side(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe<dim>(),
    		update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
    FESideValues<dim,3> fsv_rt(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe_rt<dim>(),
       		update_values);
    FEValues<dim-1,3> fv_rt(*feo->mapping<dim-1>(), *feo->q<dim-1>(), *feo->fe_rt<dim-1>(),
       		update_values);

    vector<FEValuesSpaceBase<3>*> fv_sb(2);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs();    // number of local dofs
    const unsigned int qsize = feo->q<dim-1>()->size();     // number of quadrature points
    unsigned int side_dof_indices[2*ndofs], n_dofs[2];
	vector<arma::vec3> velocity_higher, velocity_lower;
	vector<arma::vec> frac_sigma(qsize, arma::vec(Model::n_substances()));
	vector<double> csection_lower(qsize), csection_higher(qsize), por_lower(qsize), por_higher(qsize);
    PetscScalar local_matrix[4*ndofs*ndofs];
    double comm_flux[2][2];

    // index 0 = element with lower dimension,
    // index 1 = side of element with higher dimension
    fv_sb[0] = &fe_values_vb;
    fv_sb[1] = &fe_values_side;

    // assemble integral over sides
    for (unsigned int inb=0; inb<feo->dh()->n_loc_nb(); inb++)
    {
    	Neighbour *nb = &Model::mesh_->vb_neighbours_[feo->dh()->nb_index(inb)];
        // skip neighbours of different dimension
        if (nb->element()->dim() != dim-1) continue;

		typename DOFHandlerBase::CellIterator cell_sub = Model::mesh_->element.full_iter(nb->element());
		feo->dh()->get_dof_indices(cell_sub, side_dof_indices);
		fe_values_vb.reinit(cell_sub);
		n_dofs[0] = fv_sb[0]->n_dofs();

		typename DOFHandlerBase::CellIterator cell = nb->side()->element();
		feo->dh()->get_dof_indices(cell, side_dof_indices+n_dofs[0]);
		fe_values_side.reinit(cell, nb->side()->el_idx());
		n_dofs[1] = fv_sb[1]->n_dofs();

		// Element id's for testing if they belong to local partition.
		int element_id[2];
		element_id[0] = cell_sub.index();
		element_id[1] = cell.index();

		fsv_rt.reinit(cell, nb->side()->el_idx());
		fv_rt.reinit(cell_sub);
		calculate_velocity(cell, velocity_higher, fsv_rt);
		calculate_velocity(cell_sub, velocity_lower, fv_rt);
		Model::compute_advection_diffusion_coefficients(fe_values_vb.point_list(), velocity_lower, cell_sub->element_accessor(), ad_coef_edg[0], dif_coef_edg[0]);
		Model::compute_advection_diffusion_coefficients(fe_values_vb.point_list(), velocity_higher, cell->element_accessor(), ad_coef_edg[1], dif_coef_edg[1]);
		data_.cross_section.value_list(fe_values_vb.point_list(), cell_sub->element_accessor(), csection_lower);
		data_.cross_section.value_list(fe_values_vb.point_list(), cell->element_accessor(), csection_higher);
		data_.fracture_sigma.value_list(fe_values_vb.point_list(), cell_sub->element_accessor(), frac_sigma);
		data_.porosity.value_list(fe_values_vb.point_list(), cell_sub->element_accessor(), por_lower);
		data_.porosity.value_list(fe_values_vb.point_list(), cell->element_accessor(), por_higher);

		for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
		{
			for (unsigned int i=0; i<n_dofs[0]+n_dofs[1]; i++)
				for (unsigned int j=0; j<n_dofs[0]+n_dofs[1]; j++)
					local_matrix[i*(n_dofs[0]+n_dofs[1])+j] = 0;

			// set transmission conditions
			for (unsigned int k=0; k<qsize; k++)
			{
				/* The communication flux has two parts:
				 * - "diffusive" term containing sigma
				 * - "advective" term representing usual upwind
				 *
				 * The calculation differs from the reference manual, since ad_coef and dif_coef have different meaning
				 * than b and A in the manual.
				 * In calculation of sigma there appears one more csection_lower in the denominator.
				 */
				double sigma = frac_sigma[k][sbi]*arma::dot(dif_coef_edg[0][sbi][k]*fe_values_side.normal_vector(k),fe_values_side.normal_vector(k))*
				        2*csection_higher[k]*csection_higher[k]/(csection_lower[k]*csection_lower[k]);

				double transport_flux = arma::dot(ad_coef_edg[1][sbi][k], fe_values_side.normal_vector(k));
				double por_lower_over_higher = por_lower[k]/por_higher[k];

				comm_flux[0][0] =  (sigma-min(0.,transport_flux*por_lower_over_higher))*fv_sb[0]->JxW(k);
				comm_flux[0][1] = -(sigma-min(0.,transport_flux*por_lower_over_higher))*fv_sb[0]->JxW(k);
				comm_flux[1][0] = -(sigma+max(0.,transport_flux))*fv_sb[0]->JxW(k);
				comm_flux[1][1] =  (sigma+max(0.,transport_flux))*fv_sb[0]->JxW(k);

				for (int n=0; n<2; n++)
				{
					if (!feo->dh()->el_is_local(element_id[n])) continue;

					for (unsigned int i=0; i<n_dofs[n]; i++)
						for (int m=0; m<2; m++)
							for (unsigned int j=0; j<n_dofs[m]; j++)
								local_matrix[(i+n*n_dofs[0])*(n_dofs[0]+n_dofs[1]) + m*n_dofs[0] + j] +=
								        comm_flux[m][n]*fv_sb[m]->shape_value(j,k)*fv_sb[n]->shape_value(i,k);
				}
			}
			ls[sbi]->mat_set_values(n_dofs[0]+n_dofs[1], (int *)side_dof_indices, n_dofs[0]+n_dofs[1], (int *)side_dof_indices, local_matrix);
		}
    }

}






template<class Model>
void TransportDG<Model>::set_boundary_conditions()
{
  START_TIMER("assemble_bc");
    if (Model::balance_ != nullptr)
    	Model::balance_->start_flux_assembly(Model::subst_idx);
	set_boundary_conditions<1>();
	set_boundary_conditions<2>();
	set_boundary_conditions<3>();
	if (Model::balance_ != nullptr)
		Model::balance_->finish_flux_assembly(Model::subst_idx);
  END_TIMER("assemble_bc");
}


template<class Model>
template<unsigned int dim>
void TransportDG<Model>::set_boundary_conditions()
{
    FESideValues<dim,3> fe_values_side(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe<dim>(),
    		update_values | update_gradients | update_normal_vectors | update_side_JxW_values | update_quadrature_points);
    FESideValues<dim,3> fsv_rt(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe_rt<dim>(),
           		update_values);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim-1>()->size();
    vector<int> side_dof_indices(ndofs);
    unsigned int loc_b=0;
    double local_rhs[ndofs];
    vector<PetscScalar> local_flux_balance_vector(ndofs);
    PetscScalar local_flux_balance_rhs;
    vector<arma::vec> bc_values(qsize, arma::vec(Model::n_substances()));
    vector<double> bc_fluxes(qsize),
			bc_sigma(qsize),
			bc_ref_values(qsize);
    vector<double> csection(qsize);
	vector<arma::vec3> velocity;

    for (unsigned int loc_el = 0; loc_el < Model::mesh_->get_el_ds()->lsize(); loc_el++)
    {
        ElementFullIter elm = Model::mesh_->element(feo->dh()->el_index(loc_el));
        if (elm->boundary_idx_ == nullptr) continue;

        FOR_ELEMENT_SIDES(elm,si)
        {
			Edge *edg = elm->side(si)->edge();
			if (edg->n_sides > 1) continue;
			// skip edges lying not on the boundary
			if (edg->side(0)->cond() == NULL) continue;

			if (edg->side(0)->dim() != dim-1)
			{
				if (edg->side(0)->cond() != nullptr) ++loc_b;
				continue;
			}

			SideIter side = edg->side(0);
			typename DOFHandlerBase::CellIterator cell = Model::mesh_->element.full_iter(side->element());
			ElementAccessor<3> ele_acc = side->cond()->element_accessor();

			arma::uvec bc_type;
			Model::get_bc_type(ele_acc, bc_type);

			fe_values_side.reinit(cell, side->el_idx());
			fsv_rt.reinit(cell, side->el_idx());
			calculate_velocity(cell, velocity, fsv_rt);

			Model::compute_advection_diffusion_coefficients(fe_values_side.point_list(), velocity, side->element()->element_accessor(), ad_coef, dif_coef);
			data_.cross_section.value_list(fe_values_side.point_list(), side->element()->element_accessor(), csection);
			// The b.c. data are fetched for all possible b.c. types since we allow
			// different bc_type for each substance.
			data_.bc_dirichlet_value.value_list(fe_values_side.point_list(), ele_acc, bc_values);

			feo->dh()->get_dof_indices(cell, (unsigned int *)&(side_dof_indices[0]));

			for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
			{
				fill_n(local_rhs, ndofs, 0);
				local_flux_balance_vector.assign(ndofs, 0);
				local_flux_balance_rhs = 0;

				double side_flux = 0;
				for (unsigned int k=0; k<qsize; k++)
					side_flux += arma::dot(ad_coef[sbi][k], fe_values_side.normal_vector(k))*fe_values_side.JxW(k);
				double transport_flux = side_flux/side->measure();

				if (bc_type[sbi] == AdvectionDiffusionModel::abc_inflow && side_flux < 0)
				{
					for (unsigned int k=0; k<qsize; k++)
					{
						double bc_term = -transport_flux*bc_values[k][sbi]*fe_values_side.JxW(k);
						for (unsigned int i=0; i<ndofs; i++)
							local_rhs[i] += bc_term*fe_values_side.shape_value(i,k);
					}
					if (Model::balance_ != nullptr)
						for (unsigned int i=0; i<ndofs; i++)
							local_flux_balance_rhs -= local_rhs[i];
				}
				else if (bc_type[sbi] == AdvectionDiffusionModel::abc_dirichlet)
				{
					for (unsigned int k=0; k<qsize; k++)
					{
						double bc_term = gamma[sbi][side->cond_idx()]*bc_values[k][sbi]*fe_values_side.JxW(k);
						arma::vec3 bc_grad = -bc_values[k][sbi]*fe_values_side.JxW(k)*dg_variant*(arma::trans(dif_coef[sbi][k])*fe_values_side.normal_vector(k));
						for (unsigned int i=0; i<ndofs; i++)
							local_rhs[i] += bc_term*fe_values_side.shape_value(i,k)
									+ arma::dot(bc_grad,fe_values_side.shape_grad(i,k));
					}
					if (Model::balance_ != nullptr)
					{
						for (unsigned int k=0; k<qsize; k++)
						{
							for (unsigned int i=0; i<ndofs; i++)
							{
								local_flux_balance_vector[i] += (arma::dot(ad_coef[sbi][k], fe_values_side.normal_vector(k))*fe_values_side.shape_value(i,k)
										- arma::dot(dif_coef[sbi][k]*fe_values_side.shape_grad(i,k),fe_values_side.normal_vector(k))
										+ gamma[sbi][side->cond_idx()]*fe_values_side.shape_value(i,k))*fe_values_side.JxW(k);
							}
						}
						if (Model::time_->tlevel() > 0)
							for (unsigned int i=0; i<ndofs; i++)
								local_flux_balance_rhs -= local_rhs[i];
					}
				}
				else if (bc_type[sbi] == AdvectionDiffusionModel::abc_total_flux)
				{
                    Model::get_flux_bc_data(sbi, fe_values_side.point_list(), ele_acc, bc_fluxes, bc_sigma, bc_ref_values);
					for (unsigned int k=0; k<qsize; k++)
					{
						double bc_term = csection[k]*(bc_sigma[k]*bc_ref_values[k]+bc_fluxes[k])*fe_values_side.JxW(k);
						for (unsigned int i=0; i<ndofs; i++)
							local_rhs[i] += bc_term*fe_values_side.shape_value(i,k);
					}

					if (Model::balance_ != nullptr)
					{
						for (unsigned int i=0; i<ndofs; i++)
						{
							for (unsigned int k=0; k<qsize; k++)
								local_flux_balance_vector[i] += csection[k]*bc_sigma[k]*fe_values_side.JxW(k)*fe_values_side.shape_value(i,k);
							local_flux_balance_rhs -= local_rhs[i];
						}
					}
				}
				else if (bc_type[sbi] == AdvectionDiffusionModel::abc_diffusive_flux)
				{
                    Model::get_flux_bc_data(sbi, fe_values_side.point_list(), ele_acc, bc_fluxes, bc_sigma, bc_ref_values);
					for (unsigned int k=0; k<qsize; k++)
					{
						double bc_term = csection[k]*(bc_sigma[k]*bc_ref_values[k]+bc_fluxes[k])*fe_values_side.JxW(k);
						for (unsigned int i=0; i<ndofs; i++)
							local_rhs[i] += bc_term*fe_values_side.shape_value(i,k);
					}

					if (Model::balance_ != nullptr)
					{
						for (unsigned int i=0; i<ndofs; i++)
						{
							for (unsigned int k=0; k<qsize; k++)
								local_flux_balance_vector[i] += csection[k]*(arma::dot(ad_coef[sbi][k], fe_values_side.normal_vector(k)) + bc_sigma[k])*fe_values_side.JxW(k)*fe_values_side.shape_value(i,k);
							local_flux_balance_rhs -= local_rhs[i];
						}
					}
				}
				else if (bc_type[sbi] == AdvectionDiffusionModel::abc_inflow && side_flux >= 0)
				{
					if (Model::balance_ != nullptr)
					{
						for (unsigned int k=0; k<qsize; k++)
						{
							for (unsigned int i=0; i<ndofs; i++)
								local_flux_balance_vector[i] += arma::dot(ad_coef[sbi][k], fe_values_side.normal_vector(k))*fe_values_side.JxW(k)*fe_values_side.shape_value(i,k);
						}
					}
				}
				ls[sbi]->rhs_set_values(ndofs, &(side_dof_indices[0]), local_rhs);

				if (Model::balance_ != nullptr)
				{
					Model::balance_->add_flux_matrix_values(Model::subst_idx[sbi], loc_b, side_dof_indices, local_flux_balance_vector);
					Model::balance_->add_flux_vec_value(Model::subst_idx[sbi], loc_b, local_flux_balance_rhs);
				}
			}
			++loc_b;
        }
    }
}



template<class Model>
template<unsigned int dim>
void TransportDG<Model>::calculate_velocity(const typename DOFHandlerBase::CellIterator &cell, vector<arma::vec3> &velocity, FEValuesBase<dim,3> &fv)
{
	OLD_ASSERT(cell->dim() == dim, "Element dimension mismatch!");

    velocity.resize(fv.n_points());

    for (unsigned int k=0; k<fv.n_points(); k++)
    {
        velocity[k].zeros();
        for (unsigned int sid=0; sid<cell->n_sides(); sid++)
            velocity[k] += fv.shape_vector(sid,k) * Model::mh_dh->side_flux( *(cell->side(sid)) );
    }
}





template<class Model>
void TransportDG<Model>::set_DG_parameters_edge(const Edge &edg,
            const int s1,
            const int s2,
            const int K_size,
            const vector<arma::mat33> &K1,
            const vector<arma::mat33> &K2,
            const vector<double> &fluxes,
            const arma::vec3 &normal_vector,
            const double alpha1,
            const double alpha2,
            double &gamma,
            double *omega,
            double &transport_flux)
{
    double delta[2];
    double h = 0;
    double local_alpha = 0;

    OLD_ASSERT(edg.side(s1)->valid(), "Invalid side of an edge.");
    SideIter s = edg.side(s1);

    // calculate the side diameter
    if (s->dim() == 0)
    {
        h = 1;
    }
    else
    {
        for (unsigned int i=0; i<s->n_nodes(); i++)
            for (unsigned int j=i+1; j<s->n_nodes(); j++)
                h = max(h, s->node(i)->distance(*s->node(j)));
    }

    // calculate the total in- and out-flux through the edge
    double pflux = 0, nflux = 0;
    for (int i=0; i<edg.n_sides; i++)
    {
        if (fluxes[i] > 0)
            pflux += fluxes[i];
        else
            nflux += fluxes[i];
    }

    // calculate the flux from s1 to s2
    if (fluxes[s2] > 0 && fluxes[s1] < 0 && s1 < s2)
        transport_flux = fluxes[s1]*fabs(fluxes[s2]/pflux);
    else if (fluxes[s2] < 0 && fluxes[s1] > 0 && s1 < s2)
        transport_flux = fluxes[s1]*fabs(fluxes[s2]/nflux);
    else if (s1==s2)
        transport_flux = fluxes[s1];
    else
        transport_flux = 0;

    gamma = 0.5*fabs(transport_flux);


    // determine local DG penalty
    local_alpha = max(alpha1, alpha2);

    if (s1 == s2)
    {
        omega[0] = 1;

        // delta is set to the average value of Kn.n on the side
        delta[0] = 0;
        for (int k=0; k<K_size; k++)
            delta[0] += dot(K1[k]*normal_vector,normal_vector);
        delta[0] /= K_size;

        gamma += local_alpha/h*(delta[0]);
    }
    else
    {
        delta[0] = 0;
        delta[1] = 0;
        for (int k=0; k<K_size; k++)
        {
            delta[0] += dot(K1[k]*normal_vector,normal_vector);
            delta[1] += dot(K2[k]*normal_vector,normal_vector);
        }
        delta[0] /= K_size;
        delta[1] /= K_size;

        double delta_sum = delta[0] + delta[1];

        if (delta_sum > numeric_limits<double>::epsilon())
        {
            omega[0] = delta[1]/delta_sum;
            omega[1] = delta[0]/delta_sum;
            gamma += local_alpha/h*(delta[0]*delta[1]/delta_sum);
        }
        else
            for (int i=0; i<2; i++) omega[i] = 0;
    }
}






template<class Model>
void TransportDG<Model>::set_DG_parameters_boundary(const SideIter side,
		    const int K_size,
            const vector<arma::mat33> &K,
            const double flux,
            const arma::vec3 &normal_vector,
            const double alpha,
            double &gamma)
{
    double delta = 0, h = 0;

    // calculate the side diameter
    if (side->dim() == 0)
    {
        h = 1;
    }
    else
    {
        for (unsigned int i=0; i<side->n_nodes(); i++)
            for (unsigned int j=i+1; j<side->n_nodes(); j++)
                h = max(h, side->node(i)->distance( *side->node(j) ));
    }

	// delta is set to the average value of Kn.n on the side
	for (int k=0; k<K_size; k++)
		delta += dot(K[k]*normal_vector,normal_vector);
	delta /= K_size;

	gamma = 0.5*fabs(flux) + alpha/h*delta;
}





template<class Model>
void TransportDG<Model>::set_initial_condition()
{
	START_TIMER("set_init_cond");
	for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
		ls[sbi]->start_allocation();
	prepare_initial_condition<1>();
	prepare_initial_condition<2>();
	prepare_initial_condition<3>();

	for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
		ls[sbi]->start_add_assembly();
	prepare_initial_condition<1>();
	prepare_initial_condition<2>();
	prepare_initial_condition<3>();

	for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
	{
		ls[sbi]->finish_assembly();
		ls[sbi]->solve();
	}
	END_TIMER("set_init_cond");
}

template<class Model>
template<unsigned int dim>
void TransportDG<Model>::prepare_initial_condition()
{
	FEValues<dim,3> fe_values(*feo->mapping<dim>(), *feo->q<dim>(), *feo->fe<dim>(),
			update_values | update_JxW_values | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
    unsigned int dof_indices[ndofs];
    double matrix[ndofs*ndofs], rhs[ndofs];
    std::vector<arma::vec> init_values(qsize);

    for (unsigned int k=0; k<qsize; k++)
    	init_values[k].resize(Model::n_substances());

    for (unsigned int i_cell=0; i_cell<Model::mesh_->get_el_ds()->lsize(); i_cell++)
    {
    	typename DOFHandlerBase::CellIterator elem = Model::mesh_->element(feo->dh()->el_index(i_cell));
    	if (elem->dim() != dim) continue;

    	ElementAccessor<3> ele_acc = elem->element_accessor();
    	feo->dh()->get_dof_indices(elem, dof_indices);
    	fe_values.reinit(elem);

   		Model::compute_init_cond(fe_values.point_list(), ele_acc, init_values);

    	for (unsigned int sbi=0; sbi<Model::n_substances(); sbi++)
    	{
    		for (unsigned int i=0; i<ndofs; i++)
    		{
    			rhs[i] = 0;
    			for (unsigned int j=0; j<ndofs; j++)
    				matrix[i*ndofs+j] = 0;
    		}

    		for (unsigned int k=0; k<qsize; k++)
    		{
    			double rhs_term = init_values[k](sbi)*fe_values.JxW(k);

    			for (unsigned int i=0; i<ndofs; i++)
    			{
    				for (unsigned int j=0; j<ndofs; j++)
    					matrix[i*ndofs+j] += fe_values.shape_value(i,k)*fe_values.shape_value(j,k)*fe_values.JxW(k);

    				rhs[i] += fe_values.shape_value(i,k)*rhs_term;
    			}
    		}
    		ls[sbi]->set_values(ndofs, (int *)dof_indices, ndofs, (int *)dof_indices, matrix, rhs);
    	}
    }
}


template<class Model>
void TransportDG<Model>::get_par_info(int * &el_4_loc, Distribution * &el_ds)
{
	el_4_loc = Model::mesh_->get_el_4_loc();
	el_ds = Model::mesh_->get_el_ds();
}


template<class Model>
void TransportDG<Model>::update_after_reactions(bool solution_changed)
{
	if (solution_changed)
	{
		for (unsigned int i_cell=0; i_cell<Model::mesh_->get_el_ds()->lsize(); i_cell++)
		{
			typename DOFHandlerBase::CellIterator elem = Model::mesh_->element(feo->dh()->el_index(i_cell));

			unsigned int n_dofs;
			switch (elem->dim())
			{
			case 1:
				n_dofs = feo->fe<1>()->n_dofs();
				break;
			case 2:
				n_dofs = feo->fe<2>()->n_dofs();
				break;
			case 3:
				n_dofs = feo->fe<3>()->n_dofs();
				break;
			}

			unsigned int dof_indices[n_dofs];
			feo->dh()->get_dof_indices(elem, dof_indices);

			for (unsigned int sbi=0; sbi<Model::n_substances(); ++sbi)
			{
				double old_average = 0;
				for (unsigned int j=0; j<n_dofs; ++j)
					old_average += ls[sbi]->get_solution_array()[dof_indices[j]-feo->dh()->distr()->begin()];
				old_average /= n_dofs;

				for (unsigned int j=0; j<n_dofs; ++j)
					ls[sbi]->get_solution_array()[dof_indices[j]-feo->dh()->distr()->begin()] += solution_elem_[sbi][i_cell] - old_average;
			}
		}
	}
    // update mass_vec for the case that mass matrix changes in next time step
    for (unsigned int sbi=0; sbi<Model::n_substances(); ++sbi)
    	MatMult(*(ls_dt[sbi]->get_matrix()), ls[sbi]->get_solution(), mass_vec[sbi]);
}

template<class Model>
int *TransportDG<Model>::get_row_4_el()
{
	return Model::mesh_->get_row_4_el();
}








template class TransportDG<ConcentrationTransportModel>;
template class TransportDG<HeatTransferModel>;




