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
 * @file    darcy_flow_mh.cc
 * @ingroup flow
 * @brief   Setup and solve linear system of mixed-hybrid discretization of the linear
 *          porous media flow with possible preferential flow in fractures and chanels.
 */

//#include <limits>
#include <vector>
//#include <iostream>
//#include <iterator>
//#include <algorithm>
#include <armadillo>

#include "petscmat.h"
#include "petscviewer.h"
#include "petscerror.h"
#include "mpi.h"

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "input/factory.hh"

#include "mesh/mesh.h"
#include "mesh/intersection.hh"
#include "mesh/partitioning.hh"
#include "la/distribution.hh"
#include "la/linsys.hh"
#include "la/linsys_PETSC.hh"
#include "la/linsys_BDDC.hh"
#include "la/schur.hh"
//#include "la/sparse_graph.hh"
#include "la/local_to_global_map.hh"

#include "flow/darcy_flow_mh.hh"

#include "flow/darcy_flow_mh_output.hh"

/*
#include "fem/mapping_p1.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/fe_rt.hh"
#include "quadrature/quadrature_lib.hh"
*/

#include "tools/time_governor.hh"
#include "fields/field_algo_base.hh"
#include "fields/field.hh"
#include "fields/field_values.hh"
#include "fields/field_add_potential.hh"

#include "coupling/balance.hh"

#include "fields/vec_seq_double.hh"
#include "darcy_flow_assembly.hh"


FLOW123D_FORCE_LINK_IN_CHILD(darcy_flow_mh);




namespace it = Input::Type;

const it::Selection & DarcyMH::get_mh_mortar_selection() {
	return it::Selection("MH_MortarMethod")
		.add_value(NoMortar, "None", "Mortar space: P0 on elements of lower dimension.")
		.add_value(MortarP0, "P0", "Mortar space: P0 on elements of lower dimension.")
		.add_value(MortarP1, "P1", "Mortar space: P1 on intersections, using non-conforming pressures.")
		.close();
}


const it::Selection & DarcyMH::EqData::get_bc_type_selection() {
	return it::Selection("Flow_Darcy_BC_Type")
        .add_value(none, "none",
            "Homogeneous Neumann boundary condition. Zero flux")
        .add_value(dirichlet, "dirichlet",
            "Dirichlet boundary condition. "
            "Specify the pressure head through the 'bc_pressure' field "
            "or the piezometric head through the 'bc_piezo_head' field.")
        .add_value(total_flux, "total_flux", "Flux boundary condition (combines Neumann and Robin type). "
            "Water inflow equal to (($q^N + \\sigma (h^R - h)$)). "
            "Specify the water inflow by the 'bc_flux' field, the transition coefficient by 'bc_robin_sigma' "
            "and the reference pressure head or pieozmetric head through 'bc_pressure' or 'bc_piezo_head' respectively.")
        .add_value(seepage, "seepage",
            "Seepage face boundary condition. Pressure and inflow bounded from above. Boundary with potential seepage flow "
            "is described by the pair of inequalities:"
            "(($h \\le h_d^D$)) and (($ q \\le q_d^N$)), where the equality holds in at least one of them. Caution! Setting $q_d^N$ strictly negative"
            "may lead to an ill posed problem since a positive outflow is enforced."
            "Parameters (($h_d^D$)) and (($q_d^N$)) are given by fields ``bc_pressure`` (or ``bc_piezo_head``) and ``bc_flux`` respectively."
            )
        .add_value(river, "river",
            "River boundary condition. For the water level above the bedrock, (($H > H^S$)), the Robin boundary condition is used with the inflow given by: "
            "(( $q^N + \\sigma(H^D - H)$)). For the water level under the bedrock, constant infiltration is used "
            "(( $q^N + \\sigma(H^D - H^S)$)). Parameters: ``bc_pressure``, ``bc_switch_pressure``,"
            " ``bc_sigma, ``bc_flux``."
            )
        .close();
}

const it::Record & DarcyMH::type_field_descriptor() {

        const it::Record &field_descriptor =
        it::Record("Flow_Darcy_MH_Data",FieldCommon::field_descriptor_record_description("Flow_Darcy_MH_Data") )
        .copy_keys( DarcyMH::EqData().make_field_descriptor_type("Flow_Darcy_MH_Data_aux") )
            .declare_key("bc_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(),
                    "Boundary piezometric head for BC types: dirichlet, robin, and river." )
            .declare_key("bc_switch_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(),
                    "Boundary switch piezometric head for BC types: seepage, river." )
            .declare_key("init_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(),
                    "Initial condition for the pressure given as the piezometric head." )
            .close();
        return field_descriptor;
}

const it::Record & DarcyMH::get_input_type() {

    it::Record ns_rec = Input::Type::Record("NonlinearSolver", "Parameters to a non-linear solver.")
        .declare_key("linear_solver", LinSys::get_input_type(), it::Default::obligatory(),
            "Linear solver for MH problem.")
        .declare_key("tolerance", it::Double(0.0), it::Default("1E-6"),
            "Residual tolerance.")
        .declare_key("max_it", it::Integer(0), it::Default("100"),
            "Maximal number of iterations (linear solves) of the non-linear solver.")
        .declare_key("converge_on_stagnation", it::Bool(), it::Default("false"),
            "If a stagnation of the nonlinear solver is detected the solver stops. "
            "A divergence is reported by default forcing the end of the simulation. Setting this flag to 'true', the solver"
            "ends with convergence success on stagnation, but report warning about it.")
        .close();

    return it::Record("Flow_Darcy_MH", "Mixed-Hybrid  solver for STEADY saturated Darcy flow.")
		.derive_from(DarcyFlowInterface::get_input_type())
        .declare_key("gravity", it::Array(it::Double(), 3,3), it::Default("[ 0, 0, -1]"),
                "Vector of the gravitational acceleration (divided by the acceleration). Dimensionless, magnitude one for the Earth conditions.")
		.declare_key("input_fields", it::Array( type_field_descriptor() ), it::Default::obligatory(),
                "Input data for Darcy flow model.")				
        .declare_key("nonlinear_solver", ns_rec, it::Default::obligatory(),
                "Non-linear solver for MH problem.")
        .declare_key("output_stream", OutputTime::get_input_type(), it::Default::obligatory(),
                "Parameters of output stream.")

        .declare_key("output", DarcyFlowMHOutput::get_input_type(), IT::Default("{ \"fields\": [ \"pressure_p0\", \"velocity_p0\" ] }"),
                "Parameters of output from MH module.")
        .declare_key("output_specific", DarcyFlowMHOutput::get_input_type_specific(), it::Default::optional(),
                "Parameters of output form MH module.")
        .declare_key("balance", Balance::get_input_type(), it::Default("{}"),
                "Settings for computing mass balance.")
        .declare_key("time", TimeGovernor::get_input_type(),
                "Time governor setting for the unsteady Darcy flow model.")
		.declare_key("n_schurs", it::Integer(0,2), it::Default("2"),
				"Number of Schur complements to perform when solving MH system.")
		.declare_key("mortar_method", get_mh_mortar_selection(), it::Default("\"None\""),
				"Method for coupling Darcy flow between dimensions." )
		.close();
}


const int DarcyMH::registrar =
		Input::register_class< DarcyMH, Mesh &, const Input::Record >("Flow_Darcy_MH") +
		DarcyMH::get_input_type().size();



DarcyMH::EqData::EqData()
{
    ADD_FIELD(anisotropy, "Anisotropy of the conductivity tensor.", "1.0" );
    	anisotropy.units( UnitSI::dimensionless() );

    ADD_FIELD(cross_section, "Complement dimension parameter (cross section for 1D, thickness for 2D).", "1.0" );
    	cross_section.units( UnitSI().m(3).md() );

    ADD_FIELD(conductivity, "Isotropic conductivity scalar.", "1.0" );
    	conductivity.units( UnitSI().m().s(-1) );

    ADD_FIELD(sigma, "Transition coefficient between dimensions.", "1.0");
    	sigma.units( UnitSI::dimensionless() );

    ADD_FIELD(water_source_density, "Water source density.", "0.0");
    	water_source_density.units( UnitSI().s(-1) );
    
    ADD_FIELD(bc_type,"Boundary condition type, possible values:", "\"none\"" );
    	// TODO: temporary solution, we should try to get rid of pointer to the selection after having generic types
        bc_type.input_selection( get_bc_type_selection() );
        bc_type.units( UnitSI::dimensionless() );

    ADD_FIELD(bc_pressure,"Prescribed pressure value on the boundary. Used for all values of 'bc_type' save the bc_type='none'."
		"See documentation of 'bc_type' for exact meaning of 'bc_pressure' in individual boundary condition types.", "0.0");
    	bc_pressure.disable_where(bc_type, {none} );
        bc_pressure.units( UnitSI().m() );

    ADD_FIELD(bc_flux,"Incoming water boundary flux. Used for bc_types : 'none', 'total_flux', 'seepage', 'river'.", "0.0");
    	bc_flux.disable_where(bc_type, {none, dirichlet} );
        bc_flux.units( UnitSI().m(4).s(-1).md() );

    ADD_FIELD(bc_robin_sigma,"Conductivity coefficient in the 'total_flux' or the 'river' boundary condition type.", "0.0");
    	bc_robin_sigma.disable_where(bc_type, {none, dirichlet, seepage} );
        bc_robin_sigma.units( UnitSI().m(3).s(-1).md() );

    ADD_FIELD(bc_switch_pressure,
            "Critical switch pressure for 'seepage' and 'river' boundary conditions.", "0.0");
    bc_switch_pressure.disable_where(bc_type, {none, dirichlet, total_flux} );
    bc_switch_pressure.units( UnitSI().m() );

    //these are for unsteady
    ADD_FIELD(init_pressure, "Initial condition as pressure", "0.0" );
    	init_pressure.units( UnitSI().m() );

    ADD_FIELD(storativity,"Storativity.", "0.0" );
    	storativity.units( UnitSI().m(-1) );

    //time_term_fields = this->subset({"storativity"});
    //main_matrix_fields = this->subset({"anisotropy", "conductivity", "cross_section", "sigma", "bc_type", "bc_robin_sigma"});
    //rhs_fields = this->subset({"water_source_density", "bc_pressure", "bc_flux"});

}






//=============================================================================
// CREATE AND FILL GLOBAL MH MATRIX OF THE WATER MODEL
// - do it in parallel:
//   - initial distribution of elements, edges
//
/*! @brief CREATE AND FILL GLOBAL MH MATRIX OF THE WATER MODEL
 *
 * Parameters {Solver,NSchurs} number of performed Schur
 * complements (0,1,2) for water flow MH-system
 *
 */
//=============================================================================
DarcyMH::DarcyMH(Mesh &mesh_in, const Input::Record in_rec)
: DarcyFlowInterface(mesh_in, in_rec),
    solution(nullptr),
    schur0(nullptr)
{

    is_linear_=true;

    START_TIMER("Darcy constructor");
    {
        Input::Record time_record;
        if ( in_rec.opt_val("time", time_record) )
            time_ = new TimeGovernor(time_record);
        else
            time_ = new TimeGovernor();
    }

    data_ = make_shared<EqData>();
    EquationBase::eq_data_ = data_.get();
    

    size = mesh_->n_elements() + mesh_->n_sides() + mesh_->n_edges();
    n_schur_compls = in_rec.val<int>("n_schurs");
    mortar_method_= in_rec.val<MortarMethod>("mortar_method");
    if (mortar_method_ != NoMortar) {
        mesh_->make_intersec_elements();
    }
    


    //side_ds->view( std::cout );
    //el_ds->view( std::cout );
    //edge_ds->view( std::cout );
    //rows_ds->view( std::cout );
    
}



void DarcyMH::init_eq_data()
//connecting data fields with mesh
{

    START_TIMER("data init");
    data_->mesh = mesh_;
    data_->mh_dh = &mh_dh;
    data_->set_mesh(*mesh_);

    auto gravity_array = input_record_.val<Input::Array>("gravity");
    std::vector<double> gvec;
    gravity_array.copy_to(gvec);
    gvec.push_back(0.0); // zero pressure shift
    data_->gravity_ =  arma::vec(gvec);
    data_->gravity_vec_ = data_->gravity_.subvec(0,2);

    data_->bc_pressure.add_factory(
        std::make_shared<FieldAddPotential<3, FieldValue<3>::Scalar>::FieldFactory>
        (data_->gravity_, "bc_piezo_head") );
    data_->bc_switch_pressure.add_factory(
            std::make_shared<FieldAddPotential<3, FieldValue<3>::Scalar>::FieldFactory>
            (data_->gravity_, "bc_switch_piezo_head") );


    data_->set_input_list( this->input_record_.val<Input::Array>("input_fields") );
    data_->mark_input_times(*time_);
}



void DarcyMH::initialize() {

    init_eq_data();
    output_object = new DarcyFlowMHOutput(this, input_record_);

    mh_dh.reinit(mesh_);
    // Initialize bc_switch_dirichlet to size of global boundary.
    bc_switch_dirichlet.resize(mesh_->bc_elements.size(), 1);


    nonlinear_iteration_=0;
    Input::AbstractRecord rec = this->input_record_
            .val<Input::Record>("nonlinear_solver")
            .val<Input::AbstractRecord>("linear_solver");

    // auxiliary set_time call  since allocation assembly evaluates fields as well
    data_->set_time(time_->step(), LimitSide::right);
    data_->system_.local_matrix = std::make_shared<arma::mat>();
    create_linear_system(rec);



    // allocate time term vectors
    VecDuplicate(schur0->get_solution(), &previous_solution);
    VecCreateMPI(PETSC_COMM_WORLD, mh_dh.rows_ds->lsize(),PETSC_DETERMINE,&(steady_diagonal));
    VecDuplicate(steady_diagonal,& new_diagonal);
    VecZeroEntries(new_diagonal);
    VecDuplicate(steady_diagonal, &steady_rhs);


    // initialization of balance object
    balance_ = std::make_shared<Balance>("water", mesh_);
    balance_->init_from_input(input_record_.val<Input::Record>("balance"), time());
    if (balance_)
    {
        data_-> water_balance_idx_ = water_balance_idx_ = balance_->add_quantity("water_volume");
        balance_->allocate(mh_dh.rows_ds->lsize(), 1);
        balance_->units(UnitSI().m(3));
    }


    data_->system_.balance = balance_;
    data_->system_.lin_sys = schur0;


    initialize_specific();
}

void DarcyMH::initialize_specific()
{}

void DarcyMH::zero_time_step()
{

    /* TODO:
     * - Allow solution reconstruction (pressure and velocity) from initial condition on user request.
     * - Steady solution as an intitial condition may be forced by setting inti_time =-1, and set data for the steady solver in that time.
     *   Solver should be able to switch from and to steady case depending on the zero time term.
     */

    data_->set_time(time_->step(), LimitSide::right);
    // zero_time_term means steady case
    bool zero_time_term_from_right = zero_time_term();


    if (zero_time_term_from_right) {
        // steady case
        VecZeroEntries(schur0->get_solution());
        //read_initial_condition(); // Possible solution guess for steady case.
        use_steady_assembly_ = true;
        solve_nonlinear(); // with right limit data
    } else {
        VecZeroEntries(schur0->get_solution());
        VecZeroEntries(previous_solution);
        read_initial_condition();
        assembly_linear_system(); // in particular due to balance
        // TODO: reconstruction of solution in zero time.
    }
    //solution_output(T,right_limit); // data for time T in any case
    output_data();
}

//=============================================================================
// COMPOSE and SOLVE WATER MH System possibly through Schur complements
//=============================================================================
void DarcyMH::update_solution()
{
    START_TIMER("Solving MH system");

    time_->next_time();

    time_->view("DARCY"); //time governor information output
    data_->set_time(time_->step(), LimitSide::left);
    bool zero_time_term_from_left=zero_time_term();

    bool jump_time = data_->storativity.is_jump_time();
    if (! zero_time_term_from_left) {
        // time term not treated as zero
        // Unsteady solution up to the T.

        // this flag is necesssary for switching BC to avoid setting zero neumann on the whole boundary in the steady case
        use_steady_assembly_ = false;
        prepare_new_time_step(); //SWAP

        solve_nonlinear(); // with left limit data
        if (jump_time) {
        	WarningOut() << "Output of solution discontinuous in time not supported yet.\n";
            //solution_output(T, left_limit); // output use time T- delta*dt
            //output_data();
        }
    }

    if (time_->is_end()) {
        // output for unsteady case, end_time should not be the jump time
        // but rether check that
        if (! zero_time_term_from_left && ! jump_time) output_data();
        return;
    }

    data_->set_time(time_->step(), LimitSide::right);
    bool zero_time_term_from_right=zero_time_term();
    if (zero_time_term_from_right) {
        // this flag is necesssary for switching BC to avoid setting zero neumann on the whole boundary in the steady case
        use_steady_assembly_ = true;
        solve_nonlinear(); // with right limit data

    } else if (! zero_time_term_from_left && jump_time) {
    	WarningOut() << "Discontinuous time term not supported yet.\n";
        //solution_transfer(); // internally call set_time(T, left) and set_time(T,right) again
        //solve_nonlinear(); // with right limit data
    }
    //solution_output(T,right_limit); // data for time T in any case
    output_data();

}

bool DarcyMH::zero_time_term() {
    return data_->storativity.field_result(mesh_->region_db().get_region_set("BULK")) == result_zeros;
}


void DarcyMH::solve_nonlinear()
{

	DebugOut().fmt("dt: {}\n", time_->step().length());
    assembly_linear_system();
    double residual_norm = schur0->compute_residual();
    unsigned int l_it=0;
    nonlinear_iteration_ = 0;
    MessageOut().fmt("[nonlin solver] norm of initial residual: {}\n", residual_norm);

    // Reduce is_linear flag.
    int is_linear_common;
    MPI_Allreduce(&is_linear_, &is_linear_common,1, MPI_INT ,MPI_MIN,PETSC_COMM_WORLD);

    Input::Record nl_solver_rec = input_record_.val<Input::Record>("nonlinear_solver");
    this->tolerance_ = nl_solver_rec.val<double>("tolerance");
    this->max_n_it_  = nl_solver_rec.val<unsigned int>("max_it");

    if (! is_linear_common) {
        // set tolerances of the linear solver unless they are set by user.
        schur0->set_tolerances(0.1, 0.1*this->tolerance_, 10);
    }
    vector<double> convergence_history;

    Vec save_solution;
    VecDuplicate(schur0->get_solution(), &save_solution);
    while (residual_norm > this->tolerance_ &&  nonlinear_iteration_ < this->max_n_it_) {
    	OLD_ASSERT_EQUAL( convergence_history.size(), nonlinear_iteration_ );
        convergence_history.push_back(residual_norm);

        // stagnation test
        if (convergence_history.size() >= 5 &&
            convergence_history[ convergence_history.size() - 1]/convergence_history[ convergence_history.size() - 2] > 0.9 &&
            convergence_history[ convergence_history.size() - 1]/convergence_history[ convergence_history.size() - 5] > 0.8) {
            // stagnation
            if (input_record_.val<Input::Record>("nonlinear_solver").val<bool>("converge_on_stagnation")) {
            	WarningOut().fmt("Accept solution on stagnation. Its: {} Residual: {}\n", nonlinear_iteration_, residual_norm);
                break;
            } else {
                THROW(ExcSolverDiverge() << EI_Reason("Stagnation."));
            }
        }

        if (! is_linear_common)
            VecCopy( schur0->get_solution(), save_solution);
        int convergedReason = schur0->solve();
        nonlinear_iteration_++;

        // hack to make BDDC work with empty compute_residual
        if (is_linear_common) break;

        double alpha = 1; // how much of new solution
        VecAXPBY(schur0->get_solution(), (1-alpha), alpha, save_solution);



        //LogOut().fmt("Linear solver ended with reason: {} \n", convergedReason );
        //OLD_ASSERT( convergedReason >= 0, "Linear solver failed to converge. Convergence reason %d \n", convergedReason );
        assembly_linear_system();

        residual_norm = schur0->compute_residual();
        MessageOut().fmt("[nonlinear solver] it: {} lin. it:{} (reason: {}) residual: {}\n",
        		nonlinear_iteration_, l_it, convergedReason, residual_norm);
    }
    this -> postprocess();

    // adapt timestep
    if (! this->zero_time_term()) {
        double mult = 1.0;
        if (nonlinear_iteration_ < 3) mult = 1.6;
        if (nonlinear_iteration_ > 7) mult = 0.7;
        int result = time_->set_upper_constraint(time_->dt() * mult, "Darcy adaptivity.");
        //DebugOut().fmt("time adaptivity, res: {} it: {} m: {} dt: {} edt: {}\n", result, nonlinear_iteration_, mult, time_->dt(), time_->estimate_dt());
    }

    solution_changed_for_scatter=true;

}


void DarcyMH::prepare_new_time_step()
{
    //VecSwap(previous_solution, schur0->get_solution());
}

void DarcyMH::postprocess() 
{
    START_TIMER("postprocess");
    //ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);

    // modify side fluxes in parallel
    // for every local edge take time term on digonal and add it to the corresponding flux
    /*
    for (unsigned int i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {
        ele = mesh_->element(el_4_loc[i_loc]);
        FOR_ELEMENT_SIDES(ele,i) {
            side_rows[i] = side_row_4_id[ mh_dh.side_dof( ele_ac.side(i) ) ];
            values[i] = -1.0 * ele_ac.measure() *
              data.cross_section.value(ele_ac.centre(), ele_ac.element_accessor()) *
              data.water_source_density.value(ele_ac.centre(), ele_ac.element_accessor()) /
              ele_ac.n_sides();
        }
        VecSetValues(schur0->get_solution(), ele_ac.n_sides(), side_rows, values, ADD_VALUES);
    }
    VecAssemblyBegin(schur0->get_solution());
    VecAssemblyEnd(schur0->get_solution());
    */
}


void DarcyMH::output_data() {
    START_TIMER("Darcy output data");
    //time_->view("DARCY"); //time governor information output
	this->output_object->output();


	if (balance_ != nullptr)
	{
	    START_TIMER("Darcy balance output");
		if (balance_->cumulative() && time_->tlevel() > 0)
		{
			balance_->calculate_cumulative_sources(water_balance_idx_, schur0->get_solution(), time_->dt());
			balance_->calculate_cumulative_fluxes(water_balance_idx_, schur0->get_solution(), time_->dt());
		}

		if ( balance_->is_current( time().step()) )
		{
			balance_->calculate_mass(water_balance_idx_, schur0->get_solution());
			balance_->calculate_source(water_balance_idx_, schur0->get_solution());
			balance_->calculate_flux(water_balance_idx_, schur0->get_solution());
			balance_->output(time_->t());
		}
	}
}


double DarcyMH::solution_precision() const
{
	return schur0->get_solution_precision();
}



void  DarcyMH::get_solution_vector(double * &vec, unsigned int &vec_size)
{
    // TODO: make class for vectors (wrapper for PETSC or other) derived from LazyDependency
    // and use its mechanism to manage dependency between vectors
    if (solution_changed_for_scatter) {

        // scatter solution to all procs
        VecScatterBegin(par_to_all, schur0->get_solution(), sol_vec, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(  par_to_all, schur0->get_solution(), sol_vec, INSERT_VALUES, SCATTER_FORWARD);
        solution_changed_for_scatter=false;
    }

    vec = solution;
    vec_size = this->size;
    OLD_ASSERT(vec != NULL, "Requested solution is not allocated!\n");
}

void  DarcyMH::get_parallel_solution_vector(Vec &vec)
{
    vec=schur0->get_solution();
    OLD_ASSERT(vec != NULL, "Requested solution is not allocated!\n");
}


/*
void DarcyMH::local_assembly_specific() {

}
*/
// ===========================================================================================
//
//   MATRIX ASSEMBLY - we use abstract assembly routine, where  LS Mat/Vec SetValues
//   are in fact pointers to allocating or filling functions - this is governed by Linsystem roitunes
//
// =======================================================================================

void DarcyMH::assembly_mh_matrix(MultidimAssembler assembler)
{
    START_TIMER("DarcyFlowMH_Steady::assembly_steady_mh_matrix");

    LinSys *ls = schur0;

    class Boundary *bcd;
    class Neighbour *ngh;

    bool fill_matrix = assembler.size() > 0;
    int side_row, edge_row, loc_b = 0;
    int tmp_rows[100];
    double local_vb[4]; // 2x2 matrix
    int side_rows[4], edge_rows[4];

    // to make space for second schur complement, max. 10 neighbour edges of one el.
    double zeros[1000];
    for(int i=0; i<1000; i++) zeros[i]=0.0;

    double minus_ones[4] = { -1.0, -1.0, -1.0, -1.0 };
    double * loc_side_rhs = data_->system_.loc_side_rhs;

    arma::mat &local_matrix = *(data_->system_.local_matrix);


    if (balance_ != nullptr && fill_matrix)
        balance_->start_flux_assembly(water_balance_idx_);

    for (unsigned int i_loc = 0; i_loc < mh_dh.el_ds->lsize(); i_loc++) {
        auto ele_ac = mh_dh.accessor(i_loc);
        unsigned int nsides = ele_ac.n_sides();
        data_->system_.dirichlet_edge.resize(nsides);

        

        for (unsigned int i = 0; i < nsides; i++) {

/*            if (! side_ds->is_local(idx_side)) {
                cout << el_ds->myp() << " : iside: " << ele.index() << " [" << el_ds->begin() << ", " << el_ds->end() << "]" << endl;
                cout << el_ds->myp() << " : iside: " << idx_side << " [" << side_ds->begin() << ", " << side_ds->end() << "]" << endl;

            }*/

            side_rows[i] = side_row = ele_ac.side_row(i);
            edge_rows[i] = edge_row = ele_ac.edge_row(i);
            bcd=ele_ac.side(i)->cond();

            // gravity term on RHS
            //
            loc_side_rhs[i] = 0;

            // set block C and C': side-edge, edge-side
            double c_val = 1.0;
            data_->system_.dirichlet_edge[i] = 0;

            if (bcd) {
                ElementAccessor<3> b_ele = bcd->element_accessor();
                EqData::BC_Type type = (EqData::BC_Type)data_->bc_type.value(b_ele.centre(), b_ele);

                double cross_section = data_->cross_section.value(ele_ac.centre(), ele_ac.element_accessor());

                if ( type == EqData::none) {
                    // homogeneous neumann
                } else if ( type == EqData::dirichlet ) {
                    double bc_pressure = data_->bc_pressure.value(b_ele.centre(), b_ele);
                    c_val = 0.0;
                    loc_side_rhs[i] -= bc_pressure;
                    ls->rhs_set_value(edge_row, -bc_pressure);
                    ls->mat_set_value(edge_row, edge_row, -1.0);
                    data_->system_.dirichlet_edge[i] = 1;

                } else if ( type == EqData::total_flux) {
                    // internally we work with outward flux
                    double bc_flux = -data_->bc_flux.value(b_ele.centre(), b_ele);
                    double bc_pressure = data_->bc_pressure.value(b_ele.centre(), b_ele);
                            double bc_sigma = data_->bc_robin_sigma.value(b_ele.centre(), b_ele);
                            ls->mat_set_value(edge_row, edge_row, -bcd->element()->measure() * bc_sigma * cross_section );
                            ls->rhs_set_value(edge_row, (bc_flux - bc_sigma * bc_pressure) * bcd->element()->measure() * cross_section);

                } else if (type==EqData::seepage) {
                    is_linear_=false;
                    //unsigned int loc_edge_idx = edge_row - rows_ds->begin() - side_ds->lsize() - el_ds->lsize();
                    unsigned int loc_edge_idx = bcd->bc_ele_idx_;
                    char & switch_dirichlet = bc_switch_dirichlet[loc_edge_idx];
                    double bc_pressure = data_->bc_switch_pressure.value(b_ele.centre(), b_ele);
                    double bc_flux = -data_->bc_flux.value(b_ele.centre(), b_ele);
                    double side_flux=bc_flux * bcd->element()->measure() * cross_section;

                    // ** Update BC type. **
                    if (switch_dirichlet) {
                        // check and possibly switch to flux BC
                        // The switch raise error on the corresponding edge row.
                        // Magnitude of the error is abs(solution_flux - side_flux).
                        ASSERT_DBG(mh_dh.rows_ds->is_local(side_row))(side_row);
                        unsigned int loc_side_row = ele_ac.side_local_row(i);
                        double & solution_flux = ls->get_solution_array()[loc_side_row];

                        if ( solution_flux < side_flux) {
                            //DebugOut().fmt("x: {}, to neum, p: {} f: {} -> f: {}\n", b_ele.centre()[0], bc_pressure, solution_flux, side_flux);
                            solution_flux = side_flux;
                            switch_dirichlet=0;

                        }
                    } else {
                        // check and possibly switch to  pressure BC
                        // TODO: What is the appropriate DOF in not local?
                        // The switch raise error on the corresponding side row.
                        // Magnitude of the error is abs(solution_head - bc_pressure)
                        // Since usually K is very large, this error would be much
                        // higher then error caused by the inverse switch, this
                        // cause that a solution  with the flux violating the
                        // flux inequality leading may be accepted, while the error
                        // in pressure inequality is always satisfied.
                        ASSERT_DBG(mh_dh.rows_ds->is_local(edge_row))(edge_row);
                        unsigned int loc_edge_row = ele_ac.edge_local_row(i);
                        double & solution_head = ls->get_solution_array()[loc_edge_row];

                        if ( solution_head > bc_pressure) {
                            //DebugOut().fmt("x: {}, to dirich, p: {} -> p: {} f: {}\n",b_ele.centre()[0], solution_head, bc_pressure, bc_flux);
                            solution_head = bc_pressure;
                            switch_dirichlet=1;
                        }
                    }

                    // ** Apply BCUpdate BC type. **
                    // Force Dirichlet type during the first iteration of the unsteady case.
                    if (switch_dirichlet || (use_steady_assembly_ && nonlinear_iteration_ == 0) ) {
                        //DebugOut().fmt("x: {}, dirich, bcp: {}\n", b_ele.centre()[0], bc_pressure);
                        c_val = 0.0;
                        loc_side_rhs[i] -= bc_pressure;
                        ls->rhs_set_value(edge_row, -bc_pressure);
                        ls->mat_set_value(edge_row, edge_row, -1.0);
                        data_->system_.dirichlet_edge[i] = 1;
                    } else {
                        //DebugOut()("x: {}, neuman, q: {}  bcq: {}\n", b_ele.centre()[0], side_flux, bc_flux);
                        ls->rhs_set_value(edge_row, side_flux);
                    }

                } else if (type==EqData::river) {
                    is_linear_=false;
                    //unsigned int loc_edge_idx = edge_row - rows_ds->begin() - side_ds->lsize() - el_ds->lsize();
                    //unsigned int loc_edge_idx = bcd->bc_ele_idx_;
                    //char & switch_dirichlet = bc_switch_dirichlet[loc_edge_idx];

                    double bc_pressure = data_->bc_pressure.value(b_ele.centre(), b_ele);
                    double bc_switch_pressure = data_->bc_switch_pressure.value(b_ele.centre(), b_ele);
                    double bc_flux = -data_->bc_flux.value(b_ele.centre(), b_ele);
                    double bc_sigma = data_->bc_robin_sigma.value(b_ele.centre(), b_ele);
                    ASSERT_DBG(mh_dh.rows_ds->is_local(edge_row))(edge_row);
                    unsigned int loc_edge_row = ele_ac.edge_local_row(i);
                    double & solution_head = ls->get_solution_array()[loc_edge_row];


                    // Force Robin type during the first iteration of the unsteady case.
                    if (solution_head > bc_switch_pressure  || (use_steady_assembly_ && nonlinear_iteration_ ==0)) {
                        // Robin BC
                        //DebugOut().fmt("x: {}, robin, bcp: {}\n", b_ele.centre()[0], bc_pressure);
                        ls->rhs_set_value(edge_row, bcd->element()->measure() * cross_section * (bc_flux - bc_sigma * bc_pressure)  );
                        ls->mat_set_value(edge_row, edge_row, -bcd->element()->measure() * bc_sigma * cross_section );
                    } else {
                        // Neumann BC
                        //DebugOut().fmt("x: {}, neuman, q: {}  bcq: {}\n", b_ele.centre()[0], bc_switch_pressure, bc_pressure);
                        double bc_total_flux = bc_flux + bc_sigma*(bc_switch_pressure - bc_pressure);
                        ls->rhs_set_value(edge_row, bc_total_flux * bcd->element()->measure() * cross_section);
                    }
                } else {
                    xprintf(UsrErr, "BC type not supported.\n");
                }

                if (balance_ != nullptr && fill_matrix)
                {
                   /*
                    DebugOut()("add_flux: {} {} {} {}\n",
                            mh_dh.el_ds->myp(),
                            ele_ac.ele_global_idx(),
                            loc_b,
                            side_row);*/
                    balance_->add_flux_matrix_values(water_balance_idx_, loc_b, {side_row}, {1});
                }
                ++loc_b;
            }
            ls->mat_set_value(side_row, edge_row, c_val);
            ls->mat_set_value(edge_row, side_row, c_val);

        }

        if (fill_matrix) {
            assembler[ele_ac.dim()-1]->assembly_local_matrix(ele_ac);

            // assemble matrix for weights in BDDCML
            // approximation to diagonal of 
            // S = -C - B*inv(A)*B'
            // as 
            // diag(S) ~ - diag(C) - 1./diag(A)
            // the weights form a partition of unity to average a discontinuous solution from neighbouring subdomains
            // to a continuous one
            // it is important to scale the effect - if conductivity is low for one subdomain and high for the other,
            // trust more the one with low conductivity - it will be closer to the truth than an arithmetic average
            if ( typeid(*ls) == typeid(LinSys_BDDC) ) {
               for(unsigned int i=0; i < nsides; i++) {
                   double val_side =  local_matrix(i,i);
                   double val_edge =  -1./local_matrix(i,i);

                   static_cast<LinSys_BDDC*>(ls)->diagonal_weights_set_value( side_rows[i], val_side );
                   static_cast<LinSys_BDDC*>(ls)->diagonal_weights_set_value( edge_rows[i], val_edge );
               }
            }
        }

        ls->rhs_set_values(nsides, side_rows, loc_side_rhs);

        
        // set block A: side-side on one element - block diagonal matrix
        ls->mat_set_values(nsides, side_rows, nsides, side_rows, local_matrix.memptr());
        // set block B, B': element-side, side-element
        int ele_row = ele_ac.ele_row();
        ls->mat_set_values(1, &ele_row, nsides, side_rows, minus_ones);
        ls->mat_set_values(nsides, side_rows, 1, &ele_row, minus_ones);


        // D block:  diagonal: element-element

        ls->mat_set_value(ele_row, ele_row, 0.0);         // maybe this should be in virtual block for schur preallocation

        if ( typeid(*ls) == typeid(LinSys_BDDC) ) {
           double val_ele =  1.;
           static_cast<LinSys_BDDC*>(ls)->diagonal_weights_set_value( ele_row, val_ele );
        }

        // D, E',E block: compatible connections: element-edge
        
        for (unsigned int i = 0; i < ele_ac.full_iter()->n_neighs_vb; i++) {
            // every compatible connection adds a 2x2 matrix involving
            // current element pressure  and a connected edge pressure
            ngh= ele_ac.full_iter()->neigh_vb[i];
            tmp_rows[0]=ele_row;
            tmp_rows[1]=mh_dh.row_4_edge[ ngh->edge_idx() ];
            
            if (fill_matrix)
                assembler[ngh->side()->dim()]->assembly_local_vb(local_vb, ele_ac.full_iter(), ngh);
            
            ls->mat_set_values(2, tmp_rows, 2, tmp_rows, local_vb);

            // update matrix for weights in BDDCML
            if ( typeid(*ls) == typeid(LinSys_BDDC) ) {
               int ind = tmp_rows[1];
               // there is -value on diagonal in block C!
               double new_val = local_vb[0];
               static_cast<LinSys_BDDC*>(ls)->diagonal_weights_set_value( ind, new_val );
            }

            if (n_schur_compls == 2) {
                // for 2. Schur: N dim edge is conected with N dim element =>
                // there are nz between N dim edge and N-1 dim edges of the element

                ls->mat_set_values(nsides, edge_rows, 1, tmp_rows+1, zeros);
                ls->mat_set_values(1, tmp_rows+1, nsides, edge_rows, zeros);

                // save all global edge indices to higher positions
                tmp_rows[2+i] = tmp_rows[1];
            }
        }


        // add virtual values for schur complement allocation
        uint n_neigh;
        switch (n_schur_compls) {
        case 2:
            n_neigh = ele_ac.full_iter()->n_neighs_vb;
            // Connections between edges of N+1 dim. elements neighboring with actual N dim element 'ele'
        	OLD_ASSERT(n_neigh*n_neigh<1000, "Too many values in E block.");
            ls->mat_set_values(ele_ac.full_iter()->n_neighs_vb, tmp_rows+2,
                    ele_ac.full_iter()->n_neighs_vb, tmp_rows+2, zeros);

        case 1: // included also for case 2
            // -(C')*(A-)*B block and its transpose conect edge with its elements
            ls->mat_set_values(1, &ele_row, ele_ac.n_sides(), edge_rows, zeros);
            ls->mat_set_values(ele_ac.n_sides(), edge_rows, 1, &ele_row, zeros);
            // -(C')*(A-)*C block conect all edges of every element
            ls->mat_set_values(ele_ac.n_sides(), edge_rows, ele_ac.n_sides(), edge_rows, zeros);
        }
    }    
    
    if (balance_ != nullptr && fill_matrix)
        balance_->finish_flux_assembly(water_balance_idx_);



    if (mortar_method_ == MortarP0) {
        P0_CouplingAssembler(*this).assembly(*ls);
    } else if (mortar_method_ == MortarP1) {
        P1_CouplingAssembler(*this).assembly(*ls);
    }  
}


void DarcyMH::assembly_source_term()
{
    START_TIMER("assembly source term");
    if (balance_ != nullptr)
    	balance_->start_source_assembly(water_balance_idx_);

    for (unsigned int i_loc = 0; i_loc < mh_dh.el_ds->lsize(); i_loc++) {
        auto ele_ac = mh_dh.accessor(i_loc);

        double cs = data_->cross_section.value(ele_ac.centre(), ele_ac.element_accessor());

        // set sources
        double source = ele_ac.measure() * cs *
                data_->water_source_density.value(ele_ac.centre(), ele_ac.element_accessor());
        schur0->rhs_set_value(ele_ac.ele_row(), -1.0 * source );

        if (balance_ != nullptr)
        	balance_->add_source_vec_values(water_balance_idx_, ele_ac.region().bulk_idx(), {(int) ele_ac.ele_row()}, {source});
    }

    if (balance_ != nullptr)
    	balance_->finish_source_assembly(water_balance_idx_);
}


void P0_CouplingAssembler::pressure_diff(int i_ele,
		vector<int> &dofs, unsigned int &ele_type, double &delta, arma::vec &dirichlet) {

	const Element *ele;

	if (i_ele == (int)(ml_it_->size()) ) { // master element .. 1D
		ele_type = 0;
		delta = -delta_0;
		ele=master_;
	} else {
		ele_type = 1;
		const Intersection &isect=intersections_[ (*ml_it_)[i_ele] ];
		delta = isect.intersection_true_size();
		ele = isect.slave_iter();
	}

	dofs.resize(ele->n_sides());
        dirichlet.resize(ele->n_sides());
        dirichlet.zeros();

	for(unsigned int i_side=0; i_side < ele->n_sides(); i_side++ ) {
		dofs[i_side]=darcy_.mh_dh.row_4_edge[ele->side(i_side)->edge_idx()];
		Boundary * bcd = ele->side(i_side)->cond();
		if (bcd) {			
			ElementAccessor<3> b_ele = bcd->element_accessor();
			auto type = (DarcyMH::EqData::BC_Type)darcy_.data_->bc_type.value(b_ele.centre(), b_ele);
			//DebugOut().fmt("bcd id: {} sidx: {} type: {}\n", ele->id(), i_side, type);
			if (type == DarcyMH::EqData::dirichlet) {
				//DebugOut().fmt("Dirichlet: {}\n", ele->index());
				dofs[i_side] = -dofs[i_side];
				double bc_pressure = darcy_.data_->bc_pressure.value(b_ele.centre(), b_ele);
				dirichlet[i_side] = bc_pressure;
			}
		} 
	}

}

/**
 * Works well but there is large error next to the boundary.
 */
 void P0_CouplingAssembler::assembly(LinSys &ls) {
	double delta_i, delta_j;
	arma::mat product;
	arma::vec dirichlet_i, dirichlet_j;
	unsigned int ele_type_i, ele_type_j;

	unsigned int i,j;
	vector<int> dofs_i,dofs_j;

	for(ml_it_ = master_list_.begin(); ml_it_ != master_list_.end(); ++ml_it_) {

    	if (ml_it_->size() == 0) continue; // skip empty masters


		// on the intersection element we consider
		// * intersection dofs for master and slave
		//   those are dofs of the space into which we interpolate
		//   base functions from individual master and slave elements
		//   For the master dofs both are usualy eqivalent.
		// * original dofs - for master same as intersection dofs, for slave
		//   all dofs of slave elements

		// form list of intersection dofs, in this case pressures in barycenters
		// but we do not use those form MH system in order to allow second schur somplement. We rather map these
		// dofs to pressure traces, i.e. we use averages of traces as barycentric values.


		master_ = intersections_[ml_it_->front()].master_iter();
		delta_0 = master_->measure();

		double master_sigma=darcy_.data_->sigma.value( master_->centre(), master_->element_accessor());

		// rows
		double check_delta_sum=0;
		for(i = 0; i <= ml_it_->size(); ++i) {
			pressure_diff(i, dofs_i, ele_type_i, delta_i, dirichlet_i);
			check_delta_sum+=delta_i;
			//columns
			for (j = 0; j <= ml_it_->size(); ++j) {
				pressure_diff(j, dofs_j, ele_type_j, delta_j, dirichlet_j);

				double scale =  -master_sigma * delta_i * delta_j / delta_0;
				product = scale * tensor_average[ele_type_i][ele_type_j];

				arma::vec rhs(dofs_i.size());
				rhs.zeros();
				ls.set_values( dofs_i, dofs_j, product, rhs, dirichlet_i, dirichlet_j);
				auto dofs_i_cp=dofs_i;
				auto dofs_j_cp=dofs_j;
				ls.set_values( dofs_i_cp, dofs_j_cp, product, rhs, dirichlet_i, dirichlet_j);
			}
		}
		OLD_ASSERT(check_delta_sum < 1E-5*delta_0, "sum err %f > 0\n", check_delta_sum/delta_0);
    } // loop over master elements
}



 void P1_CouplingAssembler::add_sides(const Element * ele, unsigned int shift, vector<int> &dofs, vector<double> &dirichlet)
 {

		for(unsigned int i_side=0; i_side < ele->n_sides(); i_side++ ) {
			dofs[shift+i_side] =  darcy_.mh_dh.row_4_edge[ele->side(i_side)->edge_idx()];
			Boundary * bcd = ele->side(i_side)->cond();

			if (bcd) {
				ElementAccessor<3> b_ele = bcd->element_accessor();
				auto type = (DarcyMH::EqData::BC_Type)darcy_.data_->bc_type.value(b_ele.centre(), b_ele);
				//DebugOut().fmt("bcd id: {} sidx: {} type: {}\n", ele->id(), i_side, type);
				if (type == DarcyMH::EqData::dirichlet) {
					//DebugOut().fmt("Dirichlet: {}\n", ele->index());
					dofs[shift + i_side] = -dofs[shift + i_side];
					double bc_pressure = darcy_.data_->bc_pressure.value(b_ele.centre(), b_ele);
					dirichlet[shift + i_side] = bc_pressure;
				}
			}
		}
 }


/**
 * P1 connection of different dimensions
 *
 * - 20.11. 2014 - very poor convergence, big error in pressure even at internal part of the fracture
 */

void P1_CouplingAssembler::assembly(LinSys &ls) {

	for (const Intersection &intersec : intersections_) {
    	const Element * master = intersec.master_iter();
       	const Element * slave = intersec.slave_iter();

	add_sides(slave, 0, dofs, dirichlet);
       	add_sides(master, 3, dofs, dirichlet);
       	
		double master_sigma=darcy_.data_->sigma.value( master->centre(), master->element_accessor());

/*
 * Local coordinates on 1D
 *         t0
 * node 0: 0.0
 * node 1: 1.0
 *
 * base fce points
 * t0 = 0.0    on side 0 node 0
 * t0 = 1.0    on side 1 node 1
 *
 * Local coordinates on 2D
 *         t0  t1
 * node 0: 0.0 0.0
 * node 1: 1.0 0.0
 * node 2: 0.0 1.0
 *
 * base fce points
 * t0=0.5, t1=0.0        on side 0 nodes (0,1)
 * t0=0.5, t1=0.5        on side 1 nodes (1,2)
 * t0=0.0, t1=0.5        on side 2 nodes (2,0)
 */



        arma::vec point_Y(1);
        point_Y.fill(1.0);
        arma::vec point_2D_Y(intersec.map_to_slave(point_Y)); // local coordinates of  Y on slave (1, t0, t1)
        arma::vec point_1D_Y(intersec.map_to_master(point_Y)); //  local coordinates of  Y on master (1, t0)

        arma::vec point_X(1);
        point_X.fill(0.0);
        arma::vec point_2D_X(intersec.map_to_slave(point_X)); // local coordinates of  X on slave (1, t0, t1)
        arma::vec point_1D_X(intersec.map_to_master(point_X)); // local coordinates of  X on master (1, t0)

        arma::mat base_2D(3, 3);
        // basis functions are numbered as sides
        // TODO:
        // Use RT finite element to evaluate following matrices.

        // Ravirat - Thomas base functions evaluated in points (0,0), (1,0), (0,1)
        // 2D RT_i(t0, t1) = a0 + a1*t0 + a2*t1
        //         a0     a1      a2
        base_2D << 1.0 << 0.0 << -2.0 << arma::endr // RT for side 0
                << 1.0 << -2.0 << 0.0 << arma::endr // RT for side 1
                << -1.0 << 2.0 << 2.0 << arma::endr;// RT for side 2
                

        arma::mat base_1D(2, 2);
        // Ravirat - Thomas base functions evaluated in points (0,0), (1,0), (0,1)
        // 1D RT_i(t0) =   a0 + a1 * t0
        //          a0     a1
        base_1D << 1.0 << -1.0 << arma::endr // RT for side 0,
                << 0.0 << 1.0 << arma::endr; // RT for side 1,



        // Consider both 2D and 1D value are defined for the test function
        // related to the each of 5 DOFs involved in the intersection.
        // One of these values is always zero.
        // Compute difference of the 2D and 1D value for every DOF.
        // Compute value of this difference in both endpoints X,Y of the intersection.

        arma::vec difference_in_Y(5);
        arma::vec difference_in_X(5);
        // slave sides 0,1,2
        difference_in_Y.subvec(0, 2) = -base_2D * point_2D_Y;
        difference_in_X.subvec(0, 2) = -base_2D * point_2D_X;
        // master sides 3,4
        difference_in_Y.subvec(3, 4) = base_1D * point_1D_Y;
        difference_in_X.subvec(3, 4) = base_1D * point_1D_X;

        // applying the Simpson's rule
        // to the product of two linear functions f, g we get
        // (b-a)/6 * ( 3*f(a)*g(a) + 3*f(b)*g(b) + 2*f(a)*g(b) + 2*f(b)*g(a) )
        arma::mat A(5, 5);
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                A(i, j) = -master_sigma * intersec.intersection_true_size() *
                        ( difference_in_Y[i] * difference_in_Y[j]
                          + difference_in_Y[i] * difference_in_X[j]/2
                          + difference_in_X[i] * difference_in_Y[j]/2
                          + difference_in_X[i] * difference_in_X[j]
                        ) * (1.0 / 3.0);

            }
        }
        auto dofs_cp=dofs;
        ls.set_values( dofs_cp, dofs_cp, A, rhs, dirichlet, dirichlet);

    }
}



/*******************************************************************************
 * COMPOSE WATER MH MATRIX WITHOUT SCHUR COMPLEMENT
 ******************************************************************************/

void DarcyMH::create_linear_system(Input::AbstractRecord in_rec) {
  
    START_TIMER("preallocation");

    if (schur0 == NULL) { // create Linear System for MH matrix
       
    	if (in_rec.type() == LinSys_BDDC::get_input_type()) {
#ifdef FLOW123D_HAVE_BDDCML
    		WarningOut() << "For BDDC no Schur complements are used.";
            mh_dh.prepare_parallel_bddc();
            n_schur_compls = 0;
            LinSys_BDDC *ls = new LinSys_BDDC(mh_dh.global_row_4_sub_row->size(), &(*mh_dh.rows_ds),
                    3,  // 3 == la::BddcmlWrapper::SPD_VIA_SYMMETRICGENERAL
                    1,  // 1 == number of subdomains per process
                    true); // swap signs of matrix and rhs to make the matrix SPD
            ls->set_from_input(in_rec);
            ls->set_solution( NULL );
            // possible initialization particular to BDDC
            START_TIMER("BDDC set mesh data");
            set_mesh_data_for_bddc(ls);
            schur0=ls;
            END_TIMER("BDDC set mesh data");
#else
            xprintf(Err, "Flow123d was not build with BDDCML support.\n");
#endif // FLOW123D_HAVE_BDDCML
        } 
        else if (in_rec.type() == LinSys_PETSC::get_input_type()) {
        // use PETSC for serial case even when user wants BDDC
            if (n_schur_compls > 2) {
            	WarningOut() << "Invalid number of Schur Complements. Using 2.";
                n_schur_compls = 2;
            }

            LinSys_PETSC *schur1, *schur2;

            if (n_schur_compls == 0) {
                LinSys_PETSC *ls = new LinSys_PETSC( &(*mh_dh.rows_ds) );

                // temporary solution; we have to set precision also for sequantial case of BDDC
                // final solution should be probably call of direct solver for oneproc case
                if (in_rec.type() != LinSys_BDDC::get_input_type()) ls->set_from_input(in_rec);
                else {
                    ls->LinSys::set_from_input(in_rec); // get only common options
                }

                ls->set_solution( NULL );
                schur0=ls;
            } else {
                IS is;
                ISCreateStride(PETSC_COMM_WORLD, mh_dh.side_ds->lsize(), mh_dh.rows_ds->begin(), 1, &is);
                //OLD_ASSERT(err == 0,"Error in ISCreateStride.");

                SchurComplement *ls = new SchurComplement(is, &(*mh_dh.rows_ds));
                ls->set_from_input(in_rec);
                ls->set_solution( NULL );

                // make schur1
                Distribution *ds = ls->make_complement_distribution();
                if (n_schur_compls==1) {
                    schur1 = new LinSys_PETSC(ds);
                    schur1->set_positive_definite();
                } else {
                    IS is;
                    ISCreateStride(PETSC_COMM_WORLD, mh_dh.el_ds->lsize(), ls->get_distribution()->begin(), 1, &is);
                    //OLD_ASSERT(err == 0,"Error in ISCreateStride.");
                    SchurComplement *ls1 = new SchurComplement(is, ds); // is is deallocated by SchurComplement
                    ls1->set_negative_definite();

                    // make schur2
                    schur2 = new LinSys_PETSC( ls1->make_complement_distribution() );
                    schur2->set_positive_definite();
                    schur2->set_from_input(in_rec);
                    ls1->set_complement( schur2 );
                    schur1 = ls1;
                }
                ls->set_complement( schur1 );
                schur0=ls;
            }

            START_TIMER("PETSC PREALLOCATION");
            schur0->set_symmetric();
            schur0->start_allocation();
            assembly_mh_matrix(MultidimAssembler()); // preallocation
    	    VecZeroEntries(schur0->get_solution());
            END_TIMER("PETSC PREALLOCATION");
        }
        else {
            xprintf(Err, "Unknown solver type. Internal error.\n");
        }
    }

    END_TIMER("preallocation");
    make_serial_scatter();
}




void DarcyMH::assembly_linear_system() {
    START_TIMER("DarcyFlowMH_Steady::assembly_linear_system");

    is_linear_=true;
    bool is_steady = zero_time_term();
	//DebugOut() << "Assembly linear system\n";
	if (data_->changed()) {
		//DebugOut()  << "Data changed\n";
		// currently we have no optimization for cases when just time term data or RHS data are changed
	    START_TIMER("full assembly");
        if (typeid(*schur0) != typeid(LinSys_BDDC)) {
            schur0->start_add_assembly(); // finish allocation and create matrix
        }

	    auto multidim_assembler = AssemblyBase::create< AssemblyMH >(data_);

        schur0->mat_zero_entries();
        schur0->rhs_zero_entries();

        assembly_source_term();
	    assembly_mh_matrix( multidim_assembler ); // fill matrix

	    schur0->finish_assembly();
	    schur0->set_matrix_changed();
            //MatView( *const_cast<Mat*>(schur0->get_matrix()), PETSC_VIEWER_STDOUT_WORLD  );
            //VecView( *const_cast<Vec*>(schur0->get_rhs()),   PETSC_VIEWER_STDOUT_WORLD);

	    if (! is_steady) {
	        START_TIMER("fix time term");
	    	//DebugOut() << "setup time term\n";
	    	// assembly time term and rhs
	    	setup_time_term();
	    	modify_system();
	    }
	    else if (balance_ != nullptr)
	    {
	    	balance_->start_mass_assembly(water_balance_idx_);
	    	balance_->finish_mass_assembly(water_balance_idx_);
	    }
	    END_TIMER("full assembly");
	} else {
		START_TIMER("modify system");
		if (! is_steady) {
			modify_system();
		} else {
			//xprintf(PrgErr, "Planned computation time for steady solver, but data are not changed.\n");
		}
		END_TIMER("modiffy system");
	}

}



void DarcyMH::set_mesh_data_for_bddc(LinSys_BDDC * bddc_ls) {
    START_TIMER("DarcyFlowMH_Steady::set_mesh_data_for_bddc");
    // prepare mesh for BDDCML
    // initialize arrays
    // auxiliary map for creating coordinates of local dofs and global-to-local numbering
    std::map<int,arma::vec3> localDofMap;
    // connectivity for the subdomain, i.e. global dof numbers on element, stored element-by-element
    // Indices of Nodes on Elements
    std::vector<int> inet;
    // number of degrees of freedom on elements - determines elementwise chunks of INET array
    // Number of Nodes on Elements
    std::vector<int> nnet;
    // Indices of Subdomain Elements in Global Numbering - for local elements, their global indices
    std::vector<int> isegn;

    // This array is currently not used in BDDCML, it was used as an interface scaling alternative to scaling
    // by diagonal. It corresponds to the rho-scaling.
    std::vector<double> element_permeability;

    // maximal and minimal dimension of elements
    uint elDimMax = 1;
    uint elDimMin = 3;
    for ( unsigned int i_loc = 0; i_loc < mh_dh.el_ds->lsize(); i_loc++ ) {
        auto ele_ac = mh_dh.accessor(i_loc);
        // for each element, create local numbering of dofs as fluxes (sides), pressure (element centre), Lagrange multipliers (edges), compatible connections

        elDimMax = std::max( elDimMax, ele_ac.dim() );
        elDimMin = std::min( elDimMin, ele_ac.dim() );

        isegn.push_back( ele_ac.ele_global_idx() );
        int nne = 0;

        FOR_ELEMENT_SIDES(ele_ac.full_iter(), si) {
            // insert local side dof
            int side_row = ele_ac.side_row(si);
            arma::vec3 coord = ele_ac.side(si)->centre();

            localDofMap.insert( std::make_pair( side_row, coord ) );
            inet.push_back( side_row );
            nne++;
        }

        // insert local pressure dof
        int el_row  = ele_ac.ele_row();
        arma::vec3 coord = ele_ac.centre();
        localDofMap.insert( std::make_pair( el_row, coord ) );
        inet.push_back( el_row );
        nne++;

        FOR_ELEMENT_SIDES(ele_ac.full_iter(), si) {
            // insert local edge dof
            int edge_row = ele_ac.edge_row(si);
            arma::vec3 coord = ele_ac.side(si)->centre();

            localDofMap.insert( std::make_pair( edge_row, coord ) );
            inet.push_back( edge_row );
            nne++;
        }

        // insert dofs related to compatible connections
        for ( unsigned int i_neigh = 0; i_neigh < ele_ac.full_iter()->n_neighs_vb; i_neigh++) {
            int edge_row = mh_dh.row_4_edge[ ele_ac.full_iter()->neigh_vb[i_neigh]->edge_idx()  ];
            arma::vec3 coord = ele_ac.full_iter()->neigh_vb[i_neigh]->edge()->side(0)->centre();

            localDofMap.insert( std::make_pair( edge_row, coord ) );
            inet.push_back( edge_row );
            nne++;
        }

        nnet.push_back( nne );

        // version for rho scaling
        // trace computation
        arma::vec3 centre = ele_ac.centre();
        double conduct = data_->conductivity.value( centre , ele_ac.element_accessor() );
        auto aniso = data_->anisotropy.value( centre, ele_ac.element_accessor() );

        // compute mean on the diagonal
        double coef = 0.;
        for ( int i = 0; i < 3; i++) {
            coef = coef + aniso.at(i,i);
        }
        // Maybe divide by cs
        coef = conduct*coef / 3;

        OLD_ASSERT( coef > 0.,
                "Zero coefficient of hydrodynamic resistance %f . \n ", coef );
        element_permeability.push_back( 1. / coef );
    }
    //convert set of dofs to vectors
    // number of nodes (= dofs) on the subdomain
    int numNodeSub = localDofMap.size();
    OLD_ASSERT_EQUAL( (unsigned int)numNodeSub, mh_dh.global_row_4_sub_row->size() );
    // Indices of Subdomain Nodes in Global Numbering - for local nodes, their global indices
    std::vector<int> isngn( numNodeSub );
    // pseudo-coordinates of local nodes (i.e. dofs)
    // they need not be exact, they are used just for some geometrical considerations in BDDCML, 
    // such as selection of corners maximizing area of a triangle, bounding boxes fro subdomains to 
    // find candidate neighbours etc.
    std::vector<double> xyz( numNodeSub * 3 ) ;
    int ind = 0;
    std::map<int,arma::vec3>::iterator itB = localDofMap.begin();
    for ( ; itB != localDofMap.end(); ++itB ) {
        isngn[ind] = itB -> first;

        arma::vec3 coord = itB -> second;
        for ( int j = 0; j < 3; j++ ) {
            xyz[ j*numNodeSub + ind ] = coord[j];
        }

        ind++;
    }
    localDofMap.clear();

    // Number of Nodal Degrees of Freedom
    // nndf is trivially one - dofs coincide with nodes
    std::vector<int> nndf( numNodeSub, 1 );

    // prepare auxiliary map for renumbering nodes
    typedef std::map<int,int> Global2LocalMap_; //! type for storage of global to local map
    Global2LocalMap_ global2LocalNodeMap;
    for ( unsigned ind = 0; ind < isngn.size(); ++ind ) {
        global2LocalNodeMap.insert( std::make_pair( static_cast<unsigned>( isngn[ind] ), ind ) );
    }

    // renumber nodes in the inet array to locals
    int indInet = 0;
    for ( unsigned int iEle = 0; iEle < isegn.size(); iEle++ ) {
        int nne = nnet[ iEle ];
        for ( int ien = 0; ien < nne; ien++ ) {

            int indGlob = inet[indInet];
            // map it to local node
            Global2LocalMap_::iterator pos = global2LocalNodeMap.find( indGlob );
            OLD_ASSERT( pos != global2LocalNodeMap.end(),
                    "Cannot remap node index %d to local indices. \n ", indGlob );
            int indLoc = static_cast<int> ( pos -> second );

            // store the node
            inet[ indInet++ ] = indLoc;
        }
    }

    int numNodes    = size;
    int numDofsInt  = size;
    int spaceDim    = 3;    // TODO: what is the proper value here?
    int meshDim     = elDimMax;

    bddc_ls -> load_mesh( spaceDim, numNodes, numDofsInt, inet, nnet, nndf, isegn, isngn, isngn, xyz, element_permeability, meshDim );
}




//=============================================================================
// DESTROY WATER MH SYSTEM STRUCTURE
//=============================================================================
DarcyMH::~DarcyMH() {
    if (schur0 != NULL) delete schur0;

	if (solution != NULL) {
	    chkerr(VecDestroy(&sol_vec));
		delete [] solution;
	}

	delete output_object;

	VecScatterDestroy(&par_to_all);
    
}


// ================================================
// PARALLLEL PART
//


void DarcyMH::make_serial_scatter() {
    START_TIMER("prepare scatter");
    // prepare Scatter form parallel to sequantial in original numbering
    {
            IS is_loc;
            int i, *loc_idx; //, si;

            // create local solution vector
            solution = new double[size];
            VecCreateSeqWithArray(PETSC_COMM_SELF,1, size, solution,
                    &(sol_vec));

            // create seq. IS to scatter par solutin to seq. vec. in original order
            // use essentialy row_4_id arrays
            loc_idx = new int [size];
            i = 0;
            FOR_ELEMENTS(mesh_, ele) {
                FOR_ELEMENT_SIDES(ele,si) {
                    loc_idx[i++] = mh_dh.side_row_4_id[ mh_dh.side_dof( ele->side(si) ) ];
                }
            }
            FOR_ELEMENTS(mesh_, ele) {
                loc_idx[i++] = mh_dh.row_4_el[ele.index()];
            }
            for(unsigned int i_edg=0; i_edg < mesh_->n_edges(); i_edg++) {
                loc_idx[i++] = mh_dh.row_4_edge[i_edg];
            }
            OLD_ASSERT( i==size,"Size of array does not match number of fills.\n");
            //DBGPRINT_INT("loc_idx",size,loc_idx);
            ISCreateGeneral(PETSC_COMM_SELF, size, loc_idx, PETSC_COPY_VALUES, &(is_loc));
            delete [] loc_idx;
            VecScatterCreate(schur0->get_solution(), is_loc, sol_vec,
                    PETSC_NULL, &par_to_all);
            ISDestroy(&(is_loc));
    }
    solution_changed_for_scatter=true;

    END_TIMER("prepare scatter");

}


/*
void mat_count_off_proc_values(Mat m, Vec v) {
    int n, first, last;
    const PetscInt *cols;
    Distribution distr(v);

    int n_off = 0;
    int n_on = 0;
    int n_off_rows = 0;
    MatGetOwnershipRange(m, &first, &last);
    for (int row = first; row < last; row++) {
        MatGetRow(m, row, &n, &cols, PETSC_NULL);
        bool exists_off = false;
        for (int i = 0; i < n; i++)
            if (distr.get_proc(cols[i]) != distr.myp())
                n_off++, exists_off = true;
            else
                n_on++;
        if (exists_off)
            n_off_rows++;
        MatRestoreRow(m, row, &n, &cols, PETSC_NULL);
    }
}
*/












void DarcyMH::read_initial_condition()
{
	double *local_sol = schur0->get_solution_array();

	// cycle over local element rows

	DebugOut().fmt("Setup with dt: {}\n", time_->dt());
	for (unsigned int i_loc_el = 0; i_loc_el < mh_dh.el_ds->lsize(); i_loc_el++) {
		auto ele_ac = mh_dh.accessor(i_loc_el);
		// set initial condition
		local_sol[ele_ac.ele_local_row()] = data_->init_pressure.value(ele_ac.centre(),ele_ac.element_accessor());
	}

	solution_changed_for_scatter=true;

}

void DarcyMH::setup_time_term() {
    // save diagonal of steady matrix
    MatGetDiagonal(*( schur0->get_matrix() ), steady_diagonal);
    // save RHS
    VecCopy(*( schur0->get_rhs()), steady_rhs);


    PetscScalar *local_diagonal;
    VecGetArray(new_diagonal,& local_diagonal);

    DebugOut().fmt("Setup with dt: {}\n", time_->dt());

    if (balance_ != nullptr)
    	balance_->start_mass_assembly(water_balance_idx_);

    //DebugOut().fmt("time_term lsize: {} {}\n", mh_dh.el_ds->myp(), mh_dh.el_ds->lsize());
    for (unsigned int i_loc_el = 0; i_loc_el < mh_dh.el_ds->lsize(); i_loc_el++) {
        auto ele_ac = mh_dh.accessor(i_loc_el);

        // set new diagonal
        double diagonal_coeff = data_->cross_section.value(ele_ac.centre(), ele_ac.element_accessor())
        		* data_->storativity.value(ele_ac.centre(), ele_ac.element_accessor())
				* ele_ac.measure();
        local_diagonal[ele_ac.ele_local_row()]= - diagonal_coeff / time_->dt();

        //DebugOut().fmt("time_term: {} {} {} {} {}\n", mh_dh.el_ds->myp(), ele_ac.ele_global_idx(), i_loc_row, i_loc_el + mh_dh.side_ds->lsize(), diagonal_coeff);
        if (balance_ != nullptr)
        	balance_->add_mass_matrix_values(water_balance_idx_,
        	        ele_ac.region().bulk_idx(), { int(ele_ac.ele_row()) }, {diagonal_coeff});
    }
    VecRestoreArray(new_diagonal,& local_diagonal);
    MatDiagonalSet(*( schur0->get_matrix() ), new_diagonal, ADD_VALUES);

    solution_changed_for_scatter=true;
    schur0->set_matrix_changed();

    if (balance_ != nullptr)
    	balance_->finish_mass_assembly(water_balance_idx_);
}

void DarcyMH::modify_system() {
	START_TIMER("modify system");
	if (time_->is_changed_dt() && time_->step().index()>0) {
        double scale_factor=time_->step(-2).length()/time_->step().length();
        if (scale_factor != 1.0) {
            // if time step has changed and setup_time_term not called
            MatDiagonalSet(*( schur0->get_matrix() ),steady_diagonal, INSERT_VALUES);

            VecScale(new_diagonal, time_->last_dt()/time_->dt());
            MatDiagonalSet(*( schur0->get_matrix() ),new_diagonal, ADD_VALUES);
            schur0->set_matrix_changed();
        }
	}

    // modify RHS - add previous solution
    //VecPointwiseMult(*( schur0->get_rhs()), new_diagonal, previous_solution);
	VecPointwiseMult(*( schur0->get_rhs()), new_diagonal, schur0->get_solution());
    VecAXPY(*( schur0->get_rhs()), 1.0, steady_rhs);
    schur0->set_rhs_changed();

    //VecSwap(previous_solution, schur0->get_solution());
}


//-----------------------------------------------------------------------------
// vim: set cindent:
