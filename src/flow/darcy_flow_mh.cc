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

#include "mesh/side_impl.hh"
#include "mesh/long_idx.hh"
#include "mesh/mesh.h"
#include "mesh/partitioning.hh"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "la/distribution.hh"
#include "la/linsys.hh"
#include "la/linsys_PETSC.hh"
#include "la/linsys_BDDC.hh"
#include "la/schur.hh"
//#include "la/sparse_graph.hh"
#include "la/local_to_global_map.hh"

#include "flow/darcy_flow_mh.hh"

#include "flow/darcy_flow_mh_output.hh"


#include "tools/time_governor.hh"
#include "fields/field_algo_base.hh"
#include "fields/field.hh"
#include "fields/field_values.hh"
#include "fields/field_add_potential.hh"

#include "coupling/balance.hh"

#include "la/vector_mpi.hh"

#include "darcy_flow_assembly.hh"

#include "intersection/mixed_mesh_intersections.hh"
#include "intersection/intersection_local.hh"


FLOW123D_FORCE_LINK_IN_CHILD(darcy_flow_mh);




namespace it = Input::Type;

const it::Selection & DarcyMH::get_mh_mortar_selection() {
	return it::Selection("MH_MortarMethod")
		.add_value(NoMortar, "None", "No Mortar method is applied.")
		.add_value(MortarP0, "P0", "Mortar space: P0 on elements of lower dimension.")
		.add_value(MortarP1, "P1", "Mortar space: P1 on intersections, using non-conforming pressures.")
		.close();
}


const it::Selection & DarcyMH::EqData::get_bc_type_selection() {
	return it::Selection("Flow_Darcy_BC_Type")
        .add_value(none, "none",
            "Homogeneous Neumann boundary condition\n(zero normal flux over the boundary).")
        .add_value(dirichlet, "dirichlet",
            "Dirichlet boundary condition. "
            "Specify the pressure head through the ``bc_pressure`` field "
            "or the piezometric head through the ``bc_piezo_head`` field.")
        .add_value(total_flux, "total_flux", "Flux boundary condition (combines Neumann and Robin type). "
            "Water inflow equal to (($ \\delta_d(q_d^N + \\sigma_d (h_d^R - h_d) )$)). "
            "Specify the water inflow by the ``bc_flux`` field, the transition coefficient by ``bc_robin_sigma`` "
            "and the reference pressure head or piezometric head through ``bc_pressure`` or ``bc_piezo_head`` respectively.")
        .add_value(seepage, "seepage",
            "Seepage face boundary condition. Pressure and inflow bounded from above. Boundary with potential seepage flow "
            "is described by the pair of inequalities: "
            "(($h_d \\le h_d^D$)) and (($ -\\boldsymbol q_d\\cdot\\boldsymbol n \\le \\delta q_d^N$)), where the equality holds in at least one of them. "
            "Caution: setting (($q_d^N$)) strictly negative "
            "may lead to an ill posed problem since a positive outflow is enforced. "
            "Parameters (($h_d^D$)) and (($q_d^N$)) are given by the fields ``bc_switch_pressure`` (or ``bc_switch_piezo_head``) and ``bc_flux`` respectively."
            )
        .add_value(river, "river",
            "River boundary condition. For the water level above the bedrock, (($H_d > H_d^S$)), the Robin boundary condition is used with the inflow given by: "
            "(( $ \\delta_d(q_d^N + \\sigma_d(H_d^D - H_d) )$)). For the water level under the bedrock, constant infiltration is used: "
            "(( $ \\delta_d(q_d^N + \\sigma_d(H_d^D - H_d^S) )$)). Parameters: ``bc_pressure``, ``bc_switch_pressure``, "
            " ``bc_sigma``, ``bc_flux``."
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

    it::Record ns_rec = Input::Type::Record("NonlinearSolver", "Non-linear solver settings.")
        .declare_key("linear_solver", LinSys::get_input_type(), it::Default("{}"),
            "Linear solver for MH problem.")
        .declare_key("tolerance", it::Double(0.0), it::Default("1E-6"),
            "Residual tolerance.")
        .declare_key("min_it", it::Integer(0), it::Default("1"),
            "Minimum number of iterations (linear solutions) to use.\nThis is usefull if the convergence criteria "
            "does not characterize your goal well enough so it converges prematurely, possibly even without a single linear solution."
            "If greater then 'max_it' the value is set to 'max_it'.")
        .declare_key("max_it", it::Integer(0), it::Default("100"),
            "Maximum number of iterations (linear solutions) of the non-linear solver.")
        .declare_key("converge_on_stagnation", it::Bool(), it::Default("false"),
            "If a stagnation of the nonlinear solver is detected the solver stops. "
            "A divergence is reported by default, forcing the end of the simulation. By setting this flag to 'true', the solver "
            "ends with convergence success on stagnation, but it reports warning about it.")
        .close();

    return it::Record("Flow_Darcy_MH", "Mixed-Hybrid  solver for saturated Darcy flow.")
		.derive_from(DarcyFlowInterface::get_input_type())
        .declare_key("gravity", it::Array(it::Double(), 3,3), it::Default("[ 0, 0, -1]"),
                "Vector of the gravity force. Dimensionless.")
		.declare_key("input_fields", it::Array( type_field_descriptor() ), it::Default::obligatory(),
                "Input data for Darcy flow model.")				
        .declare_key("nonlinear_solver", ns_rec, it::Default("{}"),
                "Non-linear solver for MH problem.")
        .declare_key("output_stream", OutputTime::get_input_type(), it::Default("{}"),
                "Output stream settings.\n Specify file format, precision etc.")

        .declare_key("output", DarcyFlowMHOutput::get_input_type(), IT::Default("{ \"fields\": [ \"pressure_p0\", \"velocity_p0\" ] }"),
                "Specification of output fields and output times.")
        .declare_key("output_specific", DarcyFlowMHOutput::get_input_type_specific(), it::Default::optional(),
                "Output settings specific to Darcy flow model.\n"
                "Includes raw output and some experimental functionality.")
        .declare_key("balance", Balance::get_input_type(), it::Default("{}"),
                "Settings for computing mass balance.")
        .declare_key("time", TimeGovernor::get_input_type(), it::Default::optional(),
                "Time governor settings for the unsteady Darcy flow model.")
		.declare_key("n_schurs", it::Integer(0,2), it::Default("2"),
				"Number of Schur complements to perform when solving MH system.")
		.declare_key("mortar_method", get_mh_mortar_selection(), it::Default("\"None\""),
				"Method for coupling Darcy flow between dimensions on incompatible meshes. [Experimental]" )
		.close();
}


const int DarcyMH::registrar =
		Input::register_class< DarcyMH, Mesh &, const Input::Record >("Flow_Darcy_MH") +
		DarcyMH::get_input_type().size();



DarcyMH::EqData::EqData()
{
    mortar_method_=NoMortar;

    *this += anisotropy.name("anisotropy")
            .description("Anisotropy of the conductivity tensor.")
            .input_default("1.0")
            .units( UnitSI::dimensionless() );

    *this += cross_section.name("cross_section")
            .description("Complement dimension parameter (cross section for 1D, thickness for 2D).")
            .input_default("1.0")
            .units( UnitSI().m(3).md() );

    *this += conductivity.name("conductivity")
            .description("Isotropic conductivity scalar.")
            .input_default("1.0")
            .units( UnitSI().m().s(-1) )
            .set_limits(0.0);

    *this += sigma.name("sigma")
            .description("Transition coefficient between dimensions.")
            .input_default("1.0")
            .units( UnitSI::dimensionless() );

    *this += water_source_density.name("water_source_density")
            .description("Water source density.")
            .input_default("0.0")
            .units( UnitSI().s(-1) );

    *this += bc_type.name("bc_type")
            .description("Boundary condition type.")
            .input_selection( get_bc_type_selection() )
            .input_default("\"none\"")
            .units( UnitSI::dimensionless() );
    
    *this += bc_pressure
            .disable_where(bc_type, {none, seepage} )
            .name("bc_pressure")
            .description("Prescribed pressure value on the boundary. Used for all values of ``bc_type`` except ``none`` and ``seepage``. "
                "See documentation of ``bc_type`` for exact meaning of ``bc_pressure`` in individual boundary condition types.")
            .input_default("0.0")
            .units( UnitSI().m() );

    *this += bc_flux
            .disable_where(bc_type, {none, dirichlet} )
            .name("bc_flux")
            .description("Incoming water boundary flux. Used for bc_types : ``total_flux``, ``seepage``, ``river``.")
            .input_default("0.0")
            .units( UnitSI().m().s(-1) );

    *this += bc_robin_sigma
            .disable_where(bc_type, {none, dirichlet, seepage} )
            .name("bc_robin_sigma")
            .description("Conductivity coefficient in the ``total_flux`` or the ``river`` boundary condition type.")
            .input_default("0.0")
            .units( UnitSI().s(-1) );

    *this += bc_switch_pressure
            .disable_where(bc_type, {none, dirichlet, total_flux} )
            .name("bc_switch_pressure")
            .description("Critical switch pressure for ``seepage`` and ``river`` boundary conditions.")
            .input_default("0.0")
            .units( UnitSI().m() );


    //these are for unsteady
    *this += init_pressure.name("init_pressure")
            .description("Initial condition for pressure in time dependent problems.")
            .input_default("0.0")
            .units( UnitSI().m() );

    *this += storativity.name("storativity")
            .description("Storativity (in time dependent problems).")
            .input_default("0.0")
            .units( UnitSI().m(-1) );
    
    *this += extra_storativity.name("extra_storativity")
            .description("Storativity added from upstream equation.")
            .units( UnitSI().m(-1) )
            .input_default("0.0")
            .flags( input_copy );
    
    *this += extra_source.name("extra_water_source_density")
            .description("Water source density added from upstream equation.")
            .input_default("0.0")
            .units( UnitSI().s(-1) )
            .flags( input_copy );

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
DarcyMH::DarcyMH(Mesh &mesh_in, const Input::Record in_rec, TimeGovernor *tm)
: DarcyFlowInterface(mesh_in, in_rec),
    output_object(nullptr),
    solution(nullptr),
    data_changed_(false),
    schur0(nullptr),
	steady_diagonal(nullptr),
	steady_rhs(nullptr),
	new_diagonal(nullptr),
	previous_solution(nullptr)
{

    START_TIMER("Darcy constructor");
    {
        auto time_rec = in_rec.find<Input::Record>("time");
        if (tm == nullptr)
        {
            if (time_rec)
                time_ = new TimeGovernor(*time_rec);
            else
                time_ = new TimeGovernor();
        }
        else
        {
            if (time_rec)
                WarningOut() << "Time governor of DarcyFlowMH is initialized from parent class - input record will be ignored!";
            time_ = tm;
        }
    }

    data_ = make_shared<EqData>();
    EquationBase::eq_data_ = data_.get();
    
    data_->is_linear=true;

    size = mesh_->n_elements() + mesh_->n_sides() + mesh_->n_edges();
    n_schur_compls = in_rec.val<int>("n_schurs");
    data_->mortar_method_= in_rec.val<MortarMethod>("mortar_method");
    if (data_->mortar_method_ != NoMortar) {
        mesh_->mixed_intersections();
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
    data_->init_pressure.add_factory(
            std::make_shared<FieldAddPotential<3, FieldValue<3>::Scalar>::FieldFactory>
            (data_->gravity_, "init_piezo_head") );


    data_->set_input_list( this->input_record_.val<Input::Array>("input_fields"), *time_ );
    // Check that the time step was set for the transient simulation.
    if (! zero_time_term(true) && time_->is_default() ) {
        //THROW(ExcAssertMsg());
        //THROW(ExcMissingTimeGovernor() << input_record_.ei_address());
        MessageOut() << "Missing the key 'time', obligatory for the transient problems." << endl;
        ASSERT(false);
    }

    data_->mark_input_times(*time_);
}



void DarcyMH::initialize() {

    init_eq_data();
    output_object = new DarcyFlowMHOutput(this, input_record_);

    mh_dh.reinit(mesh_);
    // Initialize bc_switch_dirichlet to size of global boundary.
    data_->bc_switch_dirichlet.resize(mesh_->n_elements()+mesh_->n_elements(true), 1);


    nonlinear_iteration_=0;
    Input::AbstractRecord rec = this->input_record_
            .val<Input::Record>("nonlinear_solver")
            .val<Input::AbstractRecord>("linear_solver");

    // auxiliary set_time call  since allocation assembly evaluates fields as well
    data_changed_ = data_->set_time(time_->step(), LimitSide::right) || data_changed_;
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
    data_->water_balance_idx = balance_->add_quantity("water_volume");
    balance_->allocate(mh_dh.rows_ds->lsize(), 1);
    balance_->units(UnitSI().m(3));


    data_->balance = balance_;
    data_->lin_sys = schur0;


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

    data_changed_ = data_->set_time(time_->step(), LimitSide::right) || data_changed_;

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
//         print_matlab_matrix("matrix_zero");
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

    solve_time_step();
}


void DarcyMH::solve_time_step(bool output)
{
    data_changed_ = data_->set_time(time_->step(), LimitSide::left) || data_changed_;
    
    bool zero_time_term_from_left=zero_time_term();

    bool jump_time = data_->storativity.is_jump_time();
    if (! zero_time_term_from_left) {
        // time term not treated as zero
        // Unsteady solution up to the T.

        // this flag is necesssary for switching BC to avoid setting zero neumann on the whole boundary in the steady case
        use_steady_assembly_ = false;

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
        if (! zero_time_term_from_left && ! jump_time && output) output_data();
        return;
    }

    data_changed_ = data_->set_time(time_->step(), LimitSide::right) || data_changed_;
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
    if (output) output_data();
}


bool DarcyMH::zero_time_term(bool time_global) {
    if (time_global) {
        return (data_->storativity.input_list_size() == 0);
    } else {
        return data_->storativity.field_result(mesh_->region_db().get_region_set("BULK")) == result_zeros;
    }
}


void DarcyMH::solve_nonlinear()
{

    assembly_linear_system();
    double residual_norm = schur0->compute_residual();
    nonlinear_iteration_ = 0;
    MessageOut().fmt("[nonlinear solver] norm of initial residual: {}\n", residual_norm);

    // Reduce is_linear flag.
    int is_linear_common;
    MPI_Allreduce(&(data_->is_linear), &is_linear_common,1, MPI_INT ,MPI_MIN,PETSC_COMM_WORLD);

    Input::Record nl_solver_rec = input_record_.val<Input::Record>("nonlinear_solver");
    this->tolerance_ = nl_solver_rec.val<double>("tolerance");
    this->max_n_it_  = nl_solver_rec.val<unsigned int>("max_it");
    this->min_n_it_  = nl_solver_rec.val<unsigned int>("min_it");
    if (this->min_n_it_ > this->max_n_it_) this->min_n_it_ = this->max_n_it_;

    if (! is_linear_common) {
        // set tolerances of the linear solver unless they are set by user.
        schur0->set_tolerances(0.1*this->tolerance_, 0.01*this->tolerance_, 100);
    }
    vector<double> convergence_history;

    Vec save_solution;
    VecDuplicate(schur0->get_solution(), &save_solution);
    while (nonlinear_iteration_ < this->min_n_it_ ||
           (residual_norm > this->tolerance_ &&  nonlinear_iteration_ < this->max_n_it_ )) {
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
        LinSys::SolveInfo si = schur0->solve();
        nonlinear_iteration_++;

        // hack to make BDDC work with empty compute_residual
        if (is_linear_common){
            // we want to print this info in linear (and steady) case
            residual_norm = schur0->compute_residual();
            MessageOut().fmt("[nonlinear solver] lin. it: {}, reason: {}, residual: {}\n",
        		si.n_iterations, si.converged_reason, residual_norm);
            break;
        }
        data_changed_=true; // force reassembly for non-linear case

        double alpha = 1; // how much of new solution
        VecAXPBY(schur0->get_solution(), (1-alpha), alpha, save_solution);

        /*
        double * sol;
        unsigned int sol_size;
        get_solution_vector(sol, sol_size);
        if (mh_dh.el_ds->myp() == 0)
            VecView(sol_vec, PETSC_VIEWER_STDOUT_SELF);
        */

        //LogOut().fmt("Linear solver ended with reason: {} \n", si.converged_reason );
        //OLD_ASSERT( si.converged_reason >= 0, "Linear solver failed to converge. Convergence reason %d \n", si.converged_reason );
        assembly_linear_system();

        residual_norm = schur0->compute_residual();
        MessageOut().fmt("[nonlinear solver] it: {} lin. it: {}, reason: {}, residual: {}\n",
        		nonlinear_iteration_, si.n_iterations, si.converged_reason, residual_norm);
    }
    chkerr(VecDestroy(&save_solution));
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
    VecSwap(previous_solution, schur0->get_solution());
}

void DarcyMH::postprocess() 
{
    START_TIMER("postprocess");

    //fix velocity when mortar method is used
    if(data_->mortar_method_ != MortarMethod::NoMortar){
        auto multidim_assembler =  AssemblyBase::create< AssemblyMH >(data_);
        for (unsigned int i_loc = 0; i_loc < mh_dh.el_ds->lsize(); i_loc++) {
            auto ele_ac = mh_dh.accessor(i_loc);
            unsigned int dim = ele_ac.dim();
            multidim_assembler[dim-1]->fix_velocity(ele_ac);
        }
    }
    //ElementAccessor<3> ele;

    // modify side fluxes in parallel
    // for every local edge take time term on digonal and add it to the corresponding flux
    /*
    for (unsigned int i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {
        ele = mesh_->element_accessor(el_4_loc[i_loc]);
        for (unsigned int i=0; i<ele->n_sides(); i++) {
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


    START_TIMER("Darcy balance output");
    balance_->calculate_cumulative(data_->water_balance_idx, schur0->get_solution());
    balance_->calculate_instant(data_->water_balance_idx, schur0->get_solution());
    balance_->output();
    
    prepare_new_time_step(); //SWAP
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


// ===========================================================================================
//
//   MATRIX ASSEMBLY - we use abstract assembly routine, where  LS Mat/Vec SetValues
//   are in fact pointers to allocating or filling functions - this is governed by Linsystem roitunes
//
// =======================================================================================
void DarcyMH::assembly_mh_matrix(MultidimAssembly& assembler)
{
    START_TIMER("DarcyFlowMH_Steady::assembly_steady_mh_matrix");

    // set auxiliary flag for switchting Dirichlet like BC
    data_->force_bc_switch = use_steady_assembly_ && (nonlinear_iteration_ == 0);
    data_->n_schur_compls = n_schur_compls;
    

    balance_->start_flux_assembly(data_->water_balance_idx);

    // TODO: try to move this into balance, or have it in the generic assembler class, that should perform the cell loop
    // including various pre- and post-actions
    for (unsigned int i_loc = 0; i_loc < mh_dh.el_ds->lsize(); i_loc++) {
        auto ele_ac = mh_dh.accessor(i_loc);
        unsigned int dim = ele_ac.dim();
        assembler[dim-1]->assemble(ele_ac);
    }    
    

    balance_->finish_flux_assembly(data_->water_balance_idx);

}


void DarcyMH::allocate_mh_matrix()
{
    START_TIMER("DarcyFlowMH_Steady::allocate_mh_matrix");

    // set auxiliary flag for switchting Dirichlet like BC
    data_->n_schur_compls = n_schur_compls;
    LinSys *ls = schur0;
   


    int local_dofs[10];

    // to make space for second schur complement, max. 10 neighbour edges of one el.
    double zeros[100000];
    for(int i=0; i<100000; i++) zeros[i] = 0.0;

    std::vector<int> tmp_rows;
    tmp_rows.reserve(200);
    
    unsigned int nsides, loc_size;

    for (unsigned int i_loc = 0; i_loc < mh_dh.el_ds->lsize(); i_loc++) {
        auto ele_ac = mh_dh.accessor(i_loc);
        nsides = ele_ac.n_sides();
        
        //allocate at once matrix [sides,ele,edges]x[sides,ele,edges]
        loc_size = 1 + 2*nsides;
        unsigned int i_side = 0;
        
        for (; i_side < nsides; i_side++) {
            local_dofs[i_side] = ele_ac.side_row(i_side);
            local_dofs[i_side+nsides] = ele_ac.edge_row(i_side);
        }
        local_dofs[i_side+nsides] = ele_ac.ele_row();
        int * edge_rows = local_dofs + nsides;
        //int ele_row = local_dofs[0];
        
        // whole local MH matrix
        ls->mat_set_values(loc_size, local_dofs, loc_size, local_dofs, zeros);
        

        // compatible neighborings rows
        unsigned int n_neighs = ele_ac.element_accessor()->n_neighs_vb();
        for (unsigned int i = 0; i < n_neighs; i++) {
            // every compatible connection adds a 2x2 matrix involving
            // current element pressure  and a connected edge pressure
            Neighbour *ngh = ele_ac.element_accessor()->neigh_vb[i];
            int neigh_edge_row = mh_dh.row_4_edge[ ngh->edge_idx() ];
            tmp_rows.push_back(neigh_edge_row);
            //DebugOut() << "CC" << print_var(tmp_rows[i]);
        }

        // allocate always also for schur 2
        ls->mat_set_values(nsides+1, edge_rows, n_neighs, tmp_rows.data(), zeros); // (edges, ele)  x (neigh edges)
        ls->mat_set_values(n_neighs, tmp_rows.data(), nsides+1, edge_rows, zeros); // (neigh edges) x (edges, ele)
        ls->mat_set_values(n_neighs, tmp_rows.data(), n_neighs, tmp_rows.data(), zeros);  // (neigh edges) x (neigh edges)

        tmp_rows.clear();

        if (data_->mortar_method_ != NoMortar) {
            auto &isec_list = mesh_->mixed_intersections().element_intersections_[ele_ac.ele_global_idx()];
            for(auto &isec : isec_list ) {
                IntersectionLocalBase *local = isec.second;
                ElementAccessor<3> slave_ele = mesh_->element_accessor( local->bulk_ele_idx() );
                //DebugOut().fmt("Alloc: {} {}", ele_ac.ele_global_idx(), local->bulk_ele_idx());
                for(unsigned int i_side=0; i_side < slave_ele->n_sides(); i_side++) {
                    tmp_rows.push_back( mh_dh.row_4_edge[ slave_ele.side(i_side)->edge_idx() ] );
                    //DebugOut() << "aedge" << print_var(tmp_rows[i_rows-1]);
                }
            }
        }
        /*
        for(unsigned int i_side=0; i_side < ele_ac.element_accessor()->n_sides(); i_side++) {
            DebugOut() << "aedge:" << print_var(edge_rows[i_side]);
        }*/

        ls->mat_set_values(nsides, edge_rows, tmp_rows.size(), tmp_rows.data(), zeros);   // master edges x neigh edges
        ls->mat_set_values(tmp_rows.size(), tmp_rows.data(), nsides, edge_rows, zeros);   // neigh edges  x master edges
        ls->mat_set_values(tmp_rows.size(), tmp_rows.data(), tmp_rows.size(), tmp_rows.data(), zeros);  // neigh edges  x neigh edges

    }
/*
    // alloc edge diagonal entries
    if(rank == 0)
    for( vector<Edge>::iterator edg = mesh_->edges.begin(); edg != mesh_->edges.end(); ++edg) {
        int edg_idx = mh_dh.row_4_edge[edg->side(0)->edge_idx()];
        
//        for( vector<Edge>::iterator edg2 = mesh_->edges.begin(); edg2 != mesh_->edges.end(); ++edg2){
//            int edg_idx2 = mh_dh.row_4_edge[edg2->side(0)->edge_idx()];
//            if(edg_idx == edg_idx2){
//                 DBGCOUT(<< "P[ " << rank << " ] " << "edg alloc: " << edg_idx << "  " << edg_idx2 << "\n");
                ls->mat_set_value(edg_idx, edg_idx, 0.0);
//            }
//        }
    }
  */
    /*
    if (mortar_method_ == MortarP0) {
        P0_CouplingAssembler(*this).assembly(*ls);
    } else if (mortar_method_ == MortarP1) {
        P1_CouplingAssembler(*this).assembly(*ls);
    }*/
}

void DarcyMH::assembly_source_term()
{
    START_TIMER("assembly source term");
   	balance_->start_source_assembly(data_->water_balance_idx);

    for (unsigned int i_loc = 0; i_loc < mh_dh.el_ds->lsize(); i_loc++) {
        auto ele_ac = mh_dh.accessor(i_loc);

        double cs = data_->cross_section.value(ele_ac.centre(), ele_ac.element_accessor());

        // set sources
        double source = ele_ac.measure() * cs *
                (data_->water_source_density.value(ele_ac.centre(), ele_ac.element_accessor())
                +data_->extra_source.value(ele_ac.centre(), ele_ac.element_accessor()));
        schur0->rhs_set_value(ele_ac.ele_row(), -1.0 * source );

        balance_->add_source_values(data_->water_balance_idx, ele_ac.region().bulk_idx(), {(LongIdx) ele_ac.ele_local_row()}, {0}, {source});
    }

    balance_->finish_source_assembly(data_->water_balance_idx);
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
            ls->set_solution();
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

                ls->set_solution();
                schur0=ls;
            } else {
                IS is;
                ISCreateStride(PETSC_COMM_WORLD, mh_dh.side_ds->lsize(), mh_dh.rows_ds->begin(), 1, &is);
                //OLD_ASSERT(err == 0,"Error in ISCreateStride.");

                SchurComplement *ls = new SchurComplement(&(*mh_dh.rows_ds), is);

                // make schur1
                Distribution *ds = ls->make_complement_distribution();
                if (n_schur_compls==1) {
                    schur1 = new LinSys_PETSC(ds);
                    schur1->set_positive_definite();
                } else {
                    IS is;
                    ISCreateStride(PETSC_COMM_WORLD, mh_dh.el_ds->lsize(), ls->get_distribution()->begin(), 1, &is);
                    //OLD_ASSERT(err == 0,"Error in ISCreateStride.");
                    SchurComplement *ls1 = new SchurComplement(ds, is); // is is deallocated by SchurComplement
                    ls1->set_negative_definite();

                    // make schur2
                    schur2 = new LinSys_PETSC( ls1->make_complement_distribution() );
                    schur2->set_positive_definite();
                    ls1->set_complement( schur2 );
                    schur1 = ls1;
                }
                ls->set_complement( schur1 );
                ls->set_from_input(in_rec);
                ls->set_solution();
                schur0=ls;
            }

            START_TIMER("PETSC PREALLOCATION");
            schur0->set_symmetric();
            schur0->start_allocation();
            
            allocate_mh_matrix();
            
    	    VecZeroEntries(schur0->get_solution());
            END_TIMER("PETSC PREALLOCATION");
        }
        else {
            xprintf(Err, "Unknown solver type. Internal error.\n");
        }

        END_TIMER("preallocation");
        make_serial_scatter();

    }

}




void DarcyMH::assembly_linear_system() {
    START_TIMER("DarcyFlowMH_Steady::assembly_linear_system");

    data_->is_linear=true;
    bool is_steady = zero_time_term();
	//DebugOut() << "Assembly linear system\n";
	if (data_changed_) {
	    data_changed_ = false;
	    //DebugOut()  << "Data changed\n";
		// currently we have no optimization for cases when just time term data or RHS data are changed
	    START_TIMER("full assembly");
        if (typeid(*schur0) != typeid(LinSys_BDDC)) {
            schur0->start_add_assembly(); // finish allocation and create matrix
        }

        schur0->mat_zero_entries();
        schur0->rhs_zero_entries();

        assembly_source_term();
        
        auto multidim_assembler =  AssemblyBase::create< AssemblyMH >(data_);
	    assembly_mh_matrix( multidim_assembler ); // fill matrix

	    schur0->finish_assembly();
//         print_matlab_matrix("matrix");
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
	    else
	    {
	    	balance_->start_mass_assembly(data_->water_balance_idx);
	    	balance_->finish_mass_assembly(data_->water_balance_idx);
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


void DarcyMH::print_matlab_matrix(std::string matlab_file)
{
    std::string output_file;
    
    if ( typeid(*schur0) == typeid(LinSys_BDDC) ){
//         WarningOut() << "Can output matrix only on a single processor.";
//         output_file = FilePath(matlab_file + "_bddc.m", FilePath::output_file);
//         ofstream os( output_file );
//         auto bddc = static_cast<LinSys_BDDC*>(schur0);
//         bddc->print_matrix(os);
    }
    else {//if ( typeid(*schur0) == typeid(LinSys_PETSC) ){
        output_file = FilePath(matlab_file + ".m", FilePath::output_file);
        PetscViewer    viewer;
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, output_file.c_str(), &viewer);
        PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
        MatView( *const_cast<Mat*>(schur0->get_matrix()), viewer);
        VecView( *const_cast<Vec*>(schur0->get_rhs()), viewer);
    }
//     else{
//         WarningOut() << "No matrix output available for the current solver.";
//         return;
//     }
    
    // compute h_min for different dimensions
    double d_max = std::numeric_limits<double>::max();
    double h1 = d_max, h2 = d_max, h3 = d_max;
    double he2 = d_max, he3 = d_max;
    for (auto ele : mesh_->elements_range()) {
        switch(ele->dim()){
            case 1: h1 = std::min(h1,ele.measure()); break;
            case 2: h2 = std::min(h2,ele.measure()); break;
            case 3: h3 = std::min(h3,ele.measure()); break;
        }
        
        for (unsigned int j=0; j<ele->n_sides(); j++) {
            switch(ele->dim()){
                case 2: he2 = std::min(he2, ele.side(j)->measure()); break;
                case 3: he3 = std::min(he3, ele.side(j)->measure()); break;
            }
        }
    }
    if(h1 == d_max) h1 = 0;
    if(h2 == d_max) h2 = 0;
    if(h3 == d_max) h3 = 0;
    if(he2 == d_max) he2 = 0;
    if(he3 == d_max) he3 = 0;
    
    FILE * file;
    file = fopen(output_file.c_str(),"a");
    fprintf(file, "nA = %d;\n", mh_dh.side_ds->size());
    fprintf(file, "nB = %d;\n", mh_dh.el_ds->size());
    fprintf(file, "nBF = %d;\n", mh_dh.edge_ds->size());
    fprintf(file, "h1 = %e;\nh2 = %e;\nh3 = %e;\n", h1, h2, h3);
    fprintf(file, "he2 = %e;\nhe3 = %e;\n", he2, he3);
    fclose(file);
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

        for (unsigned int si=0; si<ele_ac.element_accessor()->n_sides(); si++) {
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

        for (unsigned int si=0; si<ele_ac.element_accessor()->n_sides(); si++) {
            // insert local edge dof
            int edge_row = ele_ac.edge_row(si);
            arma::vec3 coord = ele_ac.side(si)->centre();

            localDofMap.insert( std::make_pair( edge_row, coord ) );
            inet.push_back( edge_row );
            nne++;
        }

        // insert dofs related to compatible connections
        for ( unsigned int i_neigh = 0; i_neigh < ele_ac.element_accessor()->n_neighs_vb(); i_neigh++) {
            int edge_row = mh_dh.row_4_edge[ ele_ac.element_accessor()->neigh_vb[i_neigh]->edge_idx()  ];
            arma::vec3 coord = ele_ac.element_accessor()->neigh_vb[i_neigh]->edge()->side(0)->centre();

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
	if (previous_solution != nullptr) chkerr(VecDestroy(&previous_solution));
	if (steady_diagonal != nullptr) chkerr(VecDestroy(&steady_diagonal));
	if (new_diagonal != nullptr) chkerr(VecDestroy(&new_diagonal));
	if (steady_rhs != nullptr) chkerr(VecDestroy(&steady_rhs));
    
    
    if (schur0 != NULL) {
        delete schur0;
        chkerr(VecScatterDestroy(&par_to_all));
    }

	if (solution != NULL) {
	    chkerr(VecDestroy(&sol_vec));
		delete [] solution;
	}

	if (output_object)	delete output_object;

    if(time_ != nullptr)
        delete time_;
    
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
            for (auto ele : mesh_->elements_range()) {
            	for (unsigned int si=0; si<ele->n_sides(); si++) {
                    loc_idx[i++] = mh_dh.side_row_4_id[ mh_dh.side_dof( ele.side(si) ) ];
                }
            }
            for (auto ele : mesh_->elements_range()) {
                loc_idx[i++] = mh_dh.row_4_el[ ele.idx() ];
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
            chkerr(ISDestroy(&(is_loc)));
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

   	balance_->start_mass_assembly(data_->water_balance_idx);

    //DebugOut().fmt("time_term lsize: {} {}\n", mh_dh.el_ds->myp(), mh_dh.el_ds->lsize());
    for (unsigned int i_loc_el = 0; i_loc_el < mh_dh.el_ds->lsize(); i_loc_el++) {
        auto ele_ac = mh_dh.accessor(i_loc_el);

        // set new diagonal
        double diagonal_coeff = data_->cross_section.value(ele_ac.centre(), ele_ac.element_accessor())
        		* ( data_->storativity.value(ele_ac.centre(), ele_ac.element_accessor())
                   +data_->extra_storativity.value(ele_ac.centre(), ele_ac.element_accessor())
                  )
				* ele_ac.measure();
        local_diagonal[ele_ac.ele_local_row()]= - diagonal_coeff / time_->dt();

        //DebugOut().fmt("time_term: {} {} {} {} {}\n", mh_dh.el_ds->myp(), ele_ac.ele_global_idx(), i_loc_row, i_loc_el + mh_dh.side_ds->lsize(), diagonal_coeff);
       	balance_->add_mass_matrix_values(data_->water_balance_idx,
       	        ele_ac.region().bulk_idx(), { LongIdx(ele_ac.ele_row()) }, {diagonal_coeff});
    }
    VecRestoreArray(new_diagonal,& local_diagonal);
    MatDiagonalSet(*( schur0->get_matrix() ), new_diagonal, ADD_VALUES);

    solution_changed_for_scatter=true;
    schur0->set_matrix_changed();

    balance_->finish_mass_assembly(data_->water_balance_idx);
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
    VecPointwiseMult(*( schur0->get_rhs()), new_diagonal, previous_solution);
// 	VecPointwiseMult(*( schur0->get_rhs()), new_diagonal, schur0->get_solution());
    VecAXPY(*( schur0->get_rhs()), 1.0, steady_rhs);
    schur0->set_rhs_changed();

    //VecSwap(previous_solution, schur0->get_solution());
}


//-----------------------------------------------------------------------------
// vim: set cindent:
