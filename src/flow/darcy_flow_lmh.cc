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
 * @file    darcy_flow_lmh.cc
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
#include "system/index_types.hh"
#include "input/factory.hh"

#include "mesh/mesh.h"
#include "mesh/bc_mesh.hh"
#include "mesh/partitioning.hh"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "la/distribution.hh"
#include "la/linsys.hh"
#include "la/linsys_PETSC.hh"
// #include "la/linsys_BDDC.hh"
#include "la/schur.hh"
//#include "la/sparse_graph.hh"
#include "la/local_to_global_map.hh"
#include "la/vector_mpi.hh"

//#include "flow/assembly_lmh_old_.hh"
#include "flow/darcy_flow_lmh.hh"
#include "flow/darcy_flow_mh_output.hh"
#include "flow/assembly_lmh.hh"
#include "flow/assembly_models.hh"

#include "tools/time_governor.hh"
#include "fields/field_algo_base.hh"
#include "fields/field.hh"
#include "fields/field_values.hh"
#include "fields/field_add_potential.hh"
#include "fields/field_fe.hh"
#include "fields/field_model.hh"
#include "fields/field_constant.hh"

#include "coupling/balance.hh"

#include "intersection/mixed_mesh_intersections.hh"
#include "intersection/intersection_local.hh"

#include "fem/fe_p.hh"


FLOW123D_FORCE_LINK_IN_CHILD(darcy_flow_lmh)




namespace it = Input::Type;

const it::Selection & DarcyLMH::get_mh_mortar_selection() {
	return it::Selection("MH_MortarMethod")
		.add_value(NoMortar, "None", "No Mortar method is applied.")
		.add_value(MortarP0, "P0", "Mortar space: P0 on elements of lower dimension.")
		.add_value(MortarP1, "P1", "Mortar space: P1 on intersections, using non-conforming pressures.")
		.close();
}

const it::Record & DarcyLMH::type_field_descriptor() {

        const it::Record &field_descriptor =
        it::Record(equation_name() + "_Data",FieldCommon::field_descriptor_record_description(equation_name() + "_Data") )
        .copy_keys( DarcyLMH::EqFields().make_field_descriptor_type(equation_name() + "_Data_aux") )
            .declare_key("bc_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(),
                    "Boundary piezometric head for BC types: dirichlet, robin, and river." )
            .declare_key("bc_switch_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(),
                    "Boundary switch piezometric head for BC types: seepage, river." )
            .declare_key("init_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(),
                    "Initial condition for the pressure given as the piezometric head." )
            .close();
        return field_descriptor;
}

const it::Record & DarcyLMH::get_input_type() {
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

    DarcyLMH::EqFields eq_fields;
    
    return it::Record(equation_name(), "Lumped Mixed-Hybrid solver for saturated Darcy flow.")
		.derive_from(DarcyFlowInterface::get_input_type())
        .copy_keys(EquationBase::record_template())
		.copy_keys(EquationBase::user_fields_template(equation_name()))
        .declare_key("gravity", it::Array(it::Double(), 3,3), it::Default("[ 0, 0, -1]"),
                "Vector of the gravity force. Dimensionless.")
		.declare_key("input_fields", it::Array( type_field_descriptor() ), it::Default::obligatory(),
                "Input data for Darcy flow model.")				
        .declare_key("nonlinear_solver", ns_rec, it::Default("{}"),
                "Non-linear solver for MH problem.")
        .declare_key("output_stream", OutputTime::get_input_type(), it::Default("{}"),
                "Output stream settings.\n Specify file format, precision etc.")

        .declare_key("output", DarcyFlowMHOutput::get_input_type(eq_fields, equation_name()),
                IT::Default("{ \"fields\": [ \"pressure_p0\", \"velocity_p0\" ] }"),
                "Specification of output fields and output times.")
        .declare_key("output_specific", DarcyFlowMHOutput::get_input_type_specific(), it::Default::optional(),
                "Output settings specific to Darcy flow model.\n"
                "Includes raw output and some experimental functionality.")
        .declare_key("balance", Balance::get_input_type(), it::Default("{}"),
                "Settings for computing mass balance.")
		.declare_key("mortar_method", get_mh_mortar_selection(), it::Default("\"None\""),
				"Method for coupling Darcy flow between dimensions on incompatible meshes. [Experimental]" )
		.close();
}


const int DarcyLMH::registrar =
		Input::register_class< DarcyLMH, Mesh &, const Input::Record >(equation_name()) +
		DarcyLMH::get_input_type().size();



DarcyLMH::EqFields::EqFields()
{
    *this += field_ele_pressure.name("pressure_p0")
             .units(UnitSI().m())
             .flags(FieldFlag::equation_result)
             .description("Pressure solution - P0 interpolation.");

    *this += field_edge_pressure.name("pressure_edge")
             .units(UnitSI().m())
             .flags(FieldFlag::input_copy)
             .description("Pressure solution - Crouzeix-Raviart interpolation.");

    *this += field_ele_piezo_head.name("piezo_head_p0")
	         .units(UnitSI().m())
             .flags(FieldFlag::equation_result)
             .description("Piezo head solution - P0 interpolation.");

	*this += field_ele_velocity.name("velocity_p0")
	         .units(UnitSI().m().s(-1))
             .flags(FieldFlag::equation_result)
             .description("Velocity solution - P0 interpolation.");

    *this += flux.name("flux")
	         .units(UnitSI().m().s(-1))
             .flags(FieldFlag::equation_result)
             .description("Darcy flow flux.");

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

    *this += gravity_field.name("gravity")
            .description("Gravity vector.")
            .input_default("0.0")
            .units( UnitSI::dimensionless() );

    *this += bc_gravity.name("bc_gravity")
            .description("Boundary gravity vector.")
            .input_default("0.0")
            .units( UnitSI::dimensionless() );

    *this += init_piezo_head.name("init_piezo_head")
	         .units(UnitSI().m())
             .input_default("0.0")
             .description("Init piezo head.");

    *this += bc_piezo_head.name("bc_piezo_head")
	         .units(UnitSI().m())
             .input_default("0.0")
             .description("Boundary piezo head.");

    *this += bc_switch_piezo_head.name("bc_switch_piezo_head")
	         .units(UnitSI().m())
             .input_default("0.0")
             .description("Boundary switch piezo head.");

    *this += ref_pressure.name("ref_pressure")
             .units(UnitSI().m())
			 .input_default("0.0")
             .flags(FieldFlag::equation_result)
             .description("Precomputed pressure of l2 difference output.");

	*this += ref_velocity.name("ref_velocity")
	         .units(UnitSI().m().s(-1))
			 .input_default("0.0")
             .flags(FieldFlag::equation_result)
             .description("Precomputed velocity of l2 difference output.");

    *this += ref_divergence.name("ref_divergence")
             .units(UnitSI().m())
			 .input_default("0.0")
             .flags(FieldFlag::equation_result)
             .description("Precomputed divergence of l2 difference output.");

    this->set_default_fieldset();
    //time_term_fields = this->subset({"storativity"});
    //main_matrix_fields = this->subset({"anisotropy", "conductivity", "cross_section", "sigma", "bc_type", "bc_robin_sigma"});
    //rhs_fields = this->subset({"water_source_density", "bc_pressure", "bc_flux"});
}



const it::Selection & DarcyLMH::EqFields::get_bc_type_selection() {
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



DarcyLMH::EqData::EqData()
{
    mortar_method_=NoMortar;
}


void DarcyLMH::EqData::init()
{
    auto size = dh_p_->get_local_to_global_map().size();
    save_local_system_.resize(size);
    bc_fluxes_reconstruted.resize(size);
    loc_system_.resize(size);
    postprocess_solution_.resize(size);
}


void DarcyLMH::EqData::reset()
{
    std::fill(save_local_system_.begin(), save_local_system_.end(), false);
    std::fill(bc_fluxes_reconstruted.begin(), bc_fluxes_reconstruted.end(), false);
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
DarcyLMH::DarcyLMH(Mesh &mesh_in, const Input::Record in_rec, TimeGovernor *tm)
: DarcyFlowInterface(mesh_in, in_rec),
    output_object(nullptr),
    data_changed_(false),
	read_init_cond_assembly_(nullptr),
	mh_matrix_assembly_(nullptr),
	reconstruct_schur_assembly_(nullptr)
{

    START_TIMER("Darcy constructor");
    {
        auto time_record = input_record_.val<Input::Record>("time");
        if (tm == nullptr)
        {
            time_ = new TimeGovernor(time_record);
        }
        else
        {
            TimeGovernor tm_from_rec(time_record);
            if (!tm_from_rec.is_default()) // is_default() == false when time record is present in input file
            { 
                MessageOut() << "Duplicate key 'time', time in flow equation is already initialized from parent class!";
                ASSERT_PERMANENT(false);
            }
            time_ = tm;
        }
    }

    eq_fields_ = make_shared<EqFields>();
    eq_data_ = make_shared<EqData>();
    this->eq_fieldset_ = eq_fields_;
    
    eq_fields_->set_mesh(*mesh_);
    
    eq_data_->is_linear=true;

    size = mesh_->n_elements() + mesh_->n_sides() + mesh_->n_edges();
    eq_data_->mortar_method_= in_rec.val<MortarMethod>("mortar_method");
    if (eq_data_->mortar_method_ != NoMortar) {
        mesh_->mixed_intersections();
    }
    

    //side_ds->view( std::cout );
    //el_ds->view( std::cout );
    //edge_ds->view( std::cout );
    //rows_ds->view( std::cout );
    
}

double DarcyLMH::solved_time()
{
    // DebugOut() << "t = " << time_->t() << " step_end " << time_->step().end() << "\n";
    if(eq_data_->use_steady_assembly_)
    {
        // In steady case, the solution is computed with the data present at time t,
        // and the steady state solution is valid until another change in data,
        // which should correspond to time (t+dt).
        // "The data change appears immediatly."
        double next_t = time_->t() + time_->estimate_dt();
        // DebugOut() << "STEADY next_t = " << next_t << "\n";
        return next_t * (1 - 2*std::numeric_limits<double>::epsilon());
    }
    else
    {
        // In unsteady case, the solution is computed with the data present at time t,
        // and the solution is valid at the time t+dt.
        // "The data change does not appear immediatly, it is integrated over time interval dt."
        // DebugOut() << "UNSTEADY\n";
        return time_->t();
    }
}

void DarcyLMH::init_eq_data()
//connecting data fields with mesh
{

    START_TIMER("data init");
    eq_data_->mesh = mesh_;

    auto gravity_array = input_record_.val<Input::Array>("gravity");
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


    eq_fields_->set_input_list( this->input_record_.val<Input::Array>("input_fields"), *time_ );

    // Check that the time step was set for the transient simulation.
    if (! zero_time_term(true) && time_->is_default() ) {
        //THROW(ExcAssertMsg());
        //THROW(ExcMissingTimeGovernor() << input_record_.ei_address());
        MessageOut() << "Missing the key 'time', obligatory for the transient problems." << endl;
        ASSERT_PERMANENT(false);
    }

    eq_fields_->mark_input_times(*time_);
}

void DarcyLMH::initialize() {

    { // init DOF handler for pressure fields
// 		std::shared_ptr< FiniteElement<0> > fe0_rt = std::make_shared<FE_RT0_disc<0>>();
		std::shared_ptr< FiniteElement<1> > fe1_rt = std::make_shared<FE_RT0_disc<1>>();
		std::shared_ptr< FiniteElement<2> > fe2_rt = std::make_shared<FE_RT0_disc<2>>();
		std::shared_ptr< FiniteElement<3> > fe3_rt = std::make_shared<FE_RT0_disc<3>>();
		std::shared_ptr< FiniteElement<0> > fe0_disc = std::make_shared<FE_P_disc<0>>(0);
		std::shared_ptr< FiniteElement<1> > fe1_disc = std::make_shared<FE_P_disc<1>>(0);
		std::shared_ptr< FiniteElement<2> > fe2_disc = std::make_shared<FE_P_disc<2>>(0);
		std::shared_ptr< FiniteElement<3> > fe3_disc = std::make_shared<FE_P_disc<3>>(0);
		std::shared_ptr< FiniteElement<0> > fe0_cr = std::make_shared<FE_CR<0>>();
		std::shared_ptr< FiniteElement<1> > fe1_cr = std::make_shared<FE_CR<1>>();
		std::shared_ptr< FiniteElement<2> > fe2_cr = std::make_shared<FE_CR<2>>();
		std::shared_ptr< FiniteElement<3> > fe3_cr = std::make_shared<FE_CR<3>>();
// 	    static FiniteElement<0> fe0_sys = FE_P_disc<0>(0); //TODO fix and use solution with FESystem<0>( {fe0_rt, fe0_disc, fe0_cr} )
		FESystem<0> fe0_sys( {fe0_disc, fe0_disc, fe0_cr} );
		FESystem<1> fe1_sys( {fe1_rt, fe1_disc, fe1_cr} );
		FESystem<2> fe2_sys( {fe2_rt, fe2_disc, fe2_cr} );
		FESystem<3> fe3_sys( {fe3_rt, fe3_disc, fe3_cr} );
	    MixedPtr<FESystem> fe_sys( std::make_shared<FESystem<0>>(fe0_sys), std::make_shared<FESystem<1>>(fe1_sys),
	                                    std::make_shared<FESystem<2>>(fe2_sys), std::make_shared<FESystem<3>>(fe3_sys) );
		std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh_, fe_sys);
		eq_data_->dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
		eq_data_->dh_->distribute_dofs(ds);
    }

    init_eq_data();
    output_object = new DarcyFlowMHOutput(this, input_record_);

    eq_fields_->add_coords_field();

    { // construct pressure, velocity and piezo head fields
		uint rt_component = 0;
		eq_data_->full_solution = eq_data_->dh_->create_vector();
        auto ele_flux_ptr = create_field_fe<3, FieldValue<3>::VectorFixed>(eq_data_->dh_, &eq_data_->full_solution, rt_component);
        eq_fields_->flux.set(ele_flux_ptr, 0.0);

		eq_fields_->field_ele_velocity.set(Model<3, FieldValue<3>::VectorFixed>::create(fn_mh_velocity(), eq_fields_->flux, eq_fields_->cross_section), 0.0);

		uint p_ele_component = 1;
        auto ele_pressure_ptr = create_field_fe<3, FieldValue<3>::Scalar>(eq_data_->dh_, &eq_data_->full_solution, p_ele_component);
        eq_fields_->field_ele_pressure.set(ele_pressure_ptr, 0.0);

        uint p_edge_component = 2;
        auto edge_pressure_ptr = create_field_fe<3, FieldValue<3>::Scalar>(eq_data_->dh_, &eq_data_->full_solution, p_edge_component);
        eq_fields_->field_edge_pressure.set(edge_pressure_ptr, 0.0);

        eq_fields_->field_ele_piezo_head.set(
                Model<3, FieldValue<3>::Scalar>::create(fn_mh_piezohead(), eq_fields_->gravity_field, eq_fields_->X(), eq_fields_->field_ele_pressure),
                0.0
        );
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
//     p_edge_solution = new VectorMPI(eq_data_->dh_cr_->distr()->lsize());
//     full_solution = new VectorMPI(eq_data_->dh_->distr()->lsize());
    // this creates mpi vector from DoFHandler, including ghost values
    eq_data_->p_edge_solution = eq_data_->dh_cr_->create_vector();
    eq_data_->p_edge_solution_previous = eq_data_->dh_cr_->create_vector();
    eq_data_->p_edge_solution_previous_time = eq_data_->dh_cr_->create_vector();

    // Initialize bc_switch_dirichlet to size of global boundary.
    eq_data_->bc_switch_dirichlet.resize(mesh_->n_elements()+mesh_->bc_mesh()->n_elements(), 1);


    eq_data_->nonlinear_iteration_=0;
    Input::AbstractRecord rec = this->input_record_
            .val<Input::Record>("nonlinear_solver")
            .val<Input::AbstractRecord>("linear_solver");

    initialize_specific();
    
    // auxiliary set_time call  since allocation assembly evaluates fields as well
    data_changed_ = eq_fields_->set_time(time_->step(), LimitSide::right) || data_changed_;
    create_linear_system(rec);


    // initialization of balance object
    balance_ = std::make_shared<Balance>("water", mesh_);
    balance_->init_from_input(input_record_.val<Input::Record>("balance"), time());
    eq_data_->water_balance_idx = balance_->add_quantity("water_volume");
    balance_->allocate(eq_data_->dh_, 1);
    balance_->units(UnitSI().m(3));

    eq_data_->balance_ = this->balance_;

    this->initialize_asm();
}

void DarcyLMH::initialize_specific()
{
	//eq_data_->multidim_assembler = AssemblyFlowBase::create< AssemblyLMH >(eq_fields_, eq_data_);
}

//void DarcyLMH::read_initial_condition()
//{
//	DebugOut().fmt("Read initial condition\n");
//
//	for ( DHCellAccessor dh_cell : eq_data_->dh_->own_range() ) {
//
//        LocDofVec p_indices = dh_cell.cell_with_other_dh(eq_data_->dh_p_.get()).get_loc_dof_indices();
//        ASSERT_DBG(p_indices.n_elem == 1);
//        LocDofVec l_indices = dh_cell.cell_with_other_dh(eq_data_->dh_cr_.get()).get_loc_dof_indices();
//        ElementAccessor<3> ele = dh_cell.elm();
//
//		// set initial condition
//        double init_value = eq_fields_->init_pressure.value(ele.centre(),ele);
//        unsigned int p_idx = eq_data_->dh_p_->parent_indices()[p_indices[0]];
//        eq_data_->full_solution.set(p_idx, init_value);
//
//        for (unsigned int i=0; i<ele->n_sides(); i++) {
//             uint n_sides_of_edge =  ele.side(i)->edge().n_sides();
//             unsigned int l_idx = eq_data_->dh_cr_->parent_indices()[l_indices[i]];
//             eq_data_->full_solution.add(l_idx, init_value/n_sides_of_edge);
//
//             eq_data_->p_edge_solution.add(l_indices[i], init_value/n_sides_of_edge);
//         }
//	}
//
//    initial_condition_postprocess();
//
//	eq_data_->full_solution.ghost_to_local_begin();
//	eq_data_->full_solution.ghost_to_local_end();
//
//	eq_data_->p_edge_solution.ghost_to_local_begin();
//    eq_data_->p_edge_solution.ghost_to_local_end();
//    eq_data_->p_edge_solution_previous_time.copy_from(eq_data_->p_edge_solution);
//}
//
//void DarcyLMH::initial_condition_postprocess()
//{}

void DarcyLMH::zero_time_step()
{

    /* TODO:
     * - Allow solution reconstruction (pressure and velocity) from initial condition on user request.
     * - Steady solution as an intitial condition may be forced by setting inti_time =-1, and set data for the steady solver in that time.
     *   Solver should be able to switch from and to steady case depending on the zero time term.
     */

    data_changed_ = eq_fields_->set_time(time_->step(), LimitSide::right) || data_changed_;

    // zero_time_term means steady case
    eq_data_->use_steady_assembly_ = zero_time_term();

    eq_data_->p_edge_solution.zero_entries();
    
    if (eq_data_->use_steady_assembly_) { // steady case
        MessageOut() << "Flow zero time step - steady case\n";
        //read_initial_condition(); // Possible solution guess for steady case.
        solve_nonlinear(); // with right limit data
    } else {
        MessageOut() << "Flow zero time step - unsteady case\n";
        eq_data_->time_step_ = time_->dt();
        assembly_linear_system();
        this->read_init_cond_asm();
        accept_time_step(); // accept zero time step, i.e. initial condition
        
        
        // we reconstruct the initial solution here
        // during the reconstruction assembly:
        // - the balance objects are actually allocated
        // - the full solution vector is computed
        START_TIMER("DarcyFlowMH::reconstruct_solution_from_schur");
        this->reconstruct_schur_assembly_->assemble(eq_data_->dh_);
        END_TIMER("DarcyFlowMH::reconstruct_solution_from_schur");
    }
    //solution_output(T,right_limit); // data for time T in any case
    output_data();
}

//=============================================================================
// COMPOSE and SOLVE WATER MH System possibly through Schur complements
//=============================================================================
void DarcyLMH::update_solution()
{
    START_TIMER("Solving MH system");

    time_->next_time();

    time_->view("DARCY"); //time governor information output

    solve_time_step();
    
    eq_data_->full_solution.local_to_ghost_begin();
    eq_data_->full_solution.local_to_ghost_end();
}

void DarcyLMH::solve_time_step(bool output)
{
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
        if(output)
            accept_time_step();
        if (jump_time) {
        	WarningOut() << "Output of solution discontinuous in time not supported yet.\n";
            //solution_output(T, left_limit); // output use time T- delta*dt
            //output_data();
        }
    }

    if (time_->is_end()) {
        // output for unsteady case, end_time should not be the jump time
        // but rether check that
        if (! zero_time_term_from_left && ! jump_time && output)
            output_data();
        return;
    }

    data_changed_ = eq_fields_->set_time(time_->step(), LimitSide::right) || data_changed_;
    bool zero_time_term_from_right=zero_time_term();
    if (zero_time_term_from_right) {
        MessageOut() << "Flow time step - steady case\n";
        // this flag is necesssary for switching BC to avoid setting zero neumann on the whole boundary in the steady case
    	eq_data_->use_steady_assembly_ = true;
        solve_nonlinear(); // with right limit data
        if(output)
            accept_time_step();

    } else if (! zero_time_term_from_left && jump_time) {
    	WarningOut() << "Discontinuous time term not supported yet.\n";
        //solution_transfer(); // internally call set_time(T, left) and set_time(T,right) again
        //solve_nonlinear(); // with right limit data
    }
    //solution_output(T,right_limit); // data for time T in any case
    if (output)
        output_data();
}

bool DarcyLMH::zero_time_term(bool time_global) {
    if (time_global) {
        return (eq_fields_->storativity.input_list_size() == 0);
    } else {
        return eq_fields_->storativity.field_result(mesh_->region_db().get_region_set("BULK")) == result_zeros;
    }
}


void DarcyLMH::solve_nonlinear()
{

    assembly_linear_system();
    double residual_norm = lin_sys_schur().compute_residual();
    eq_data_->nonlinear_iteration_ = 0;
    MessageOut().fmt("[nonlinear solver] norm of initial residual: {}\n", residual_norm);

    // Reduce is_linear flag.
    int is_linear_common;
    MPI_Allreduce(&(eq_data_->is_linear), &is_linear_common,1, MPI_INT ,MPI_MIN,PETSC_COMM_WORLD);

    Input::Record nl_solver_rec = input_record_.val<Input::Record>("nonlinear_solver");
    this->tolerance_ = nl_solver_rec.val<double>("tolerance");
    this->max_n_it_  = nl_solver_rec.val<unsigned int>("max_it");
    this->min_n_it_  = nl_solver_rec.val<unsigned int>("min_it");
    if (this->min_n_it_ > this->max_n_it_) this->min_n_it_ = this->max_n_it_;

    if (! is_linear_common) {
        // set tolerances of the linear solver unless they are set by user.
        lin_sys_schur().set_tolerances(0.1*this->tolerance_, 0.01*this->tolerance_, 10000, 100);
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
                THROW(ExcSolverDiverge() << EI_Reason("Stagnation."));
            }
        }

        if (! is_linear_common){
        	eq_data_->p_edge_solution_previous.copy_from(eq_data_->p_edge_solution);
        	eq_data_->p_edge_solution_previous.local_to_ghost_begin();
        	eq_data_->p_edge_solution_previous.local_to_ghost_end();
        }

        LinSys::SolveInfo si = lin_sys_schur().solve();
        MessageOut().fmt("[schur solver] lin. it: {}, reason: {}, residual: {}\n",
        		si.n_iterations, si.converged_reason, lin_sys_schur().compute_residual());
        
        eq_data_->nonlinear_iteration_++;

        // hack to make BDDC work with empty compute_residual
        if (is_linear_common){
            // we want to print this info in linear (and steady) case
            residual_norm = lin_sys_schur().compute_residual();
            MessageOut().fmt("[nonlinear solver] lin. it: {}, reason: {}, residual: {}\n",
        		si.n_iterations, si.converged_reason, residual_norm);
            break;
        }
        data_changed_=true; // force reassembly for non-linear case

        double alpha = 1; // how much of new solution
        VecAXPBY(eq_data_->p_edge_solution.petsc_vec(), (1-alpha), alpha, eq_data_->p_edge_solution_previous.petsc_vec());

        //LogOut().fmt("Linear solver ended with reason: {} \n", si.converged_reason );
        //ASSERT_PERMANENT_GE( si.converged_reason, 0).error("Linear solver failed to converge.\n");
        assembly_linear_system();

        residual_norm = lin_sys_schur().compute_residual();
        MessageOut().fmt("[nonlinear solver] it: {} lin. it: {}, reason: {}, residual: {}\n",
                eq_data_->nonlinear_iteration_, si.n_iterations, si.converged_reason, residual_norm);
    }
    
//    reconstruct_solution_from_schur(eq_data_->multidim_assembler);
    START_TIMER("DarcyFlowMH::reconstruct_solution_from_schur");
    this->reconstruct_schur_assembly_->assemble(eq_data_->dh_);
    END_TIMER("DarcyFlowMH::reconstruct_solution_from_schur");

    // adapt timestep
    if (! this->zero_time_term()) {
        double mult = 1.0;
        if (eq_data_->nonlinear_iteration_ < 3) mult = 1.6;
        if (eq_data_->nonlinear_iteration_ > 7) mult = 0.7;
        time_->set_upper_constraint(time_->dt() * mult, "Darcy adaptivity.");
        // int result = time_->set_upper_constraint(time_->dt() * mult, "Darcy adaptivity.");
        //DebugOut().fmt("time adaptivity, res: {} it: {} m: {} dt: {} edt: {}\n", result, nonlinear_iteration_, mult, time_->dt(), time_->estimate_dt());
    }
}


void DarcyLMH::accept_time_step()
{
	eq_data_->p_edge_solution_previous_time.copy_from(eq_data_->p_edge_solution);
	eq_data_->p_edge_solution_previous_time.local_to_ghost_begin();
	eq_data_->p_edge_solution_previous_time.local_to_ghost_end();
}


void DarcyLMH::output_data() {
    START_TIMER("Darcy output data");
    
    // print_matlab_matrix("matrix_" + std::to_string(time_->step().index()));
    
    //time_->view("DARCY"); //time governor information output
	this->output_object->output();


    START_TIMER("Darcy balance output");
    balance_->calculate_cumulative(eq_data_->water_balance_idx, eq_data_->full_solution.petsc_vec());
    balance_->calculate_instant(eq_data_->water_balance_idx, eq_data_->full_solution.petsc_vec());
    balance_->output();
}


//double DarcyLMH::solution_precision() const
//{
//    return eq_data_->lin_sys_schur->get_solution_precision();
//}


// ===========================================================================================
//
//   MATRIX ASSEMBLY - we use abstract assembly routine, where  LS Mat/Vec SetValues
//   are in fact pointers to allocating or filling functions - this is governed by Linsystem roitunes
//
// =======================================================================================
//void DarcyLMH::assembly_mh_matrix(FMT_UNUSED MultidimAssembly& assembler)
//{
//    START_TIMER("DarcyLMH::assembly_steady_mh_matrix");
//
//    // DebugOut() << "assembly_mh_matrix \n";
//    // set auxiliary flag for switchting Dirichlet like BC
//    eq_data_->force_no_neumann_bc = eq_data_->use_steady_assembly_ && (eq_data_->nonlinear_iteration_ == 0);
//
//    balance_->start_flux_assembly(eq_data_->water_balance_idx);
//    balance_->start_source_assembly(eq_data_->water_balance_idx);
//    balance_->start_mass_assembly(eq_data_->water_balance_idx);
//
//    // TODO: try to move this into balance, or have it in the generic assembler class, that should perform the cell loop
//    // including various pre- and post-actions
////    for ( DHCellAccessor dh_cell : eq_data_->dh_->own_range() ) {
////        unsigned int dim = dh_cell.dim();
////        assembler[dim-1]->assemble(dh_cell);
////    }
//    this->mh_matrix_assembly_->assemble(eq_data_->dh_);
//
//
//    balance_->finish_mass_assembly(eq_data_->water_balance_idx);
//    balance_->finish_source_assembly(eq_data_->water_balance_idx);
//    balance_->finish_flux_assembly(eq_data_->water_balance_idx);
//
//}


void DarcyLMH::allocate_mh_matrix()
{
    START_TIMER("DarcyLMH::allocate_mh_matrix");

    // to make space for second schur complement, max. 10 neighbour edges of one el.
    double zeros[100000];
    for(int i=0; i<100000; i++) zeros[i] = 0.0;

    std::vector<LongIdx> tmp_rows;
    tmp_rows.reserve(200);

    std::vector<LongIdx> dofs, dofs_ngh;
    dofs.reserve(eq_data_->dh_cr_->max_elem_dofs());
    dofs_ngh.reserve(eq_data_->dh_cr_->max_elem_dofs());

    // DebugOut() << "Allocate new schur\n";
    for ( DHCellAccessor dh_cell : eq_data_->dh_cr_->own_range() ) {
        ElementAccessor<3> ele = dh_cell.elm(); 

        const uint ndofs = dh_cell.n_dofs();
        dofs.resize(dh_cell.n_dofs());
        dh_cell.get_dof_indices(dofs);

        int* dofs_ptr = dofs.data();
        lin_sys_schur().mat_set_values(ndofs, dofs_ptr, ndofs, dofs_ptr, zeros);
        
        tmp_rows.clear();
        
        // compatible neighborings rows
        unsigned int n_neighs = ele->n_neighs_vb();
        for ( DHCellSide neighb_side : dh_cell.neighb_sides() ) {
            // every compatible connection adds a 2x2 matrix involving
            // current element pressure  and a connected edge pressure

            // read neighbor dofs (dh_cr dofhandler)
            // neighbor cell owning neighb_side
            DHCellAccessor dh_neighb_cell = neighb_side.cell();
            
            const uint ndofs_ngh = dh_neighb_cell.n_dofs();
            dofs_ngh.resize(ndofs_ngh);
            dh_neighb_cell.get_dof_indices(dofs_ngh);

            // local index of pedge dof on neighboring cell
            tmp_rows.push_back(dofs_ngh[neighb_side.side().side_idx()]);
        }
        
        lin_sys_schur().mat_set_values(ndofs, dofs_ptr, n_neighs, tmp_rows.data(), zeros); // (edges)  x (neigh edges)
        lin_sys_schur().mat_set_values(n_neighs, tmp_rows.data(), ndofs, dofs_ptr, zeros); // (neigh edges) x (edges)
        lin_sys_schur().mat_set_values(n_neighs, tmp_rows.data(), n_neighs, tmp_rows.data(), zeros);  // (neigh edges) x (neigh edges)

        tmp_rows.clear();
//        if (eq_data_->mortar_method_ != NoMortar) {
//            auto &isec_list = mesh_->mixed_intersections().element_intersections_[ele.idx()];
//            for(auto &isec : isec_list ) {
//                IntersectionLocalBase *local = isec.second;
//                DHCellAccessor dh_cell_slave = eq_data_->dh_cr_->cell_accessor_from_element(local->bulk_ele_idx());
//
//                const uint ndofs_slave = dh_cell_slave.n_dofs();
//                dofs_ngh.resize(ndofs_slave);
//                dh_cell_slave.get_dof_indices(dofs_ngh);
//
//                //DebugOut().fmt("Alloc: {} {}", ele.idx(), local->bulk_ele_idx());
//                for(unsigned int i_side=0; i_side < dh_cell_slave.elm()->n_sides(); i_side++) {
//                    tmp_rows.push_back( dofs_ngh[i_side] );
//                    //DebugOut() << "aedge" << print_var(tmp_rows[tmp_rows.size()-1]);
//                }
//            }
//        }

        lin_sys_schur().mat_set_values(ndofs, dofs_ptr, tmp_rows.size(), tmp_rows.data(), zeros);   // master edges x slave edges
        lin_sys_schur().mat_set_values(tmp_rows.size(), tmp_rows.data(), ndofs, dofs_ptr, zeros);   // slave edges  x master edges
        lin_sys_schur().mat_set_values(tmp_rows.size(), tmp_rows.data(), tmp_rows.size(), tmp_rows.data(), zeros);  // slave edges  x slave edges
    }
    // DebugOut() << "end Allocate new schur\n";
    
    // int local_dofs[10];
    // unsigned int nsides;
    // for ( DHCellAccessor dh_cell : eq_data_->dh_->own_range() ) {
    //     LocalElementAccessorBase<3> ele_ac(dh_cell);
    //     nsides = ele_ac.n_sides();
        
    //     //allocate at once matrix [sides,ele,edges]x[sides,ele,edges]
    //     loc_size = 1 + 2*nsides;
    //     unsigned int i_side = 0;
        
    //     for (; i_side < nsides; i_side++) {
    //         local_dofs[i_side] = ele_ac.side_row(i_side);
    //         local_dofs[i_side+nsides] = ele_ac.edge_row(i_side);
    //     }
    //     local_dofs[i_side+nsides] = ele_ac.ele_row();
    //     int * edge_rows = local_dofs + nsides;
    //     //int ele_row = local_dofs[0];
        
    //     // whole local MH matrix
    //     ls->mat_set_values(loc_size, local_dofs, loc_size, local_dofs, zeros);
        

    //     // compatible neighborings rows
    //     unsigned int n_neighs = ele_ac.element_accessor()->n_neighs_vb();
    //     unsigned int i=0;
    //     for ( DHCellSide neighb_side : dh_cell.neighb_sides() ) {
    //     //for (unsigned int i = 0; i < n_neighs; i++) {
    //         // every compatible connection adds a 2x2 matrix involving
    //         // current element pressure  and a connected edge pressure
    //         Neighbour *ngh = ele_ac.element_accessor()->neigh_vb[i];
    //         DHCellAccessor cell_higher_dim = eq_data_->dh_->cell_accessor_from_element(neighb_side.elem_idx());
    //         LocalElementAccessorBase<3> acc_higher_dim( cell_higher_dim );
    //         for (unsigned int j = 0; j < neighb_side.element().dim()+1; j++)
    //         	if (neighb_side.element()->edge_idx(j) == ngh->edge_idx()) {
    //         		int neigh_edge_row = acc_higher_dim.edge_row(j);
    //         		tmp_rows.push_back(neigh_edge_row);
    //         		break;
    //         	}
    //         //DebugOut() << "CC" << print_var(tmp_rows[i]);
    //         ++i;
    //     }

    //     // allocate always also for schur 2
    //     ls->mat_set_values(nsides+1, edge_rows, n_neighs, tmp_rows.data(), zeros); // (edges, ele)  x (neigh edges)
    //     ls->mat_set_values(n_neighs, tmp_rows.data(), nsides+1, edge_rows, zeros); // (neigh edges) x (edges, ele)
    //     ls->mat_set_values(n_neighs, tmp_rows.data(), n_neighs, tmp_rows.data(), zeros);  // (neigh edges) x (neigh edges)

    //     tmp_rows.clear();

    //     if (eq_data_->mortar_method_ != NoMortar) {
    //         auto &isec_list = mesh_->mixed_intersections().element_intersections_[ele_ac.ele_global_idx()];
    //         for(auto &isec : isec_list ) {
    //             IntersectionLocalBase *local = isec.second;
    //             LocalElementAccessorBase<3> slave_acc( eq_data_->dh_->cell_accessor_from_element(local->bulk_ele_idx()) );
    //             //DebugOut().fmt("Alloc: {} {}", ele_ac.ele_global_idx(), local->bulk_ele_idx());
    //             for(unsigned int i_side=0; i_side < slave_acc.dim()+1; i_side++) {
    //                 tmp_rows.push_back( slave_acc.edge_row(i_side) );
    //                 //DebugOut() << "aedge" << print_var(tmp_rows[tmp_rows.size()-1]);
    //             }
    //         }
    //     }
    //     /*
    //     for(unsigned int i_side=0; i_side < ele_ac.element_accessor()->n_sides(); i_side++) {
    //         DebugOut() << "aedge:" << print_var(edge_rows[i_side]);
    //     }*/

    //     ls->mat_set_values(nsides, edge_rows, tmp_rows.size(), tmp_rows.data(), zeros);   // master edges x neigh edges
    //     ls->mat_set_values(tmp_rows.size(), tmp_rows.data(), nsides, edge_rows, zeros);   // neigh edges  x master edges
    //     ls->mat_set_values(tmp_rows.size(), tmp_rows.data(), tmp_rows.size(), tmp_rows.data(), zeros);  // neigh edges  x neigh edges

    // }
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



/*******************************************************************************
 * COMPOSE WATER MH MATRIX WITHOUT SCHUR COMPLEMENT
 ******************************************************************************/

void DarcyLMH::create_linear_system(Input::AbstractRecord in_rec) {
  
    START_TIMER("preallocation");

    // if (schur0 == NULL) { // create Linear System for MH matrix
       
//     	if (in_rec.type() == LinSys_BDDC::get_input_type()) {
// #ifdef FLOW123D_HAVE_BDDCML
//     		WarningOut() << "For BDDC no Schur complements are used.";
//             n_schur_compls = 0;
//             LinSys_BDDC *ls = new LinSys_BDDC(&(*eq_data_->dh_->distr()),
//                     true); // swap signs of matrix and rhs to make the matrix SPD
//             ls->set_from_input(in_rec);
//             ls->set_solution( eq_data_->full_solution.petsc_vec() );
//             // possible initialization particular to BDDC
//             START_TIMER("BDDC set mesh data");
//             set_mesh_data_for_bddc(ls);
//             schur0=ls;
//             END_TIMER("BDDC set mesh data");
// #else
//             Exception
//             THROW( ExcBddcmlNotSupported() );
// #endif // FLOW123D_HAVE_BDDCML
//         } 
//         else
        if (in_rec.type() == LinSys_PETSC::get_input_type()) {
        // use PETSC for serial case even when user wants BDDC

        	eq_data_->lin_sys_schur = std::make_shared<LinSys_PETSC>( &(*eq_data_->dh_cr_->distr()) );
            lin_sys_schur().set_from_input(in_rec);
            lin_sys_schur().set_positive_definite();
            lin_sys_schur().set_solution( eq_data_->p_edge_solution.petsc_vec() );
            lin_sys_schur().set_symmetric();
            ((LinSys_PETSC *)&lin_sys_schur())->set_initial_guess_nonzero(true);
            
//             LinSys_PETSC *schur1, *schur2;

//             if (n_schur_compls == 0) {
//                 LinSys_PETSC *ls = new LinSys_PETSC( &(*eq_data_->dh_->distr()) );

//                 // temporary solution; we have to set precision also for sequantial case of BDDC
//                 // final solution should be probably call of direct solver for oneproc case
// //                 if (in_rec.type() != LinSys_BDDC::get_input_type()) ls->set_from_input(in_rec);
// //                 else {
// //                     ls->LinSys::set_from_input(in_rec); // get only common options
// //                 }
//                 ls->set_from_input(in_rec);

// //                 ls->set_solution( eq_data_->full_solution.petsc_vec() );
//                 schur0=ls;
//             } else {
//                 IS is;
//                 auto side_dofs_vec = get_component_indices_vec(0);

//                 ISCreateGeneral(PETSC_COMM_SELF, side_dofs_vec.size(), &(side_dofs_vec[0]), PETSC_COPY_VALUES, &is);
//                 //ISView(is, PETSC_VIEWER_STDOUT_SELF);
//                 //ASSERT_PERMANENT(err == 0).error("Error in ISCreateStride.");

//                 SchurComplement *ls = new SchurComplement(&(*eq_data_->dh_->distr()), is);

//                 // make schur1
//                 Distribution *ds = ls->make_complement_distribution();
//                 if (n_schur_compls==1) {
//                     schur1 = new LinSys_PETSC(ds);
//                     schur1->set_positive_definite();
//                 } else {
//                     IS is;
//                     auto elem_dofs_vec = get_component_indices_vec(1);

//                     const PetscInt *b_indices;
//                     ISGetIndices(ls->IsB, &b_indices);
//                     uint b_size = ls->loc_size_B;
//                     for(uint i_b=0, i_bb=0; i_b < b_size && i_bb < elem_dofs_vec.size(); i_b++) {
//                         if (b_indices[i_b] == elem_dofs_vec[i_bb])
//                             elem_dofs_vec[i_bb++] = i_b + ds->begin();
//                     }
//                     ISRestoreIndices(ls->IsB, &b_indices);


//                     ISCreateGeneral(PETSC_COMM_SELF, elem_dofs_vec.size(), &(elem_dofs_vec[0]), PETSC_COPY_VALUES, &is);
//                     //ISView(is, PETSC_VIEWER_STDOUT_SELF);
//                     //ASSERT_PERMANENT(err == 0).error("Error in ISCreateStride.");
//                     SchurComplement *ls1 = new SchurComplement(ds, is); // is is deallocated by SchurComplement
//                     ls1->set_negative_definite();

//                     // make schur2
//                     schur2 = new LinSys_PETSC( ls1->make_complement_distribution() );
//                     schur2->set_positive_definite();
//                     ls1->set_complement( schur2 );
//                     schur1 = ls1;
//                 }
//                 ls->set_complement( schur1 );
//                 ls->set_from_input(in_rec);
// //                 ls->set_solution( eq_data_->full_solution.petsc_vec() );
//                 schur0=ls;
            // }

            START_TIMER("PETSC PREALLOCATION");
            lin_sys_schur().start_allocation();
            
            allocate_mh_matrix();
            
            eq_data_->full_solution.zero_entries();
            eq_data_->p_edge_solution.zero_entries();
            END_TIMER("PETSC PREALLOCATION");
        }
        else {
            THROW( ExcUnknownSolver() );
        }

        END_TIMER("preallocation");
}

void DarcyLMH::postprocess()
{}

//void DarcyLMH::reconstruct_solution_from_schur(MultidimAssembly& assembler)
//{
//    START_TIMER("DarcyFlowMH::reconstruct_solution_from_schur");
//
//    eq_data_->full_solution.zero_entries();
//    eq_data_->p_edge_solution.local_to_ghost_begin();
//    eq_data_->p_edge_solution.local_to_ghost_end();
//
//    balance_->start_flux_assembly(eq_data_->water_balance_idx);
//    balance_->start_source_assembly(eq_data_->water_balance_idx);
//    balance_->start_mass_assembly(eq_data_->water_balance_idx);
//
//    for ( DHCellAccessor dh_cell : eq_data_->dh_->own_range() ) {
//        unsigned int dim = dh_cell.dim();
//        assembler[dim-1]->assemble_reconstruct(dh_cell);
//    }
//
//    eq_data_->full_solution.local_to_ghost_begin();
//    eq_data_->full_solution.local_to_ghost_end();
//
//    balance_->finish_mass_assembly(eq_data_->water_balance_idx);
//    balance_->finish_source_assembly(eq_data_->water_balance_idx);
//    balance_->finish_flux_assembly(eq_data_->water_balance_idx);
//}

void DarcyLMH::assembly_linear_system() {
    START_TIMER("DarcyFlowMH::assembly_linear_system");
//     DebugOut() << "DarcyLMH::assembly_linear_system\n";

    eq_data_->p_edge_solution.local_to_ghost_begin();
    eq_data_->p_edge_solution.local_to_ghost_end();

    eq_data_->is_linear=true;
    //DebugOut() << "Assembly linear system\n";
//  if (data_changed_) {
//      data_changed_ = false;
    {
        //DebugOut()  << "Data changed\n";
        // currently we have no optimization for cases when just time term data or RHS data are changed
        START_TIMER("full assembly");
//         if (typeid(*schur0) != typeid(LinSys_BDDC)) {
//             schur0->start_add_assembly(); // finish allocation and create matrix
//             schur_compl->start_add_assembly();
//         }
        
        lin_sys_schur().start_add_assembly();
        
        lin_sys_schur().mat_zero_entries();
        lin_sys_schur().rhs_zero_entries();
        
        eq_data_->time_step_ = time_->dt();

        START_TIMER("DarcyLMH::assembly_steady_mh_matrix");
        this->mh_matrix_assembly_->assemble(eq_data_->dh_);; // fill matrix
        END_TIMER("DarcyLMH::assembly_steady_mh_matrix");
//        assembly_mh_matrix( eq_data_->multidim_assembler ); // fill matrix

        lin_sys_schur().finish_assembly();
        lin_sys_schur().set_matrix_changed();

        // print_matlab_matrix("matrix");
    }
}


void DarcyLMH::print_matlab_matrix(std::string matlab_file)
{
    std::string output_file;
    
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
    fprintf(file, "nA = %d;\n", eq_data_->dh_cr_disc_->distr()->size());
    fprintf(file, "nB = %d;\n", eq_data_->dh_->mesh()->get_el_ds()->size());
    fprintf(file, "nBF = %d;\n", eq_data_->dh_cr_->distr()->size());
    fprintf(file, "h1 = %e;\nh2 = %e;\nh3 = %e;\n", h1, h2, h3);
    fprintf(file, "he2 = %e;\nhe3 = %e;\n", he2, he3);
    fclose(file);

    {
        output_file = FilePath(matlab_file + "_sch_new.m", FilePath::output_file);
        PetscViewer    viewer;
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, output_file.c_str(), &viewer);
        PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
        MatView( *const_cast<Mat*>(lin_sys_schur().get_matrix()), viewer);
        VecView( *const_cast<Vec*>(lin_sys_schur().get_rhs()), viewer);
        VecView( *const_cast<Vec*>(&(lin_sys_schur().get_solution())), viewer);
        VecView( *const_cast<Vec*>(&(eq_data_->full_solution.petsc_vec())), viewer);
    }
}


//template <int dim>
//std::vector<arma::vec3> dof_points(DHCellAccessor cell, const Mapping<dim, 3> &mapping) {
//
//
//    vector<arma::vec::fixed<dim+1>> bary_dof_points = cell->fe()->dof_points();
//
//    std::vector<arma::vec3> points(20);
//    points.resize(0);
//
//}
//

// void DarcyLMH::set_mesh_data_for_bddc(LinSys_BDDC * bddc_ls) {
//     START_TIMER("DarcyFlowMH_Steady::set_mesh_data_for_bddc");
//     // prepare mesh for BDDCML
//     // initialize arrays
//     // auxiliary map for creating coordinates of local dofs and global-to-local numbering
//     std::map<int, arma::vec3> localDofMap;
//     // connectivity for the subdomain, i.e. global dof numbers on element, stored element-by-element
//     // Indices of Nodes on Elements
//     std::vector<int> inet;
//     // number of degrees of freedom on elements - determines elementwise chunks of INET array
//     // Number of Nodes on Elements
//     std::vector<int> nnet;
//     // Indices of Subdomain Elements in Global Numbering - for local elements, their global indices
//     std::vector<int> isegn;
// 
//     // This array is currently not used in BDDCML, it was used as an interface scaling alternative to scaling
//     // by diagonal. It corresponds to the rho-scaling.
//     std::vector<double> element_permeability;
// 
//     // maximal and minimal dimension of elements
//     uint elDimMax = 1;
//     uint elDimMin = 3;
//     std::vector<LongIdx> cell_dofs_global(10);
// 
// 
// 
//     for ( DHCellAccessor dh_cell : eq_data_->dh_->own_range() ) {
//         // LocalElementAccessorBase<3> ele_ac(dh_cell);
//         // for each element, create local numbering of dofs as fluxes (sides), pressure (element centre), Lagrange multipliers (edges), compatible connections
// 
//         dh_cell.get_dof_indices(cell_dofs_global);
// 
//         inet.insert(inet.end(), cell_dofs_global.begin(), cell_dofs_global.end());
//         uint n_inet = cell_dofs_global.size();
// 
// 
//         uint dim = dh_cell.elm().dim();
//         elDimMax = std::max( elDimMax, dim );
//         elDimMin = std::min( elDimMin, dim );
// 
//         // TODO: this is consistent with previous implementation, but may be wrong as it use global element numbering
//         // used in sequential mesh, do global numbering of distributed elements.
//         isegn.push_back( dh_cell.elm_idx());
// 
//         // TODO: use FiniteElement::dof_points
//         for (unsigned int si=0; si<dh_cell.elm()->n_sides(); si++) {
//             arma::vec3 coord = dh_cell.elm().side(si)->centre();
//             // flux dof points
//             localDofMap.insert( std::make_pair( cell_dofs_global[si], coord ) );
//             // pressure trace dof points
//             localDofMap.insert( std::make_pair( cell_dofs_global[si+dim+2], coord ) );
//         }
//         // pressure dof points
//         arma::vec3 elm_centre = dh_cell.elm().centre();
//         localDofMap.insert( std::make_pair( cell_dofs_global[dim+1], elm_centre ) );
// 
//         // insert dofs related to compatible connections
//         //const Element *ele = dh_cell.elm().element();
//         for(DHCellSide side : dh_cell.neighb_sides()) {
//             uint neigh_dim = side.cell().elm().dim();
//             side.cell().get_dof_indices(cell_dofs_global);
//             int edge_row = cell_dofs_global[neigh_dim+2+side.side_idx()];
//             localDofMap.insert( std::make_pair( edge_row, side.centre() ) );
//             inet.push_back( edge_row );
//             n_inet++;
//         }
//         nnet.push_back(n_inet);
// 
// 
//         // version for rho scaling
//         // trace computation
//         double conduct = eq_fields_->conductivity.value( elm_centre , dh_cell.elm() );
//         auto aniso = eq_fields_->anisotropy.value( elm_centre , dh_cell.elm() );
// 
//         // compute mean on the diagonal
//         double coef = 0.;
//         for ( int i = 0; i < 3; i++) {
//             coef = coef + aniso.at(i,i);
//         }
//         // Maybe divide by cs
//         coef = conduct*coef / 3;
// 
//         ASSERT_PERMANENT_GT(coef, 0).error("Zero coefficient of hydrodynamic resistance.\n");
//         element_permeability.push_back( 1. / coef );
//     }
// //    uint i_inet = 0;
// //    for(int n_dofs : nnet) {
// //        DebugOut() << "nnet: " << n_dofs;
// //        for(int j=0; j < n_dofs; j++, i_inet++) {
// //            DebugOut() << "inet: " << inet[i_inet];
// //        }
// //    }
// 
//     auto distr = eq_data_->dh_->distr();
// //    for(auto pair : localDofMap) {
// //        DebugOut().every_proc() << "r: " << distr->myp() << " gi: " << pair.first << "xyz: " << pair.second[0];
// //
// //    }
// 
// 
//     //convert set of dofs to vectors
//     // number of nodes (= dofs) on the subdomain
//     int numNodeSub = localDofMap.size();
//     //ASSERT_PERMANENT_EQ( (unsigned int)numNodeSub, eq_data_->dh_->lsize() );
//     // Indices of Subdomain Nodes in Global Numbering - for local nodes, their global indices
//     std::vector<int> isngn( numNodeSub );
//     // pseudo-coordinates of local nodes (i.e. dofs)
//     // they need not be exact, they are used just for some geometrical considerations in BDDCML, 
//     // such as selection of corners maximizing area of a triangle, bounding boxes fro subdomains to 
//     // find candidate neighbours etc.
//     std::vector<double> xyz( numNodeSub * 3 ) ;
//     int ind = 0;
//     std::map<int,arma::vec3>::iterator itB = localDofMap.begin();
//     for ( ; itB != localDofMap.end(); ++itB ) {
//         isngn[ind] = itB -> first;
// 
//         arma::vec3 coord = itB -> second;
//         for ( int j = 0; j < 3; j++ ) {
//             xyz[ j*numNodeSub + ind ] = coord[j];
//         }
// 
//         ind++;
//     }
//     localDofMap.clear();
// 
//     // Number of Nodal Degrees of Freedom
//     // nndf is trivially one - dofs coincide with nodes
//     std::vector<int> nndf( numNodeSub, 1 );
// 
//     // prepare auxiliary map for renumbering nodes
//     typedef std::map<int,int> Global2LocalMap_; //! type for storage of global to local map
//     Global2LocalMap_ global2LocalNodeMap;
//     for ( unsigned ind = 0; ind < isngn.size(); ++ind ) {
//         global2LocalNodeMap.insert( std::make_pair( static_cast<unsigned>( isngn[ind] ), ind ) );
//     }
// 
//     // renumber nodes in the inet array to locals
//     int indInet = 0;
//     for ( unsigned int iEle = 0; iEle < isegn.size(); iEle++ ) {
//         int nne = nnet[ iEle ];
//         for ( int ien = 0; ien < nne; ien++ ) {
// 
//             int indGlob = inet[indInet];
//             // map it to local node
//             Global2LocalMap_::iterator pos = global2LocalNodeMap.find( indGlob );
//             ASSERT_PERMANENT( pos != global2LocalNodeMap.end())(indGlob).error("Cannot remap node index to local indices. \n " );
//             int indLoc = static_cast<int> ( pos -> second );
// 
//             // store the node
//             inet[ indInet++ ] = indLoc;
//         }
//     }
// 
//     int numNodes    = size;
//     int numDofsInt  = size;
//     int spaceDim    = 3;    // TODO: what is the proper value here?
//     int meshDim     = elDimMax;
// 
//     /**
//      * We need:
//      * - local to global element map (possibly mesh->el_4_loc
//      * - inet, nnet - local dof numbers per element, local numbering of only those dofs that are on owned elements
//      *   1. collect DH local dof indices  on elements, manage map from DH local indices to BDDC local dof indices
//      *   2. map collected DH indices to BDDC indices using the map
//      * - local BDDC dofs to global dofs, use DH to BDDC map with DH local to global map
//      * - XYZ - permuted, collect in main loop into array of size of all DH local dofs, compress and rearrange latter
//      * - element_permeability - in main loop
//      */
//     bddc_ls -> load_mesh( LinSys_BDDC::BDDCMatrixType::SPD_VIA_SYMMETRICGENERAL, spaceDim, numNodes, numDofsInt, inet, nnet, nndf, isegn, isngn, isngn, xyz, element_permeability, meshDim );
// }




//=============================================================================
// DESTROY WATER MH SYSTEM STRUCTURE
//=============================================================================
DarcyLMH::~DarcyLMH() {
	if (output_object)	delete output_object;

    if(time_ != nullptr)
        delete time_;
    
    if (read_init_cond_assembly_!=nullptr) {
        delete read_init_cond_assembly_;
        read_init_cond_assembly_ = nullptr;
    }
    if (mh_matrix_assembly_!=nullptr) {
        delete mh_matrix_assembly_;
        mh_matrix_assembly_ = nullptr;
    }
    if (reconstruct_schur_assembly_!=nullptr) {
        delete reconstruct_schur_assembly_;
        reconstruct_schur_assembly_ = nullptr;
    }
}


/// Helper method fills range (min and max) of given component
void dofs_range(unsigned int n_dofs, unsigned int &min, unsigned int &max, unsigned int component) {
    if (component==0) {
        min = 0;
        max = n_dofs/2;
    } else if (component==1) {
        min = n_dofs/2;
        max = (n_dofs+1)/2;
    } else {
        min = (n_dofs+1)/2;
        max = n_dofs;
    }
}


std::vector<int> DarcyLMH::get_component_indices_vec(unsigned int component) const {
	ASSERT_LT(component, 3).error("Invalid component!");
	unsigned int i, n_dofs, min, max;
    std::vector<int> dof_vec;
    std::vector<LongIdx> dof_indices(eq_data_->dh_->max_elem_dofs());
	for ( DHCellAccessor dh_cell : eq_data_->dh_->own_range() ) {
        n_dofs = dh_cell.get_dof_indices(dof_indices);
        dofs_range(n_dofs, min, max, component);
        for (i=min; i<max; ++i) dof_vec.push_back(dof_indices[i]);
    }
	return dof_vec;
}


void DarcyLMH::initialize_asm() {
    this->read_init_cond_assembly_ = new GenericAssembly< ReadInitCondAssemblyLMH >(eq_fields_.get(), eq_data_.get());
    this->mh_matrix_assembly_ = new GenericAssembly< MHMatrixAssemblyLMH >(eq_fields_.get(), eq_data_.get());
    this->reconstruct_schur_assembly_ = new GenericAssembly< ReconstructSchurAssemblyLMH >(eq_fields_.get(), eq_data_.get());
}


void DarcyLMH::read_init_cond_asm() {
    this->read_init_cond_assembly_->assemble(eq_data_->dh_cr_);
}


//-----------------------------------------------------------------------------
// vim: set cindent:
