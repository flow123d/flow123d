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
 * @file    hm_iterative.cc
 * @brief   
 * @author  Jan Stebel
 */

#include "hm_iterative.hh"
#include "system/sys_profiler.hh"
#include "input/input_type.hh"
#include "flow/richards_lmh.hh"
#include "fields/field_fe.hh"         // for create_field()


FLOW123D_FORCE_LINK_IN_CHILD(coupling_iterative)


namespace it = Input::Type;



/** Create elementwise FieldFE with parallel VectorMPI */
template <int spacedim, class Value>
std::shared_ptr<FieldFE<spacedim, Value> > create_field(Mesh & mesh, int n_comp)
{
	FiniteElement<0> *fe0; // Finite element objects (allow to create DOF handler)
	FiniteElement<1> *fe1;
	FiniteElement<2> *fe2;
	FiniteElement<3> *fe3;

	switch (n_comp) { // prepare FEM objects for DOF handler by number of components
		case 1: { // scalar
			fe0 = new FE_P_disc<0>(0);
			fe1 = new FE_P_disc<1>(0);
			fe2 = new FE_P_disc<2>(0);
			fe3 = new FE_P_disc<3>(0);
			break;
		}
		case -1: { // scalar with dof on sides
			fe0 = new FE_CR<0>;
			fe1 = new FE_CR<1>;
			fe2 = new FE_CR<2>;
			fe3 = new FE_CR<3>;
			break;
		}
		case 3: { // vector
			std::shared_ptr< FiniteElement<0> > fe0_ptr = std::make_shared< FE_P_disc<0> >(0);
			std::shared_ptr< FiniteElement<1> > fe1_ptr = std::make_shared< FE_P_disc<1> >(0);
			std::shared_ptr< FiniteElement<2> > fe2_ptr = std::make_shared< FE_P_disc<2> >(0);
			std::shared_ptr< FiniteElement<3> > fe3_ptr = std::make_shared< FE_P_disc<3> >(0);
			fe0 = new FESystem<0>(fe0_ptr, FEType::FEVector, 3);
			fe1 = new FESystem<1>(fe1_ptr, FEType::FEVector, 3);
			fe2 = new FESystem<2>(fe2_ptr, FEType::FEVector, 3);
			fe3 = new FESystem<3>(fe3_ptr, FEType::FEVector, 3);
			break;
		}
		case 9: { // tensor
			std::shared_ptr< FiniteElement<0> > fe0_ptr = std::make_shared< FE_P_disc<0> >(0);
			std::shared_ptr< FiniteElement<1> > fe1_ptr = std::make_shared< FE_P_disc<1> >(0);
			std::shared_ptr< FiniteElement<2> > fe2_ptr = std::make_shared< FE_P_disc<2> >(0);
			std::shared_ptr< FiniteElement<3> > fe3_ptr = std::make_shared< FE_P_disc<3> >(0);
			fe0 = new FESystem<0>(fe0_ptr, FEType::FETensor, 9);
			fe1 = new FESystem<1>(fe1_ptr, FEType::FETensor, 9);
			fe2 = new FESystem<2>(fe2_ptr, FEType::FETensor, 9);
			fe3 = new FESystem<3>(fe3_ptr, FEType::FETensor, 9);
			break;
		}
		default:
			ASSERT(false).error("Should not happen!\n");
	}

	// Prepare DOF handler
	std::shared_ptr<DOFHandlerMultiDim> dh_par = std::make_shared<DOFHandlerMultiDim>(mesh);
	std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( &mesh, fe0, fe1, fe2, fe3);
	dh_par->distribute_dofs(ds);

	// Construct FieldFE
	std::shared_ptr< FieldFE<spacedim, Value> > field_ptr = std::make_shared< FieldFE<spacedim, Value> >();
	field_ptr->set_fe_data( dh_par, 0, dh_par->create_vector() );
	return field_ptr;
}



const it::Record & HM_Iterative::get_input_type() {
    return it::Record("Coupling_Iterative",
            "Record with data for iterative coupling of flow and mechanics.\n")
        .derive_from( DarcyFlowInterface::get_input_type() )
		.declare_key("flow_equation", RichardsLMH::get_input_type(),
		        it::Default::obligatory(),
				"Flow equation, provides the velocity field as a result.")
		.declare_key("mechanics_equation", Elasticity::get_input_type(),
				"Mechanics, provides the displacement field.")
        .declare_key("time", TimeGovernor::get_input_type(), it::Default::obligatory(),
                    "Time governor setting for the HM coupling.")
        .declare_key("input_fields", it::Array(
		        HM_Iterative::EqData()
		            .make_field_descriptor_type("Coupling_Iterative")),
		        IT::Default::obligatory(),
		        "Input fields of the HM coupling.")
        .declare_key( "iteration_parameter", it::Double(), it::Default("1"),
                "Tuning parameter for iterative splitting. Its default value"
                "corresponds to a theoretically optimal value with fastest convergence." )
        .declare_key( "max_it", it::Integer(0), it::Default("100"),
                "Maximal count of HM iterations." )
        .declare_key( "min_it", it::Integer(0), it::Default("1"),
                "Minimal count of HM iterations." )
        .declare_key( "a_tol", it::Double(0), it::Default("0"),
                "Absolute tolerance for difference in HM iteration." )
        .declare_key( "r_tol", it::Double(0), it::Default("1e-7"),
                "Relative tolerance for difference in HM iteration." )
		.close();
}


const int HM_Iterative::registrar = Input::register_class< HM_Iterative, Mesh &, const Input::Record >("Coupling_Iterative")
                                    + HM_Iterative::get_input_type().size();


HM_Iterative::EqData::EqData()
{
    *this += alpha.name("biot_alpha")
                     .units(UnitSI().dimensionless())
                     .input_default("0.0")
                     .flags_add(FieldFlag::in_rhs);
    
    *this += density.name("fluid_density")
                     .units(UnitSI().kg().m(-3))
                     .input_default("0.0")
                     .flags_add(FieldFlag::in_rhs);
    
    *this += gravity.name("gravity")
                     .units(UnitSI().m().s(-2))
                     .input_default("9.81")
                     .flags_add(FieldFlag::in_rhs);
    
    *this += beta.name("relaxation_beta")
                     .units(UnitSI().dimensionless())
                     .flags(FieldFlag::equation_external_output);
    
    *this += pressure_potential.name("pressure_potential")
                     .units(UnitSI().m())
                     .flags(FieldFlag::equation_result);
    
    *this += flow_source.name("extra_flow_source")
                     .units(UnitSI().s(-1))
                     .flags(FieldFlag::equation_result);
}


void HM_Iterative::EqData::initialize(Mesh &mesh)
{
    // initialize coupling fields with FieldFE
    set_mesh(mesh);
    
    potential_ptr_ = create_field<3, FieldValue<3>::Scalar>(mesh, -1);
    pressure_potential.set_field(mesh.region_db().get_region_set("ALL"), potential_ptr_);
    
    beta_ptr_ = create_field<3, FieldValue<3>::Scalar>(mesh, 1);
    beta.set_field(mesh.region_db().get_region_set("ALL"), beta_ptr_);
    
    flow_source_ptr_ = create_field<3, FieldValue<3>::Scalar>(mesh, 1);
    flow_source.set_field(mesh.region_db().get_region_set("ALL"), flow_source_ptr_);
    
    old_pressure_ptr_ = create_field<3, FieldValue<3>::Scalar>(mesh, 1);
    old_iter_pressure_ptr_ = create_field<3, FieldValue<3>::Scalar>(mesh, 1);
    div_u_ptr_ = create_field<3, FieldValue<3>::Scalar>(mesh, 1);
    old_div_u_ptr_ = create_field<3, FieldValue<3>::Scalar>(mesh, 1);
}

                                    

HM_Iterative::HM_Iterative(Mesh &mesh, Input::Record in_record)
: DarcyFlowInterface(mesh, in_record)
{
	START_TIMER("HM constructor");
    using namespace Input;

    time_ = new TimeGovernor(in_record.val<Record>("time"));
    
    // setup flow equation
    Record flow_rec = in_record.val<Record>("flow_equation");
    // Need explicit template types here, since reference is used (automatically passing by value)
    flow_ = std::make_shared<RichardsLMH>(*mesh_, flow_rec, time_);
    flow_->initialize();
    std::stringstream ss; // print warning message with table of uninitialized fields
    if ( FieldCommon::print_message_table(ss, "flow") )
        WarningOut() << ss.str();
    
    // setup mechanics
    Record mech_rec = in_record.val<Record>("mechanics_equation");
    mechanics_ = std::make_shared<Elasticity>(*mesh_, mech_rec, this->time_);
    mechanics_->data()["cross_section"].copy_from(flow_->data()["cross_section"]);
    mechanics_->initialize();
    
    // read parameters controlling the iteration
    beta_ = in_record.val<double>("iteration_parameter");
    min_it_ = in_record.val<unsigned int>("min_it");
    max_it_ = in_record.val<unsigned int>("max_it");
    a_tol_ = in_record.val<double>("a_tol");
    r_tol_ = in_record.val<double>("r_tol");

    this->eq_data_ = &data_;
    
    // setup input fields
    data_.set_input_list( in_record.val<Input::Array>("input_fields"), time() );

    data_.initialize(*mesh_);
    mechanics_->set_potential_load(data_.pressure_potential);
}


void HM_Iterative::initialize()
{
}


template<int dim, class Value>
void copy_field(const Field<dim, Value> &from_field, FieldFE<dim, Value> &to_field)
{
    auto dh = to_field.get_dofhandler();
    auto vec = to_field.get_data_vec();
    
    for ( auto cell : dh->own_range() )
        vec[cell.local_idx()] = from_field.value(cell.elm().centre(), cell.elm());
    
//     vec.local_to_ghost_begin();
//     vec.local_to_ghost_end();
}


template<int dim, class Value>
void update_field_from_mh_dofhandler(const MH_DofHandler &mh_dh, FieldFE<dim, Value> &field)
{
    auto dh = field.get_dofhandler();
    auto vec = field.get_data_vec();
    
    for ( auto cell : dh->own_range() )
    {
        auto elm = cell.elm();
        vec[cell.local_idx()] = mh_dh.element_scalar(elm);
    }
}




void HM_Iterative::zero_time_step()
{
    data_.set_time(time_->step(), LimitSide::right);
    std::stringstream ss;
    if ( FieldCommon::print_message_table(ss, "coupling_iterative") )
        WarningOut() << ss.str();
    
    flow_->zero_time_step();
    update_potential();
    mechanics_->zero_time_step();
    
    update_field_from_mh_dofhandler(flow_->get_mh_dofhandler(), *data_.old_pressure_ptr_);
    update_field_from_mh_dofhandler(flow_->get_mh_dofhandler(), *data_.old_iter_pressure_ptr_);
    copy_field(mechanics_->data().output_divergence, *data_.div_u_ptr_);
}


void HM_Iterative::update_solution()
{
    unsigned it = 0;
    double difference = 0;
    double norm = 1;
    
    time_->next_time();
    time_->view("HM");
    data_.set_time(time_->step(), LimitSide::right);
    
    while (it < min_it_ || 
           (it < max_it_ && difference > a_tol_ && difference/norm > r_tol_)
          )
    {
        it++;

        // pass displacement (divergence) to flow
        // and solve flow problem
        update_flow_fields();
        flow_->solve_time_step(false);
        
        // pass pressure to mechanics and solve mechanics
        update_potential();
        mechanics_->solve_linear_system();
        mechanics_->output_vector_gather();
        
        // update displacement divergence
        copy_field(mechanics_->data().output_divergence, *data_.div_u_ptr_);
        
        // TODO: compute difference of iterates
        compute_iteration_error(difference, norm);
        
        MessageOut().fmt("HM Iteration {} abs. difference: {}  rel. difference: {}\n"
                         "--------------------------------------------------------",
                         it, difference, difference/norm);
        update_field_from_mh_dofhandler(flow_->get_mh_dofhandler(), *data_.old_iter_pressure_ptr_);
    }
    
    flow_->output_data();
    mechanics_->output_data();
    
    update_field_from_mh_dofhandler(flow_->get_mh_dofhandler(), *data_.old_pressure_ptr_);
    copy_field(mechanics_->data().output_divergence, *data_.old_div_u_ptr_);
}


void HM_Iterative::update_potential()
{
    auto potential_vec_ = data_.potential_ptr_->get_data_vec();
    auto dh = data_.potential_ptr_->get_dofhandler();
    double difference2 = 0, norm2 = 0;
    std::vector<int> dof_indices(dh->max_elem_dofs());
    for ( auto ele : dh->local_range() )
    {
        auto elm = ele.elm();
        ele.get_loc_dof_indices(dof_indices);
        for ( auto side : ele.side_range() )
        {
            double alpha = data_.alpha.value(side.centre(), elm);
            double density = data_.density.value(side.centre(), elm);
            double gravity = data_.gravity.value(side.centre(), elm);
            double pressure = flow_->get_mh_dofhandler().side_scalar(side.side());
            double potential = -alpha*density*gravity*pressure;
        
            if (ele.is_own())
            {
                difference2 += pow(potential_vec_[dof_indices[side.side_idx()]] - potential,2);
                norm2 += pow(potential,2);
            }
        
            potential_vec_[dof_indices[side.side_idx()]] = potential;
        }
    }
    
    double send_data[] = { difference2, norm2 };
    double recv_data[2];
    MPI_Allreduce(&send_data, &recv_data, 2, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    difference2 = recv_data[0];
    norm2       = recv_data[1];

    double dif2norm;
    if (norm2 == 0)
        dif2norm = (difference2 == 0)?0:numeric_limits<double>::max();
    else 
        dif2norm = sqrt(difference2/norm2);
    DebugOut() << "Relative potential difference: " << dif2norm << endl;
    
    if (dif2norm > numeric_limits<double>::epsilon())
    {
        data_.pressure_potential.set_time_result_changed();
        mechanics_->set_potential_load(data_.pressure_potential);
    }
}


void HM_Iterative::update_flow_fields()
{
    auto beta_vec = data_.beta_ptr_->get_data_vec();
    auto src_vec = data_.flow_source_ptr_->get_data_vec();
    auto dh = data_.beta_ptr_->get_dofhandler();
    double beta_diff2 = 0, beta_norm2 = 0, src_diff2 = 0, src_norm2 = 0;
    for ( auto ele : dh->local_range() )
    {
        auto elm = ele.elm();
        
        double alpha = data_.alpha.value(elm.centre(), elm);
        double young = mechanics_->data().young_modulus.value(elm.centre(), elm);
        double poisson = mechanics_->data().poisson_ratio.value(elm.centre(), elm);
        double beta = beta_ * 0.5*alpha*alpha/(2*lame_mu(young, poisson)/elm.dim() + lame_lambda(young, poisson));
        
        double old_p = data_.old_pressure_ptr_->value(elm.centre(), elm);
        double p = flow_->get_mh_dofhandler().element_scalar(elm);
        double div_u = data_.div_u_ptr_->value(elm.centre(), elm);
        double old_div_u = data_.old_div_u_ptr_->value(elm.centre(), elm);
        double src = (beta*(p-old_p) + alpha*(old_div_u - div_u)) / time_->dt();
        
        if (ele.is_own())
        {
            beta_diff2 += pow(beta_vec[ele.local_idx()] - beta,2);
            beta_norm2 += pow(beta,2);
            src_diff2 += pow(src_vec[ele.local_idx()] - src,2);
            src_norm2 += pow(src,2);
        }
        
        beta_vec[ele.local_idx()] = beta;
        src_vec[ele.local_idx()] = src;
    }
    
    double send_data[] = { beta_diff2, beta_norm2, src_diff2, src_norm2 };
    double recv_data[4];
    MPI_Allreduce(&send_data, &recv_data, 4, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    beta_diff2 = recv_data[0];
    beta_norm2 = recv_data[1];
    src_diff2  = recv_data[2];
    src_norm2  = recv_data[3];
    
    double beta_dif2norm, src_dif2norm;
    if (beta_norm2 == 0)
        beta_dif2norm = (beta_diff2 == 0)?0:numeric_limits<double>::max();
    else 
        beta_dif2norm = sqrt(beta_diff2/beta_norm2);
    if (src_norm2 == 0)
        src_dif2norm = (src_diff2 == 0)?0:numeric_limits<double>::max();
    else 
        src_dif2norm = sqrt(src_diff2/src_norm2);
    DebugOut() << "Relative difference in beta: " << beta_dif2norm << endl;
    DebugOut() << "Relative difference in extra_source: " << src_dif2norm << endl;
    
    if (beta_dif2norm > numeric_limits<double>::epsilon())
    {
        data_.beta.set_time_result_changed();
        flow_->set_extra_storativity(data_.beta);
    }
    if (src_dif2norm > numeric_limits<double>::epsilon())
    {
        data_.flow_source.set_time_result_changed();
        flow_->set_extra_source(data_.flow_source);
    }
}


void HM_Iterative::compute_iteration_error(double& difference, double& norm)
{
    auto dh = data_.beta_ptr_->get_dofhandler();
    double p_dif2 = 0, p_norm2 = 0;
    for (auto cell : dh->own_range())
    {
        auto elm = cell.elm();
        double new_p = flow_->get_mh_dofhandler().element_scalar(elm);
        double old_p = data_.old_iter_pressure_ptr_->value(elm.centre(), elm);
        p_dif2 += pow(new_p - old_p, 2)*elm.measure();
        p_norm2 += pow(old_p, 2)*elm.measure();
    }
    
    double send_data[] = { p_dif2, p_norm2 };
    double recv_data[2];
    MPI_Allreduce(&send_data, &recv_data, 2, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    difference = sqrt(recv_data[0]);
    norm = sqrt(recv_data[1]);
}



const MH_DofHandler & HM_Iterative::get_mh_dofhandler()
{ 
    return flow_->get_mh_dofhandler(); 
}


HM_Iterative::~HM_Iterative() {
	flow_.reset();
    mechanics_.reset();
}



