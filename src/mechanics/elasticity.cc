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
#include "mechanics/elasticity.hh"

#include "io/output_time.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/fe_values.hh"
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_system.hh"
#include "fields/field_fe.hh"
#include "la/linsys_PETSC.hh"
#include "coupling/balance.hh"
#include "mesh/neighbours.h"

#include "fields/multi_field.hh"
#include "fields/generic_field.hh"
#include "input/factory.hh"




using namespace Input::Type;


const Record & Elasticity::get_input_type() {
    std::string equation_name = std::string(Elasticity::EqData::name()) + "_FE";
	return IT::Record(
                std::string(equation_name),
                "FEM for linear elasticity.")
           .copy_keys(EquationBase::record_template())
           .declare_key("balance", Balance::get_input_type(), Default("{}"),
                    "Settings for computing balance.")
           .declare_key("output_stream", OutputTime::get_input_type(), Default::obligatory(),
                    "Parameters of output stream.")
           .declare_key("solver", LinSys_PETSC::get_input_type(), Default::obligatory(),
				"Linear solver for elasticity.")
		   .declare_key("input_fields", Array(
		        Elasticity::EqData()
		            .make_field_descriptor_type(equation_name)),
		        IT::Default::obligatory(),
		        "Input fields of the equation.")
           .declare_key("output",
                EqData().output_fields.make_output_type(equation_name, ""),
                IT::Default("{ \"fields\": [ "+Elasticity::EqData::default_output_field()+" ] }"),
                "Setting of the field output.")
		   .close();
}

const int Elasticity::registrar =
		Input::register_class< Elasticity, Mesh &, const Input::Record>(std::string(Elasticity::EqData::name()) + "_FE") +
		Elasticity::get_input_type().size();



namespace Mechanics {

FEObjects::FEObjects(Mesh *mesh_, unsigned int fe_order)
: q_(QGauss::make_array(2))
{
	switch (fe_order)
	{
	case 1: {
        MixedPtr<FE_P> fe_p(1);
        fe_ = mixed_fe_system(fe_p, FEVector, 3);
		break;
    }
	default:
		xprintf(PrgErr, "Unsupported polynomial order %d for finite elements in Elasticity", fe_order);
		break;
	}


    ds_ = std::make_shared<EqualOrderDiscreteSpace>(mesh_, fe_);
	dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);

	dh_->distribute_dofs(ds_);
    
    
    MixedPtr<FE_P_disc> fe_p(0);
    dh_scalar_ = make_shared<DOFHandlerMultiDim>(*mesh_);
	std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh_, fe_p);
	dh_scalar_->distribute_dofs(ds);
    

    MixedPtr<FiniteElement> fe_t = mixed_fe_system(MixedPtr<FE_P_disc>(0), FEType::FETensor, 9);
    dh_tensor_ = make_shared<DOFHandlerMultiDim>(*mesh_);
	std::shared_ptr<DiscreteSpace> dst = std::make_shared<EqualOrderDiscreteSpace>( mesh_, fe_t);
	dh_tensor_->distribute_dofs(dst);
}


FEObjects::~FEObjects()
{
}

template<> std::shared_ptr<FiniteElement<0>> FEObjects::fe<0>() { return fe_[0_d]; }
template<> std::shared_ptr<FiniteElement<1>> FEObjects::fe<1>() { return fe_[1_d]; }
template<> std::shared_ptr<FiniteElement<2>> FEObjects::fe<2>() { return fe_[2_d]; }
template<> std::shared_ptr<FiniteElement<3>> FEObjects::fe<3>() { return fe_[3_d]; }

std::shared_ptr<DOFHandlerMultiDim> FEObjects::dh() { return dh_; }
std::shared_ptr<DOFHandlerMultiDim> FEObjects::dh_scalar() { return dh_scalar_; }
std::shared_ptr<DOFHandlerMultiDim> FEObjects::dh_tensor() { return dh_tensor_; }



} // namespace Mechanics


double lame_mu(double young, double poisson)
{
    return young*0.5/(poisson+1.);
}


double lame_lambda(double young, double poisson)
{
    return young*poisson/((poisson+1.)*(1.-2.*poisson));
}






const Selection & Elasticity::EqData::get_bc_type_selection() {
    return Selection("Elasticity_BC_Type", "Types of boundary conditions for heat transfer model.")
            .add_value(bc_type_displacement, "displacement",
                  "Prescribed displacement.")
            .add_value(bc_type_displacement_normal, "displacement_n",
                  "Prescribed displacement in the normal direction to the boundary.")
            .add_value(bc_type_traction, "traction",
                  "Prescribed traction.")
            .close();
}


Elasticity::EqData::EqData()
{
    *this+=bc_type
        .name("bc_type")
        .description(
        "Type of boundary condition.")
        .units( UnitSI::dimensionless() )
        .input_default("\"traction\"")
        .input_selection( get_bc_type_selection() )
        .flags_add(FieldFlag::in_rhs & FieldFlag::in_main_matrix);
        
    *this+=bc_displacement
        .name("bc_displacement")
        .description("Prescribed displacement.")
        .units( UnitSI().m() )
        .input_default("0.0")
        .flags_add(in_rhs);
        
    *this+=bc_traction
        .name("bc_traction")
        .description("Prescribed traction.")
        .units( UnitSI().Pa() )
        .input_default("0.0")
        .flags_add(in_rhs);
        
    *this+=load
        .name("load")
        .description("Prescribed load.")
        .units( UnitSI().N().m(-3) )
        .input_default("0.0")
        .flags_add(in_rhs);
    
    *this+=young_modulus
        .name("young_modulus")
        .description("Young modulus.")
        .units( UnitSI().Pa() )
         .input_default("0.0")
        .flags_add(in_main_matrix & in_rhs);
    
    *this+=poisson_ratio
        .name("poisson_ratio")
        .description("Poisson ratio.")
        .units( UnitSI().dimensionless() )
         .input_default("0.0")
        .flags_add(in_main_matrix & in_rhs);
        
    *this+=fracture_sigma
            .name("fracture_sigma")
            .description(
            "Coefficient of diffusive transfer through fractures (for each substance).")
            .units( UnitSI::dimensionless() )
            .input_default("1.0")
            .flags_add(in_main_matrix & in_rhs);

    *this += region_id.name("region_id")
    	        .units( UnitSI::dimensionless())
    	        .flags(FieldFlag::equation_external_output);
                
    *this += subdomain.name("subdomain")
      .units( UnitSI::dimensionless() )
      .flags(FieldFlag::equation_external_output);
      
    *this+=cross_section
      .name("cross_section")
      .units( UnitSI().m(3).md() )
      .flags(input_copy & in_time_term & in_main_matrix & in_rhs);
      
    *this+=potential_load
      .name("potential_load")
      .units( UnitSI().m() )
      .flags(input_copy & in_rhs);

    *this+=output_field
            .name("displacement")
            .units( UnitSI().m() )
            .flags(equation_result);
    
    *this += output_stress
            .name("stress")
            .units( UnitSI().Pa() )
            .flags(equation_result);
    
    *this += output_von_mises_stress
            .name("von_mises_stress")
            .units( UnitSI().Pa() )
            .flags(equation_result);
    
    *this += output_cross_section
            .name("cross_section_updated")
            .units( UnitSI().m() )
            .flags(equation_result);
            
    *this += output_divergence
            .name("displacement_divergence")
            .units( UnitSI().dimensionless() )
            .flags(equation_result);

    // add all input fields to the output list
    output_fields += *this;

}

Elasticity::Elasticity(Mesh & init_mesh, const Input::Record in_rec, TimeGovernor *tm)
        : EquationBase(init_mesh, in_rec),
		  input_rec(in_rec),
		  allocation_done(false)
{
	// Can not use name() + "constructor" here, since START_TIMER only accepts const char *
	// due to constexpr optimization.
	START_TIMER(EqData::name());

	this->eq_data_ = &data_;
    
    auto time_rec = in_rec.val<Input::Record>("time");
    if (tm == nullptr)
    {
        time_ = new TimeGovernor(time_rec);
    }
    else
    {
        TimeGovernor time_from_rec(time_rec);
        ASSERT( time_from_rec.is_default() ).error("Duplicate key 'time', time in elasticity is already initialized from parent class!");
        time_ = tm;
    }


    // Set up physical parameters.
    data_.set_mesh(init_mesh);
    data_.region_id = GenericField<3>::region_id(*mesh_);
    data_.subdomain = GenericField<3>::subdomain(*mesh_);
    
    // create finite element structures and distribute DOFs
    feo = new Mechanics::FEObjects(mesh_, 1);
    DebugOut().fmt("Mechanics: solution size {}\n", feo->dh()->n_global_dofs());
    
}


void Elasticity::initialize()
{
    output_stream_ = OutputTime::create_output_stream("mechanics", input_rec.val<Input::Record>("output_stream"), time().get_unit_string());
    
    data_.set_components({"displacement"});
    data_.set_input_list( input_rec.val<Input::Array>("input_fields"), time() );
    
//     balance_ = std::make_shared<Balance>("mechanics", mesh_);
//     balance_->init_from_input(input_rec.val<Input::Record>("balance"), *time_);
    // initialization of balance object
//     data_.balance_idx_ = {
//         balance_->add_quantity("force_x"),
//         balance_->add_quantity("force_y"),
//         balance_->add_quantity("force_z")
//     };
//     balance_->units(UnitSI().kg().m().s(-2));

    // create shared pointer to a FieldFE, pass FE data and push this FieldFE to output_field on all regions
    data_.output_field_ptr = create_field_fe<3, FieldValue<3>::VectorFixed>(feo->dh());
    data_.output_field.set_field(mesh_->region_db().get_region_set("ALL"), data_.output_field_ptr, 0.);
    
    // setup output stress
    data_.output_stress_ptr = create_field_fe<3, FieldValue<3>::TensorFixed>(feo->dh_tensor());
    data_.output_stress.set_field(mesh_->region_db().get_region_set("ALL"), data_.output_stress_ptr);
    
    // setup output von Mises stress
    data_.output_von_mises_stress_ptr = create_field_fe<3, FieldValue<3>::Scalar>(feo->dh_scalar());
    data_.output_von_mises_stress.set_field(mesh_->region_db().get_region_set("ALL"), data_.output_von_mises_stress_ptr);
    
    // setup output cross-section
    data_.output_cross_section_ptr = create_field_fe<3, FieldValue<3>::Scalar>(feo->dh_scalar());
    data_.output_cross_section.set_field(mesh_->region_db().get_region_set("ALL"), data_.output_cross_section_ptr);
    
    // setup output divergence
    data_.output_div_ptr = create_field_fe<3, FieldValue<3>::Scalar>(feo->dh_scalar());
    data_.output_divergence.set_field(mesh_->region_db().get_region_set("ALL"), data_.output_div_ptr);
    
    data_.output_field.output_type(OutputTime::CORNER_DATA);

    // set time marks for writing the output
    data_.output_fields.initialize(output_stream_, mesh_, input_rec.val<Input::Record>("output"), this->time());

    // equation default PETSc solver options
    std::string petsc_default_opts;
    petsc_default_opts = "-ksp_type cg -pc_type hypre -pc_hypre_type boomeramg";
    
    // allocate matrix and vector structures
    ls = new LinSys_PETSC(feo->dh()->distr().get(), petsc_default_opts);
    ( (LinSys_PETSC *)ls )->set_from_input( input_rec.val<Input::Record>("solver") );
    ls->set_solution(data_.output_field_ptr->get_data_vec().petsc_vec());

    // initialization of balance object
//     balance_->allocate(feo->dh()->distr()->lsize(),
//             max(feo->fe<1>()->n_dofs(), max(feo->fe<2>()->n_dofs(), feo->fe<3>()->n_dofs())));

}


Elasticity::~Elasticity()
{
//     delete time_;

    MatDestroy(&stiffness_matrix);
    VecDestroy(&rhs);
    delete ls;
    delete feo;

}



void Elasticity::update_output_fields()
{
    // update ghost values of solution vector
    data_.output_field_ptr->get_data_vec().local_to_ghost_begin();
    data_.output_field_ptr->get_data_vec().local_to_ghost_end();

    // compute new output fields depending on solution (stress, divergence etc.)
    data_.output_stress_ptr->get_data_vec().zero_entries();
    data_.output_cross_section_ptr->get_data_vec().zero_entries();
    data_.output_div_ptr->get_data_vec().zero_entries();
    compute_output_fields<1>();
    compute_output_fields<2>();
    compute_output_fields<3>();

    // update ghost values of computed fields
    data_.output_stress_ptr->get_data_vec().local_to_ghost_begin();
    data_.output_stress_ptr->get_data_vec().local_to_ghost_end();
    data_.output_von_mises_stress_ptr->get_data_vec().local_to_ghost_begin();
    data_.output_von_mises_stress_ptr->get_data_vec().local_to_ghost_end();
    data_.output_cross_section_ptr->get_data_vec().local_to_ghost_begin();
    data_.output_cross_section_ptr->get_data_vec().local_to_ghost_end();
    data_.output_div_ptr->get_data_vec().local_to_ghost_begin();
    data_.output_div_ptr->get_data_vec().local_to_ghost_end();
}




void Elasticity::zero_time_step()
{
	START_TIMER(EqData::name());
	data_.mark_input_times( *time_ );
	data_.set_time(time_->step(), LimitSide::right);
	std::stringstream ss; // print warning message with table of uninitialized fields
	if ( FieldCommon::print_message_table(ss, "mechanics") ) {
		WarningOut() << ss.str();
	}

    ( (LinSys_PETSC *)ls )->set_initial_guess_nonzero();

    // check first time assembly - needs preallocation
    if (!allocation_done) preallocate();
    

    // after preallocation we assemble the matrices and vectors required for balance of forces
//     for (auto subst_idx : data_.balance_idx_)
//         balance_->calculate_instant(subst_idx, ls->get_solution());
    
//     update_solution();
    ls->start_add_assembly();
    MatSetOption(*ls->get_matrix(), MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    ls->mat_zero_entries();
    ls->rhs_zero_entries();
    assemble_stiffness_matrix();
    assemble_rhs();
    ls->finish_assembly();
    ls->solve();
    output_data();
}



void Elasticity::preallocate()
{
    // preallocate system matrix
    ls->start_allocation();
    stiffness_matrix = NULL;
    rhs = NULL;

	assemble_stiffness_matrix();
    assemble_rhs();

	allocation_done = true;
}




void Elasticity::update_solution()
{
	START_TIMER("DG-ONE STEP");

    next_time();
	solve_linear_system();

    calculate_cumulative_balance();
    
    output_data();

    END_TIMER("DG-ONE STEP");
}

void Elasticity::next_time()
{
    time_->next_time();
    time_->view("MECH");
    
}



void Elasticity::solve_linear_system()
{
    START_TIMER("data reinit");
    data_.set_time(time_->step(), LimitSide::right);
    END_TIMER("data reinit");
    
    // assemble stiffness matrix
    if (stiffness_matrix == NULL
        || data_.subset(FieldFlag::in_main_matrix).changed())
    {
        DebugOut() << "Mechanics: Assembling matrix.\n";
        ls->start_add_assembly();
        ls->mat_zero_entries();
        assemble_stiffness_matrix();
        ls->finish_assembly();

        if (stiffness_matrix == NULL)
            MatConvert(*( ls->get_matrix() ), MATSAME, MAT_INITIAL_MATRIX, &stiffness_matrix);
        else
            MatCopy(*( ls->get_matrix() ), stiffness_matrix, DIFFERENT_NONZERO_PATTERN);
    }

    // assemble right hand side (due to sources and boundary conditions)
    if (rhs == NULL
        || data_.subset(FieldFlag::in_rhs).changed())
    {
        DebugOut() << "Mechanics: Assembling right hand side.\n";
        ls->start_add_assembly();
        ls->rhs_zero_entries();
    	assemble_rhs();
        ls->finish_assembly();

        if (rhs == nullptr) VecDuplicate(*( ls->get_rhs() ), &rhs);
        VecCopy(*( ls->get_rhs() ), rhs);
    }

    START_TIMER("solve");
    ls->solve();
    END_TIMER("solve");
}






void Elasticity::output_data()
{
    START_TIMER("MECH-OUTPUT");

    // gather the solution from all processors
    data_.output_fields.set_time( this->time().step(), LimitSide::left);
    //if (data_.output_fields.is_field_output_time(data_.output_field, this->time().step()) )
    update_output_fields();
    data_.output_fields.output(this->time().step());
    output_stream_->write_time_frame();

//     START_TIMER("MECH-balance");
//     balance_->calculate_instant(subst_idx, ls->get_solution());
//     balance_->output();
//     END_TIMER("MECH-balance");

    END_TIMER("MECH-OUTPUT");
}


template<unsigned int dim>
void Elasticity::compute_output_fields()
{
    QGauss q(dim, 0), q_sub(dim-1, 0);
    FEValues<3> fv(q, *feo->fe<dim>(),
    		update_values | update_gradients | update_quadrature_points);
    FEValues<3> fsv(q_sub, *feo->fe<dim>(),
    		update_values | update_normal_vectors | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs();
    auto vec = fv.vector_view(0);
    auto vec_side = fsv.vector_view(0);
    VectorMPI output_vec = data_.output_field_ptr->get_data_vec();
    VectorMPI output_stress_vec = data_.output_stress_ptr->get_data_vec();
    VectorMPI output_von_mises_stress_vec = data_.output_von_mises_stress_ptr->get_data_vec();
    VectorMPI output_cross_sec_vec = data_.output_cross_section_ptr->get_data_vec();
    VectorMPI output_div_vec = data_.output_div_ptr->get_data_vec();
    
    DHCellAccessor cell_tensor = *feo->dh_tensor()->own_range().begin();
    DHCellAccessor cell_scalar = *feo->dh_scalar()->own_range().begin();
    for (auto cell : feo->dh()->own_range())
    {
        if (cell.dim() == dim)
        {
            auto elm = cell.elm();
            
            double poisson = data_.poisson_ratio.value(elm.centre(), elm);
            double young = data_.young_modulus.value(elm.centre(), elm);
            double mu = lame_mu(young, poisson);
            double lambda = lame_lambda(young, poisson);
            
            fv.reinit(elm);
            LocDofVec dof_indices        = cell.get_loc_dof_indices();
            LocDofVec dof_indices_scalar = cell_scalar.get_loc_dof_indices();
            LocDofVec dof_indices_tensor = cell_tensor.get_loc_dof_indices();
            
            arma::mat33 stress = arma::zeros(3,3);
            double div = 0;
            for (unsigned int i=0; i<ndofs; i++)
            {
                stress += (2*mu*vec.sym_grad(i,0) + lambda*vec.divergence(i,0)*arma::eye(3,3))*output_vec[dof_indices[i]];
                div += vec.divergence(i,0)*output_vec[dof_indices[i]];
            }
            
            arma::mat33 stress_dev = stress - arma::trace(stress)/3*arma::eye(3,3);
            double von_mises_stress = sqrt(1.5*arma::dot(stress_dev, stress_dev));
            output_div_vec[dof_indices_scalar[0]] += div;
            
            for (unsigned int i=0; i<3; i++)
                for (unsigned int j=0; j<3; j++)
                    output_stress_vec[dof_indices_tensor[i*3+j]] += stress(i,j);
            output_von_mises_stress_vec[dof_indices_scalar[0]] = von_mises_stress;
            
            output_cross_sec_vec[dof_indices_scalar[0]] += data_.cross_section.value(fv.point(0), elm);
        } 
        else if (cell.dim() == dim-1)
        {
            auto elm = cell.elm();
            double normal_displacement = 0;
            double csection = data_.cross_section.value(fsv.point(0), elm);
            arma::mat33 normal_stress = arma::zeros(3,3);

            double poisson = data_.poisson_ratio.value(elm.centre(), elm);
            double young = data_.young_modulus.value(elm.centre(), elm);
            double mu = lame_mu(young, poisson);
            double lambda = lame_lambda(young, poisson);

            for (unsigned int inb=0; inb<elm->n_neighs_vb(); inb++)
            {
                auto side = elm->neigh_vb[inb]->side();
                auto cell_side = side->element();
                fsv.reinit(*side);
                LocDofVec side_dof_indices =
                    feo->dh()->cell_accessor_from_element(cell_side.idx()).get_loc_dof_indices();
                
                for (unsigned int i=0; i<ndofs; i++)
                {
                    normal_displacement -= arma::dot(vec_side.value(i,0)*output_vec[side_dof_indices[i]], fsv.normal_vector(0));
                    arma::mat33 grad = -arma::kron(vec_side.value(i,0)*output_vec[side_dof_indices[i]], fsv.normal_vector(0).t()) / csection;
                    normal_stress += mu*(grad+grad.t()) + lambda*arma::trace(grad)*arma::eye(3,3);
                }
            }
            LocDofVec dof_indices_scalar = cell_scalar.get_loc_dof_indices();
            LocDofVec dof_indices_tensor = cell_tensor.get_loc_dof_indices();
            for (unsigned int i=0; i<3; i++)
                for (unsigned int j=0; j<3; j++)
                    output_stress_vec[dof_indices_tensor[i*3+j]] += normal_stress(i,j);
            output_cross_sec_vec[dof_indices_scalar[0]] += normal_displacement;
            output_div_vec[dof_indices_scalar[0]] += normal_displacement / csection;
        }
        cell_scalar.inc();
        cell_tensor.inc();
    }
}





void Elasticity::calculate_cumulative_balance()
{
//     if (balance_->cumulative())
//     {
//         balance_->calculate_cumulative(subst_idx, ls->get_solution());
//     }
}








void Elasticity::assemble_stiffness_matrix()
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

   START_TIMER("assemble_matrix_elem_side");
    assemble_matrix_element_side<1>();
    assemble_matrix_element_side<2>();
    assemble_matrix_element_side<3>();
   END_TIMER("assemble_matrix_elem_side");
  END_TIMER("assemble_stiffness");
}



template<unsigned int dim>
void Elasticity::assemble_volume_integrals()
{
    FEValues<3> fe_values(*feo->q<dim>(), *feo->fe<dim>(),
    		update_values | update_gradients | update_JxW_values | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
    vector<int> dof_indices(ndofs);
    vector<arma::vec3> velocity(qsize);
    vector<double> young(qsize), poisson(qsize), csection(qsize);
    PetscScalar local_matrix[ndofs*ndofs];
    auto vec = fe_values.vector_view(0);

	// assemble integral over elements
    for (auto cell : feo->dh()->own_range())
    {
        if (cell.dim() != dim) continue;
        ElementAccessor<3> elm_acc = cell.elm();

        fe_values.reinit(elm_acc);
        cell.get_dof_indices(dof_indices);

        data_.cross_section.value_list(fe_values.point_list(), elm_acc, csection);
        data_.young_modulus.value_list(fe_values.point_list(), elm_acc, young);
        data_.poisson_ratio.value_list(fe_values.point_list(), elm_acc, poisson);
        
        // assemble the local stiffness matrix
        for (unsigned int i=0; i<ndofs; i++)
            for (unsigned int j=0; j<ndofs; j++)
                local_matrix[i*ndofs+j] = 0;

        for (unsigned int k=0; k<qsize; k++)
        {
            double mu = lame_mu(young[k], poisson[k]);
            double lambda = lame_lambda(young[k], poisson[k]);
            for (unsigned int i=0; i<ndofs; i++)
            {
                for (unsigned int j=0; j<ndofs; j++)
                    local_matrix[i*ndofs+j] += csection[k]*(
                                                2*mu*arma::dot(vec.sym_grad(j,k), vec.sym_grad(i,k))
                                                + lambda*vec.divergence(j,k)*vec.divergence(i,k)
                                               )*fe_values.JxW(k);
            }
        }
        ls->mat_set_values(ndofs, dof_indices.data(), ndofs, dof_indices.data(), local_matrix);
    }
}



void Elasticity::assemble_rhs()
{
  START_TIMER("assemble_rhs");
//     balance_->start_source_assembly(subst_idx);
//     balance_->start_flux_assembly(subst_idx);
	assemble_sources<1>();
	assemble_sources<2>();
	assemble_sources<3>();
    assemble_rhs_element_side<1>();
    assemble_rhs_element_side<2>();
    assemble_rhs_element_side<3>();
    assemble_boundary_conditions<1>();
    assemble_boundary_conditions<2>();
    assemble_boundary_conditions<3>();
// 	balance_->finish_flux_assembly(subst_idx);
// 	balance_->finish_source_assembly(subst_idx);
  END_TIMER("assemble_rhs");
}


template<unsigned int dim>
void Elasticity::assemble_sources()
{
    FEValues<3> fe_values(*feo->q<dim>(), *feo->fe<dim>(),
    		update_values | update_gradients | update_JxW_values | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
    vector<arma::vec3> load(qsize);
    vector<double> csection(qsize), potential(qsize);
    vector<int> dof_indices(ndofs);
    PetscScalar local_rhs[ndofs];
    vector<PetscScalar> local_source_balance_vector(ndofs), local_source_balance_rhs(ndofs);
    auto vec = fe_values.vector_view(0);

	// assemble integral over elements
    for (auto cell : feo->dh()->own_range())
    {
        if (cell.dim() != dim) continue;
        ElementAccessor<3> elm_acc = cell.elm();

        fe_values.reinit(elm_acc);
        cell.get_dof_indices(dof_indices);

        // assemble the local stiffness matrix
        fill_n(local_rhs, ndofs, 0);
        local_source_balance_vector.assign(ndofs, 0);
        local_source_balance_rhs.assign(ndofs, 0);
        
        data_.cross_section.value_list(fe_values.point_list(), elm_acc, csection);
        data_.load.value_list(fe_values.point_list(), elm_acc, load);
        data_.potential_load.value_list(fe_values.point_list(), elm_acc, potential);

        // compute sources
        for (unsigned int k=0; k<qsize; k++)
        {
            for (unsigned int i=0; i<ndofs; i++)
                local_rhs[i] += (
                                 arma::dot(load[k], vec.value(i,k))
                                 -potential[k]*vec.divergence(i,k)
                                )*csection[k]*fe_values.JxW(k);
        }
        ls->rhs_set_values(ndofs, dof_indices.data(), local_rhs);

//         for (unsigned int i=0; i<ndofs; i++)
//         {
//             for (unsigned int k=0; k<qsize; k++)
//                 local_source_balance_vector[i] -= 0;//sources_sigma[k]*fe_values[vec].value(i,k)*fe_values.JxW(k);
// 
//             local_source_balance_rhs[i] += local_rhs[i];
//         }
//         balance_->add_source_matrix_values(subst_idx, elm_acc.region().bulk_idx(), dof_indices, local_source_balance_vector);
//         balance_->add_source_vec_values(subst_idx, elm_acc.region().bulk_idx(), dof_indices, local_source_balance_rhs);
    }
}



double Elasticity::dirichlet_penalty(SideIter side)
{
    double young = data_.young_modulus.value(side->centre(), side->element());
    double poisson = data_.poisson_ratio.value(side->centre(), side->element());
    return 1e3*(2*lame_mu(young, poisson) + lame_lambda(young, poisson)) / side->measure();
}



template<unsigned int dim>
void Elasticity::assemble_fluxes_boundary()
{
    FEValues<3> fe_values_side(*feo->q<dim-1>(), *feo->fe<dim>(),
    		update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim-1>()->size();
    vector<int> side_dof_indices(ndofs);
    PetscScalar local_matrix[ndofs*ndofs];
    auto vec = fe_values_side.vector_view(0);

    // assemble boundary integral
    for (unsigned int iedg=0; iedg<feo->dh()->n_loc_edges(); iedg++)
    {
    	Edge edg = mesh_->edge(feo->dh()->edge_index(iedg));
    	if (edg.n_sides() > 1) continue;
    	// check spatial dimension
    	if (edg.side(0)->dim() != dim-1) continue;
    	// skip edges lying not on the boundary
    	if (! edg.side(0)->is_boundary()) continue;

    	SideIter side = edg.side(0);
        ElementAccessor<3> cell = side->element();
        DHCellAccessor dh_cell = feo->dh()->cell_accessor_from_element(cell.idx());
        DHCellSide dh_side(dh_cell, side->side_idx());
        dh_cell.get_dof_indices(side_dof_indices);
        unsigned int bc_type = data_.bc_type.value(side->centre(), side->cond().element_accessor());
        fe_values_side.reinit(*side);

        for (unsigned int i=0; i<ndofs; i++)
            for (unsigned int j=0; j<ndofs; j++)
                local_matrix[i*ndofs+j] = 0;
        
        if (bc_type == EqData::bc_type_displacement)
        {
            for (unsigned int k=0; k<qsize; k++)
                for (unsigned int i=0; i<ndofs; i++)
                    for (unsigned int j=0; j<ndofs; j++)
                        local_matrix[i*ndofs+j] += dirichlet_penalty(side)*arma::dot(vec.value(i,k),vec.value(j,k))*fe_values_side.JxW(k);
        }
        else if (bc_type == EqData::bc_type_displacement_normal)
        {
            for (unsigned int k=0; k<qsize; k++)
                for (unsigned int i=0; i<ndofs; i++)
                    for (unsigned int j=0; j<ndofs; j++)
                        local_matrix[i*ndofs+j] += dirichlet_penalty(side)*arma::dot(vec.value(i,k),fe_values_side.normal_vector(k))*arma::dot(vec.value(j,k),fe_values_side.normal_vector(k))*fe_values_side.JxW(k);
        }
        
        ls->mat_set_values(ndofs, side_dof_indices.data(), ndofs, side_dof_indices.data(), local_matrix);
    }
}


arma::mat33 mat_t(const arma::mat33 &m, const arma::vec3 &n)
{
  arma::mat33 mt = m - m*arma::kron(n,n.t());
  return mt;
}


template<unsigned int dim>
void Elasticity::assemble_matrix_element_side()
{
	if (dim == 1) return;
    FEValues<3> fe_values_sub(*feo->q<dim-1>(), *feo->fe<dim-1>(),
    		update_values | update_gradients | update_JxW_values | update_quadrature_points);
    FEValues<3> fe_values_side(*feo->q<dim-1>(), *feo->fe<dim>(),
    		update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
 
    const unsigned int ndofs_side = feo->fe<dim>()->n_dofs();    // number of local dofs
    const unsigned int ndofs_sub  = feo->fe<dim-1>()->n_dofs();
    const unsigned int qsize = feo->q<dim-1>()->size();     // number of quadrature points
    vector<vector<int> > side_dof_indices(2);
    vector<unsigned int> n_dofs = { ndofs_sub, ndofs_side };
	vector<double> frac_sigma(qsize);
	vector<double> csection_lower(qsize), csection_higher(qsize), young(qsize), poisson(qsize), alpha(qsize);
    PetscScalar local_matrix[2][2][(ndofs_side)*(ndofs_side)];
    auto vec_side = fe_values_side.vector_view(0);
    auto vec_sub = fe_values_sub.vector_view(0);

    // index 0 = element with lower dimension,
    // index 1 = side of element with higher dimension
    side_dof_indices[0].resize(ndofs_sub);
    side_dof_indices[1].resize(ndofs_side);

    // assemble integral over sides
    for (unsigned int inb=0; inb<feo->dh()->n_loc_nb(); inb++)
    {
    	Neighbour *nb = &mesh_->vb_neighbours_[feo->dh()->nb_index(inb)];
        // skip neighbours of different dimension
        if (nb->element()->dim() != dim-1) continue;

		ElementAccessor<3> cell_sub = nb->element();
        DHCellAccessor dh_cell_sub = feo->dh()->cell_accessor_from_element(cell_sub.idx());
        dh_cell_sub.get_dof_indices(side_dof_indices[0]);
		fe_values_sub.reinit(cell_sub);

		ElementAccessor<3> cell = nb->side()->element();
		DHCellAccessor dh_cell = feo->dh()->cell_accessor_from_element(cell.idx());
        DHCellSide dh_side(dh_cell, nb->side()->side_idx());
        dh_cell.get_dof_indices(side_dof_indices[1]);
		fe_values_side.reinit(dh_side.side());

		// Element id's for testing if they belong to local partition.
		bool own_element_id[2];
		own_element_id[0] = feo->dh()->cell_accessor_from_element(cell_sub.idx()).is_own();
		own_element_id[1] = feo->dh()->cell_accessor_from_element(cell.idx()).is_own();
        
		data_.cross_section.value_list(fe_values_sub.point_list(), cell_sub, csection_lower);
		data_.cross_section.value_list(fe_values_sub.point_list(), cell, csection_higher);
 		data_.fracture_sigma.value_list(fe_values_sub.point_list(), cell_sub, frac_sigma);
        data_.young_modulus.value_list(fe_values_sub.point_list(), cell_sub, young);
        data_.poisson_ratio.value_list(fe_values_sub.point_list(), cell_sub, poisson);
        
        for (unsigned int n=0; n<2; ++n)
            for (unsigned int i=0; i<ndofs_side; i++)
                for (unsigned int m=0; m<2; ++m)
                    for (unsigned int j=0; j<ndofs_side; j++)
                        local_matrix[n][m][i*(ndofs_side)+j] = 0;

        // set transmission conditions
        for (unsigned int k=0; k<qsize; k++)
        {
            arma::vec3 nv = fe_values_side.normal_vector(k);
            double mu = lame_mu(young[k], poisson[k]);
            double lambda = lame_lambda(young[k], poisson[k]);
            
            for (int n=0; n<2; n++)
            {
                if (!own_element_id[n]) continue;

                for (unsigned int i=0; i<n_dofs[n]; i++)
                {
                    arma::vec3 vi = (n==0)?arma::zeros(3):vec_side.value(i,k);
                    arma::vec3 vf = (n==1)?arma::zeros(3):vec_sub.value(i,k);
                    arma::mat33 gvft = (n==0)?vec_sub.grad(i,k):arma::zeros(3,3);
                    
                    for (int m=0; m<2; m++)
                        for (unsigned int j=0; j<n_dofs[m]; j++) {
                            arma::vec3 ui = (m==0)?arma::zeros(3):vec_side.value(j,k);
                            arma::vec3 uf = (m==1)?arma::zeros(3):vec_sub.value(j,k);
                            arma::mat33 guit = (m==1)?mat_t(vec_side.grad(j,k),nv):arma::zeros(3,3);
                            double divuit = (m==1)?arma::trace(guit):0;
                            
                            local_matrix[n][m][i*n_dofs[m] + j] +=
                                    frac_sigma[k]*(
                                     arma::dot(vf-vi,
                                      2/csection_lower[k]*(mu*(uf-ui)+(mu+lambda)*(arma::dot(uf-ui,nv)*nv))
                                      + mu*arma::trans(guit)*nv
                                      + lambda*divuit*nv
                                     )
                                     - arma::dot(gvft, mu*arma::kron(nv,ui.t()) + lambda*arma::dot(ui,nv)*arma::eye(3,3))
                                    )*fe_values_sub.JxW(k);
                        }
                    
                }
            }
        }
        
        for (unsigned int n=0; n<2; ++n)
            for (unsigned int m=0; m<2; ++m)
                ls->mat_set_values(n_dofs[n], side_dof_indices[n].data(), n_dofs[m], side_dof_indices[m].data(), local_matrix[n][m]);
    }
}



template<unsigned int dim>
void Elasticity::assemble_rhs_element_side()
{
	if (dim == 1) return;
    FEValues<3> fe_values_sub(*feo->q<dim-1>(), *feo->fe<dim-1>(),
    		update_values | update_JxW_values | update_quadrature_points);
    FEValues<3> fe_values_side(*feo->q<dim-1>(), *feo->fe<dim>(),
    		update_values | update_normal_vectors);
 
    const unsigned int ndofs_side = feo->fe<dim>()->n_dofs();    // number of local dofs
    const unsigned int ndofs_sub  = feo->fe<dim-1>()->n_dofs();
    const unsigned int qsize = feo->q<dim-1>()->size();     // number of quadrature points
    vector<vector<int> > side_dof_indices(2);
    vector<unsigned int> n_dofs = { ndofs_sub, ndofs_side };
	vector<double> frac_sigma(qsize), potential(qsize);
    PetscScalar local_rhs[2][ndofs_side];
    auto vec_side = fe_values_side.vector_view(0);
    auto vec_sub = fe_values_sub.vector_view(0);

    // index 0 = element with lower dimension,
    // index 1 = side of element with higher dimension
    side_dof_indices[0].resize(ndofs_sub);
    side_dof_indices[1].resize(ndofs_side);

    // assemble integral over sides
    for (unsigned int inb=0; inb<feo->dh()->n_loc_nb(); inb++)
    {
    	Neighbour *nb = &mesh_->vb_neighbours_[feo->dh()->nb_index(inb)];
        // skip neighbours of different dimension
        if (nb->element()->dim() != dim-1) continue;

		ElementAccessor<3> cell_sub = nb->element();
		feo->dh()->cell_accessor_from_element(cell_sub.idx()).get_dof_indices(side_dof_indices[0]);
		fe_values_sub.reinit(cell_sub);

		ElementAccessor<3> cell = nb->side()->element();
		feo->dh()->cell_accessor_from_element(cell.idx()).get_dof_indices(side_dof_indices[1]);
		fe_values_side.reinit(*nb->side());

		// Element id's for testing if they belong to local partition.
		bool own_element_id[2];
		own_element_id[0] = feo->dh()->cell_accessor_from_element(cell_sub.idx()).is_own();
		own_element_id[1] = feo->dh()->cell_accessor_from_element(cell.idx()).is_own();
        
 		data_.fracture_sigma.value_list(fe_values_sub.point_list(), cell_sub, frac_sigma);
        data_.potential_load.value_list(fe_values_sub.point_list(), cell, potential);
        
        for (unsigned int n=0; n<2; ++n)
            for (unsigned int i=0; i<ndofs_side; i++)
                local_rhs[n][i] = 0;

        // set transmission conditions
        for (unsigned int k=0; k<qsize; k++)
        {
            arma::vec3 nv = fe_values_side.normal_vector(k);
            
            for (int n=0; n<2; n++)
            {
                if (!own_element_id[n]) continue;

                for (unsigned int i=0; i<n_dofs[n]; i++)
                {
                    arma::vec3 vi = (n==0)?arma::zeros(3):vec_side.value(i,k);
                    arma::vec3 vf = (n==1)?arma::zeros(3):vec_sub.value(i,k);
                    
                    local_rhs[n][i] -= frac_sigma[k]*arma::dot(vf-vi,potential[k]*nv)*fe_values_sub.JxW(k);
                }
            }
        }
        
        for (unsigned int n=0; n<2; ++n)
            ls->rhs_set_values(n_dofs[n], side_dof_indices[n].data(), local_rhs[n]);
    }
}





template<unsigned int dim>
void Elasticity::assemble_boundary_conditions()
{
    FEValues<3> fe_values_side(*feo->q<dim-1>(), *feo->fe<dim>(),
    		update_values | update_normal_vectors | update_side_JxW_values | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim-1>()->size();
    vector<int> side_dof_indices(ndofs);
    // unsigned int loc_b=0;
    double local_rhs[ndofs];
    vector<PetscScalar> local_flux_balance_vector(ndofs);
    // PetscScalar local_flux_balance_rhs;
    vector<arma::vec3> bc_values(qsize), bc_traction(qsize);
    vector<double> csection(qsize), bc_potential(qsize);
    auto vec = fe_values_side.vector_view(0);

    for (DHCellAccessor cell : feo->dh()->own_range())
    {
        ElementAccessor<3> elm = cell.elm();
        if (elm->boundary_idx_ == nullptr) continue;

        for (unsigned int si=0; si<elm->n_sides(); ++si)
        {
			Edge edg = elm.side(si)->edge();
            Side side = *cell.elm().side(si);

			if (edg.n_sides() > 1) continue;
			// skip edges lying not on the boundary
			if ( ! side.is_boundary()) continue;

			if (side.dim() != dim-1)
			{
				// if (edg.side(0)->cond() != nullptr) ++loc_b;
				continue;
			}

			ElementAccessor<3> bc_cell = side.cond().element_accessor();
 			unsigned int bc_type = data_.bc_type.value(side.centre(), bc_cell);

			fe_values_side.reinit(side);

			data_.cross_section.value_list(fe_values_side.point_list(), elm, csection);
			// The b.c. data are fetched for all possible b.c. types since we allow
			// different bc_type for each substance.
			data_.bc_displacement.value_list(fe_values_side.point_list(), bc_cell, bc_values);
            data_.bc_traction.value_list(fe_values_side.point_list(), bc_cell, bc_traction);
            data_.potential_load.value_list(fe_values_side.point_list(), elm, bc_potential);

			cell.get_dof_indices(side_dof_indices);

            fill_n(local_rhs, ndofs, 0);
            local_flux_balance_vector.assign(ndofs, 0);
            // local_flux_balance_rhs = 0;

            if (bc_type == EqData::bc_type_displacement)
            {
                for (unsigned int k=0; k<qsize; k++)
                    for (unsigned int i=0; i<ndofs; i++)
                        local_rhs[i] += dirichlet_penalty(side)*arma::dot(bc_values[k],vec.value(i,k))*fe_values_side.JxW(k);
            }
            else if (bc_type == EqData::bc_type_displacement_normal)
            {
                for (unsigned int k=0; k<qsize; k++)
                    for (unsigned int i=0; i<ndofs; i++)
                        local_rhs[i] += dirichlet_penalty(side)*arma::dot(bc_values[k],fe_values_side.normal_vector(k))*arma::dot(vec.value(i,k),fe_values_side.normal_vector(k))*fe_values_side.JxW(k);
            }
            else if (bc_type == EqData::bc_type_traction)
            {
              for (unsigned int k=0; k<qsize; k++)
              {
                for (unsigned int i=0; i<ndofs; i++)
                  local_rhs[i] += csection[k]*arma::dot(vec.value(i,k),bc_traction[k] + bc_potential[k]*fe_values_side.normal_vector(k))*fe_values_side.JxW(k);
              }
            }
            ls->rhs_set_values(ndofs, side_dof_indices.data(), local_rhs);


            
//             balance_->add_flux_matrix_values(subst_idx, loc_b, side_dof_indices, local_flux_balance_vector);
//             balance_->add_flux_vec_value(subst_idx, loc_b, local_flux_balance_rhs);
			// ++loc_b;
        }
    }
}




