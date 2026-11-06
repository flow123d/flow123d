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
#include "mechanics/assembly_elasticity.hh"

#include "io/output_time.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_system.hh"
#include "fields/field_fe.hh"
#include "la/linsys_PETSC.hh"
#include "la/linsys_PERMON.hh"
#include "coupling/balance.hh"
#include "mesh/neighbours.h"
#include "coupling/generic_assembly.hh"

#include "fields/multi_field.hh"
#include "fields/generic_field.hh"
#include "fields/field_model.hh"
#include "input/factory.hh"
#include <algorithm>

namespace {

#ifdef FLOW123D_HAVE_PERMON
/**
 * Import requested original global constraint rows into an owned local copy.
 * All source-vector ranks participate, including ranks requesting no rows.
 * Keep the source borrowed; destroy the sequential buffer and scatter locally.
 */
std::vector<PetscScalar> import_constraint_rows(Vec source, const std::vector<PetscInt> &rows)
{
    if (!source) return {};
    struct Import {
        Vec local = nullptr;
        IS indices = nullptr;
        VecScatter scatter = nullptr;
        ~Import() {
            VecScatterDestroy(&scatter);
            ISDestroy(&indices);
            VecDestroy(&local);
        }
    } work;
    chkerr(VecCreateSeq(PETSC_COMM_SELF, rows.size(), &work.local));
    chkerr(ISCreateGeneral(PETSC_COMM_SELF, rows.size(), rows.data(), PETSC_COPY_VALUES, &work.indices));
    chkerr(VecScatterCreate(source, work.indices, work.local, nullptr, &work.scatter));
    chkerr(VecScatterBegin(work.scatter, source, work.local, INSERT_VALUES, SCATTER_FORWARD));
    chkerr(VecScatterEnd(work.scatter, source, work.local, INSERT_VALUES, SCATTER_FORWARD));
    const PetscScalar *values;
    chkerr(VecGetArrayRead(work.local, &values));
    std::vector<PetscScalar> result(rows.size());
    if (!result.empty()) std::copy_n(values, result.size(), result.begin());
    chkerr(VecRestoreArrayRead(work.local, &values));
    return result;
}
#endif

} // namespace

using namespace Input::Type;



const Record & Elasticity::get_input_type() {
    std::string equation_name = std::string(name_) + "_FE";
	return IT::Record(
                std::string(equation_name),
                "FEM for linear elasticity.")
           .copy_keys(EquationBase::record_template())
		   .copy_keys(EquationBase::user_fields_template(equation_name))
           .declare_key("balance", Balance::get_input_type(), Default("{}"),
                    "Settings for computing balance.")
           .declare_key("output_stream", OutputTime::get_input_type(), Default::obligatory(),
                    "Parameters of output stream.")
           .declare_key("solver", LinSys_PETSC::get_input_type(), Default::obligatory(),
				"Linear solver for elasticity.")
		   .declare_key("input_fields", Array(
		        Elasticity::EqFields()
		            .make_field_descriptor_type(equation_name)),
		        IT::Default::obligatory(),
		        "Input fields of the equation.")
           .declare_key("output",
                EqFields().output_fields.make_output_type(equation_name, ""),
                IT::Default("{ \"fields\": [ \"displacement\" ] }"),
                "Setting of the field output.")
           .declare_key("contact", Bool(), IT::Default("false"), "Indicates the use of contact conditions on fractures.")
           .declare_key("tresca_friction", Bool(), IT::Default("false"),
                    "Enable incremental Tresca friction on two-sided line fractures. Requires contact=true. "
                    "The prescribed friction_threshold is independent of normal traction, including on open fractures. "
                    "Positive thresholds on triangular fracture surfaces are not yet supported.")
		   .close();
}

const int Elasticity::registrar =
		Input::register_class< Elasticity, Mesh &, const Input::Record>(std::string(name_) + "_FE") +
		Elasticity::get_input_type().size();



double lame_mu(double young, double poisson)
{
    return young*0.5/(poisson+1.);
}


double lame_lambda(double young, double poisson)
{
    return young*poisson/((poisson+1.)*(1.-2.*poisson));
}

// Functor computing lame_mu
struct fn_lame_mu {
	inline double operator() (double young, double poisson) {
        return young * 0.5 / (poisson+1.);
    }
};

// Functor computing lame_lambda
struct fn_lame_lambda {
	inline double operator() (double young, double poisson) {
        return young * poisson / ((poisson+1.)*(1.-2.*poisson));
    }
};

// Functor computing base of dirichlet_penalty (without dividing by side meassure)
struct fn_dirichlet_penalty {
	inline double operator() (double lame_mu, double lame_lambda) {
        return 1e3 * (2 * lame_mu + lame_lambda);
    }
};


arma::mat33 Elasticity::EqFields::stress_tensor(BulkPoint &p, const arma::mat33 &strain_tensor)
{
    return 2*lame_mu(p)*strain_tensor + lame_lambda(p)*arma::trace(strain_tensor)*arma::eye(3,3);
}




const Selection & Elasticity::EqFields::get_bc_type_selection() {
    return Selection("Elasticity_BC_Type", "Types of boundary conditions for mechanics.")
            .add_value(bc_type_displacement, "displacement",
                  "Prescribed displacement.")
            .add_value(bc_type_displacement_normal, "displacement_n",
                  "Prescribed displacement in the normal direction to the boundary.")
            .add_value(bc_type_traction, "traction",
                  "Prescribed traction.")
            .add_value(bc_type_stress, "stress",
                  "Prescribed stress tensor.")
            .close();
}


Elasticity::EqFields::EqFields()
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
        .description("Prescribed displacement on boundary.")
        .units( UnitSI().m() )
        .input_default("0.0")
        .flags_add(in_rhs);
        
    *this+=bc_traction
        .name("bc_traction")
        .description("Prescribed traction on boundary.")
        .units( UnitSI().Pa() )
        .input_default("0.0")
        .flags_add(in_rhs);
    
    *this+=bc_stress
        .name("bc_stress")
        .description("Prescribed stress on boundary.")
        .units( UnitSI().Pa() )
        .input_default("0.0")
        .flags_add(in_rhs);
        
    *this+=load
        .name("load")
        .description("Prescribed bulk load.")
        .units( UnitSI().N().m(-3) )
        .input_default("0.0")
        .flags_add(in_rhs);
    
    *this+=young_modulus
        .name("young_modulus")
        .description("Young's modulus.")
        .units( UnitSI().Pa() )
         .input_default("0.0")
        .flags_add(in_main_matrix & in_rhs);
    
    *this+=poisson_ratio
        .name("poisson_ratio")
        .description("Poisson's ratio.")
        .units( UnitSI().dimensionless() )
         .input_default("0.0")
        .flags_add(in_main_matrix & in_rhs);
        
    *this+=fracture_sigma
            .name("fracture_sigma")
            .description(
            "Coefficient of transfer of forces through fractures.")
            .units( UnitSI::dimensionless() )
            .input_default("1.0")
            .flags_add(in_main_matrix & in_rhs);

    *this+=initial_stress
        .name("initial_stress")
        .description("Initial stress tensor.")
        .units( UnitSI().Pa() )
        .input_default("0.0")
        .flags_add(in_rhs);

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
    
    *this+=cross_section_min
      .name("cross_section_min")
      .description("Minimal cross-section of fractures.")
      .units( UnitSI().m(3).md() )
      .input_default("0.0");
      
    *this+=potential_load
      .name("potential_load")
      .units( UnitSI().m() )
      .flags(input_copy & in_rhs);

    *this += friction_threshold.name("friction_threshold")
      .description("Prescribed nonnegative Tresca shear traction threshold on two-sided line fractures. "
                   "Used only with tresca_friction=true; zero gives frictionless tangential motion.")
      .units(UnitSI().Pa())
      .input_default("0.0");

    *this+=output_field
            .name("displacement")
            .description("Displacement vector field output.")
            .units( UnitSI().m() )
            .flags(equation_result);
    
    *this += output_stress
            .name("stress")
            .description("Stress tensor output.")
            .units( UnitSI().Pa() )
            .flags(equation_result);
    
    *this += output_von_mises_stress
            .name("von_mises_stress")
            .description("von Mises stress output.")
            .units( UnitSI().Pa() )
            .flags(equation_result);
    
    *this += output_mean_stress
            .name("mean_stress")
            .description("mean stress output.")
            .units( UnitSI().Pa() )
            .flags(equation_result);

    *this += output_cross_section
            .name("cross_section_updated")
            .description("Cross-section after deformation - output.")
            .units( UnitSI().m() )
            .flags(equation_result);
            
    *this += output_divergence
            .name("displacement_divergence")
            .description("Displacement divergence output.")
            .units( UnitSI().dimensionless() )
            .flags(equation_result);

    *this += output_contact_normal_stress
            .name("contact_normal_stress")
            .description("Element-average normal contact stress; compression is negative. "
                         "Zero on cells without contact constraints.")
            .units(UnitSI().Pa())
            .flags(equation_result);

    *this += output_contact_tangential_stress
            .name("contact_tangential_stress")
            .description("Element-average Cartesian Tresca reaction on the bulk neighbour with the lowest input ID. "
                         "Its magnitude is independent of side selection; its sign follows that neighbour. "
                         "Zero on cells without Tresca friction.")
            .units(UnitSI().Pa())
            .flags(equation_result);

    *this += lame_mu.name("lame_mu")
            .description("Field lame_mu.")
            .input_default("0.0")
            .units( UnitSI().Pa() );

    *this += lame_lambda.name("lame_lambda")
            .description("Field lame_lambda.")
            .input_default("0.0")
            .units( UnitSI().Pa() );

    *this += dirichlet_penalty.name("dirichlet_penalty")
            .description("Field dirichlet_penalty.")
            .input_default("0.0")
            .units( UnitSI().Pa() );

    // add all input fields to the output list
    output_fields += *this;

    this->add_coords_field();
    this->set_default_fieldset();

}

void Elasticity::EqData::create_dh(Mesh * mesh)
{
	ASSERT_EQ(quad_order(), 1).error("Unsupported polynomial order for finite elements in Elasticity");
    MixedPtr<FE_P> fe_p( quad_order() );
    MixedPtr<FiniteElement> fe = mixed_fe_system(fe_p, FEVector, 3);

    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe);
	dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh);

	dh_->distribute_dofs(ds);
}

void Elasticity::EqData::allocate_tangential_jump_matrix()
{
    ASSERT_PERMANENT(tangential_jump_matrix == nullptr).error("Tangential jump matrix is already allocated.");
    PetscInt n_local_rows = 0;
    for (auto cell : dh_->own_range())
        if (cell.dim() > 0 && cell.dim() < 3 && cell.elm()->n_neighs_vb() == 2)
            n_local_rows += cell.dim();

    PetscInt max_nnz = 2 * dh_->max_elem_dofs();
    PetscInt n_local_cols = dh_->lsize();
    PetscInt n_global_cols = dh_->n_global_dofs();
    chkerr(MatCreateAIJ(PETSC_COMM_WORLD, n_local_rows, n_local_cols, PETSC_DECIDE,
                       n_global_cols, std::min(max_nnz, n_local_cols), nullptr,
                       std::min(max_nnz, n_global_cols - n_local_cols), nullptr,
                       &tangential_jump_matrix));

    PetscInt first_row, end_row;
    chkerr(MatGetOwnershipRange(tangential_jump_matrix, &first_row, &end_row));
    tangential_jump_idx.clear();
    for (auto cell : dh_->own_range()) {
        if (cell.dim() > 0 && cell.dim() < 3 && cell.elm()->n_neighs_vb() == 2) {
            tangential_jump_idx[cell.elm_idx()] = first_row;
            first_row += cell.dim();
        }
    }
    ASSERT_EQ(first_row, end_row);
}


arma::mat Elasticity::EqData::tangential_frame(DHCellAccessor cell)
{
    ASSERT_PERMANENT(cell.dim() == 1 || cell.dim() == 2)
        .error("Tangential frames require a line or triangular fracture.");
    auto fracture = cell.elm();
    arma::mat frame(3, cell.dim());
    frame.col(0) = *fracture.node(1) - *fracture.node(0);
    frame.col(0) /= arma::norm(frame.col(0), 2);
    if (cell.dim() == 2) {
        frame.col(1) = *fracture.node(2) - *fracture.node(0);
        frame.col(1) -= arma::dot(frame.col(1), frame.col(0)) * frame.col(0);
        frame.col(1) /= arma::norm(frame.col(1), 2);
    }
    return frame;
}


void Elasticity::OutputEqData::create_dh(Mesh * mesh)
{
	ASSERT_EQ(quad_order(), 0).error("Unsupported polynomial order for output in Elasticity");

	MixedPtr<FE_P_disc> fe_p_disc( quad_order() );
    dh_scalar_ = make_shared<DOFHandlerMultiDim>(*mesh);
	std::shared_ptr<DiscreteSpace> ds_scalar = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe_p_disc);
	dh_scalar_->distribute_dofs(ds_scalar);


    MixedPtr<FiniteElement> fe_t = mixed_fe_system(MixedPtr<FE_P_disc>( quad_order() ), FEType::FETensor, 9);
    dh_tensor_ = make_shared<DOFHandlerMultiDim>(*mesh);
	std::shared_ptr<DiscreteSpace> dst = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe_t);
	dh_tensor_->distribute_dofs(dst);

    MixedPtr<FiniteElement> fe_v = mixed_fe_system(MixedPtr<FE_P_disc>(quad_order()), FEVector, 3);
    dh_vector_ = make_shared<DOFHandlerMultiDim>(*mesh);
    std::shared_ptr<DiscreteSpace> ds_vector = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe_v);
    dh_vector_->distribute_dofs(ds_vector);
}


Elasticity::Elasticity(Mesh & init_mesh, const Input::Record in_rec, TimeGovernor *tm)
        : EquationBase(init_mesh, in_rec),
		  input_rec(in_rec),
		  stiffness_assembly_(nullptr),
		  rhs_assembly_(nullptr),
          constraint_assembly_(nullptr),
          tangential_jump_assembly_(nullptr),
          tresca_threshold_assembly_(nullptr),
		  output_fields_assembly_(nullptr)
{
	// Can not use name() + "constructor" here, since START_TIMER only accepts const char *
	// due to constexpr optimization.
	START_TIMER("Mechanics constructor");

    eq_fields_ = std::make_shared<EqFields>();
    eq_data_ = std::make_shared<EqData>(eq_fields_);
    output_eq_data_ = std::make_shared<OutputEqData>(eq_fields_);
    this->eq_fieldset_ = eq_fields_;
    
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
    eq_fields_->set_mesh(init_mesh);
    eq_fields_->region_id = GenericField<3>::region_id(*mesh_);
    eq_fields_->subdomain = GenericField<3>::subdomain(*mesh_);
    eq_data_->balance_ = this->balance();
    
    // create finite element structures and distribute DOFs
    eq_data_->create_dh(mesh_);
    output_eq_data_->create_dh(mesh_);
    DebugOut().fmt("Mechanics: solution size {}\n", eq_data_->dh_->n_global_dofs());
    END_TIMER("Mechanics constructor");
}


void Elasticity::initialize()
{
    output_stream_ = OutputTime::create_output_stream("mechanics", input_rec.val<Input::Record>("output_stream"), time().get_unit_conversion());
    
    eq_fields_->set_components({"displacement"});
    eq_fields_->set_input_list( input_rec.val<Input::Array>("input_fields"), time() );
    
//     balance_ = std::make_shared<Balance>("mechanics", mesh_);
//     balance_->init_from_input(input_rec.val<Input::Record>("balance"), *time_);
    // initialization of balance object
//     eq_data_->balance_idx_ = {
//         balance_->add_quantity("force_x"),
//         balance_->add_quantity("force_y"),
//         balance_->add_quantity("force_z")
//     };
//     balance_->units(UnitSI().kg().m().s(-2));

    // create shared pointer to a FieldFE, pass FE data and push this FieldFE to output_field on all regions
    eq_fields_->output_field_ptr = create_field_fe<3, FieldValue<3>::VectorFixed>(eq_data_->dh_);
    eq_fields_->output_field.set(eq_fields_->output_field_ptr, 0.);
    
    // setup output stress
    eq_fields_->output_stress_ptr = create_field_fe<3, FieldValue<3>::TensorFixed>(output_eq_data_->dh_tensor_);
    eq_fields_->output_stress.set(eq_fields_->output_stress_ptr, 0.);
    
    // setup output von Mises stress
    eq_fields_->output_von_mises_stress_ptr = create_field_fe<3, FieldValue<3>::Scalar>(output_eq_data_->dh_scalar_);
    eq_fields_->output_von_mises_stress.set(eq_fields_->output_von_mises_stress_ptr, 0.);

    // setup output mean stress
    eq_fields_->output_mean_stress_ptr = create_field_fe<3, FieldValue<3>::Scalar>(output_eq_data_->dh_scalar_);
    eq_fields_->output_mean_stress.set(eq_fields_->output_mean_stress_ptr, 0.);
    
    // setup output cross-section
    eq_fields_->output_cross_section_ptr = create_field_fe<3, FieldValue<3>::Scalar>(output_eq_data_->dh_scalar_);
    eq_fields_->output_cross_section.set(eq_fields_->output_cross_section_ptr, 0.);
    
    // setup output divergence
    eq_fields_->output_div_ptr = create_field_fe<3, FieldValue<3>::Scalar>(output_eq_data_->dh_scalar_);
    eq_fields_->output_divergence.set(eq_fields_->output_div_ptr, 0.);

    // Allocate independent P0 storage for multiplier-derived fracture stresses.
    eq_fields_->output_contact_normal_stress_ptr =
        create_field_fe<3, FieldValue<3>::Scalar>(output_eq_data_->dh_scalar_);
    eq_fields_->output_contact_normal_stress.set(eq_fields_->output_contact_normal_stress_ptr, 0.);
    eq_fields_->output_contact_normal_stress.output_type(OutputTime::ELEM_DATA);
    eq_fields_->output_contact_tangential_stress_ptr =
        create_field_fe<3, FieldValue<3>::VectorFixed>(output_eq_data_->dh_vector_);
    eq_fields_->output_contact_tangential_stress.set(eq_fields_->output_contact_tangential_stress_ptr, 0.);
    eq_fields_->output_contact_tangential_stress.output_type(OutputTime::ELEM_DATA);

    // read optional user fields
    Input::Array user_fields_arr;
    if (input_rec.opt_val("user_fields", user_fields_arr)) {
       	this->init_user_fields(user_fields_arr, eq_fields_->output_fields);
    }
    
    eq_fields_->output_fields.set_mesh(*mesh_);
    eq_fields_->output_field.output_type(OutputTime::CORNER_DATA);

    // set time marks for writing the output
    eq_fields_->output_fields.initialize(output_stream_, mesh_, input_rec.val<Input::Record>("output"), this->time());

    // set instances of FieldModel
    eq_fields_->lame_mu.set(Model<3, FieldValue<3>::Scalar>::create(fn_lame_mu(), eq_fields_->young_modulus, eq_fields_->poisson_ratio), 0.0);
    eq_fields_->lame_lambda.set(Model<3, FieldValue<3>::Scalar>::create(fn_lame_lambda(), eq_fields_->young_modulus, eq_fields_->poisson_ratio), 0.0);
    eq_fields_->dirichlet_penalty.set(Model<3, FieldValue<3>::Scalar>::create(fn_dirichlet_penalty(), eq_fields_->lame_mu, eq_fields_->lame_lambda), 0.0);

    // equation default PETSc solver options
    std::string petsc_default_opts;
    petsc_default_opts = "-ksp_type cg -pc_type hypre -pc_hypre_type boomeramg";
    
    // allocate matrix and vector structures
    LinSys *ls;
    has_contact_ = input_rec.val<bool>("contact");
    has_tresca_ = input_rec.val<bool>("tresca_friction");
    ASSERT_PERMANENT(!has_tresca_ || has_contact_).error("Tresca friction requires contact=true.");
    if (has_contact_) {
#ifndef FLOW123D_HAVE_PERMON
        ASSERT(false).error("Flow123d was not built with PERMON library, therefore contact conditions are unsupported.");
#endif //FLOW123D_HAVE_PERMON
        ls = new LinSys_PERMON(eq_data_->dh_->distr().get(), petsc_default_opts);

        // allocate constraint matrix and vector
        unsigned int n_own_constraints = 0; // count locally owned cells with neighbours
        for (auto cell : eq_data_->dh_->own_range())
            if (cell.elm()->n_neighs_vb() > 0)
                n_own_constraints++;
        unsigned int n_constraints = 0; // count all cells with neighbours
        for (auto elm : mesh_->elements_range())
            if (elm->n_neighs_vb() > 0)
                eq_data_->constraint_idx[elm.idx()] = n_constraints++;
        unsigned int nnz = eq_data_->dh_->ds()->fe()[1_d]->n_dofs()*mesh_->max_edge_sides(1) +
                        eq_data_->dh_->ds()->fe()[2_d]->n_dofs()*mesh_->max_edge_sides(2) +
                        eq_data_->dh_->ds()->fe()[3_d]->n_dofs()*mesh_->max_edge_sides(3);
        MatCreateAIJ(PETSC_COMM_WORLD, n_own_constraints, eq_data_->dh_->lsize(), PETSC_DECIDE, PETSC_DECIDE, nnz, 0, nnz, 0, &eq_data_->constraint_matrix);
        VecCreateMPI(PETSC_COMM_WORLD, n_own_constraints, PETSC_DECIDE, &eq_data_->constraint_vec);
        ((LinSys_PERMON*)ls)->set_inequality(eq_data_->constraint_matrix,eq_data_->constraint_vec);

        constraint_assembly_ = new GenericAssembly< ConstraintAssemblyDim >(eq_data_.get(), eq_data_->dh_.get());
        eq_data_->allocate_tangential_jump_matrix();
        tangential_jump_assembly_ = new GenericAssembly< TangentialJumpAssemblyDim >(eq_data_.get(), eq_data_->dh_.get());
        if (has_tresca_) {
            chkerr(MatCreateVecs(eq_data_->tangential_jump_matrix, nullptr, &eq_data_->tresca_force_bound));
            chkerr(VecDuplicate(eq_data_->tresca_force_bound, &eq_data_->tresca_force_lower_bound));
            chkerr(VecDuplicate(eq_data_->tresca_force_bound, &eq_data_->tresca_reference_jump));
            chkerr(VecSet(eq_data_->tresca_reference_jump, 0.0));
            friction_reference_time_ = time_->t();
            static_cast<LinSys_PERMON *>(ls)->set_dual_box(eq_data_->tangential_jump_matrix,
                    eq_data_->tresca_reference_jump, eq_data_->tresca_force_lower_bound, eq_data_->tresca_force_bound);
            tresca_threshold_assembly_ = new GenericAssembly<TrescaThresholdAssemblyDim>(eq_data_.get(), eq_data_->dh_.get());
        }
    } else {
        ls = new LinSys_PETSC(eq_data_->dh_->distr().get(), petsc_default_opts);
        ((LinSys_PETSC*)ls)->set_initial_guess_nonzero();
    }
    ls->set_from_input( input_rec.val<Input::Record>("solver") );
    ls->set_solution(eq_fields_->output_field_ptr->vec().petsc_vec());
    eq_data_->ls = ls;

    stiffness_assembly_ = new GenericAssembly< StiffnessAssemblyDim >(eq_data_.get(), eq_data_->dh_.get());
    rhs_assembly_ = new GenericAssembly< RhsAssemblyDim >(eq_data_.get(), eq_data_->dh_.get());
    output_fields_assembly_ = new GenericAssembly< OutpuFieldsAssemblyDim >(output_eq_data_.get(), eq_data_->dh_.get());

    // initialization of balance object
//     balance_->allocate(eq_data_->dh_->distr()->lsize(),
//             max(feo->fe<1>()->n_dofs(), max(feo->fe<2>()->n_dofs(), feo->fe<3>()->n_dofs())));

}


Elasticity::~Elasticity()
{
//     delete time_;

    if (stiffness_assembly_!=nullptr) delete stiffness_assembly_;
    if (rhs_assembly_!=nullptr) delete rhs_assembly_;
    if (constraint_assembly_ != nullptr) delete constraint_assembly_;
    if (tangential_jump_assembly_ != nullptr) delete tangential_jump_assembly_;
    if (tresca_threshold_assembly_ != nullptr) delete tresca_threshold_assembly_;
    if (output_fields_assembly_!=nullptr) delete output_fields_assembly_;

    eq_data_.reset();
    eq_fields_.reset();
}



void Elasticity::update_output_fields()
{
    eq_fields_->set_time(time_->step(), LimitSide::right);
    
    // update ghost values of solution vector and prepare dependent fields
	eq_fields_->output_field_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_stress_ptr->vec().zero_entries();
    eq_fields_->output_von_mises_stress_ptr->vec().zero_entries();
    eq_fields_->output_mean_stress_ptr->vec().zero_entries();
	eq_fields_->output_cross_section_ptr->vec().zero_entries();
	eq_fields_->output_div_ptr->vec().zero_entries();
	eq_fields_->output_field_ptr->vec().local_to_ghost_end();

    // compute new output fields depending on solution (stress, divergence etc.)
    output_fields_assembly_->assemble(eq_data_->dh_);
    update_contact_stress_fields();

    // finish assembly of vectors
    eq_fields_->output_stress_ptr->vec().assembly_begin();
    eq_fields_->output_stress_ptr->vec().assembly_end();
    eq_fields_->output_von_mises_stress_ptr->vec().assembly_begin();
    eq_fields_->output_von_mises_stress_ptr->vec().assembly_end();
    eq_fields_->output_mean_stress_ptr->vec().assembly_begin();
    eq_fields_->output_mean_stress_ptr->vec().assembly_end();
    eq_fields_->output_cross_section_ptr->vec().assembly_begin();
    eq_fields_->output_cross_section_ptr->vec().assembly_end();
    eq_fields_->output_div_ptr->vec().assembly_begin();
    eq_fields_->output_div_ptr->vec().assembly_end();
    eq_fields_->output_contact_normal_stress_ptr->vec().assembly_begin();
    eq_fields_->output_contact_normal_stress_ptr->vec().assembly_end();
    eq_fields_->output_contact_tangential_stress_ptr->vec().assembly_begin();
    eq_fields_->output_contact_tangential_stress_ptr->vec().assembly_end();

    // update ghost values of computed fields
    eq_fields_->output_stress_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_von_mises_stress_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_mean_stress_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_cross_section_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_div_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_contact_normal_stress_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_contact_tangential_stress_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_stress_ptr->vec().local_to_ghost_end();
    eq_fields_->output_von_mises_stress_ptr->vec().local_to_ghost_end();
    eq_fields_->output_mean_stress_ptr->vec().local_to_ghost_end();
    eq_fields_->output_cross_section_ptr->vec().local_to_ghost_end();
    eq_fields_->output_div_ptr->vec().local_to_ghost_end();
    eq_fields_->output_contact_normal_stress_ptr->vec().local_to_ghost_end();
    eq_fields_->output_contact_tangential_stress_ptr->vec().local_to_ghost_end();
}


void Elasticity::update_contact_stress_fields()
{
    // Clear every cell so disabled friction and newly inactive reactions cannot leave stale values.
    auto &normal_output = eq_fields_->output_contact_normal_stress_ptr->vec();
    auto &tangential_output = eq_fields_->output_contact_tangential_stress_ptr->vec();
    normal_output.zero_entries();
    tangential_output.zero_entries();
    if (!has_contact_) return;

#ifdef FLOW123D_HAVE_PERMON
    // Borrow the original-row solver results and retain zero output during HM initialization,
    // which requests fields before the first mechanics solve.
    auto *solver = static_cast<LinSys_PERMON *>(eq_data_->ls);
    Vec normal = solver->dual_ineq_multipliers();
    Vec tangential = has_tresca_ ? solver->dual_box_multipliers() : nullptr;
    if (!normal && !tangential) return;
    PetscInt n_normal, n_tangential = 0;
    chkerr(MatGetSize(eq_data_->constraint_matrix, &n_normal, nullptr));
    if (has_tresca_) chkerr(MatGetSize(eq_data_->tangential_jump_matrix, &n_tangential, nullptr));
    ASSERT_PERMANENT(n_normal == 0 || normal != nullptr)
        .error("Normal contact stress requires solved inequality multipliers.");
    ASSERT_PERMANENT(n_tangential == 0 || tangential != nullptr)
        .error("Tangential contact stress requires solved Tresca multipliers.");

    // Request rows for owned fracture cells, independently of the multiplier vectors' ownership.
    std::vector<PetscInt> normal_rows, tangential_rows;
    for (auto cell : eq_data_->dh_->own_range()) {
        auto normal_row = eq_data_->constraint_idx.find(cell.elm_idx());
        if (normal_row != eq_data_->constraint_idx.end()) normal_rows.push_back(normal_row->second);
        auto tangential_row = eq_data_->tangential_jump_idx.find(cell.elm_idx());
        if (has_tresca_ && tangential_row != eq_data_->tangential_jump_idx.end())
            for (unsigned int k = 0; k < cell.dim(); ++k)
                tangential_rows.push_back(tangential_row->second + k);
    }
    auto normal_values = import_constraint_rows(normal, normal_rows);
    auto tangential_values = import_constraint_rows(tangential, tangential_rows);

    // Consume imported rows in the same owned-cell order and normalize average-jump reactions.
    std::size_t normal_pos = 0, tangential_pos = 0;
    for (auto cell : eq_data_->dh_->own_range()) {
        auto normal_row = eq_data_->constraint_idx.find(cell.elm_idx());
        auto tangential_row = eq_data_->tangential_jump_idx.find(cell.elm_idx());
        const bool has_normal = normal_row != eq_data_->constraint_idx.end();
        const bool has_tangential = has_tresca_ && tangential_row != eq_data_->tangential_jump_idx.end();
        if (!has_normal && !has_tangential) continue;
        double measure = cell.elm().measure();
        ASSERT_PERMANENT_GT(measure, 0.0).error("Contact stress requires positive fracture measure.");
        if (has_normal) {
            auto output_cell = cell.cell_with_other_dh(output_eq_data_->dh_scalar_.get());
            normal_output.set(output_cell.get_loc_dof_indices()[0], -normal_values[normal_pos++] / measure);
        }
        if (has_tangential) {
            // Reconstruct the Cartesian reaction on the lowest-ID neighbour, without an aperture factor.
            arma::mat frame = EqData::tangential_frame(cell);
            arma::vec3 traction;
            traction.zeros();
            for (unsigned int k = 0; k < cell.dim(); ++k)
                traction -= frame.col(k) * tangential_values[tangential_pos++] / measure;
            auto output_cell = cell.cell_with_other_dh(output_eq_data_->dh_vector_.get());
            auto dofs = output_cell.get_loc_dof_indices();
            for (unsigned int i = 0; i < dofs.size(); ++i)
                tangential_output.set(dofs[i], arma::dot(output_cell.cell_dof(i).coefs, traction));
        }
    }
#endif
}




void Elasticity::zero_time_step()
{
	START_TIMER("Mechanics zero time step");
	eq_fields_->mark_input_times( *time_ );
	eq_fields_->set_time(time_->step(), LimitSide::right);
	std::stringstream ss; // print warning message with table of uninitialized fields
	if ( FieldCommon::print_message_table(ss, "mechanics") ) {
		WarningOut() << ss.str();
	}

    preallocate();
    

    // after preallocation we assemble the matrices and vectors required for balance of forces
//     for (auto subst_idx : eq_data_->balance_idx_)
//         balance_->calculate_instant(subst_idx, eq_data_->ls->get_solution());
    
//     update_solution();
    eq_data_->ls->start_add_assembly();
    MatSetOption(*eq_data_->ls->get_matrix(), MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    eq_data_->ls->mat_zero_entries();
    eq_data_->ls->rhs_zero_entries();
    stiffness_assembly_->assemble(eq_data_->dh_);
    rhs_assembly_->assemble(eq_data_->dh_);
    eq_data_->ls->finish_assembly();
    if (has_tresca_) prepare_tresca();
    LinSys::SolveInfo si = eq_data_->ls->solve();
    if (has_tresca_)
        ASSERT_PERMANENT_GT(si.converged_reason, 0).error("Tresca solve did not converge; load step is not accepted.");
    MessageOut().fmt("[mech solver] lin. it: {}, reason: {}, residual: {}\n",
        		si.n_iterations, si.converged_reason, eq_data_->ls->compute_residual());
    output_data();
    END_TIMER("Mechanics zero time step");
}



void Elasticity::preallocate()
{
    // preallocate system matrix
	eq_data_->ls->start_allocation();
    stiffness_assembly_->assemble(eq_data_->dh_);
    rhs_assembly_->assemble(eq_data_->dh_);

    if (has_contact_) {
        assemble_constraint_matrix();
        assemble_tangential_jump_matrix();
    }
}



void Elasticity::next_time()
{
    time_->next_time();
    time_->view("MECH");
    
}



void Elasticity::solve_linear_system()
{
    START_TIMER("Mechanics step");
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

    START_TIMER("solve");
    if (has_tresca_) prepare_tresca();
    LinSys::SolveInfo si = eq_data_->ls->solve();
    if (has_tresca_)
        ASSERT_PERMANENT_GT(si.converged_reason, 0).error("Tresca solve did not converge; load step is not accepted.");
    MessageOut().fmt("[mech solver] lin. it: {}, reason: {}, residual: {}\n",
        		si.n_iterations, si.converged_reason, eq_data_->ls->compute_residual());
    END_TIMER("solve");
    END_TIMER("Mechanics step");
}




void Elasticity::output_data()
{
    START_TIMER("Mechanics output");

    // gather the solution from all processors
    eq_fields_->output_fields.set_time( this->time().step(), LimitSide::left);
    //if (eq_fields_->output_fields.is_field_output_time(eq_fields_->output_field, this->time().step()) )
    update_output_fields();
    eq_fields_->output_fields.output(this->time().step());

//     START_TIMER("Mechanics-balance");
//     balance_->calculate_instant(subst_idx, eq_data_->ls->get_solution());
//     balance_->output();
//     END_TIMER("Mechanics-balance");

    END_TIMER("Mechanics output");
}




void Elasticity::calculate_cumulative_balance()
{
//     if (balance_->cumulative())
//     {
//         balance_->calculate_cumulative(subst_idx, eq_data_->ls->get_solution());
//     }
}


void Elasticity::assemble_constraint_matrix()
{
    MatZeroEntries(eq_data_->constraint_matrix);
    VecZeroEntries(eq_data_->constraint_vec);
    constraint_assembly_->assemble(eq_data_->dh_);
    MatAssemblyBegin(eq_data_->constraint_matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(eq_data_->constraint_matrix, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(eq_data_->constraint_vec);
    VecAssemblyEnd(eq_data_->constraint_vec);
}


void Elasticity::assemble_tangential_jump_matrix()
{
    chkerr(MatZeroEntries(eq_data_->tangential_jump_matrix));
    tangential_jump_assembly_->assemble(eq_data_->dh_);
    chkerr(MatAssemblyBegin(eq_data_->tangential_jump_matrix, MAT_FINAL_ASSEMBLY));
    chkerr(MatAssemblyEnd(eq_data_->tangential_jump_matrix, MAT_FINAL_ASSEMBLY));
}


void Elasticity::prepare_tresca()
{
    // HM shares the time governor and may call solve_linear_system repeatedly.
    // Capture the last converged displacement only when physical time advances.
    if (time_->t() != friction_reference_time_) {
        chkerr(MatMult(eq_data_->tangential_jump_matrix, eq_data_->ls->get_solution(),
                       eq_data_->tresca_reference_jump));
        friction_reference_time_ = time_->t();
    }
    assemble_constraint_matrix();
    tresca_threshold_assembly_->assemble(eq_data_->dh_);
    chkerr(VecCopy(eq_data_->tresca_force_bound, eq_data_->tresca_force_lower_bound));
    chkerr(VecScale(eq_data_->tresca_force_lower_bound, -1.0));
}


