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
 * @file    conc_dispersion_model.cc
 * @brief   Discontinuous Galerkin method for equation of transport with dispersion.
 * @author  Jan Brezina
 *
 * Transport model with:
 * - full dispersion tensor prescribed by 6 tensor fields (UGLY).
 * - support for switching both time and convection term to zero (additional multiplicative factor)
 * - meant for experimental homogenisation
 */

#include "input/input_type.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"

#include "transport/transport_operator_splitting.hh"
#include "conc_dispersion_model.hh"
#include "tools/unit_si.hh"
#include "coupling/balance.hh"
#include "fields/field_model.hh"



using namespace std;
using namespace Input::Type;


/*******************************************************************************
 * Functors of FieldModels
 */
using Sclr = double;
using Vect = arma::vec3;
using Tens = arma::mat33;

// Functor computing velocity norm
struct fn_conc_v_norm {
    inline Sclr operator() (Vect vel) {
        return arma::norm(vel, 2);
    }
};

// Functor computing mass matrix coefficients (cross_section * water_content)
struct fn_conc_mass_matrix {
	inline Sclr operator() (Sclr csec, Sclr wcont) {
        return csec * wcont;
    }
};

struct fn_conc_mass_matrix_static {
	inline Sclr operator() () {
        return 0;
    }
};

// Functor computing retardation coefficients:
// (1-porosity) * rock_density * sorption_coefficient * rock_density
struct fn_conc_retardation {
    inline Sclr operator() (Sclr csec, Sclr por_m, Sclr rho_s, Sclr sorp_mult) {
        return (1.-por_m)*rho_s*sorp_mult*csec;
    }
};

struct fn_conc_retardation_static {
    inline Sclr operator() () {
        return 0;
    }
};

// Functor computing sources density output (cross_section * sources_density)
struct fn_conc_sources_dens {
    inline Sclr operator() (Sclr csec, Sclr sdens) {
        return csec * sdens;
    }
};

// Functor computing sources sigma output (cross_section * sources_sigma)
struct fn_conc_sources_sigma {
    inline Sclr operator() (Sclr csec, Sclr ssigma) {
        return csec * ssigma;
    }
};

// Functor computing sources concentration output (sources_conc)
struct fn_conc_sources_conc {
    inline Sclr operator() (Sclr sconc) {
        return sconc;
    }
};

// Functor computing advection coefficient (velocity)
struct fn_conc_ad_coef {
    inline Vect operator() (Vect velocity) {
        return velocity;
    }
};

struct fn_conc_ad_coef_static {
    inline Vect operator() () {
        return arma::zeros(3);
    }
};

// Functor computing diffusion coefficient (see notes in function)

// Fourth order tnesor in partially Voigt notation.
// In particular we use is to specify full dispersion tensor.

struct fn_conc_diff_coef_full {
    inline Tens operator() (Tens duxx, Tens duyy, Tens duzz, Tens duyz, Tens duxz, Tens duxy, Vect velocity, Sclr v_norm) {

        // result
        auto v = velocity;
        Tens K =   (duxx * v[0] * v[0] + duyy * v[1] * v[1] + duzz * v[2] * v[2]
                 + 2 * ( duyz * v[1] * v[2] + duxz * v[0] * v[2] + duxy * v[0] * v[1])) / v_norm;
        return K;
    }
};





ConcDispersionModel::ModelEqFields::ModelEqFields()
: TransportEqFields()
{
    *this+=bc_type
            .name("bc_type")
            .description(
            "Type of boundary condition.")
            .units( UnitSI::dimensionless() )
            .input_default("\"inflow\"")
            .input_selection( get_bc_type_selection() )
            .flags_add(FieldFlag::in_rhs & FieldFlag::in_main_matrix);
    *this+=bc_dirichlet_value
            .name("bc_conc")
            .units( UnitSI().kg().m(-3) )
            .description("Dirichlet boundary condition (for each substance).")
            .input_default("0.0")
            .flags_add( in_rhs );
	*this+=bc_flux
	        .disable_where(bc_type, { bc_inflow, bc_dirichlet })
			.name("bc_flux")
			.description("Flux in Neumann boundary condition.")
			.units( UnitSI().kg().m().s(-1).md() )
			.input_default("0.0")
			.flags_add(FieldFlag::in_rhs);
	*this+=bc_robin_sigma
	        .disable_where(bc_type, { bc_inflow, bc_dirichlet })
			.name("bc_robin_sigma")
			.description("Conductivity coefficient in Robin boundary condition.")
			.units( UnitSI().m(4).s(-1).md() )
			.input_default("0.0")
			.flags_add(FieldFlag::in_rhs & FieldFlag::in_main_matrix);
    *this+=init_condition
            .name("init_conc")
            .units( UnitSI().kg().m(-3) )
            .description("Initial values for concentration of substances.")
            .input_default("0.0");
    *this+=disp_uxx
            .name("dispersion_uxx")
            .description("Dispersion 4-th order tensor component ux D ux. Common for all substances right now.")
            .units( UnitSI().m() )  // ??
            .input_default("1.0")
            .flags_add( in_main_matrix );
    *this+=disp_uyy
            .name("dispersion_uyy")
            .description("Dispersion 4-th order tensor component uy D uy. Common for all substances right now.")
            .units( UnitSI().m() )  // ??
            .input_default("1.0")
            .flags_add( in_main_matrix );
    *this+=disp_uzz
            .name("dispersion_uzz")
            .description("Dispersion 4-th order tensor component uz D uz. Common for all substances right now.")
            .units( UnitSI().m() )  // ??
            .input_default("1.0")
            .flags_add( in_main_matrix );
    *this+=disp_uyz
            .name("dispersion_uyz")
            .description("Dispersion 4-th order tensor component uy D uz. Common for all substances right now.")
            .units( UnitSI().m() )  // ??
            .input_default("0.0")
            .flags_add( in_main_matrix );
    *this+=disp_uxz
            .name("dispersion_uxz")
            .description("Dispersion 4-th order tensor component ux D uz. Common for all substances right now.")
            .units( UnitSI().m() )  // ??
            .input_default("1.0")
            .flags_add( in_main_matrix );
    *this+=disp_uxy
            .name("dispersion_uxy")
            .description("Dispersion 4-th order tensor component ux D uy. Common for all substances right now.")
            .units( UnitSI().m() )  // ??
            .input_default("1.0")
            .flags_add( in_main_matrix );

    *this+=rock_density
    		.name("rock_density")
			.description("Rock matrix density.")
			.units(UnitSI().kg().m(-3))
			.input_default("0.0")
			.flags_add( in_time_term );
    *this+=sorption_coefficient
    		.name("sorption_coefficient")
			.description("Coefficient of linear sorption.")
			.units(UnitSI().m(3).kg(-1))
			.input_default("0.0")
			.flags_add( in_time_term );

	*this+=output_field
	        .name("conc")
            .description("Concentration solution.")
	        .units( UnitSI().kg().m(-3) )
	        .flags( equation_result );


	// initiaization of FieldModels
    *this += v_norm.name("v_norm")
            .description("Velocity norm field.")
            .input_default("0.0")
            .units( UnitSI().m().s(-1) );

    *this += mass_matrix_coef.name("mass_matrix_coef")
            .description("Matrix coefficients computed by model in mass assemblation.")
            .input_default("0.0")
            .units( UnitSI().m(3).md() );

    *this += retardation_coef.name("retardation_coef")
            .description("Retardation coefficients computed by model in mass assemblation.")
            .input_default("0.0")
            .units( UnitSI().m(3).md() );
    *this += sources_density_out.name("sources_density_out")
            .description("Concentration sources output - density of substance source, only positive part is used..")
            .input_default("0.0")
            .units( UnitSI().kg().s(-1).md() );

    *this += sources_sigma_out.name("sources_sigma_out")
            .description("Concentration sources - Robin type, in_flux = sources_sigma * (sources_conc - mobile_conc).")
            .input_default("0.0")
            .units( UnitSI().s(-1).m(3).md() );

    *this += sources_conc_out.name("sources_conc_out")
            .description("Concentration sources output.")
            .input_default("0.0")
            .units( UnitSI().kg().m(-3) );

    *this += advection_coef.name("advection_coef")
            .description("Advection coefficients model.")
            .input_default("0.0")
            .units( UnitSI().m().s(-1) );

    *this += diffusion_coef.name("diffusion_coef")
            .description("Diffusion coefficients model.")
            .input_default("0.0")
            .units( UnitSI().m(2).s(-1) );

}


void ConcDispersionModel::ModelEqFields::initialize(FMT_UNUSED Input::Record transport_rec)
{
	bool static_flag = transport_rec.val<bool>("static_model");

	// initialize multifield components
	sorption_coefficient.setup_components();
    sources_conc.setup_components();
    sources_density.setup_components();
    sources_sigma.setup_components();

    // create FieldModels
    v_norm.set(Model<3, FieldValue<3>::Scalar>::create(
    		fn_conc_v_norm(), flow_flux), 0.0);
    sources_density_out.set(Model<3, FieldValue<3>::Scalar>::create_multi(
    		fn_conc_sources_dens(), cross_section, sources_density), 0.0);
    sources_sigma_out.set(Model<3, FieldValue<3>::Scalar>::create_multi(
    		fn_conc_sources_sigma(), cross_section, sources_sigma), 0.0);
    sources_conc_out.set(Model<3, FieldValue<3>::Scalar>::create_multi(
    		fn_conc_sources_conc(), sources_conc), 0.0);

    if (static_flag) {
    	advection_coef.set(Model<3, FieldValue<3>::Vector>::create_multi(
    		fn_conc_ad_coef_static()), 0.0);
    	mass_matrix_coef.set(Model<3, FieldValue<3>::Scalar>::create(
    		fn_conc_mass_matrix_static()), 0.0);
        retardation_coef.set(Model<3, FieldValue<3>::Scalar>::create_multi(
        		fn_conc_retardation_static()), 0.0);

    } else {
    	advection_coef.set(Model<3, FieldValue<3>::Vector>::create_multi(
    		fn_conc_ad_coef(), flow_flux), 0.0);
    	mass_matrix_coef.set(Model<3, FieldValue<3>::Scalar>::create(
    		fn_conc_mass_matrix(), cross_section, water_content), 0.0);
        retardation_coef.set(Model<3, FieldValue<3>::Scalar>::create_multi(
        		fn_conc_retardation(), cross_section, porosity, rock_density, sorption_coefficient), 0.0);
    }
    diffusion_coef.set(Model<3, FieldValue<3>::Tensor>::create_multi(
    		fn_conc_diff_coef_full(), disp_uxx, disp_uyy, disp_uzz, disp_uyz, disp_uxz, disp_uxy,
			flow_flux, v_norm), 0.0);
}



const Selection & ConcDispersionModel::ModelEqFields::get_bc_type_selection() {
	return Selection("Solute_AdvectionDiffusion_BC_Type", "Types of boundary conditions for advection-diffusion solute transport model.")
              .add_value(bc_inflow, "inflow",
            		  "Default transport boundary condition.\n"
            		  "On water inflow (($(q_w \\le 0)$)), total flux is given by the reference concentration 'bc_conc'. "
            		  "On water outflow we prescribe zero diffusive flux, "
            		  "i.e. the mass flows out only due to advection.")
              .add_value(bc_dirichlet, "dirichlet",
            		  "Dirichlet boundary condition (($ c = c_D $)).\n"
            		  "The prescribed concentration (($c_D$)) is specified by the field 'bc_conc'.")
              .add_value(bc_total_flux, "total_flux",
            		  "Total mass flux boundary condition.\n"
            		  "The prescribed total incoming flux can have the general form (($\\delta(f_N+\\sigma_R(c_R-c) )$)), "
            		  "where the absolute flux (($f_N$)) is specified by the field 'bc_flux', "
            		  "the transition parameter (($\\sigma_R$)) by 'bc_robin_sigma', "
            		  "and the reference concentration (($c_R$)) by 'bc_conc'.")
              .add_value(bc_diffusive_flux, "diffusive_flux",
            		  "Diffusive flux boundary condition.\n"
            		  "The prescribed incoming mass flux due to diffusion can have the general form (($\\delta(f_N+\\sigma_R(c_R-c) )$)), "
            		  "where the absolute flux (($f_N$)) is specified by the field 'bc_flux', "
            		  "the transition parameter (($\\sigma_R$)) by 'bc_robin_sigma', "
            		  "and the reference concentration (($c_R$)) by 'bc_conc'.")
			  .close();
}




IT::Selection ConcDispersionModel::ModelEqData::get_output_selection()
{
    // Return empty selection just to provide model specific selection name and description.
    // The fields are added by TransportDG using an auxiliary selection.
	return IT::Selection(
				std::string(ModelEqData::name()) + "_DG_output_fields",
				"Selection of output fields for Diffusive Solute Transport DG model.");
}


IT::Record ConcDispersionModel::get_input_type(const string &implementation, const string &description)
{
	return IT::Record(
				std::string(ModelEqData::name()) + "_" + implementation,
				description + " for solute transport.")
			.derive_from(ConcentrationTransportBase::get_input_type())
			.declare_key("static_model", IT::Bool(), IT::Default("false"),
					"Forces zero time and advection term.")
			.declare_key("solvent_density", IT::Double(0), IT::Default("1.0"),
					"Density of the solvent [ (($kg.m^{-3}$)) ].");
}


ConcDispersionModel::ConcDispersionModel(Mesh &mesh, const Input::Record &in_rec) :
		ConcentrationTransportBase(mesh, in_rec)
{}


void ConcDispersionModel::init_from_input(const Input::Record &in_rec)
{
	solvent_density_ = in_rec.val<double>("solvent_density");
}


ConcDispersionModel::~ConcDispersionModel()
{}


void ConcDispersionModel::set_balance_object(std::shared_ptr<Balance> balance)
{
	balance_ = balance;
	eq_data().subst_idx_ = balance_->add_quantities(eq_data().substances_.names());
}


void ConcDispersionModel::init_balance(FMT_UNUSED const Input::Record &in_rec) {}



