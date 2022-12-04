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
 * @file    concentration_model.cc
 * @brief   Discontinuous Galerkin method for equation of transport with dispersion.
 * @author  Jan Stebel
 */

#include "input/input_type.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"

#include "transport/transport_operator_splitting.hh"
#include "concentration_model.hh"
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

// Functor computing retardation coefficients:
// (1-porosity) * rock_density * sorption_coefficient * rock_density
struct fn_conc_retardation {
    inline Sclr operator() (Sclr csec, Sclr por_m, Sclr rho_s, Sclr sorp_mult) {
        return (1.-por_m)*rho_s*sorp_mult*csec;
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

// Functor computing diffusion coefficient (see notes in function)
struct fn_conc_diff_coef {
    inline Tens operator() (Tens diff_m, Vect velocity, Sclr v_norm, Sclr alphaL, Sclr alphaT,
    		FMT_UNUSED Sclr water_content, FMT_UNUSED  Sclr porosity, Sclr cross_sec) {

        // used tortuosity model dues to Millington and Quirk(1961) (should it be with power 10/3 ?)
        // for an overview of other models see: Chou, Wu, Zeng, Chang (2011)

    	// double tortuosity = pow(water_content, 7.0 / 3.0)/ (porosity * porosity);

        // result
        Tens K;

        // Note that the velocity vector is in fact the Darcian flux,
        // so we need not to multiply vnorm by water_content and cross_section.
	    //K = ((alphaL-alphaT) / vnorm) * K + (alphaT*vnorm + Dm*tortuosity*cross_cut*water_content) * arma::eye(3,3);

        if (fabs(v_norm) > 0) {
            /*
            for (int i=0; i<3; i++)
                for (int j=0; j<3; j++)
                    K(i,j) = (velocity[i]/vnorm)*(velocity[j]);
            */
            K = ((alphaL - alphaT) / v_norm) * arma::kron(velocity.t(), velocity);

            //arma::mat33 abs_diff_mat = arma::abs(K -  kk);
            //double diff = arma::min( arma::min(abs_diff_mat) );
            //ASSERT_PERMANENT(  diff < 1e-12 )(diff)(K)(kk);
        } else
            K.zeros();

        // Note that the velocity vector is in fact the Darcian flux,
        // so to obtain |v| we have to divide vnorm by porosity and cross_section.
        //K += alphaT*v_norm*arma::eye(3,3) + diff_m*(tortuosity*cross_sec*water_content);
        // We should name it: diffusivity_effective
        K += alphaT*v_norm*arma::eye(3,3) + diff_m*cross_sec;

        return K;
    }
};





ConcentrationTransportModel::ModelEqFields::ModelEqFields()
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
    *this+=disp_l
            .name("disp_l")
            .description("Longitudinal dispersivity in the liquid (for each substance).")
            .units( UnitSI().m() )
            .input_default("0.0")
            .flags_add( in_main_matrix & in_rhs );
    *this+=disp_t
            .name("disp_t")
            .description("Transverse dispersivity in the liquid (for each substance).")
            .units( UnitSI().m() )
            .input_default("0.0")
            .flags_add( in_main_matrix & in_rhs );
    *this+=diff_m
            .name("diff_m")
            .description("Molecular diffusivity in the liquid (for each substance).")
            .units( UnitSI().m(2).s(-1) )
            .input_default("0.0")
            .flags_add( in_main_matrix & in_rhs );
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


void ConcentrationTransportModel::ModelEqFields::initialize(FMT_UNUSED Input::Record transport_rec)
{
    // initialize multifield components
	sorption_coefficient.setup_components();
    sources_conc.setup_components();
    sources_density.setup_components();
    sources_sigma.setup_components();
    diff_m.setup_components();
    disp_l.setup_components();
    disp_t.setup_components();

    // create FieldModels
    v_norm.set(Model<3, FieldValue<3>::Scalar>::create(fn_conc_v_norm(), flow_flux), 0.0);
    mass_matrix_coef.set(Model<3, FieldValue<3>::Scalar>::create(fn_conc_mass_matrix(), cross_section, water_content), 0.0);
    retardation_coef.set(
        Model<3, FieldValue<3>::Scalar>::create_multi(fn_conc_retardation(), cross_section, porosity, rock_density, sorption_coefficient),
		0.0
    );
    sources_density_out.set(Model<3, FieldValue<3>::Scalar>::create_multi(fn_conc_sources_dens(), cross_section, sources_density), 0.0);
    sources_sigma_out.set(Model<3, FieldValue<3>::Scalar>::create_multi(fn_conc_sources_sigma(), cross_section, sources_sigma), 0.0);
    sources_conc_out.set(Model<3, FieldValue<3>::Scalar>::create_multi(fn_conc_sources_conc(), sources_conc), 0.0);
    std::vector<typename Field<3, FieldValue<3>::Vector>::FieldBasePtr> ad_coef_ptr_vec;
    for (unsigned int sbi=0; sbi<sorption_coefficient.size(); sbi++)
        ad_coef_ptr_vec.push_back( Model<3, FieldValue<3>::Vector>::create(fn_conc_ad_coef(), flow_flux) );
    advection_coef.set(ad_coef_ptr_vec, 0.0);
    diffusion_coef.set(
        Model<3, FieldValue<3>::Tensor>::create_multi(
            fn_conc_diff_coef(), diff_m, flow_flux, v_norm, disp_l, disp_t, water_content, porosity, cross_section
        ),
        0.0
    );
}



const Selection & ConcentrationTransportModel::ModelEqFields::get_bc_type_selection() {
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




IT::Selection ConcentrationTransportModel::ModelEqData::get_output_selection()
{
    // Return empty selection just to provide model specific selection name and description.
    // The fields are added by TransportDG using an auxiliary selection.
	return IT::Selection(
				std::string(ModelEqData::name()) + "_DG_output_fields",
				"Selection of output fields for Diffusive Solute Transport DG model.");
}


IT::Record ConcentrationTransportModel::get_input_type(const string &implementation, const string &description)
{
	return IT::Record(
				std::string(ModelEqData::name()) + "_" + implementation,
				description + " for solute transport.")
			.derive_from(ConcentrationTransportBase::get_input_type())
			.declare_key("solvent_density", IT::Double(0), IT::Default("1.0"),
					"Density of the solvent [ (($kg.m^{-3}$)) ].");
}


ConcentrationTransportModel::ConcentrationTransportModel(Mesh &mesh, const Input::Record &in_rec) :
		ConcentrationTransportBase(mesh, in_rec)
{}


void ConcentrationTransportModel::init_from_input(const Input::Record &in_rec)
{
	solvent_density_ = in_rec.val<double>("solvent_density");
}


ConcentrationTransportModel::~ConcentrationTransportModel()
{}


void ConcentrationTransportModel::set_balance_object(std::shared_ptr<Balance> balance)
{
	balance_ = balance;
	eq_data().subst_idx_ = balance_->add_quantities(eq_data().substances_.names());
}


void ConcentrationTransportModel::init_balance(FMT_UNUSED const Input::Record &in_rec) {}



