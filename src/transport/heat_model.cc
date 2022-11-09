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
 * @file    heat_model.cc
 * @brief   Discontinuous Galerkin method for equation of transport with dispersion.
 * @author  Jan Stebel
 */

#include "input/input_type.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
//#include "transport/transport_operator_splitting.hh"
#include "heat_model.hh"
#include "tools/unit_si.hh"
#include "coupling/balance.hh"
#include "fields/field_model.hh"
#include "fields/field_constant.hh"



using namespace std;
using namespace Input::Type;


/*******************************************************************************
 * Functors of FieldModels
 */
using Sclr = double;
using Vect = arma::vec3;
using Tens = arma::mat33;

// Functor computing velocity norm
struct fn_heat_v_norm {
    inline Sclr operator() (Vect vel) {
        return arma::norm(vel, 2);
    }
};

/**
 * Functor computing mass matrix coefficients:
 * cross_section * (porosity*fluid_density*fluid_heat_capacity + (1.-porosity)*solid_density*solid_heat_capacity)
 */
struct fn_heat_mass_matrix {
    inline Sclr operator() (Sclr csec, Sclr por, Sclr f_rho, Sclr f_c, Sclr s_rho, Sclr s_c) {
        return csec * (por*f_rho*f_c + (1.-por)*s_rho*s_c);
    }
};

/**
 * Functor computing sources density output:
 * cross_section * (porosity*fluid_thermal_source + (1-porosity)*solid_thermal_source)
 */
struct fn_heat_sources_dens {
    inline Sclr operator() (Sclr csec, Sclr por, Sclr f_source, Sclr s_source) {
        return csec * (por * f_source + (1.-por) * s_source);
    }
};

/**
 * Functor computing sources sigma output:
 * cross_section * (porosity*fluid_density*fluid_heat_capacity*fluid_heat_exchange_rate + (1-porosity)*solid_density*solid_heat_capacity*solid_thermal_source)
 */
struct fn_heat_sources_sigma {
    inline Sclr operator() (Sclr csec, Sclr por, Sclr f_rho, Sclr f_cap, Sclr f_sigma, Sclr s_rho, Sclr s_cap, Sclr s_sigma) {
        return csec * (por * f_rho * f_cap * f_sigma + (1.-por) * s_rho * s_cap * s_sigma);
    }
};

/**
 * Functor computing sources concentration output for positive heat_sources_sigma (otherwise return 0):
 * cross_section * (porosity*fluid_density*fluid_heat_capacity*fluid_heat_exchange_rate*fluid_ref_temperature + (1-porosity)*solid_density*solid_heat_capacity*solid_heat_exchange_rate*solid_ref_temperature)
 */
struct fn_heat_sources_conc {
    inline Sclr operator() (Sclr csec, Sclr por, Sclr f_rho, Sclr f_cap, Sclr f_sigma, Sclr f_temp, Sclr s_rho, Sclr s_cap, Sclr s_sigma, Sclr s_temp, Sclr sigma) {
        if (fabs(sigma) > numeric_limits<double>::epsilon())
            return csec * (por * f_rho * f_cap * f_sigma * f_temp + (1.-por) * s_rho * s_cap * s_sigma * s_temp);
        else {
            return 0;
        }
    }
};

/**
 * Functor computing advection coefficient
 * velocity * fluid_density * fluid_heat_capacity
 */
struct fn_heat_ad_coef {
    inline Vect operator() (Sclr f_rho, Sclr f_cap, Vect velocity) {
        return velocity * f_rho * f_cap;
    }
};

/**
 * Functor computing diffusion coefficient (see notes in function)
 */
struct fn_heat_diff_coef {
    inline Tens operator() (Vect velocity, Sclr v_norm, Sclr f_rho, Sclr disp_l, Sclr disp_t, Sclr f_cond, Sclr s_cond, Sclr c_sec, Sclr por) {
        // result
        Tens dif_coef;

        // dispersive part of thermal diffusion
        // Note that the velocity vector is in fact the Darcian flux,
        // so to obtain |v| we have to divide vnorm by porosity and cross_section.
        if ( fabs(v_norm) > 0 )
            for (int i=0; i<3; i++)
                for (int j=0; j<3; j++)
                    dif_coef(i,j) = ( (velocity(i)/v_norm) * (velocity(j)/v_norm) * (disp_l-disp_t) + disp_t*(i==j?1:0))*v_norm*f_rho*f_cond;
        else
            dif_coef.zeros();

        // conductive part of thermal diffusion
        dif_coef += c_sec * (por*f_cond + (1.-por)*s_cond) * arma::eye(3,3);
        return dif_coef;
    }
};








HeatTransferModel::ModelEqFields::ModelEqFields()
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
            .name("bc_temperature")
            .description("Boundary value of temperature.")
            .units( UnitSI().K() )
            .input_default("0.0")
            .flags_add(in_rhs);
	*this+=bc_flux
	        .disable_where(bc_type, { bc_dirichlet, bc_inflow })
			.name("bc_flux")
			.description("Flux in Neumann boundary condition.")
			.units( UnitSI().kg().m().s(-1).md() )
			.input_default("0.0")
			.flags_add(FieldFlag::in_rhs);
	*this+=bc_robin_sigma
	        .disable_where(bc_type, { bc_dirichlet, bc_inflow })
			.name("bc_robin_sigma")
			.description("Conductivity coefficient in Robin boundary condition.")
			.units( UnitSI().m(4).s(-1).md() )
			.input_default("0.0")
			.flags_add(FieldFlag::in_rhs & FieldFlag::in_main_matrix);

    *this+=init_condition
            .name("init_temperature")
            .description("Initial temperature.")
            .units( UnitSI().K() )
            .input_default("0.0");

    *this+=porosity
            .name("porosity")
            .description("Porosity.")
            .units( UnitSI::dimensionless() )
            .input_default("1.0")
            .flags_add(in_main_matrix & in_time_term)
			.set_limits(0.0);

    *this+=water_content
            .name("water_content")
            .units( UnitSI::dimensionless() )
            .input_default("1.0")
            .flags_add(input_copy & in_main_matrix & in_time_term);

    *this += flow_flux.name("flow_flux")
               .flags( FieldFlag::input_copy )
               .flags_add(in_time_term & in_main_matrix & in_rhs);

    *this+=fluid_density
            .name("fluid_density")
            .description("Density of fluid.")
            .units( UnitSI().kg().m(-3) )
            .input_default("1000")
            .flags_add(in_main_matrix & in_time_term);

    *this+=fluid_heat_capacity
            .name("fluid_heat_capacity")
            .description("Heat capacity of fluid.")
            .units( UnitSI::J() * UnitSI().kg(-1).K(-1) )
            .flags_add(in_main_matrix & in_time_term);

    *this+=fluid_heat_conductivity
            .name("fluid_heat_conductivity")
            .description("Heat conductivity of fluid.")
            .units( UnitSI::W() * UnitSI().m(-1).K(-1) )
            .flags_add(in_main_matrix)
			.set_limits(0.0);


    *this+=solid_density
            .name("solid_density")
            .description("Density of solid (rock).")
            .units( UnitSI().kg().m(-3) )
            .flags_add(in_time_term);

    *this+=solid_heat_capacity
            .name("solid_heat_capacity")
            .description("Heat capacity of solid (rock).")
            .units( UnitSI::J() * UnitSI().kg(-1).K(-1) )
            .flags_add(in_time_term);

    *this+=solid_heat_conductivity
            .name("solid_heat_conductivity")
            .description("Heat conductivity of solid (rock).")
            .units( UnitSI::W() * UnitSI().m(-1).K(-1) )
            .flags_add(in_main_matrix)
			.set_limits(0.0);

    *this+=disp_l
            .name("disp_l")
            .description("Longitudinal heat dispersivity in fluid.")
            .units( UnitSI().m() )
            .input_default("0.0")
            .flags_add(in_main_matrix);

    *this+=disp_t
            .name("disp_t")
            .description("Transverse heat dispersivity in fluid.")
            .units( UnitSI().m() )
            .input_default("0.0")
            .flags_add(in_main_matrix);

    *this+=fluid_thermal_source
            .name("fluid_thermal_source")
            .description("Density of thermal source in fluid.")
            .units( UnitSI::W() * UnitSI().m(-3) )
            .input_default("0.0")
            .flags_add(in_rhs);

    *this+=solid_thermal_source
            .name("solid_thermal_source")
            .description("Density of thermal source in solid.")
            .units( UnitSI::W() * UnitSI().m(-3) )
            .input_default("0.0")
            .flags_add(in_rhs);

    *this+=fluid_heat_exchange_rate
            .name("fluid_heat_exchange_rate")
            .description("Heat exchange rate of source in fluid.")
            .units( UnitSI().s(-1) )
            .input_default("0.0")
            .flags_add(in_rhs);

    *this+=solid_heat_exchange_rate
            .name("solid_heat_exchange_rate")
            .description("Heat exchange rate of source in solid.")
            .units( UnitSI().s(-1) )
            .input_default("0.0")
            .flags_add(in_rhs);

    *this+=fluid_ref_temperature
            .name("fluid_ref_temperature")
            .description("Reference temperature of source in fluid.")
            .units( UnitSI().K() )
            .input_default("0.0")
            .flags_add(in_rhs);

    *this+=solid_ref_temperature
            .name("solid_ref_temperature")
            .description("Reference temperature in solid.")
            .units( UnitSI().K() )
            .input_default("0.0")
            .flags_add(in_rhs);

    *this+=cross_section
            .name("cross_section")
            .units( UnitSI().m(3).md() )
            .flags(input_copy & in_time_term & in_main_matrix);

    *this+=output_field
            .name("temperature")
            .description("Temperature solution.")
            .units( UnitSI().K() )
            .flags(equation_result);


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

    *this += sources_conc_out.name("sources_conc_out")
            .description("Concentration sources output.")
            .input_default("0.0")
            .units( UnitSI().kg().m(-3) );

    *this += sources_density_out.name("sources_density_out")
            .description("Concentration sources output - density of substance source, only positive part is used..")
            .input_default("0.0")
            .units( UnitSI().kg().s(-1).md() );

    *this += sources_sigma_out.name("sources_sigma_out")
            .description("Concentration sources - Robin type, in_flux = sources_sigma * (sources_conc - mobile_conc).")
            .input_default("0.0")
            .units( UnitSI().s(-1).m(3).md() );

    *this += advection_coef.name("advection_coef")
            .description("Advection coefficients model.")
            .input_default("0.0")
            .units( UnitSI().m().s(-1) );

    *this += diffusion_coef.name("diffusion_coef")
            .description("Diffusion coefficients model.")
            .input_default("0.0")
            .units( UnitSI().m(2).s(-1) );
}


void HeatTransferModel::ModelEqFields::initialize()
{
    // create FieldModels
    v_norm.set(Model<3, FieldValue<3>::Scalar>::create(fn_heat_v_norm(), flow_flux), 0.0);
    mass_matrix_coef.set(
        Model<3, FieldValue<3>::Scalar>::create(
            fn_heat_mass_matrix(), cross_section, porosity, fluid_density, fluid_heat_capacity, solid_density, solid_heat_capacity
	    ),
	    0.0
	);
    retardation_coef.set(std::make_shared< FieldConstant<3, FieldValue<3>::Scalar> >(), 0.0);
    sources_density_out.set(
        Model<3, FieldValue<3>::Scalar>::create_multi(fn_heat_sources_dens(), cross_section, porosity, fluid_thermal_source, solid_thermal_source),
        0.0
    );
    sources_sigma_out.set(
        Model<3, FieldValue<3>::Scalar>::create_multi(
            fn_heat_sources_sigma(), cross_section, porosity, fluid_density, fluid_heat_capacity, fluid_heat_exchange_rate, solid_density,
            solid_heat_capacity, solid_heat_exchange_rate
        ),
        0.0
    );
    sources_conc_out.set(
        Model<3, FieldValue<3>::Scalar>::create_multi(
            fn_heat_sources_conc(), cross_section, porosity, fluid_density, fluid_heat_capacity, fluid_heat_exchange_rate, fluid_ref_temperature,
            solid_density, solid_heat_capacity, solid_heat_exchange_rate, solid_ref_temperature, sources_sigma_out
        ),
        0.0
    );
    advection_coef.set(Model<3, FieldValue<3>::Vector>::create(fn_heat_ad_coef(), fluid_density, fluid_heat_capacity, flow_flux), 0.0);
    diffusion_coef.set(
        Model<3, FieldValue<3>::Tensor>::create(
            fn_heat_diff_coef(), flow_flux, v_norm, fluid_density, disp_l, disp_t, fluid_heat_conductivity, solid_heat_conductivity, cross_section, porosity
        ),
        0.0
    );
}


const Selection & HeatTransferModel::ModelEqFields::get_bc_type_selection() {
	return Selection("Heat_BC_Type", "Types of boundary conditions for heat transfer model.")
            .add_value(bc_inflow, "inflow",
          		  "Default heat transfer boundary condition.\n"
          		  "On water inflow (($(q_w \\le 0)$)), total energy flux is given by the reference temperature 'bc_temperature'. "
          		  "On water outflow we prescribe zero diffusive flux, "
          		  "i.e. the energy flows out only due to advection.")
            .add_value(bc_dirichlet, "dirichlet",
          		  "Dirichlet boundary condition (($T = T_D $)).\n"
          		  "The prescribed temperature (($T_D$)) is specified by the field 'bc_temperature'.")
            .add_value(bc_total_flux, "total_flux",
          		  "Total energy flux boundary condition.\n"
          		  "The prescribed incoming total flux can have the general form (($\\delta(f_N+\\sigma_R(T_R-T) )$)), "
          		  "where the absolute flux (($f_N$)) is specified by the field 'bc_flux', "
          		  "the transition parameter (($\\sigma_R$)) by 'bc_robin_sigma', "
          		  "and the reference temperature (($T_R$)) by 'bc_temperature'.")
            .add_value(bc_diffusive_flux, "diffusive_flux",
          		  "Diffusive flux boundary condition.\n"
          		  "The prescribed incoming energy flux due to diffusion can have the general form (($\\delta(f_N+\\sigma_R(T_R-T) )$)), "
          		  "where the absolute flux (($f_N$)) is specified by the field 'bc_flux', "
          		  "the transition parameter (($\\sigma_R$)) by 'bc_robin_sigma', "
          		  "and the reference temperature (($T_R$)) by 'bc_temperature'.")
			  .close();
}


IT::Selection HeatTransferModel::ModelEqData::get_output_selection()
{
    // Return empty selection just to provide model specific selection name and description.
    // The fields are added by TransportDG using an auxiliary selection.
	return IT::Selection(
				std::string(ModelEqData::name()) + "_DG_output_fields",
				"Selection of output fields for Heat Transfer DG model.");
}



IT::Record HeatTransferModel::get_input_type(const string &implementation, const string &description)
{
	return IT::Record(
				std::string(ModelEqData::name()) + "_" + implementation,
				description + " for heat transfer.")
			.derive_from(AdvectionProcessBase::get_input_type())
			.copy_keys(EquationBase::record_template())
			.declare_key("balance", Balance::get_input_type(), Default("{}"),
					"Settings for computing balance.")
			.declare_key("output_stream", OutputTime::get_input_type(), Default("{}"),
					"Parameters of output stream.");
}


HeatTransferModel::HeatTransferModel(Mesh &mesh, const Input::Record in_rec) :
		AdvectionProcessBase(mesh, in_rec)
{
	time_ = new TimeGovernor(in_rec.val<Input::Record>("time"));
	ASSERT( time_->is_default() == false ).error("Missing key 'time' in Heat_AdvectionDiffusion_DG.");

    output_stream_ = OutputTime::create_output_stream("heat", in_rec.val<Input::Record>("output_stream"), time().get_unit_conversion());
    //output_stream_->add_admissible_field_names(in_rec.val<Input::Array>("output_fields"));
}


void HeatTransferModel::init_balance(const Input::Record &in_rec)
{
    balance_ = std::make_shared<Balance>("energy", mesh_);
    balance_->init_from_input(in_rec.val<Input::Record>("balance"), *time_);
    // initialization of balance object
    eq_data().subst_idx_ = {balance_->add_quantity("energy")};
    balance_->units(UnitSI().m(2).kg().s(-2));
}


HeatTransferModel::~HeatTransferModel()
{}




