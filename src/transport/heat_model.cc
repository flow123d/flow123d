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
#include "transport/transport_operator_splitting.hh"
#include "heat_model.hh"
#include "fields/unit_si.hh"
#include "coupling/balance.hh"



using namespace std;
using namespace Input::Type;








const Selection & HeatTransferModel::ModelEqData::get_bc_type_selection() {
	return Selection("HeatTransfer_BC_Type", "Types of boundary conditions for heat transfer model.")
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


HeatTransferModel::ModelEqData::ModelEqData()
{
    *this+=bc_type
            .name("bc_type")
            .description(
            "Type of boundary condition.")
            .units( UnitSI::dimensionless() )
            .input_default("\"inflow\"")
            .input_selection( &get_bc_type_selection() )
            .flags_add(FieldFlag::in_rhs & FieldFlag::in_main_matrix);
    *this+=bc_dirichlet_value
            .name("bc_temperature")
            .description("Boundary value of temperature.")
            .units( UnitSI().K() )
            .input_default("0.0")
            .flags_add(in_rhs);
	*this+=bc_flux
			.name("bc_flux")
			.description("Flux in Neumann boundary condition.")
			.units( UnitSI().kg().m().s(-1).md() )
			.input_default("0.0")
			.flags_add(FieldFlag::in_rhs);
	*this+=bc_robin_sigma
			.name("bc_robin_sigma")
			.description("Conductivity coefficient in Robin boundary condition.")
			.units( UnitSI().m(4).s(-1).md() )
			.input_default("0.0")
			.flags_add(FieldFlag::in_rhs & FieldFlag::in_main_matrix);

    *this+=init_temperature
            .name("init_temperature")
            .description("Initial temperature.")
            .units( UnitSI().K() )
            .input_default("0.0");

    *this+=porosity
            .name("porosity")
            .description("Porosity.")
            .units( UnitSI::dimensionless() )
            .input_default("1.0")
            .flags_add(in_main_matrix & in_time_term);

    *this+=fluid_density
            .name("fluid_density")
            .description("Density of fluid.")
            .units( UnitSI().kg().m(-3) )
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
            .flags_add(in_main_matrix);


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
            .flags_add(in_main_matrix);

    *this+=disp_l
            .name("disp_l")
            .description("Longitudal heat dispersivity in fluid.")
            .units( UnitSI().m() )
            .input_default("0.0")
            .flags_add(in_main_matrix);

    *this+=disp_t
            .name("disp_t")
            .description("Transversal heat dispersivity in fluid.")
            .units( UnitSI().m() )
            .input_default("0.0")
            .flags_add(in_main_matrix);

    *this+=fluid_thermal_source
            .name("fluid_thermal_source")
            .description("Thermal source density in fluid.")
            .units( UnitSI::W() * UnitSI().m(-3) )
            .input_default("0.0")
            .flags_add(in_rhs);

    *this+=solid_thermal_source
            .name("solid_thermal_source")
            .description("Thermal source density in solid.")
            .units( UnitSI::W() * UnitSI().m(-3) )
            .input_default("0.0")
            .flags_add(in_rhs);

    *this+=fluid_heat_exchange_rate
            .name("fluid_heat_exchange_rate")
            .description("Heat exchange rate in fluid.")
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
            .units( UnitSI().K() )
            .flags(equation_result);
}



IT::Record HeatTransferModel::get_input_type(const string &implementation, const string &description)
{
	return IT::Record(
				std::string(ModelEqData::name()) + "_" + implementation,
				description + " for heat transfer.")
			.derive_from(AdvectionProcessBase::get_input_type())
			.declare_key("time", TimeGovernor::get_input_type(), Default::obligatory(),
					"Time governor setting for the secondary equation.")
			.declare_key("balance", Balance::get_input_type(), Default::obligatory(),
					"Settings for computing balance.")
			.declare_key("output_stream", OutputTime::get_input_type(), Default::obligatory(),
					"Parameters of output stream.");
}


IT::Selection HeatTransferModel::ModelEqData::get_output_selection()
{
    // Return empty selection just to provide model specific selection name and description.
    // The fields are added by TransportDG using an auxiliary selection.
	return IT::Selection(
				std::string(ModelEqData::name()) + "_DG_output_fields",
				"Selection of output fields for Heat Transfer DG model.");
}


HeatTransferModel::HeatTransferModel(Mesh &mesh, const Input::Record in_rec) :
		AdvectionProcessBase(mesh, in_rec),
		flux_changed(true),
		mh_dh(nullptr)
{
	time_ = new TimeGovernor(in_rec.val<Input::Record>("time"));
	substances_.initialize({""});

    output_stream_ = OutputTime::create_output_stream(in_rec.val<Input::Record>("output_stream"));
    output_stream_->add_admissible_field_names(in_rec.val<Input::Array>("output_fields"));

    // initialization of balance object
    Input::Iterator<Input::Record> it = in_rec.find<Input::Record>("balance");
    if (it->val<bool>("balance_on"))
    {
    	balance_ = boost::make_shared<Balance>("energy", mesh_, *it);
    	subst_idx = {balance_->add_quantity("energy")};
    	balance_->units(UnitSI().m(2).kg().s(-2));
    }
}


void HeatTransferModel::set_components(SubstanceList &substances, const Input::Record &in_rec)
{
	substances.initialize({""});
}

void HeatTransferModel::output_data()
{
	output_stream_->write_time_frame();
	if (balance_ != nullptr)
	{
		calculate_instant_balance();
		balance_->output(time_->t());
	}
}


void HeatTransferModel::compute_mass_matrix_coefficient(const std::vector<arma::vec3 > &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector<double> &mm_coef)
{
	vector<double> elem_csec(point_list.size()),
			por(point_list.size()),
			f_rho(point_list.size()),
			s_rho(point_list.size()),
			f_c(point_list.size()),
			s_c(point_list.size());

	data().cross_section.value_list(point_list, ele_acc, elem_csec);
	data().porosity.value_list(point_list, ele_acc, por);
	data().fluid_density.value_list(point_list, ele_acc, f_rho);
	data().fluid_heat_capacity.value_list(point_list, ele_acc, f_c);
	data().solid_density.value_list(point_list, ele_acc, s_rho);
	data().solid_heat_capacity.value_list(point_list, ele_acc, s_c);

	for (unsigned int i=0; i<point_list.size(); i++)
		mm_coef[i] = elem_csec[i]*(por[i]*f_rho[i]*f_c[i] + (1.-por[i])*s_rho[i]*s_c[i]);
}


void HeatTransferModel::compute_advection_diffusion_coefficients(const std::vector<arma::vec3 > &point_list,
		const std::vector<arma::vec3> &velocity,
		const ElementAccessor<3> &ele_acc,
		std::vector<std::vector<arma::vec3> > &ad_coef,
		std::vector<std::vector<arma::mat33> > &dif_coef)
{
	const unsigned int qsize = point_list.size();
	std::vector<double> f_rho(qsize), f_cap(qsize), f_cond(qsize),
			s_cond(qsize), por(qsize), csection(qsize), disp_l(qsize), disp_t(qsize);

	data().fluid_density.value_list(point_list, ele_acc, f_rho);
	data().fluid_heat_capacity.value_list(point_list, ele_acc, f_cap);
	data().fluid_heat_conductivity.value_list(point_list, ele_acc, f_cond);
	data().solid_heat_conductivity.value_list(point_list, ele_acc, s_cond);
	data().disp_l.value_list(point_list, ele_acc, disp_l);
	data().disp_t.value_list(point_list, ele_acc, disp_t);
	data().porosity.value_list(point_list, ele_acc, por);
	data().cross_section.value_list(point_list, ele_acc, csection);

	for (unsigned int k=0; k<qsize; k++) {
		ad_coef[0][k] = velocity[k]*f_rho[k]*f_cap[k];

		// dispersive part of thermal diffusion
		// Note that the velocity vector is in fact the Darcian flux,
		// so to obtain |v| we have to divide vnorm by porosity and cross_section.
		double vnorm = arma::norm(velocity[k], 2);
		if (fabs(vnorm) > 0)
			for (int i=0; i<3; i++)
				for (int j=0; j<3; j++)
					dif_coef[0][k](i,j) = ((velocity[k][i]/vnorm)*(velocity[k][j]/vnorm)*(disp_l[k]-disp_t[k]) + disp_t[k]*(i==j?1:0))
											*vnorm*f_rho[k]*f_cond[k];
		else
			dif_coef[0][k].zeros();

		// conductive part of thermal diffusion
		dif_coef[0][k] += csection[k]*(por[k]*f_cond[k] + (1.-por[k])*s_cond[k])*arma::eye(3,3);
	}
}


void HeatTransferModel::compute_init_cond(const std::vector<arma::vec3> &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector< arma::vec > &init_values)
{
	vector<double> init_value(point_list.size());
	data().init_temperature.value_list(point_list, ele_acc, init_value);
	for (unsigned int i=0; i<point_list.size(); i++)
		init_values[i] = init_value[i];
}


void HeatTransferModel::get_bc_type(const ElementAccessor<3> &ele_acc,
			arma::uvec &bc_types)
{
	// Currently the bc types for HeatTransfer are numbered in the same way as in TransportDG.
	// In general we should use some map here.
	bc_types = { data().bc_type.value(ele_acc.centre(), ele_acc) };
}


void HeatTransferModel::get_flux_bc_data(const std::vector<arma::vec3> &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector< arma::vec > &bc_flux,
		std::vector< arma::vec > &bc_sigma,
		std::vector< arma::vec > &bc_ref_value)
{
	data().bc_flux.value_list(point_list, ele_acc, bc_flux);
	data().bc_robin_sigma.value_list(point_list, ele_acc, bc_sigma);
	data().bc_dirichlet_value.value_list(point_list, ele_acc, bc_ref_value);
	
	// Change sign in bc_flux since internally we work with outgoing fluxes.
	for (auto f : bc_flux) f = -f;
}

void HeatTransferModel::get_flux_bc_sigma(const std::vector<arma::vec3> &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector< arma::vec > &bc_sigma)
{
	data().bc_robin_sigma.value_list(point_list, ele_acc, bc_sigma);
}


void HeatTransferModel::compute_source_coefficients(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<arma::vec> &sources_value,
			std::vector<arma::vec> &sources_density,
			std::vector<arma::vec> &sources_sigma)
{
	const unsigned int qsize = point_list.size();
	std::vector<double> por(qsize), csection(qsize), f_rho(qsize), s_rho(qsize), f_cap(qsize), s_cap(qsize),
			f_source(qsize), s_source(qsize), f_sigma(qsize), s_sigma(qsize), f_temp(qsize), s_temp(qsize);
	data().porosity.value_list(point_list, ele_acc, por);
	data().cross_section.value_list(point_list, ele_acc, csection);
	data().fluid_density.value_list(point_list, ele_acc, f_rho);
	data().solid_density.value_list(point_list, ele_acc, s_rho);
	data().fluid_heat_capacity.value_list(point_list, ele_acc, f_cap);
	data().solid_heat_capacity.value_list(point_list, ele_acc, s_cap);
	data().fluid_thermal_source.value_list(point_list, ele_acc, f_source);
	data().solid_thermal_source.value_list(point_list, ele_acc, s_source);
	data().fluid_heat_exchange_rate.value_list(point_list, ele_acc, f_sigma);
	data().solid_heat_exchange_rate.value_list(point_list, ele_acc, s_sigma);
	data().fluid_ref_temperature.value_list(point_list, ele_acc, f_temp);
	data().solid_ref_temperature.value_list(point_list, ele_acc, s_temp);

	for (unsigned int k=0; k<point_list.size(); k++)
	{
		sources_density[k].resize(1);
		sources_sigma[k].resize(1);
		sources_value[k].resize(1);

		sources_density[k][0] = csection[k]*(por[k]*f_source[k] + (1.-por[k])*s_source[k]);
		sources_sigma[k][0] = csection[k]*(por[k]*f_rho[k]*f_cap[k]*f_sigma[k] + (1.-por[k])*s_rho[k]*s_cap[k]*s_sigma[k]);
		if (fabs(sources_sigma[k][0]) > numeric_limits<double>::epsilon())
			sources_value[k][0] = csection[k]*(por[k]*f_rho[k]*f_cap[k]*f_sigma[k]*f_temp[k]
		                   + (1.-por[k])*s_rho[k]*s_cap[k]*s_sigma[k]*s_temp[k])/sources_sigma[k][0];
		else
			sources_value[k][0] = 0;
	}
}


void HeatTransferModel::compute_sources_sigma(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<arma::vec> &sources_sigma)
{
	const unsigned int qsize = point_list.size();
	std::vector<double> por(qsize), csection(qsize), f_rho(qsize), s_rho(qsize), f_cap(qsize), s_cap(qsize),
			f_source(qsize), s_source(qsize), f_sigma(qsize), s_sigma(qsize), f_temp(qsize), s_temp(qsize);
	data().porosity.value_list(point_list, ele_acc, por);
	data().cross_section.value_list(point_list, ele_acc, csection);
	data().fluid_density.value_list(point_list, ele_acc, f_rho);
	data().solid_density.value_list(point_list, ele_acc, s_rho);
	data().fluid_heat_capacity.value_list(point_list, ele_acc, f_cap);
	data().solid_heat_capacity.value_list(point_list, ele_acc, s_cap);
	data().fluid_heat_exchange_rate.value_list(point_list, ele_acc, f_sigma);
	data().solid_heat_exchange_rate.value_list(point_list, ele_acc, s_sigma);
	for (unsigned int k=0; k<point_list.size(); k++)
	{
		sources_sigma[k].resize(1);
		sources_sigma[k][0] = csection[k]*(por[k]*f_rho[k]*f_cap[k]*f_sigma[k] + (1.-por[k])*s_rho[k]*s_cap[k]*s_sigma[k]);
	}
}


HeatTransferModel::~HeatTransferModel()
{}




