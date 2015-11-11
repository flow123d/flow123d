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
//#include "flow/darcy_flow_mh.hh"

#include "transport/transport_operator_splitting.hh"
#include "concentration_model.hh"
#include "fields/unit_si.hh"
#include "coupling/balance.hh"



using namespace std;
using namespace Input::Type;



const Selection & ConcentrationTransportModel::ModelEqData::get_bc_type_selection() {
	return Selection("SoluteTransport_BC_Type", "Types of boundary conditions for solute transport model.")
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
            		  "The prescribed total flux can have the general form (($\\delta(f_N+\\sigma_R(c-c_R) )+q_wc_A$)), "
            		  "where the absolute flux (($f_N$)) is specified by the field 'bc_flux', "
            		  "the advected concentration (($c_A$)) by 'bc_ad_conc', "
            		  "the transition parameter (($\\sigma_R$)) by 'bc_robin_sigma', "
            		  "and the reference concentration (($c_R$)) by 'bc_conc'.")
              .add_value(bc_diffusive_flux, "diffusive_flux",
            		  "Diffusive flux boundary condition.\n"
            		  "The prescribed mass flux due to diffusion can have the general form (($\\delta(f_N+\\sigma_R(c-c_R) )$)), "
            		  "where the absolute flux (($f_N$)) is specified by the field 'bc_flux', "
            		  "the transition parameter (($\\sigma_R$)) by 'bc_robin_sigma', "
            		  "and the reference concentration (($c_R$)) by 'bc_conc'.")
			  .close();
}




ConcentrationTransportModel::ModelEqData::ModelEqData()
: TransportBase::TransportEqData()
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
            .name("bc_conc")
            .units( UnitSI().kg().m(-3) )
            .description("Dirichlet boundary condition (for each substance).")
            .input_default("0.0")
            .flags_add( in_rhs );
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
    *this+=init_conc
            .name("init_conc")
            .units( UnitSI().kg().m(-3) )
            .description("Initial concentrations.")
            .input_default("0.0");
    *this+=disp_l
            .name("disp_l")
            .description("Longitudal dispersivity (for each substance).")
            .units( UnitSI().m() )
            .input_default("0.0")
            .flags_add( in_main_matrix & in_rhs );
    *this+=disp_t
            .name("disp_t")
            .description("Transversal dispersivity (for each substance).")
            .units( UnitSI().m() )
            .input_default("0.0")
            .flags_add( in_main_matrix & in_rhs );
    *this+=diff_m
            .name("diff_m")
            .description("Molecular diffusivity (for each substance).")
            .units( UnitSI().m(2).s(-1) )
            .input_default("0.0")
            .flags_add( in_main_matrix & in_rhs );

	*this+=output_field
	        .name("conc")
	        .units( UnitSI().kg().m(-3) )
	        .flags( equation_result );
}



UnitSI ConcentrationTransportModel::balance_units()
{
	return UnitSI().kg();
}


IT::Record ConcentrationTransportModel::get_input_type(const string &implementation, const string &description)
{
	return IT::Record(
				std::string(ModelEqData::name()) + "_" + implementation,
				description + " for solute transport.")
			.derive_from(AdvectionProcessBase::get_input_type())
			.declare_key("time", TimeGovernor::get_input_type(), Default::obligatory(),
					"Time governor setting for the secondary equation.")
			.declare_key("balance", Balance::get_input_type(), Default::obligatory(),
					"Settings for computing balance.")
			.declare_key("output_stream", OutputTime::get_input_type(), Default::obligatory(),
					"Parameters of output stream.")
			.declare_key("substances", IT::Array( Substance::get_input_type() ), IT::Default::obligatory(),
					"Names of transported substances.");
}

IT::Selection ConcentrationTransportModel::ModelEqData::get_output_selection_input_type(const string &implementation, const string &description)
{
	return IT::Selection(
				std::string(ModelEqData::name()) + "_" + implementation + "_Output",
				"Output record for " + description + " for solute transport.");
}


ConcentrationTransportModel::ConcentrationTransportModel() :
		flux_changed(true)
{}


void ConcentrationTransportModel::set_components(SubstanceList &substances, const Input::Record &in_rec)
{
	substances.initialize(in_rec.val<Input::Array>("substances"));
}


void ConcentrationTransportModel::compute_mass_matrix_coefficient(const std::vector<arma::vec3 > &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector<double> &mm_coef)
{
	vector<double> elem_csec(point_list.size()), por_m(point_list.size());

	data().cross_section.value_list(point_list, ele_acc, elem_csec);
	data().porosity.value_list(point_list, ele_acc, por_m);

	for (unsigned int i=0; i<point_list.size(); i++)
		mm_coef[i] = elem_csec[i]*por_m[i];
}


void ConcentrationTransportModel::calculate_dispersivity_tensor(const arma::vec3 &velocity,
		double Dm, double alphaL, double alphaT, double porosity, double cross_cut, arma::mat33 &K)
{
    double vnorm = arma::norm(velocity, 2);

	if (fabs(vnorm) > 0)
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				K(i,j) = (velocity[i]/vnorm)*(velocity[j]/vnorm)*(alphaL-alphaT) + alphaT*(i==j?1:0);
	else
		K.zeros();

	// Note that the velocity vector is in fact the Darcian flux,
	// so to obtain |v| we have to divide vnorm by porosity and cross_section.
	K = (vnorm*K + (Dm*pow(porosity, 1./3)*porosity*cross_cut)*arma::eye(3,3));
}




void ConcentrationTransportModel::compute_advection_diffusion_coefficients(const std::vector<arma::vec3 > &point_list,
		const std::vector<arma::vec3> &velocity,
		const ElementAccessor<3> &ele_acc,
		std::vector<std::vector<arma::vec3> > &ad_coef,
		std::vector<std::vector<arma::mat33> > &dif_coef)
{
	const unsigned int qsize = point_list.size();
	const unsigned int n_subst = dif_coef.size();
	std::vector<arma::vec> Dm(qsize, arma::vec(n_subst) ), alphaL(qsize, arma::vec(n_subst) ), alphaT(qsize, arma::vec(n_subst) );
	std::vector<double> por_m(qsize), csection(qsize);

	data().diff_m.value_list(point_list, ele_acc, Dm);
	data().disp_l.value_list(point_list, ele_acc, alphaL);
	data().disp_t.value_list(point_list, ele_acc, alphaT);
	data().porosity.value_list(point_list, ele_acc, por_m);
	data().cross_section.value_list(point_list, ele_acc, csection);

	for (unsigned int i=0; i<qsize; i++) {
		for (unsigned int sbi=0; sbi<n_subst; sbi++) {
			ad_coef[sbi][i] = velocity[i];
			calculate_dispersivity_tensor(velocity[i],
					Dm[i][sbi], alphaL[i][sbi], alphaT[i][sbi], por_m[i], csection[i],
					dif_coef[sbi][i]);
		}
	}
}


void ConcentrationTransportModel::compute_init_cond(const std::vector<arma::vec3> &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector< arma::vec > &init_values)
{
	data().init_conc.value_list(point_list, ele_acc, init_values);
}

void ConcentrationTransportModel::get_bc_type(const ElementAccessor<3> &ele_acc,
			arma::uvec &bc_types)
{
	// Currently the bc types for ConcentrationTransport are numbered in the same way as in TransportDG.
	// In general we should use some map here.
	bc_types = data().bc_type.value(ele_acc.centre(), ele_acc);
}


void ConcentrationTransportModel::get_flux_bc_data(const std::vector<arma::vec3> &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector< arma::vec > &bc_flux,
		std::vector< arma::vec > &bc_sigma,
		std::vector< arma::vec > &bc_ref_value)
{
	data().bc_flux.value_list(point_list, ele_acc, bc_flux);
	data().bc_robin_sigma.value_list(point_list, ele_acc, bc_sigma);
	data().bc_dirichlet_value.value_list(point_list, ele_acc, bc_ref_value);
}

void ConcentrationTransportModel::get_flux_bc_sigma(const std::vector<arma::vec3> &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector< arma::vec > &bc_sigma)
{
	data().bc_robin_sigma.value_list(point_list, ele_acc, bc_sigma);
}


void ConcentrationTransportModel::compute_source_coefficients(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<arma::vec> &sources_value,
			std::vector<arma::vec> &sources_density,
			std::vector<arma::vec> &sources_sigma)
{
	const unsigned int qsize = point_list.size();
	vector<double> csection(qsize);
	data().cross_section.value_list(point_list, ele_acc, csection);
	data().sources_conc.value_list(point_list, ele_acc, sources_value);
	data().sources_density.value_list(point_list, ele_acc, sources_density);
	data().sources_sigma.value_list(point_list, ele_acc, sources_sigma);

	for (unsigned int k=0; k<qsize; k++)
	{
		sources_density[k] *= csection[k];
		sources_sigma[k] *= csection[k];
	}
}


void ConcentrationTransportModel::compute_sources_sigma(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<arma::vec> &sources_sigma)
{
	const unsigned int qsize = point_list.size();
	vector<double> csection(qsize);
	data().cross_section.value_list(point_list, ele_acc, csection);
	data().sources_sigma.value_list(point_list, ele_acc, sources_sigma);

	for (unsigned int k=0; k<qsize; k++)
		sources_sigma[k] *= csection[k];
}


ConcentrationTransportModel::~ConcentrationTransportModel()
{}



