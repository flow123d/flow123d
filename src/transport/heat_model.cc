/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief Discontinuous Galerkin method for equation of transport with dispersion.
 *  @author Jan Stebel
 */

#include "input/input_type.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "flow/darcy_flow_mh.hh"
#include "transport/transport_operator_splitting.hh"
#include "heat_model.hh"



using namespace std;
using namespace Input::Type;










HeatTransferModel::ModelEqData::ModelEqData() : EqDataBase("HeatTransfer")
{
	ADD_FIELD(bc_temperature, "Boundary value of temperature.", "0.0");

	ADD_FIELD(init_temperature, "Initial temperature.", "0.0");
	ADD_FIELD(porosity, "Porosity.", "1.0" );
	ADD_FIELD(fluid_density, "Density of fluid.");
	ADD_FIELD(fluid_heat_capacity, "Heat capacity of fluid.");
	ADD_FIELD(fluid_heat_conductivity, "Heat conductivity of fluid.");
	ADD_FIELD(solid_density, "Density of solid (rock).");
	ADD_FIELD(solid_heat_capacity, "Heat capacity of solid (rock).");
	ADD_FIELD(solid_heat_conductivity, "Heat conductivity of solid (rock).");
	ADD_FIELD(heat_dispersivity, "Heat dispersivity", "0.0" );
}


IT::Record &HeatTransferModel::get_input_type(const string &implementation, const string &description)
{
	static IT::Record input_type = IT::Record("HeatTransfer_" + implementation, description + " for heat transfer.")
			.derive_from(AdvectionProcessBase::input_type);

	return input_type;

}


HeatTransferModel::HeatTransferModel() :
		flux_changed(true)
{}


void HeatTransferModel::init_data(unsigned int n_subst_)
{}


void HeatTransferModel::set_cross_section_field(Field< 3, FieldValue<3>::Scalar >* cross_section)
{
	data().cross_section = cross_section;
}


void HeatTransferModel::set_component_names(std::vector<string> &names, const Input::Record &in_rec)
{
	names.clear();
	names.push_back("temperature");
}


bool HeatTransferModel::mass_matrix_changed()
{
	return (data().porosity.changed() ||
			data().fluid_density.changed() ||
			data().fluid_heat_capacity.changed() ||
			data().solid_density.changed() ||
			data().solid_heat_capacity.changed() ||
			data().cross_section->changed());
}


bool HeatTransferModel::stiffness_matrix_changed()
{
	return (flux_changed ||
			data().fluid_density.changed() ||
			data().fluid_heat_capacity.changed() ||
			data().porosity.changed() ||
			data().fluid_heat_conductivity.changed() ||
			data().solid_heat_conductivity.changed() ||
			data().heat_dispersivity.changed() ||
			data().cross_section->changed());
}


bool HeatTransferModel::rhs_changed()
{
	return (flux_changed ||
			data().bc_temperature.changed());
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

	data().cross_section->value_list(point_list, ele_acc, elem_csec);
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
			s_cond(qsize), por(qsize), csection(qsize), disp(qsize);

	data().fluid_density.value_list(point_list, ele_acc, f_rho);
	data().fluid_heat_capacity.value_list(point_list, ele_acc, f_cap);
	data().fluid_heat_conductivity.value_list(point_list, ele_acc, f_cond);
	data().solid_heat_conductivity.value_list(point_list, ele_acc, s_cond);
	data().heat_dispersivity.value_list(point_list, ele_acc, disp);
	data().porosity.value_list(point_list, ele_acc, por);
	data().cross_section->value_list(point_list, ele_acc, csection);

	for (unsigned int i=0; i<qsize; i++) {
		ad_coef[0][i] = velocity[i]*f_rho[i]*f_cap[i];
		dif_coef[0][i] = csection[i]*(disp[i] + por[i]*f_cond[i] + (1.-por[i])*s_cond[i])*arma::eye(3,3);
	}
}


void HeatTransferModel::compute_init_cond(const std::vector<arma::vec3> &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector< arma::vec > &init_values)
{
	vector<double> init_value(point_list.size());
	data().init_temperature.value_list(point_list, ele_acc, init_value);
	for (int i=0; i<point_list.size(); i++)
		init_values[i] = init_value[i];
}


void HeatTransferModel::compute_dirichlet_bc(const std::vector<arma::vec3> &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector< arma::vec > &bc_values)
{
	vector<double> bc_value(point_list.size());
	data().bc_temperature.value_list(point_list, ele_acc, bc_value);
	for (int i=0; i<point_list.size(); i++)
		bc_values[i] = bc_value[i];
}


void HeatTransferModel::compute_source_coefficients(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<arma::vec> &sources_conc,
			std::vector<arma::vec> &sources_density,
			std::vector<arma::vec> &sources_sigma)
{
//	data().sources_conc.value_list(point_list, ele_acc, sources_conc);
//	data().sources_density.value_list(point_list, ele_acc, sources_density);
//	data().sources_sigma.value_list(point_list, ele_acc, sources_sigma);
	for (int k=0; k<point_list.size(); k++)
	{
		sources_conc[k].resize(1);
		sources_density[k].resize(1);
		sources_sigma[k].resize(1);
	}
}


void HeatTransferModel::compute_sources_sigma(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<arma::vec> &sources_sigma)
{
//	data().sources_sigma.value_list(point_list, ele_acc, sources_sigma);
	for (int k=0; k<point_list.size(); k++)
		sources_sigma[k].resize(1);
}


HeatTransferModel::~HeatTransferModel()
{}

