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
#include "concentration_model.hh"



using namespace std;
using namespace Input::Type;







ConcentrationTransportModel::ModelEqData::ModelEqData()
: TransportBase::TransportEqData()
{
    *this+=bc_conc
            .name("bc_conc")
            .description("Dirichlet boundary condition (for each substance).")
            .input_default("0.0")
            .flags_add( in_rhs );
    *this+=init_conc
            .name("init_conc")
            .description("Initial concentrations.")
            .input_default("0.0");
    *this+=disp_l
            .name("disp_l")
            .description("Longitudal dispersivity (for each substance).")
            .input_default("0.0")
            .flags_add( in_main_matrix );
    *this+=disp_t
            .name("disp_t")
            .description("Transversal dispersivity (for each substance).")
            .input_default("0.0")
            .flags_add( in_main_matrix );
    *this+=diff_m
            .name("diff_m")
            .description("Molecular diffusivity (for each substance).")
            .input_default("0.0")
            .flags_add( in_main_matrix );

	*this+=output_field
	        .name("conc")
	        .units("M/L^3")
	        .flags( equation_result );
}





IT::Record &ConcentrationTransportModel::get_input_type(const string &implementation, const string &description)
{
	static IT::Record rec = IT::Record(ModelEqData::name() + "_" + implementation, description + " for solute transport.")
			.derive_from(AdvectionProcessBase::input_type)
			.declare_key("substances", IT::Array(IT::String()), IT::Default::obligatory(),
					"Names of transported substances.");

	return rec;
}

IT::Selection &ConcentrationTransportModel::ModelEqData::get_output_selection_input_type(const string &implementation, const string &description)
{
	static IT::Selection sel = IT::Selection(ModelEqData::name() + "_" + implementation + "_Output", "Output record for " + description + " for solute transport.");

	return sel;
}


ConcentrationTransportModel::ConcentrationTransportModel() :
		flux_changed(true)
{}


void ConcentrationTransportModel::set_component_names(std::vector<string> &names, const Input::Record &in_rec)
{
	in_rec.val<Input::Array>("substances").copy_to(names);
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

	if (fabs(vnorm) > sqrt(numeric_limits<double>::epsilon()))
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				K(i,j) = velocity[i]*velocity[j]/(vnorm*vnorm)*(alphaL-alphaT) + alphaT*(i==j?1:0);
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


void ConcentrationTransportModel::compute_dirichlet_bc(const std::vector<arma::vec3> &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector< arma::vec > &bc_values)
{
	data().bc_conc.value_list(point_list, ele_acc, bc_values);
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



