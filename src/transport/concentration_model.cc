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
#include "tools/unit_si.hh"
#include "coupling/balance.hh"



using namespace std;
using namespace Input::Type;



const Selection & ConcentrationTransportModel::ModelEqData::get_bc_type_selection() {
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






ConcentrationTransportModel::ModelEqData::ModelEqData()
: TransportEqData()
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
    *this+=init_conc
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

IT::Selection ConcentrationTransportModel::ModelEqData::get_output_selection()
{
    // Return empty selection just to provide model specific selection name and description.
    // The fields are added by TransportDG using an auxiliary selection.
	return IT::Selection(
				std::string(ModelEqData::name()) + "_DG_output_fields",
				"Selection of output fields for Diffusive Solute Transport DG model.");
}


ConcentrationTransportModel::ConcentrationTransportModel(Mesh &mesh, const Input::Record &in_rec) :
		ConcentrationTransportBase(mesh, in_rec),
		flux_changed(true)
{}


void ConcentrationTransportModel::init_from_input(const Input::Record &in_rec)
{
	solvent_density_ = in_rec.val<double>("solvent_density");
}



void ConcentrationTransportModel::compute_mass_matrix_coefficient(const Armor::array &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector<double> &mm_coef)
{
	vector<double> elem_csec(point_list.size()), wc(point_list.size()); //por_m(point_list.size()),

	data().cross_section.value_list(point_list, ele_acc, elem_csec);
	//data().porosity.value_list(point_list, ele_acc, por_m);
	data().water_content.value_list(point_list, ele_acc, wc);

	for (unsigned int i=0; i<point_list.size(); i++)
		mm_coef[i] = elem_csec[i]*wc[i];
}


void ConcentrationTransportModel::compute_retardation_coefficient(const Armor::array &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector<std::vector<double> > &ret_coef)
{
	vector<double> elem_csec(point_list.size()),
			por_m(point_list.size()),
			rho_l(point_list.size()),
			rho_s(point_list.size()),
			sorp_mult(point_list.size());

	data().cross_section.value_list(point_list, ele_acc, elem_csec);
	data().porosity.value_list(point_list, ele_acc, por_m);
	data().rock_density.value_list(point_list, ele_acc, rho_s);

	// Note: Noe effective water content here, since sorption happen only in the rock (1-porosity).
	for (unsigned int sbi=0; sbi<substances_.size(); sbi++)
	{
        data().sorption_coefficient[sbi].value_list(point_list, ele_acc, sorp_mult);
		for (unsigned int i=0; i<point_list.size(); i++)
		{
			ret_coef[sbi][i] = (1.-por_m[i])*rho_s[i]*sorp_mult[i]*elem_csec[i];
		}
	}
}


void ConcentrationTransportModel::calculate_dispersivity_tensor(const arma::vec3 &velocity,
		const arma::mat33 &Dm, double alphaL, double alphaT, double water_content, double porosity, double cross_cut, arma::mat33 &K)
{
    double vnorm = arma::norm(velocity, 2);


	// used tortuosity model dues to Millington and Quirk(1961) (should it be with power 10/3 ?)
	// for an overview of other models see: Chou, Wu, Zeng, Chang (2011)
	double tortuosity = pow(water_content, 7.0 / 3.0)/ (porosity * porosity);

    // Note that the velocity vector is in fact the Darcian flux,
    // so we need not to multiply vnorm by water_content and cross_section.
	//K = ((alphaL-alphaT) / vnorm) * K + (alphaT*vnorm + Dm*tortuosity*cross_cut*water_content) * arma::eye(3,3);

    if (fabs(vnorm) > 0) {
        /*
        for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
               K(i,j) = (velocity[i]/vnorm)*(velocity[j]);
        */
        K = ((alphaL - alphaT) / vnorm) * arma::kron(velocity.t(), velocity);

        //arma::mat33 abs_diff_mat = arma::abs(K -  kk);
        //double diff = arma::min( arma::min(abs_diff_mat) );
        //ASSERT(  diff < 1e-12 )(diff)(K)(kk);
    } else
        K.zeros();

   // Note that the velocity vector is in fact the Darcian flux,
   // so to obtain |v| we have to divide vnorm by porosity and cross_section.
   K += alphaT*vnorm*arma::eye(3,3) + Dm*(tortuosity*cross_cut*water_content);

}




void ConcentrationTransportModel::compute_advection_diffusion_coefficients(const Armor::array &point_list,
		const std::vector<arma::vec3> &velocity,
		const ElementAccessor<3> &ele_acc,
		std::vector<std::vector<arma::vec3> > &ad_coef,
		std::vector<std::vector<arma::mat33> > &dif_coef)
{
	const unsigned int qsize = point_list.size();
	const unsigned int n_subst = dif_coef.size();
    std::vector<arma::mat33> Dm(qsize);
	std::vector<double> alphaL(qsize), alphaT(qsize), por_m(qsize), csection(qsize), wc(qsize);

	data().porosity.value_list(point_list, ele_acc, por_m);
	data().water_content.value_list(point_list, ele_acc, wc);
	data().cross_section.value_list(point_list, ele_acc, csection);

    for (unsigned int sbi=0; sbi<n_subst; sbi++)
    {
        data().diff_m[sbi].value_list(point_list, ele_acc, Dm);
        data().disp_l[sbi].value_list(point_list, ele_acc, alphaL);
        data().disp_t[sbi].value_list(point_list, ele_acc, alphaT);
		for (unsigned int i=0; i<qsize; i++)
        {
			ad_coef[sbi][i] = velocity[i];
			calculate_dispersivity_tensor(velocity[i],
					Dm[i], alphaL[i], alphaT[i], wc[i], por_m[i], csection[i],
					dif_coef[sbi][i]);
		}
	}
}


void ConcentrationTransportModel::compute_init_cond(const Armor::array &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector<std::vector<double> > &init_values)
{
    for (unsigned int sbi=0; sbi<n_substances(); sbi++)
        data().init_conc[sbi].value_list(point_list, ele_acc, init_values[sbi]);
}

void ConcentrationTransportModel::get_bc_type(const ElementAccessor<3> &ele_acc,
			arma::uvec &bc_types)
{
	// Currently the bc types for ConcentrationTransport are numbered in the same way as in TransportDG.
	// In general we should use some map here.
    bc_types.resize(n_substances());
    for (unsigned int sbi=0; sbi<n_substances(); sbi++)
        bc_types[sbi] = data().bc_type[sbi].value(ele_acc.centre(), ele_acc);
}


void ConcentrationTransportModel::get_flux_bc_data(unsigned int index,
        const Armor::array &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector< double > &bc_flux,
		std::vector< double > &bc_sigma,
		std::vector< double > &bc_ref_value)
{
	data().bc_flux[index].value_list(point_list, ele_acc, bc_flux);
	data().bc_robin_sigma[index].value_list(point_list, ele_acc, bc_sigma);
	data().bc_dirichlet_value[index].value_list(point_list, ele_acc, bc_ref_value);
	
	// Change sign in bc_flux since internally we work with outgoing fluxes.
	for (auto f : bc_flux) f = -f;
}

void ConcentrationTransportModel::get_flux_bc_sigma(unsigned int index,
        const Armor::array &point_list,
		const ElementAccessor<3> &ele_acc,
		std::vector< double > &bc_sigma)
{
	data().bc_robin_sigma[index].value_list(point_list, ele_acc, bc_sigma);
}


void ConcentrationTransportModel::compute_source_coefficients(const Armor::array &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<std::vector<double> > &sources_value,
			std::vector<std::vector<double> > &sources_density,
			std::vector<std::vector<double> > &sources_sigma)
{
	const unsigned int qsize = point_list.size();
	vector<double> csection(qsize);
	data().cross_section.value_list(point_list, ele_acc, csection);
    for (unsigned int sbi=0; sbi<n_substances(); sbi++)
    {
      data().sources_conc[sbi].value_list(point_list, ele_acc, sources_value[sbi]);
      data().sources_density[sbi].value_list(point_list, ele_acc, sources_density[sbi]);
      data().sources_sigma[sbi].value_list(point_list, ele_acc, sources_sigma[sbi]);
      
      for (unsigned int k=0; k<qsize; k++)
      {
          sources_density[sbi][k] *= csection[k];
          sources_sigma[sbi][k] *= csection[k];
      }
    }
}


void ConcentrationTransportModel::compute_sources_sigma(const Armor::array &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<std::vector<double> > &sources_sigma)
{
	const unsigned int qsize = point_list.size();
	vector<double> csection(qsize);
	data().cross_section.value_list(point_list, ele_acc, csection);
    for (unsigned int sbi=0; sbi<n_substances(); sbi++)
    {
        data().sources_sigma[sbi].value_list(point_list, ele_acc, sources_sigma[sbi]);

        for (unsigned int k=0; k<qsize; k++)
            sources_sigma[sbi][k] *= csection[k];
    }
}


ConcentrationTransportModel::~ConcentrationTransportModel()
{}


void ConcentrationTransportModel::set_balance_object(std::shared_ptr<Balance> balance)
{
	balance_ = balance;
	subst_idx = balance_->add_quantities(substances_.names());
}



