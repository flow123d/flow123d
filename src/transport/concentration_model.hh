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
 * @file    concentration_model.hh
 * @brief   Discontinuous Galerkin method for equation of transport with dispersion.
 * @author  Jan Stebel
 */

#ifndef CONC_TRANS_MODEL_HH_
#define CONC_TRANS_MODEL_HH_

#include "advection_diffusion_model.hh"

#include "fields/bc_field.hh"
#include "fields/field.hh"
#include "fields/multi_field.hh"




class ConcentrationTransportModel : public AdvectionDiffusionModel {
public:

	class ModelEqData : public TransportBase::TransportEqData {
	public:

		enum Concentration_bc_types {
			bc_inflow,
			bc_dirichlet,
			bc_total_flux,
			bc_diffusive_flux
		};

		/// Type of boundary condition (see also BC_Type)
        BCField<3, FieldValue<3>::EnumVector > bc_type;
		/// Prescribed concentration for Dirichlet/reference concentration for flux b.c.
		BCField<3, FieldValue<3>::Vector> bc_dirichlet_value;
		/// Flux value in total/diffusive flux b.c.
		BCField<3, FieldValue<3>::Vector > bc_flux;
		/// Transition coefficient in total/diffusive flux b.c.
		BCField<3, FieldValue<3>::Vector > bc_robin_sigma;
		/// Initial concentrations.
		Field<3, FieldValue<3>::Vector> init_conc;
		/// Longitudal dispersivity (for each substance).
		Field<3, FieldValue<3>::Vector> disp_l;
		/// Transversal dispersivity (for each substance).
		Field<3, FieldValue<3>::Vector> disp_t;
		/// Molecular diffusivity (for each substance).
		Field<3, FieldValue<3>::Vector> diff_m;


		MultiField<3, FieldValue<3>::Scalar> output_field;



		ModelEqData();

		static constexpr const char * name() { return "SoluteTransport"; }

		static string default_output_field() { return "\"conc\""; }

        static const Input::Type::Selection & get_bc_type_selection();

		static IT::Selection get_output_selection();

	};

protected:

	/// Derived class should implement getter for ModelEqData instance.
	virtual ModelEqData &data() = 0;

	/**
	 * Create input type that can be passed to the derived class.
	 * @param implementation String characterizing the numerical method, e.g. DG, FEM, FVM.
	 * @param description    Comment used to describe the record key.
	 * @return
	 */
	static IT::Record get_input_type(const string &implementation, const string &description);

	/**
	 * Formula to calculate the dispersivity tensor.
	 * @param velocity  Fluid velocity.
	 * @param Dm        Molecular diffusivity.
	 * @param alphaL    Longitudal dispersivity.
	 * @param alphaT    Transversal dispersivity.
	 * @param porosity  Porosity.
	 * @param cross_cut Cross-section.
	 * @param K         Dispersivity tensor (output).
	 */
	void calculate_dispersivity_tensor(const arma::vec3 &velocity,
			double Dm,
			double alphaL,
			double alphaT,
			double porosity,
			double cross_cut,
			arma::mat33 &K);

	/// Indicator of change in advection vector field.
	bool flux_changed;


public:

	ConcentrationTransportModel();

	static string balance_prefix() { return "mass"; }

	UnitSI balance_units();

	void set_components(SubstanceList &substances, const Input::Record &in_rec) override;

	void compute_mass_matrix_coefficient(const std::vector<arma::vec3 > &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<double> &mm_coef) override;



	void compute_advection_diffusion_coefficients(const std::vector<arma::vec3 > &point_list,
			const std::vector<arma::vec3> &velocity,
			const ElementAccessor<3> &ele_acc,
			std::vector<std::vector<arma::vec3> > &ad_coef,
			std::vector<std::vector<arma::mat33> > &dif_coef) override;

	void compute_init_cond(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector< arma::vec > &init_values) override;

	void get_bc_type(const ElementAccessor<3> &ele_acc,
				arma::uvec &bc_types) override;

	void get_flux_bc_data(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector< arma::vec > &bc_flux,
			std::vector< arma::vec > &bc_sigma,
			std::vector< arma::vec > &bc_ref_value) override;

	void get_flux_bc_sigma(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector< arma::vec > &bc_sigma) override;

	void compute_source_coefficients(const std::vector<arma::vec3> &point_list,
				const ElementAccessor<3> &ele_acc,
				std::vector<arma::vec> &sources_conc,
				std::vector<arma::vec> &sources_density,
				std::vector<arma::vec> &sources_sigma) override;

	void compute_sources_sigma(const std::vector<arma::vec3> &point_list,
				const ElementAccessor<3> &ele_acc,
				std::vector<arma::vec> &sources_sigma) override;

	~ConcentrationTransportModel() override;
};








#endif /* CONC_TRANS_MODEL_HH_ */
