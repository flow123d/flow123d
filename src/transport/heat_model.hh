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

#ifndef HEAT_MODEL_HH_
#define HEAT_MODEL_HH_

#include "advection_diffusion_model.hh"




class HeatTransferModel : public AdvectionDiffusionModel {
public:

	class ModelEqData : public EqDataBase {
	public:

		/// Dirichlet boundary condition for temperature.
		BCField<3, FieldValue<3>::Scalar> bc_temperature;
		/// Initial temperature.
		Field<3, FieldValue<3>::Scalar> init_temperature;
		/// Porosity of solid.
		Field<3, FieldValue<3>::Scalar> porosity;
		/// Density of fluid.
		Field<3, FieldValue<3>::Scalar> fluid_density;
		/// Heat capacity of fluid.
		Field<3, FieldValue<3>::Scalar> fluid_heat_capacity;
		/// Heat conductivity of fluid.
		Field<3, FieldValue<3>::Scalar> fluid_heat_conductivity;
		/// Density of solid.
		Field<3, FieldValue<3>::Scalar> solid_density;
		/// Heat capacity of solid.
		Field<3, FieldValue<3>::Scalar> solid_heat_capacity;
		/// Heat conductivity of solid.
		Field<3, FieldValue<3>::Scalar> solid_heat_conductivity;
		/// Heat dispersivity.
		Field<3, FieldValue<3>::Scalar> heat_dispersivity;

		/// Pointer to DarcyFlow field cross_section
		Field<3, FieldValue<3>::Scalar > *cross_section;


		ModelEqData();

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
	static IT::Record &get_input_type(const string &implementation, const string &description);

	/// Indicator of change in advection vector field.
	bool flux_changed;


public:

	HeatTransferModel();

	void init_data(unsigned int n_subst_) override;

	void set_eq_data(DarcyFlowMH::EqData &water_data) override;

	void set_component_names(std::vector<string> &names, const Input::Record &in_rec) override;

	bool mass_matrix_changed() override;

	bool stiffness_matrix_changed() override;

	bool rhs_changed() override;

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

	void compute_dirichlet_bc(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector< arma::vec > &bc_values) override;

	void compute_source_coefficients(const std::vector<arma::vec3> &point_list,
				const ElementAccessor<3> &ele_acc,
				std::vector<arma::vec> &sources_conc,
				std::vector<arma::vec> &sources_density,
				std::vector<arma::vec> &sources_sigma) override;

	void compute_sources_sigma(const std::vector<arma::vec3> &point_list,
				const ElementAccessor<3> &ele_acc,
				std::vector<arma::vec> &sources_sigma) override;

	~HeatTransferModel() override;

};



#endif /* HEAT_MODEL_HH_ */
