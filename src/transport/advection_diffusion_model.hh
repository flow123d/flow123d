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

#ifndef AD_MODEL_HH_
#define AD_MODEL_HH_

#include <armadillo>
#include <vector>

namespace IT = Input::Type;

class AdvectionDiffusionModel {
public:

	virtual void init_data(unsigned int n_subst_) = 0;

	virtual void set_cross_section_field(Field< 3, FieldValue<3>::Scalar >* cross_section) = 0;

	virtual bool mass_matrix_changed() = 0;

	virtual bool stiffness_matrix_changed() = 0;

	virtual bool rhs_changed() = 0;

	virtual void compute_mass_matrix_coefficient(const std::vector<arma::vec3 > &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<double> &mm_coef) = 0;

	virtual void compute_advection_diffusion_coefficients(const std::vector<arma::vec3> &point_list,
			const std::vector<arma::vec3> &velocity,
			const ElementAccessor<3> &ele_acc,
			std::vector<std::vector<arma::vec3> > &ad_coef,
			std::vector<std::vector<arma::mat33> > &dif_coef) = 0;

	virtual void compute_init_cond(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector< arma::vec > &init_values) = 0;

	virtual void compute_dirichlet_bc(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector< arma::vec > &bc_values) = 0;

	virtual void compute_source_coefficients(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<arma::vec> &sources_conc,
			std::vector<arma::vec> &sources_density,
			std::vector<arma::vec> &sources_sigma) = 0;

	virtual void compute_sources_sigma(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<arma::vec> &sources_sigma) = 0;

	virtual ~AdvectionDiffusionModel() {};



};


class ConcentrationTransportModel : public AdvectionDiffusionModel {
public:

	class ModelEqData : public TransportBase::TransportEqData {
	public:

		/**
		 * Boundary conditions (Dirichlet) for concentrations.
		 * They are applied only on water inflow part of the boundary.
		 */
		BCField<3, FieldValue<3>::Vector> bc_conc;

		/// Initial concentrations.
		Field<3, FieldValue<3>::Vector> init_conc;
		Field<3, FieldValue<3>::Vector> disp_l;     ///< Longitudal dispersivity (for each substance).
		Field<3, FieldValue<3>::Vector> disp_t;     ///< Transversal dispersivity (for each substance).
		Field<3, FieldValue<3>::Vector> diff_m;     ///< Molecular diffusivity (for each substance).


		ModelEqData();

        RegionSet read_boundary_list_item(Input::Record rec);

	};

protected:

	virtual ModelEqData &data() = 0;

	static IT::Record &get_input_type(const string &implementation, const string &description) {
		return IT::Record("ConcentrationTransport_" + implementation, description + " for solute transport.")
				.derive_from(AdvectionProcessBase::input_type)
				.declare_key("substances", IT::Array(IT::String()), IT::Default::obligatory(),
						"Names of transported substances.")
						// input data
				.declare_key("sorption_enable", IT::Bool(), IT::Default("false"),
						"Model of sorption.")
				.declare_key("dual_porosity", IT::Bool(), IT::Default("false"),
						"Dual porosity model.")
				.declare_key("output", TransportBase::input_type_output_record, IT::Default::obligatory(),
						"Parameters of output stream.");
	}

	bool flux_changed;

	static const char *input_key_name_;


public:

	ConcentrationTransportModel();

	void init_data(unsigned int n_subst_);

	void set_cross_section_field(Field< 3, FieldValue<3>::Scalar >* cross_section);

	static string input_key_name() { return input_key_name_; }

	bool mass_matrix_changed();

	bool stiffness_matrix_changed();

	bool rhs_changed();

	void compute_mass_matrix_coefficient(const std::vector<arma::vec3 > &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<double> &mm_coef);

	void calculate_dispersivity_tensor(arma::mat33 &K, const arma::vec3 &velocity,
			double Dm, double alphaL, double alphaT, double porosity, double cross_cut);

	void compute_advection_diffusion_coefficients(const std::vector<arma::vec3 > &point_list,
			const std::vector<arma::vec3> &velocity,
			const ElementAccessor<3> &ele_acc,
			std::vector<std::vector<arma::vec3> > &ad_coef,
			std::vector<std::vector<arma::mat33> > &dif_coef);

	void compute_init_cond(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector< arma::vec > &init_values);

	void compute_dirichlet_bc(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector< arma::vec > &bc_values);

	void compute_source_coefficients(const std::vector<arma::vec3> &point_list,
				const ElementAccessor<3> &ele_acc,
				std::vector<arma::vec> &sources_conc,
				std::vector<arma::vec> &sources_density,
				std::vector<arma::vec> &sources_sigma);

	void compute_sources_sigma(const std::vector<arma::vec3> &point_list,
				const ElementAccessor<3> &ele_acc,
				std::vector<arma::vec> &sources_sigma);

	~ConcentrationTransportModel();
};





class HeatTransferModel : public AdvectionDiffusionModel {
public:

	class ModelEqData : public EqDataBase {
	public:

		/// Dirichlet boundary condition for temperature.
		BCField<3, FieldValue<3>::Scalar> bc_conc;

		/// Initial temperature.
		Field<3, FieldValue<3>::Vector> init_temperature;
		Field<3, FieldValue<3>::Scalar> porosity;
		Field<3, FieldValue<3>::Scalar> fluid_density;
		Field<3, FieldValue<3>::Scalar> fluid_heat_capacity;
		Field<3, FieldValue<3>::Scalar> fluid_heat_conductivity;
		Field<3, FieldValue<3>::Scalar> solid_density;
		Field<3, FieldValue<3>::Scalar> solid_heat_capacity;
		Field<3, FieldValue<3>::Scalar> solid_heat_conductivity;
		Field<3, FieldValue<3>::Scalar> heat_dispersivity;


		ModelEqData() : EqDataBase("HeatTransferDG") {};

	};

protected:

	virtual ModelEqData &data() = 0;

	static IT::Record &get_input_type(const string &implementation, const string &description) {
		static IT::Record input_type = IT::Record("HeatTransfer_" + implementation, description + " for heat transfer.")
				.derive_from(AdvectionProcessBase::input_type);

		return input_type;

	}

	bool flux_changed;

	static const char *input_key_name_;


public:

	HeatTransferModel() {};

	void init_data(unsigned int n_subst_) {};

	void set_cross_section_field(Field< 3, FieldValue<3>::Scalar >* cross_section) {};

	static string input_key_name() { return input_key_name_; }

	bool mass_matrix_changed() {};

	bool stiffness_matrix_changed() {};

	bool rhs_changed() {};

	void compute_mass_matrix_coefficient(const std::vector<arma::vec3 > &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<double> &mm_coef) {};

	void calculate_dispersivity_tensor(arma::mat33 &K, const arma::vec3 &velocity,
			double Dm, double alphaL, double alphaT, double porosity, double cross_cut) {};

	void compute_advection_diffusion_coefficients(const std::vector<arma::vec3 > &point_list,
			const std::vector<arma::vec3> &velocity,
			const ElementAccessor<3> &ele_acc,
			std::vector<std::vector<arma::vec3> > &ad_coef,
			std::vector<std::vector<arma::mat33> > &dif_coef) {};

	void compute_init_cond(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector< arma::vec > &init_values) {};

	void compute_dirichlet_bc(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector< arma::vec > &bc_values) {};

	void compute_source_coefficients(const std::vector<arma::vec3> &point_list,
				const ElementAccessor<3> &ele_acc,
				std::vector<arma::vec> &sources_conc,
				std::vector<arma::vec> &sources_density,
				std::vector<arma::vec> &sources_sigma) {};

	void compute_sources_sigma(const std::vector<arma::vec3> &point_list,
				const ElementAccessor<3> &ele_acc,
				std::vector<arma::vec> &sources_sigma) {};

	~HeatTransferModel() {};

	static Input::Type::AbstractRecord input_type;
};



#endif /* AD_MODEL_HH_ */
