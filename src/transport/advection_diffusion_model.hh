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

/**
 * AdvectionDiffusionModel is a base class for description of a physical process described
 * by the advection-diffusion partial differential equation (PDE). The derived classes define input parameters
 * and implement methods that calculate coefficients of the PDE. These methods are then used by a template
 * class for numerical solution, whose specialization derives from the model class.
 */
class AdvectionDiffusionModel {
public:

	/// Initialize model data. E.g. set vector field dimensions.
	virtual void init_data(unsigned int n_subst_) = 0;

<<<<<<< HEAD
	virtual void set_cross_section_field(Field< 3, FieldValue<3>::Scalar >* cross_section) = 0;
=======
	/// Temporary solution for sharing data from other equations.
	virtual void set_eq_data(DarcyFlowMH::EqData &water_data) = 0;
>>>>>>> 15478b6... Documentation of class AdvectionDiffusionModel (ADM).

	/// Read or set names of solution components.
	virtual void set_component_names(std::vector<string> &names, const Input::Record &in_rec) = 0;

	/// Check if mass matrix coefficients have changed.
	virtual bool mass_matrix_changed() = 0;

	/// Check if stiffness matrix coefficients have changed.
	virtual bool stiffness_matrix_changed() = 0;

	/// Check if right hand side coefficients have changed.
	virtual bool rhs_changed() = 0;

	/**
	 * Compute coefficients of mass matrix.
	 * @param point_list   Points at which to evaluate.
	 * @param ele_acc      Element accessor.
	 * @param mm_coef      Coefficient vector (output).
	 */
	virtual void compute_mass_matrix_coefficient(const std::vector<arma::vec3 > &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<double> &mm_coef) = 0;

	/**
	 * Compute coefficients of stiffness matrix.
	 * @param point_list  Points at which to evaluate.
	 * @param velocity    Velocity field (input). Temporary solution before we can pass data from other
	 *                    equations.
	 * @param ele_acc     Element accessor.
	 * @param ad_coef     Coefficients of advection (output).
	 * @param dif_coef    Coefficients of diffusion (output).
	 */
	virtual void compute_advection_diffusion_coefficients(const std::vector<arma::vec3> &point_list,
			const std::vector<arma::vec3> &velocity,
			const ElementAccessor<3> &ele_acc,
			std::vector<std::vector<arma::vec3> > &ad_coef,
			std::vector<std::vector<arma::mat33> > &dif_coef) = 0;

	/**
	 * Compute initial conditions.
	 * @param point_list   Points at which to evaluate.
	 * @param ele_acc      Element accessor.
	 * @param init_values  Vector of intial values (output).
	 */
	virtual void compute_init_cond(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector< arma::vec > &init_values) = 0;

	/**
	 * Computes the Dirichlet boundary condition values.
	 * @param point_list   Points at which to evaluate.
	 * @param ele_acc      Element accessor.
	 * @param bc_values    Vector of b.c. values (output).
	 */
	virtual void compute_dirichlet_bc(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector< arma::vec > &bc_values) = 0;

	/**
	 * Compute coefficients of volume sources.
	 * @param point_list      Points at which to evaluate.
	 * @param ele_acc         Element accessor.
	 * @param sources_conc    Source concentrations (output).
	 * @param sources_density Source densities (output).
	 * @param sources_sigma   Source sigmas (output).
	 */
	virtual void compute_source_coefficients(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<arma::vec> &sources_conc,
			std::vector<arma::vec> &sources_density,
			std::vector<arma::vec> &sources_sigma) = 0;

	/**
	 * Compute coefficients of volume sources.
	 * @param point_list      Points at which to evaluate.
	 * @param ele_acc         Element accessor.
	 * @param sources_sigma   Source sigmas (output).
	 */
	virtual void compute_sources_sigma(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<arma::vec> &sources_sigma) = 0;

	/// Destructor.
	virtual ~AdvectionDiffusionModel() {};



};


class ConcentrationTransportModel : public AdvectionDiffusionModel {
public:

	class ModelEqData : public TransportBase::TransportEqData {
	public:

		/// Boundary conditions (Dirichlet) for concentrations.
		BCField<3, FieldValue<3>::Vector> bc_conc;
		/// Initial concentrations.
		Field<3, FieldValue<3>::Vector> init_conc;
		/// Longitudal dispersivity (for each substance).
		Field<3, FieldValue<3>::Vector> disp_l;
		/// Transversal dispersivity (for each substance).
		Field<3, FieldValue<3>::Vector> disp_t;
		/// Molecular diffusivity (for each substance).
		Field<3, FieldValue<3>::Vector> diff_m;


		ModelEqData();

		/// Reads boundary conditions in old format.
        RegionSet read_boundary_list_item(Input::Record rec);

	};

protected:

	/// Derived class should implement getter for ModelEqData instance.
	virtual ModelEqData &data() = 0;

	/// Create input type that can be passed to the derived class.
	static IT::Record &get_input_type(const string &implementation, const string &description);

	/// Indicator of change in advection vector field.
	bool flux_changed;

	/// Name of the physical process (used in definition of input record).
	static const char *input_key_name_;


public:

	ConcentrationTransportModel();

	void init_data(unsigned int n_subst_);

	void set_cross_section_field(Field< 3, FieldValue<3>::Scalar >* cross_section);

	void set_component_names(std::vector<string> &names, const Input::Record &in_rec);

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

	/// Create input type that can be passed to the derived class.
	static IT::Record &get_input_type(const string &implementation, const string &description);

	/// Indicator of change in advection vector field.
	bool flux_changed;

	/// Name of the physical process (used in definition of input record).
	static const char *input_key_name_;


public:

	HeatTransferModel();

	void init_data(unsigned int n_subst_);

	void set_cross_section_field(Field< 3, FieldValue<3>::Scalar >* cross_section);

	void set_component_names(std::vector<string> &names, const Input::Record &in_rec);

	static string input_key_name() { return input_key_name_; }

	bool mass_matrix_changed();

	bool stiffness_matrix_changed();

	bool rhs_changed();

	void compute_mass_matrix_coefficient(const std::vector<arma::vec3 > &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<double> &mm_coef);

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

	~HeatTransferModel();

};



#endif /* AD_MODEL_HH_ */
