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
 * @file    advection_diffusion_model.hh
 * @brief   Discontinuous Galerkin method for equation of transport with dispersion.
 * @author  Jan Stebel
 */

#ifndef AD_MODEL_HH_
#define AD_MODEL_HH_

#include <armadillo>
#include <vector>

class SubstanceList;

namespace IT = Input::Type;

/**
 * AdvectionDiffusionModel is a base class for description of a physical process described
 * by the advection-diffusion partial differential equation (PDE). The derived classes define input parameters
 * and implement methods that calculate coefficients of the PDE. These methods are then used by a template
 * class for numerical solution, whose specialization derives from the model class.
 */
class AdvectionDiffusionModel {
public:

    enum Abstract_bc_types {
//        abc_none,
        abc_inflow,
        abc_dirichlet,
        abc_total_flux,
        abc_diffusive_flux
    };

	/// Read necessary data from input record.
	virtual void init_from_input(const Input::Record &in_rec) = 0;

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
	 * Compute retardation coefficients due to sorption.
	 * @param point_list  Points at which to evaluate.
	 * @param ele_acc	  Element accessor.
	 * @param ret_coef    Coefficient vector (output).
	 */
	virtual void compute_retardation_coefficient(const std::vector<arma::vec3 > &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<std::vector<double> > &ret_coef) = 0;

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
	 * Return types of boundary conditions for each solution component.
	 * @param ele_acc  Element accessor.
	 * @param bc_types Vector of bc. types (output, see BC_Type)
	 */
	virtual void get_bc_type(const ElementAccessor<3> &ele_acc,
			arma::uvec &bc_types) = 0;

	/**
	 * \brief Return data for diffusive or total flux b.c.
	 *
	 * The flux can in general take the form
	 *
	 *   cross_section*(flux + sigma*(solution - ref_value))
	 *
     * @param index        Component index.
	 * @param point_list   Points at which to evaluate.
	 * @param ele_acc      Element accessor.
	 * @param bc_flux      Neumann flux (output).
	 * @param bc_sigma     Transition parameter (output).
	 * @param bc_ref_value Reference value (output).
	 */
	virtual void get_flux_bc_data(unsigned int index,
            const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector< double > &bc_flux,
			std::vector< double > &bc_sigma,
			std::vector< double > &bc_ref_value) = 0;

	/**
	 * \brief Return transition coefficient for flux b.c.
	 *
	 * In assembly of system matrix one does not teed all data for total/diffusive flux b.c.
	 * This method therefore returns only the sigma coefficient.
	 *
     * @param index        Component index.
	 * @param point_list   Points at which to evaluate.
	 * @param ele_acc      Element accessor.
	 * @param bc_sigma     Transition parameter (output).
	 */
	virtual void get_flux_bc_sigma(unsigned int index,
            const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector< double > &bc_sigma) = 0;

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






#endif /* AD_MODEL_HH_ */
