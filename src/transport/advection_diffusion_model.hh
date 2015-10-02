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

	/// Read or set names of solution components.
	virtual void set_components(SubstanceList &substances, const Input::Record &in_rec) = 0;

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






#endif /* AD_MODEL_HH_ */
