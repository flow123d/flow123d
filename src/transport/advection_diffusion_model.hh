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

namespace Input {
    class Record;
}


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

	/// Destructor.
	virtual ~AdvectionDiffusionModel() {};

};






#endif /* AD_MODEL_HH_ */
