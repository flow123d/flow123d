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
 * @file    linear_ode_solver.cc
 * @brief   
 */

#include "reaction/linear_ode_solver.hh"

#include "armadillo"
#include "input/accessors.hh"

using namespace Input::Type;

Abstract & LinearODESolverBase::get_input_type() {
	return Abstract("LinearODESolver",
			"Solver of a linear system of ODEs.")
			.close();
}
    
LinearODESolverBase::LinearODESolverBase()
:step_(0), step_changed_(true)
{
}

LinearODESolverBase::~LinearODESolverBase()
{
}

void LinearODESolverBase::set_system_matrix(const arma::mat& matrix)
{
    system_matrix_ = matrix;
    step_changed_ = true;
}

void LinearODESolverBase::set_step(double step)
{
    step_ = step;
    step_changed_ = true;
}
