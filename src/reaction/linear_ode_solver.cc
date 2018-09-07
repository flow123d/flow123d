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

    
LinearODESolver::LinearODESolver()
: step_(0), step_changed_(true),
  system_matrix_changed_(false)
{
}

LinearODESolver::~LinearODESolver()
{
}

void LinearODESolver::set_system_matrix(const arma::mat& matrix)
{
    system_matrix_ = matrix;
    system_matrix_changed_ = true;
}

void LinearODESolver::set_step(double step)
{
    step_ = step;
    step_changed_ = true;
}

void LinearODESolver::update_solution(arma::vec& init_vector, arma::vec& output_vec)
{
    if(step_changed_ || system_matrix_changed_)
    {
        solution_matrix_ = arma::expmat(system_matrix_*step_);    //coefficients multiplied by time
        step_changed_ = false;
        system_matrix_changed_ = false;
    }
    
    output_vec = solution_matrix_ * init_vector;
}


