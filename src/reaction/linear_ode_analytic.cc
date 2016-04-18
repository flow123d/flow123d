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
 * @file    linear_ode_analytic.cc
 * @brief   
 */

#include "reaction/linear_ode_analytic.hh"
#include "input/accessors.hh"
#include "input/factory.hh"
#include "system/sys_profiler.hh"

FLOW123D_FORCE_LINK_IN_CHILD(linearODEAnalytic)


using namespace Input::Type;

const Record & LinearODEAnalytic::get_input_type() {
    return Record("LinearODEAnalytic", "Evaluate analytic solution of the system of ODEs.")
    		.derive_from(LinearODESolverBase::get_input_type())
			.close();
}
    
const int LinearODEAnalytic::registrar =
		Input::register_class< LinearODEAnalytic, Input::Record >("LinearODEAnalytic") +
		LinearODEAnalytic::get_input_type().size();

LinearODEAnalytic::LinearODEAnalytic(Input::Record in_rec)
{
}

LinearODEAnalytic::~LinearODEAnalytic()
{
}

void LinearODEAnalytic::update_solution(arma::vec& init_vector, arma::vec& output_vec)
{
    if(step_changed_)
    {
        compute_matrix();
        step_changed_ = false;
    }
    
    output_vec = solution_matrix_ * init_vector;
}

void LinearODEAnalytic::compute_matrix()
{
    START_TIMER("ODEAnalytic::compute_matrix");

    ASSERT(system_matrix_.n_cols == system_matrix_.n_rows, "Matrix is not square.");
    solution_matrix_.copy_size(system_matrix_);
    
    double exponential, temp;
    for(unsigned int i = 0; i < solution_matrix_.n_rows; i++)
    {
        exponential = std::exp(system_matrix_(i,i) * step_);
        for(unsigned int j = 0; j < solution_matrix_.n_cols; j++)
        {
            temp = 1.0;
            if( (i != j) && (system_matrix_(i,i) != 0.0) )
                temp = -system_matrix_(j,i)/system_matrix_(i,i);
            
            solution_matrix_(j,i) = (1-exponential)*temp;
        }
        solution_matrix_(i,i) = exponential;
    }
}

bool LinearODEAnalytic::evaluate_time_constraint(double &time_constraint)
{
    if (!system_matrix_changed_) return false;

    system_matrix_changed_ = false;
    
    // set time constraint to 1/max(diag(reaction_matrix_))
    time_constraint = 0;
    for (unsigned int k=0; k<system_matrix_.n_rows; k++)
      time_constraint = std::max(time_constraint, -system_matrix_(k,k));
    
    time_constraint = 1 / time_constraint;
    
    return true;
}
