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
 * @file    linear_ode_solver.hh
 * @brief   
 */

#ifndef LINEAR_ODE_SOLVER_H_
#define LINEAR_ODE_SOLVER_H_

#include <boost/exception/detail/error_info_impl.hpp>  // for error_info
#include <boost/exception/info.hpp>                    // for operator<<
#include <iosfwd>                                      // for stringstream
#include <string>                                      // for string, basic_...
#include <vector>                                      // for vector
#include "armadillo"
#include "system/exc_common.hh"                        // for ExcAssertMsg
#include "system/exceptions.hh"                        // for ExcAssertMsg::...
#include "system/global_defs.h"                        // for msg, rank, OLD...
#include "system/asserts.hh"


/// @brief Class for linear ODE solver.
/** This class represents the solver of a system of linear ordinary differential 
 *  equations with constant coefficients which uses matrix exponential to compute
 *  the solution at given times.
 */
class LinearODESolver
{
public:
    
    LinearODESolver();
    ~LinearODESolver();
    
    void set_system_matrix(const arma::mat &matrix);  ///< Sets the matrix of ODE system.
    void set_step(double step);                 ///< Sets the step of the numerical method.
    
    /// Updates solution of the ODEs system.
    /**
     * @param init_vec is the column initial vector
     * @param output_vec is the column output vector containing the result
     */
    void update_solution(arma::vec &init_vec, arma::vec &output_vec);
    
    /// Estimate upper bound for time step. Return true if constraint was set.
     virtual bool evaluate_time_constraint(FMT_UNUSED double &time_constraint) { return false; }
                                 
protected:
    arma::mat system_matrix_;     ///< the square matrix of ODE system
    arma::mat solution_matrix_;   ///< the square solution matrix (exponential of system matrix)
    arma::vec rhs_;               ///< the column vector of RHS values (not used currently)
    double step_;           ///< the step of the numerical method
    bool step_changed_;     ///< flag is true if the step has been changed
    bool system_matrix_changed_; ///< Indicates that the system_matrix_ was recently updated.
};



#endif // LINEAR_ODE_SOLVER_H_
