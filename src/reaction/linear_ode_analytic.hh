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
 * @file    linear_ode_analytic.hh
 * @brief   
 */

#ifndef LINEAR_ODE_ANALYTIC_H_
#define LINEAR_ODE_ANALYTIC_H_

#include "reaction/linear_ode_solver.hh"

/** @brief This class implements the analytic solution of a system of linear ODEs with constant matrix.
 *
 * The analytic solution can be obtained in the special case due to a physical nature of the problem.
 * The problem is then solved only by a matrix multiplication.
 * 
 * The assumption is made that the equations are independent. Each quantity is decreased (supposing 
 * negative diagonal) to \f$ e^{a_{ii} t} \f$. The decrement \f$ \left( 1-e^{a_{ii} t} \right) \f$
 * is then distributed among other quantities according to the given fraction.
 * 
 * In case of the decays and first order reactions the elements of the solution matrix are:
 * \f{eqnarray*}{
 *      a_{ii} &=& e^{-\lambda_i t} \\
 *      a_{ji} &=& \left( 1-e^{-\lambda_i t} \right) b_{ji} \frac{M_j}{M_i}
 * \f}
 * where \f$ b_{ji} \f$ is the branching ratio of \f$ i \f$-th reactant and \f$ \frac{M_j}{M_i} \f$ is 
 * the fraction of molar masses.
 * 
 * The fractions \f$ b_{ji} \frac{M_j}{M_i} \f$ are then obtained from the system matrix by dividing 
 * \f$ -\frac{a_{ji}}{a_{ii}} \f$.
 * 
 * <B>Drawback:</B>  These assumptions (equation independence) are adequate when very small time step is 
 * applied. This will lead to huge amount of evaluations of the exponential functions which can be expensive,
 * so other numerical methods might be more appropriate. 
 * When the time step is large then the assumption is quite inadequate.
 * 
 */
class LinearODEAnalytic : public LinearODESolver<LinearODEAnalytic>
{
public:
	typedef LinearODESolverBase FactoryBaseType;

    /**
     * Input record for class LinearODE_analytic.
     */
    static const Input::Type::Record & get_input_type();
    
    ///Default constructor is possible because the input record is not needed.
    LinearODEAnalytic(){};
    
    /// Constructor from the input data.
    LinearODEAnalytic(Input::Record in_rec);

    /// Destructor.
    ~LinearODEAnalytic(void);
    
    void update_solution(arma::vec &init_vector, arma::vec &output_vec) override;
    
    bool evaluate_time_constraint(double &time_constraint) override;
    
protected:
    /**
     *   Computes the standard fundamental matrix.
     */
    void compute_matrix();
    
    /// The solution is computed only by a matrix multiplication (standard fundamental matrix).
    arma::mat solution_matrix_;
    
private:
    /// Registrar of class to factory
    static const int registrar;

};

#endif // LINEAR_ODE_ANALYTIC_H_
