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
 * @file    pade_approximant.cc
 * @brief   
 */

#include <armadillo>

#include "system/global_defs.h"
#include "input/accessors.hh"
#include "input/factory.hh"
#include "system/sys_profiler.hh"
#include "reaction/linear_ode_solver.hh"
#include "reaction/pade_approximant.hh"


// FLOW123D_FORCE_LINK_IN_CHILD(padeApproximant)


using namespace Input::Type;
    
    
const Record & PadeApproximant::get_input_type() {
    return Record("PadeApproximant", 
                  "Record with an information about pade approximant parameters."
                  "Note that stable method is guaranteed only if d-n=1 or d-n=2, "
                  "where d=degree of denominator and n=degree of nominator. "
                  "In those cases the Pade approximant corresponds to an implicit "
                  "Runge-Kutta method which is both A- and L-stable. "
                  "The default values n=2, d=3 yield relatively good precision "
                  "while keeping the order moderately low.")
//     	.derive_from(LinearODESolverBase::get_input_type())
		.declare_key("pade_nominator_degree", Integer(1), Default("1"),
                "Polynomial degree of the nominator of Pade approximant.")
		.declare_key("pade_denominator_degree", Integer(1), Default("3"),
                "Polynomial degree of the denominator of Pade approximant")
		.close();
}

const int PadeApproximant::registrar =
		Input::register_class< PadeApproximant, Input::Record >("PadeApproximant") +
		PadeApproximant::get_input_type().size();

PadeApproximant::PadeApproximant(Input::Record in_rec)
{
    nominator_degree_ = in_rec.val<int>("pade_nominator_degree");
    denominator_degree_ = in_rec.val<int>("pade_denominator_degree");
    if (nominator_degree_+1 != denominator_degree_ &&
        nominator_degree_+2 != denominator_degree_)
    	WarningOut() << "Pade approximation can be unstable since (denominator_degree-nominator_degree) is not 1 or 2.\n";
}

PadeApproximant::PadeApproximant(unsigned int nominator_degree, unsigned int denominator_degree)
:   nominator_degree_(nominator_degree), denominator_degree_(denominator_degree)
{
}

PadeApproximant::~PadeApproximant()
{
}

void PadeApproximant::update_solution(arma::vec& init_vector, arma::vec& output_vec)
{
    if(step_changed_)
    {
        solution_matrix_ = system_matrix_*step_;    //coefficients multiplied by time
        approximate_matrix(solution_matrix_);
        step_changed_ = false;
    }
    
    output_vec = solution_matrix_ * init_vector;
}


void PadeApproximant::approximate_matrix(arma::mat &matrix)
{
    START_TIMER("ODEAnalytic::compute_matrix");

    OLD_ASSERT(matrix.n_rows == matrix.n_cols, "Matrix is not square.");
    
    unsigned int size = matrix.n_rows;
    
    //compute Pade Approximant
    arma::mat nominator_matrix(size, size),
              denominator_matrix(size, size);
        
    nominator_matrix.fill(0);
    denominator_matrix.fill(0);

    std::vector<double> nominator_coefs(nominator_degree_+1),
                        denominator_coefs(denominator_degree_+1);
    
    // compute Pade approximant polynomials for the function e^x
    compute_exp_coefs(nominator_degree_, denominator_degree_, nominator_coefs, denominator_coefs);  
    // evaluation of polynomials of Pade approximant where x = -kt = R
    evaluate_matrix_polynomial(nominator_matrix, matrix, nominator_coefs);
    evaluate_matrix_polynomial(denominator_matrix, matrix, denominator_coefs);
    // compute P(R(t)) / Q(R(t))
    matrix = nominator_matrix * inv(denominator_matrix);
}

void PadeApproximant::compute_exp_coefs(unsigned int nominator_degree, 
                                        unsigned int denominator_degree, 
                                        std::vector< double >& nominator_coefs, 
                                        std::vector< double >& denominator_coefs)
{
    // compute factorials in forward
    std::vector<unsigned int> factorials(nominator_degree+denominator_degree+1);
    factorials[0] = 1;
    for(unsigned int i = 1; i < factorials.size(); i++)
        factorials[i] = factorials[i-1]*i;
    
    int sign;   // variable for denominator sign alternation
    
    for(int j = nominator_degree; j >= 0; j--)
    {
        nominator_coefs[j] = 
            (double)(factorials[nominator_degree + denominator_degree - j] * factorials[nominator_degree]) 
            / (factorials[nominator_degree + denominator_degree] * factorials[j] * factorials[nominator_degree - j]);
    }

    for(int i = denominator_degree; i >= 0; i--)
    {
        if(i % 2 == 0) sign = 1; else sign = -1;
        denominator_coefs[i] = sign * 
            (double)(factorials[nominator_degree + denominator_degree - i] * factorials[denominator_degree])
            / (factorials[nominator_degree + denominator_degree] * factorials[i] * factorials[denominator_degree - i]);
    } 
}

void PadeApproximant::evaluate_matrix_polynomial(arma::mat& polynomial_matrix, 
                                                 const arma::mat& input_matrix, 
                                                 const std::vector< double >& coefs)
{
    arma::mat identity = arma::eye(input_matrix.n_rows, input_matrix.n_cols);

    ///Horner scheme for evaluating polynomial a0 + [a1 + [a2 + [a3 +...]*R(t)]*R(t)]*R(t)
    for(int i = coefs.size()-1; i >= 0; i--)
    {
        polynomial_matrix = coefs[i] * identity + (polynomial_matrix * input_matrix);
    }
    //polynomial_matrix.print();
}
