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
 * @file    pade_approximant.hh
 * @brief   
 */

#ifndef PADE_APPROXIMANT_H_
#define PADE_APPROXIMANT_H_

#include "input/accessors_forward.hh"
#include "reaction/linear_ode_solver.hh"
#include "armadillo"

/** @brief This class implements the Pade approximation of exponential function. 
 *
 * The exponential function is considered in the form \f$ e^{At} \f$ where \f$ A \f$ is a constant matrix.
 * It is then approximated by a fraction of polynomials
 * \f[ e^{At} = \frac{P(t)}{Q(t)},\f]
 * where the degrees of polynomials in nominator and denominator, \f$ P(t) \f$ and \f$ Q(t) \f$, are
 * set from the input record.
 * 
 */
class PadeApproximant : public LinearODESolver<PadeApproximant>
{
public:
	typedef LinearODESolverBase FactoryBaseType;

    /**
     * Input record for class PadeApproximant.
     */
    static const Input::Type::Record & get_input_type();
    
    /// Constructor from input record.
    PadeApproximant(Input::Record in_rec);
    
    /// Constructor.
    PadeApproximant(unsigned int nominator_degree, unsigned int denominator_degree);

    /// Destructor.
    ~PadeApproximant(void);
    
    void update_solution(arma::vec &init_vector, arma::vec &output_vec) override;
    
    bool evaluate_time_constraint(double &time_constraint) override { return false; }
    
protected:
    ///Hide default constructor.
    PadeApproximant(){};
    
    /**
     *   Approximate the matrix function.
     */
    void approximate_matrix(arma::mat &matrix);
    
    /// Evaluates nominator and denominator coeficients of PadeApproximant for exponencial function.
    /** @param nominator_degree is the degree of polynomial in the nominator
     * @param denominator_degree is the degree of polynomial in the denominator
     * @param nominator_coefs is the vector of coeficients of the polynomial in the nominator
     * @param denominator_coefs is the vector of coeficients of the polynomial in the denominator
     */
    void compute_exp_coefs(unsigned int nominator_degree, unsigned int denominator_degree,
                           std::vector<double> &nominator_coefs, std::vector<double> &denominator_coefs);
    
    /// Evaluates the matrix polynomial by Horner scheme.
    /** @param polynomial_matrix is the output matrix
     * @param input_matrix is the input matrix (with elements -kt)
     * @param coefs is the vector of coeficients of the polynomial
     */
    void evaluate_matrix_polynomial(arma::mat &polynomial_matrix, 
                                    const arma::mat &input_matrix, 
                                    const std::vector<double> &coefs);
    
    int nominator_degree_;      ///< Degree of the polynomial in the nominator.
    int denominator_degree_;    ///< Degree of the polynomial in the denominator.
    
    arma::mat solution_matrix_;       ///< Solution matrix \f$ e^{At} \f$.

private:
    /// Registrar of class to factory
    static const int registrar;
};

#endif // PADE_APPROXIMANT_H_
