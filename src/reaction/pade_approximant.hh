/*
 * pade_approximant.h
 *
 *  Created on: Apr 2, 2012
 *      Author: lukas
 */

#ifndef PADE_APPROXIMANT_H_
#define PADE_APPROXIMANT_H_

#include <input/input_type.hh>
#include "armadillo"

using namespace arma;

class Mesh;
class LinearReaction;

class PadeApproximant: public LinearReaction
{
public:
    /**
     * Input record for class PadeApproximant.
     */
    static Input::Type::Record input_type;
    /**
     * Input record which defines particular decay step.
     */
    static Input::Type::Record input_type_one_decay_substep;
   
    /// Constructor.
    PadeApproximant(Mesh &mesh, Input::Record in_rec);

    /// Destructor.
    ~PadeApproximant(void);

    void initialize() override;
    void zero_time_step() override;

protected:
    /**
    *   Evaluates Pade approximant from Reaction_matrix.
    */
    void modify_reaction_matrix(void) override;
    
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
     * @param reaction_matrix is the reaction matrix (with elements -kt)
     * @param coefs is the vector of coeficients of the polynomial
     */
    void evaluate_matrix_polynomial(mat &polynomial_matrix, 
                                    const mat &reaction_matrix, 
                                    const std::vector<double> &coefs);
    
    int nominator_degree_;      ///< Degree of the polynomial in the nominator.
    int denominator_degree_;    ///< Degree of the polynomial in the denominator.
};

#endif // PADE_APPROXIMANT_H_
