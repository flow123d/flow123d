/*
 * pade_approximant.h
 *
 * Author: PavelExner
 */

#ifndef PADE_APPROXIMANT_H_
#define PADE_APPROXIMANT_H_

#include "input/accessors.hh"
#include "armadillo"

using namespace arma;

/// Base class for the numerical method used in reaction.
class NumericalMethod
{
public:
    /**
     * Abstract record for numerical method.
     */
    static Input::Type::AbstractRecord input_type;
    
    /// Type of numerical method.
    typedef enum { analytic,            ///< Analytic solution is computed.
                   pade_approximant,    ///< Pade approximation.
    } Type;
};

/** @brief This class implements the Pade approximation of exponential function. 
 *
 * The exponential function is considered in the form \f$ e^{At} \f$ where \f$ A \f$ is a constant matrix.
 * It is then approximated by a fraction of polynomials
 * \f[ e^{At} = \frac{P(t)}{Q(t)},\f]
 * where the degrees of polynomials in nominator and denominator, \f$ P(t) \f$ and \f$ Q(t) \f$, are
 * set from the input record.
 * 
 */
class PadeApproximant : NumericalMethod
{
public:
    /**
     * Input record for class PadeApproximant.
     */
    static Input::Type::Record input_type;
    
    /// Constructor.
    PadeApproximant(Input::Record in_rec);

    /// Destructor.
    ~PadeApproximant(void);
    
    /**
    *   Evaluates Pade approximant from Reaction_matrix.
    */
    void approximate_matrix(mat &matrix);
    
protected:
    
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
    
    /// Computes factorial of @p k.
    //unsigned int factorial(int k);
    
    int nominator_degree_;      ///< Degree of the polynomial in the nominator.
    int denominator_degree_;    ///< Degree of the polynomial in the denominator.
};

#endif // PADE_APPROXIMANT_H_
