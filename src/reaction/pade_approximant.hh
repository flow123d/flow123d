/*
 * pade_approximant.h
 *
 * Author: PavelExner
 */

#ifndef PADE_APPROXIMANT_H_
#define PADE_APPROXIMANT_H_

#include "input/accessors.hh"
#include "reaction/linear_ode_solver.hh"
#include "armadillo"

using namespace arma;

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
    /**
     * Input record for class PadeApproximant.
     */
    static Input::Type::Record input_type;
    
    /// Constructor from input record.
    PadeApproximant(Input::Record in_rec);
    
    /// Constructor.
    PadeApproximant(unsigned int nominator_degree, unsigned int denominator_degree);

    /// Destructor.
    ~PadeApproximant(void);
    
    void update_solution(vec &init_vector, vec &output_vec) override;
//     void update_solution(mat &init_vectors, mat &output_vecs, 
//                          const std::vector<unsigned int> &mask = std::vector<unsigned int>(0)) override;
    
protected:
    ///Hide default constructor.
    PadeApproximant(){};
    
    /**
     *   Approximate the matrix function.
     */
    void approximate_matrix(mat &matrix);
    
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
    void evaluate_matrix_polynomial(mat &polynomial_matrix, 
                                    const mat &input_matrix, 
                                    const std::vector<double> &coefs);
    
    int nominator_degree_;      ///< Degree of the polynomial in the nominator.
    int denominator_degree_;    ///< Degree of the polynomial in the denominator.
    
    mat solution_matrix_;       ///< Solution matrix \f$ e^{At} \f$.
};

#endif // PADE_APPROXIMANT_H_
