#include "reaction/pade_approximant.hh"
#include "reaction/linear_reaction_base.hh"
#include "system/global_defs.h"

#include <input/accessors.hh>

#include "armadillo"

using namespace arma;
using namespace Input::Type;

AbstractRecord NumericalMethod::input_type
    = AbstractRecord("NumericalMethod", "Numerical method used in reaction computation.");
    
    
Record PadeApproximant::input_type
    = Record("PadeApproximant", "Record with an information about pade approximant parameters.")
    .derive_from( NumericalMethod::input_type)
    .declare_key("nominator_degree", Integer(), Default("2"),
                "Polynomial degree of the nominator of Pade approximant.")
    .declare_key("denominator_degree", Integer(), Default("2"),
                "Polynomial degree of the nominator of Pade approximant");

PadeApproximant::PadeApproximant(Input::Record in_rec)
{
    //DBGMSG("PadeApproximant constructor.\n");
    nominator_degree_ = in_rec.val<int>("nominator_degree");
    denominator_degree_ = in_rec.val<int>("denominator_degree");
    if(nominator_degree_ < 0) xprintf(UsrErr,"Wrong nominator degree in PadeApproximant.");
    if(denominator_degree_ < 0) xprintf(UsrErr,"Wrong denominator degree in PadeApproximant.");
}

PadeApproximant::~PadeApproximant()
{
}

void PadeApproximant::approximate_matrix(mat &matrix)
{
    ASSERT(matrix.n_rows == matrix.n_cols, "Matrix is not square.");
    
    unsigned int size = matrix.n_rows;
    
    //compute Pade Approximant
    mat nominator_matrix(size, size),
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

void PadeApproximant::evaluate_matrix_polynomial(mat& polynomial_matrix, 
                                                 const mat& reaction_matrix, 
                                                 const std::vector< double >& coefs)
{
    mat identity = eye(reaction_matrix.n_rows, reaction_matrix.n_cols);

    ///Horner scheme for evaluating polynomial a0 + [a1 + [a2 + [a3 +...]*R(t)]*R(t)]*R(t)
    for(int i = coefs.size()-1; i >= 0; i--)
    {
        polynomial_matrix = coefs[i] * identity + (polynomial_matrix * reaction_matrix);
    }
    //polynomial_matrix.print();
}
