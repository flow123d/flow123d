
#include <armadillo>

#include "system/global_defs.h"
#include "input/accessors.hh"
#include "input/factory.hh"
#include "system/sys_profiler.hh"
#include "reaction/linear_ode_solver.hh"
#include "reaction/pade_approximant.hh"


FLOW123D_FORCE_LINK_IN_CHILD(padeApproximant)


using namespace Input::Type;
    
    
Record & PadeApproximant::get_input_type() {
    static Record type = Record("PadeApproximant", "Record with an information about pade approximant parameters.")
    	.derive_from(LinearODESolverBase::get_input_type())
		.declare_key("nominator_degree", Integer(1), Default("2"),
                "Polynomial degree of the nominator of Pade approximant.")
		.declare_key("denominator_degree", Integer(1), Default("2"),
                "Polynomial degree of the nominator of Pade approximant")
		.close();

    return type;
}

const int PadeApproximant::registrar =
		( Input::register_class< PadeApproximant, Input::Record >("PadeApproximant"),
		LinearODESolverBase::get_input_type().add_child(PadeApproximant::get_input_type()) );

PadeApproximant::PadeApproximant(Input::Record in_rec)
{
    //DBGMSG("PadeApproximant constructor.\n");
    nominator_degree_ = in_rec.val<int>("nominator_degree");
    denominator_degree_ = in_rec.val<int>("denominator_degree");
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

    ASSERT(matrix.n_rows == matrix.n_cols, "Matrix is not square.");
    
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
