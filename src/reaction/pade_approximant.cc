#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "system/global_defs.h"

#include "mesh/mesh.h"
#include "armadillo"

using namespace arma;
using namespace std;
using namespace Input::Type;

Record PadeApproximant::input_type_one_decay_substep
	= Record("Substep", "Equation for reading information about radioactive decays.")
	.declare_key("parent", String(), Default::obligatory(),
				"Identifier of an isotope.")
    .declare_key("half_life", Double(), Default::optional(),
                "Half life of the parent substance.")
    .declare_key("kinetic", Double(), Default::optional(),
                "Kinetic constants describing first order reactions.")
	.declare_key("products", Array(String()), Default::obligatory(),
				"Identifies isotopes which decays parental atom to.")
	.declare_key("branch_ratios", Array(Double()), Default("1.0"),  //default is one product, with ratio = 1.0
				"Decay chain branching percentage.");


Record PadeApproximant::input_type
	= Record("PadeApproximant", "Abstract record with an information about pade approximant parameters.")
	.derive_from( ReactionTerm::input_type )
    .declare_key("decays", Array( PadeApproximant::input_type_one_decay_substep ), Default::obligatory(),
                "Description of particular decay chain substeps.")
	.declare_key("nom_pol_deg", Integer(), Default("2"),
				"Polynomial degree of the nominator of Pade approximant.")
	.declare_key("den_pol_deg", Integer(), Default("2"),
				"Polynomial degree of the nominator of Pade approximant");


PadeApproximant::PadeApproximant(Mesh &init_mesh, Input::Record in_rec)
      : LinearReaction(init_mesh, in_rec)
{
}

PadeApproximant::~PadeApproximant()
{
}

void PadeApproximant::initialize()
{
    LinearReaction::initialize();
    
    // init from input
    nominator_degree_ = input_record_.val<int>("nom_pol_deg");
    denominator_degree_ = input_record_.val<int>("den_pol_deg");
    if((nominator_degree_ + denominator_degree_) < 0){
        xprintf(UsrErr, "You did not specify Pade approximant required polynomial degrees.");
    }
}

void PadeApproximant::zero_time_step()
{
    LinearReaction::zero_time_step();
}

void PadeApproximant::modify_reaction_matrix(void )
{
    // create decay matrix
    mat r_reaction_matrix_ = zeros(n_substances_, n_substances_);
    unsigned int reactant_index, product_index; //global indices of the substances
    double exponent;    //temporary variable
    for (unsigned int i_decay = 0; i_decay < half_lives_.size(); i_decay++) {
        reactant_index = substance_ids_[i_decay][0];
        exponent = log(2) * time_->dt() / half_lives_[i_decay];
        r_reaction_matrix_(reactant_index, reactant_index) = -exponent;
        
        for (unsigned int i_product = 1; i_product < substance_ids_[i_decay].size(); ++i_product){
            product_index = substance_ids_[i_decay][i_product];
            r_reaction_matrix_(product_index, reactant_index) = exponent * bifurcation_[i_decay][i_product-1];
        }
    }
    //DBGMSG("reactions_matrix_created\n");
    //r_reaction_matrix_.print();
    
    //compute Pade Approximant
    mat nominator_matrix(n_substances_, n_substances_),
        denominator_matrix(n_substances_, n_substances_),
        pade_approximant_matrix(n_substances_, n_substances_);
        
    nominator_matrix.fill(0);
    denominator_matrix.fill(0);
    pade_approximant_matrix.fill(0);

    std::vector<double> nominator_coefs(nominator_degree_+1),
                        denominator_coefs(denominator_degree_+1);
    
    // compute Pade approximant polynomials for the function e^x
    compute_exp_coefs(nominator_degree_, denominator_degree_, nominator_coefs, denominator_coefs);  
    // evaluation of polynomials of Pade approximant where x = -kt = R
    evaluate_matrix_polynomial(nominator_matrix, r_reaction_matrix_, nominator_coefs);
    evaluate_matrix_polynomial(denominator_matrix, r_reaction_matrix_, denominator_coefs);
    // compute P(R(t)) / Q(R(t))
    pade_approximant_matrix = nominator_matrix * inv(denominator_matrix);
    //pade_approximant_matrix.print();
    
    // write matrix to reaction matrix
    unsigned int rows, cols;
    for(rows = 0; rows < n_substances_; rows++)
    {
        for(cols = 0; cols < n_substances_ ; cols++)
        {
            reaction_matrix_[rows][cols] = pade_approximant_matrix(rows,cols);
        }
    }
    //print_reaction_matrix();
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
        //DBGMSG("p(%d)=%f\n",j,nominator_coefs[j]);
    }

    for(int i = denominator_degree; i >= 0; i--)
    {
        if(i % 2 == 0) sign = 1; else sign = -1;
        denominator_coefs[i] = sign * 
            (double)(factorials[nominator_degree + denominator_degree - i] * factorials[denominator_degree])
            / (factorials[nominator_degree + denominator_degree] * factorials[i] * factorials[denominator_degree - i]);
        //DBGMSG("q(%d)=%f\n",i,denominator_coefs[i]);
    } 
}

void PadeApproximant::evaluate_matrix_polynomial(mat& polynomial_matrix, 
                                                 const mat& reaction_matrix, 
                                                 const std::vector< double >& coefs)
{
    //DBGMSG("evaluate_matrix_polynomial\n");
    mat identity = eye(n_substances_, n_substances_);

    ///Horner scheme for evaluating polynomial a0 + [a1 + [a2 + [a3 +...]*R(t)]*R(t)]*R(t)
    for(int i = coefs.size()-1; i >= 0; i--)
    {
        polynomial_matrix = coefs[i] * identity + (polynomial_matrix * reaction_matrix);
    }
    //polynomial_matrix.print();
}

// unsigned int PadeApproximant::factorial(int k)
// {
//     ASSERT(k >= 0, "Cannot compute factorial of negative number.");
//     
//     unsigned int fact = 1;
//     while(k > 1)
//     {
//             fact *= k;
//             k--;
//     }
//     return fact;
// }
