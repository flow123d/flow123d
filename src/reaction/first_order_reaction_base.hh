#ifndef FIRST_ORDER_REACTION_BASE_H_
#define FIRST_ORDER_REACTION_BASE_H_

#include <vector>

#include "reaction/reaction_term.hh"
#include "input/accessors.hh"

#include "armadillo"

class Mesh;
class LinearODESolverBase;

/** @brief Base class for linear reactions and decay chain.
 *
 * The class implements common interface for linear reactions and decay chains.
 * One step of the linear reaction or the decay is represented as a product of a reaction matrix and
 * a vector of concentrations of transported substances on a single element.
 * 
 * It uses armadillo to compute the reaction matrix which then multiplies to concetration vector.
 * This class also resolves the choice of the numerical method which is used to compute the reaction matrix.
 */
class FirstOrderReactionBase: public ReactionTerm
{
public:
    /// Constructor.
    FirstOrderReactionBase(Mesh &init_mesh, Input::Record in_rec);

    /// Destructor.
    ~FirstOrderReactionBase(void);
                
    /// Prepares the object to usage.
    /**
     * Allocating memory, reading input, initialization of fields.
     */
    void initialize() override;
  
    /// Moves the model to zero time. 
    /** The assembly of the system matrix is called here.
     */
    void zero_time_step() override;
                
    /// Updates the solution. 
    /**
     * Goes through local distribution of elements and calls @p compute_reaction.
     */
    void update_solution(void) override;
    
protected:
    /// Assembles the matrix of the ODEs.
    /**
     * We solve the system of \f$N\f$ equations
     * \f[ 
     *  \frac{\textrm{d} c_i}{\textrm{d}t}=-\sum\limits_{j=1}^N \frac{M_i}{M_j} \lambda_{i} b_j c_j, \qquad \textrm i=1,\ldots,N
     * \f]
     * where \f$M_i, M_j\f$ are the molar masses of the parent substances and products, \f$\lambda_i\f$ are the 
     * reaction rate constants (in case of decays converted from half_lives) and \f$b_j\f$ are the branching ratios.
     * The constant coefficients \f$\mathbf{A}_{ij}=\frac{M_i}{M_j} \lambda_{i} b_j\f$ are the elements of the system matrix.
     * 
     * It is pure virtual and must be implemented in descendants.
     */
    virtual void assemble_ode_matrix(void) = 0;
    
    /// Computes the reaction on a specified element.
    virtual double **compute_reaction(double **concentrations, int loc_el) override;
            
    /// Initializes private members of sorption from the input record.
    virtual void initialize_from_input() = 0;
    
    /** Help function to create mapping of substance indices. 
     * Finds a position of a string in specified array.
     */
    unsigned int find_subst_name(const std::string &name);
    
    /**
    *   Sequence of integers describing an order of isotopes.
    *   substance_ids_[reactant][local_product_idx] = global_substance_idx
    */
    std::vector< std::vector <unsigned int> >substance_ids_;
    
    /**
    *   Two dimensional array contains mass percentage of every single decay bifurcation on every single row.
    */
    std::vector<std::vector<double> > bifurcation_;
    arma::mat bifurcation_matrix_;
    
    /// Number of all transported substances. It is the dimension of the reaction matrix.
    unsigned int n_substances_;
    
    arma::mat reaction_matrix_;   ///< Reaction matrix.
    arma::vec prev_conc_;      ///< Column vector storing previous concetrations on an element.
    
    arma::mat molar_matrix_;      ///< Diagonal matrix with molar masses of substances.
    arma::mat molar_mat_inverse_; ///< Inverse of @p molar_matrix_.

    LinearODESolverBase *linear_ode_solver_;
};

#endif  // FIRST_ORDER_REACTION_BASE_H_
