#ifndef LINEAR_REACTION_BASE_H
#define LINEAR_REACTION_BASE_H

#include <vector>
#include <ostream>

#include "reaction/reaction.hh"
#include "reaction/pade_approximant.hh"
#include "input/accessors.hh"

#include "armadillo"

using namespace arma;

class Mesh;
class ReactionTerm;
class PadeApproximant;

/** @brief Base class for linear reactions and decay chain.
 *
 * The class implements common interface for linear reactions and decay chains.
 * One step of the linear reaction or the decay is represented as a product of a reaction matrix and
 * a vector of concentrations of transported substances on a single element.
 * 
 * It uses armadillo to compute the reaction matrix which then multiplies to concetration vector.
 * This class also resolves the choice of the numerical method which is used to compute the reaction matrix.
 */
class LinearReactionBase: public ReactionTerm
{
public:
    /// Constructor.
    LinearReactionBase(Mesh &init_mesh, Input::Record in_rec);

    /// Destructor.
    ~LinearReactionBase(void);
                
    /// Prepares the object to usage.
    /**
     * Allocating memory, reading input, initialization of fields.
     */
    void initialize() override;
  
    void zero_time_step() override;
                
    /// Updates the solution. 
    /**
     * Goes through local distribution of elements and calls @p compute_reaction.
     */
    void update_solution(void) override;
    
protected:
    /// Initializes and prepares the reaction matrix.
    /**
     * It is pure virtual and must be implemented in descendants.
     */
    virtual void prepare_reaction_matrix(void) = 0;
    
    /// Computes the reaction matrix analyticaly.
    /**
     * Evaluation can be expensive.
     * It is pure virtual and must be implemented in descendants.
     */
    virtual void prepare_reaction_matrix_analytic(void) = 0;
    
    /// Chooses the numerical method and computes the reaction matrix.
    void compute_reaction_matrix(void);
    
    /// Computes the reaction on a specified element.
    virtual double **compute_reaction(double **concentrations, int loc_el) override;
            
    /// Initializes private members of sorption from the input record.
    virtual void initialize_from_input();
            
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
    
    /// Number of all transported substances. It is the dimension of the reaction matrix.
    unsigned int n_substances_;
    
    mat reaction_matrix_;   ///< Reaction matrix.
    colvec prev_conc_;      ///< Column vector storing previous concetrations on an element.
    
    NumericalMethod::Type numerical_method_;    ///< Numerical method selection.
    PadeApproximant *pade_approximant_;         ///< Pade approximant object.
};

#endif  // LINEAR_REACTION_H
