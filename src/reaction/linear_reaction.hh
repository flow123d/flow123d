#ifndef LINEAR_REACTION_H
#define LINEAR_REACTION_H

#include <vector>

#include "reaction/linear_reaction_base.hh"
#include "input/accessors.hh"

class Mesh;
class Distribution;

/** @brief Class implements the linear reactions.
 * 
 * This class implements the user interface for linear reactions and prepares the reaction matrix.
 * Common features are inherited from the @p LinearReactionBase class.
 */
class LinearReaction: public LinearReactionBase
{
public:
    /**
     * Static variable for new input data types input
     */
    static Input::Type::Record input_type;
    /**
     * Static variable gets information about particular decay step
     */
    static Input::Type::Record input_type_single_reaction;

    /// Constructor.
    LinearReaction(Mesh &init_mesh, Input::Record in_rec);

    /// Destructor.
    ~LinearReaction(void);
    
protected:

    /// Implements reaction matrix initialization and preparation.
    void prepare_reaction_matrix(void) override;
    
    /// Implements reaction matrix analytical computation.
    void prepare_reaction_matrix_analytic(void) override;
    
    /// Initializes private members of sorption from the input record.
    void initialize_from_input();

    /// Prints the half-lives on the screen.
    void print_half_lives();
    
    std::vector<double> half_lives_;    ///< Vector of half-lives of the transported substances.
};

#endif  // LINEAR_REACTION_H
