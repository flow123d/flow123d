#ifndef DECAY_CHAIN_H_
#define DECAY_CHAIN_H_

#include <vector>

#include "reaction/linear_reaction_base.hh"
#include "input/accessors.hh"

class Mesh;
class LinearReactionBase;

/** @brief Class implements the radioactive decay chain.
 * 
 * This class implements the user interface for radioactive decay chain and prepares the reaction matrix,
 * The decay behaves like the linear reaction thus everything else is inherited from 
 * the @p LinearReactionBase class.
 * 
 * TODO: fix the mass balance - we do not take in account the emitted particles.
 */
class DecayChain: public LinearReactionBase
{
public:
    /**
     * Input record for class DecayChain.
     */
    static Input::Type::Record input_type;
    /**
     * Input record which defines particular decay step.
     */
    static Input::Type::Record input_type_single_decay;
   
    /// Constructor.
    DecayChain(Mesh &mesh, Input::Record in_rec);

    /// Destructor.
    ~DecayChain(void);


protected:
    /// Initializes private members of sorption from the input record.
    void initialize_from_input() override;
    
    /// Implements reaction matrix initialization and preparation.
    void prepare_reaction_matrix(void) override;
    
    /// Implements reaction matrix analytical computation.
    void prepare_reaction_matrix_analytic(void) override;
    
    std::vector<double> half_lives_;    ///< Half-lives of the substances.
};

#endif // DECAY_CHAIN_H_
