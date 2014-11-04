#ifndef DECAY_CHAIN_H_
#define DECAY_CHAIN_H_

#include <vector>

#include "reaction/linear_reaction_base.hh"
#include "input/accessors.hh"

class Mesh;

/** @brief Class implements the radioactive decay chain.
 * 
 * This class implements the user interface for radioactive decay chain and prepares the reaction matrix,
 * The decay behaves like the linear reaction thus everything else is inherited from 
 * the @p FirstOrderReactionBase class.
 * 
 * TODO: fix the mass balance - we do not take in account the emitted particles.
 */
class RadioactiveDecay: public FirstOrderReactionBase
{
public:
    /**
     * Input record for class RadioactiveDecay.
     */
    static Input::Type::Record input_type;
    /**
     * Input record which defines particular decay step.
     */
    static Input::Type::Record input_type_single_decay;
   
    /// Constructor.
    RadioactiveDecay(Mesh &mesh, Input::Record in_rec);

    /// Destructor.
    ~RadioactiveDecay(void);


protected:
    /// Initializes private members of sorption from the input record.
    void initialize_from_input() override;
    
    /// Implements the assembly of the system matrix of the ODEs.
    void assemble_ode_matrix(void) override;
    
    std::vector<double> half_lives_;    ///< Half-lives of the substances.
};

#endif // DECAY_CHAIN_H_
