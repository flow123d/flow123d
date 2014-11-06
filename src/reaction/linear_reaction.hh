/** @brief class Linear_reaction is used to enable simulation of simple chemical reactions
 *
 * Class in this file makes it possible to realize  simulation of reaction of the first order by simple matrix multiplication.
 * One step of the linear reaction is represented as a product of a matrix containing concentrations of observed speciesin elements in rows multiplied by so called
 * reaction_matrix. Through this way radioactive decay can bee also realized and that was exactly what we did at the begining of journey. :-)
 * Matrix containing concentrations has a dimension Nxn, where N is a number of elements in mesh and n denotes a number of transported chemical species.
 * The reaction_matrix is a square matrix and it has a dimension nxn.
 *
 */
#ifndef LINEAR_REACTION_H
#define LINEAR_REACTION_H

#include <vector>
#include <input/input_type.hh>

class Mesh;
class Distribution;
class ReactionTerm;

class LinearReaction: public ReactionTerm
{
public:
    /**
     * Static variable for new input data types input
     */
    static Input::Type::Record input_type;
    /**
     * Static variable gets information about particular decay step
     */
    static Input::Type::Record input_type_one_decay_substep;

    /// Constructor.
    LinearReaction(Mesh &init_mesh, Input::Record in_rec);

    /// Destructor.
    ~LinearReaction(void);
                
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

    /**
    *   This method modificates reaction matrix as described in ini-file a single section [Decay_i] or [FoReact_i]. It is used when bifurcation is switched off.
    */
    virtual void modify_reaction_matrix(void);
    
    /**
     *   For simulation of chemical reaction in just one element either inside of MOBILE or IMMOBILE pores.
     */
    virtual double **compute_reaction(double **concentrations, int loc_el) override;
        
    
    /// Resets reaction matrix as eye matrix.
    void reset_reaction_matrix();
            
    /// Initializes private members of sorption from the input record.
    void initialize_from_input();
	/**
	*	For control printing of a matrix describing simple chemical raections.
	*/
	void print_reaction_matrix(void);
	/**
	*	For printing nr_of_isotopes identifies of isotopes in a current decay chain.
	*/
	void print_indices(int dec_nr, int n_subst);

	/**
	*	For printing (nr_of_isotopes - 1) doubles containing half-lives belonging to particular isotopes on screen.
	*/
	void print_half_lives();
            
    /**
    *       Finds a position of a string in specified array.
    */
    unsigned int find_subst_name(const std::string &name);
    
	/**
	*	Small (nr_of_species x nr_of_species) square matrix for realization of radioactive decay and first order reactions simulation.
	*/
	std::vector<std::vector<double> > reaction_matrix_;
    /**
     *       Pointer to reference previous concentration array used in compute_reaction().
     */
    std::vector<double> prev_conc_;
	/**
	*	Sequence of (nr_of_isotopes - 1) doubles containing half-lives belonging to particular isotopes.
	*/
	std::vector<double> half_lives_;
	/**
	*	Sequence of integers describing an order of isotopes in decay chain or first order reaction.
	*/
	std::vector< std::vector <unsigned int> >substance_ids_;
	/**
	*	Two dimensional array contains mass percentage of every single decay bifurcation on every single row.
	*/
	std::vector<std::vector<double> > bifurcation_;
    
    unsigned int n_substances_;
};

#endif  // LINEAR_REACTION_H
