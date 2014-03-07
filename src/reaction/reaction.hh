/** @brief class Linear_reaction is used to enable simulation of simple chemical reactions
 *
 * Class in this file makes it possible to realize  simulation of reaction of the first order by simple matrix multiplication.
 * One step of the linear reaction is represented as a product of a matrix containing concentrations of observed speciesin elements in rows multiplied by so called
 * reaction_matrix. Through this way radioactive decay can bee also realized and that was exactly what we did at the begining of journey. :-)
 * Matrix containing concentrations has a dimension Nxn, where N is a number of elements in mesh and n denotes a number of transported chemical species.
 * The reaction_matrix is a square matrix and it has a dimension nxn.
 *
 */
#ifndef REACT
#define REACT

#include "input/accessors.hh"
#include "coupling/equation.hh"
#include "mesh/elements.h"

class Mesh;
class Distribution;


class Reaction: public EquationBase
{
public:
  enum Reaction_type {No_reaction, Linear_react, Linear_react_Pade, General_react_Semch, Lim_Sorp}; 
  
  /**
   * Static variable for new input data types input
   */
  static Input::Type::AbstractRecord input_type;

  /**
   *  Constructor with parameter for initialization of a new declared class member
   *  TODO: parameter description
   */
  Reaction(Mesh &init_mesh, Input::Record in_rec, const std::vector<string> &names);
  /**
   * Destructor.
   */
  ~Reaction(void);
  
  /// Initialize from record in input file.
  virtual void init_from_input(Input::Record in_rec) = 0;
  /**
   * For simulation of chemical raection in just one element either inside of MOBILE or IMMOBILE pores.
   */
  virtual double **compute_reaction(double **concentrations, int loc_el) = 0;
  /**
   * Virtual method that is reimplemented in ascendants. Computes new solution of the reaction model.
   */
  virtual void update_solution(void) = 0;

  /**
   * Returns number of substances involved in reactions. This should be same as number of substances in transport.
   */
  inline unsigned int n_substances()
    { return names_.size(); }
  
  ///@name Setters
  //@{
  /**
   * Sets the whole concentration matrix for both mobile and immobile phase, all substances and on all elements.
   */
  void set_concentration_matrix(double ***ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc);
  
  /// Set mesh used by the model.
  void set_mesh(Mesh &mesh);
  
  /// Set names of substances.
  void set_names(const std::vector<string> &names);
  
  /**
   * TODO: implement in ascendants
   */
  virtual void set_concentration_vector(Vec &vec){};
  
  //@}
  
  /// @name Inherited and not used.
  //@{
  virtual void choose_next_time(void);
  virtual void set_time_step_constrain(double dt);
  virtual void get_parallel_solution_vector(Vec &vc);
  virtual void get_solution_vector(double* &vector, unsigned int &size);
  //@}
                
  Element * get_element_for_dof_index(unsigned int idx);
                
protected:
		/**
		*	This method disables to use constructor without parameters.
		*/
		Reaction();
		/**
		*	Finds a position of a string in specified array.
		*/
		unsigned int find_subst_name(const std::string &name);

  /**
   * Pointer to threedimensional array[mobile/immobile][species][elements] containing concentrations.
   */
  double ***concentration_matrix;
  /**
   * Distribution of elements between processors?
   */
  int *el_4_loc;
  /**
   * Pointer to reference to distribution of elements between processors.
   */
  Distribution *distribution;
  /**
   * Names belonging to substances. Should be same as in the transport.
   */
  vector<string> names_;
};

#endif
