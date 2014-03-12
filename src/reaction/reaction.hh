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

#include "coupling/equation.hh"

class Mesh;
class Element;
class Distribution;


class Reaction: public EquationBase
{
public:
  
  /**
   * Static variable for definition of common input record in reactions.
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
  
  /** Some of the ascendants need to do some job after setting time, mesh, distribution etc.
   * This method can be overriden just for that purpose.
   */
  virtual void initialize(void) {};
  /**
   * For simulation of chemical raection in just one element either inside of MOBILE or IMMOBILE pores.
   */
  virtual double **compute_reaction(double **concentrations, int loc_el);
  /**
   * Virtual method that is reimplemented in ascendants. Computes new solution of the reaction model.
   */
  virtual void update_solution(void) = 0;
  
  ///@name Setters
  //@{
  /**
   * Sets the concentration matrix for the mobile zone, all substances and on all elements.
   */
  virtual void set_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc);
  
  /// Set mesh used by the model.
  inline void set_mesh(Mesh &mesh) 
    { mesh_ = &mesh; };
  
  /** Set names of substances.
   * TODO: this will not work with substance index mapping (initialized from input record)
   * Is this not redundant?
   */
  inline void set_names(const std::vector<string> &names)
    { names_ = names; };
  
  /**
   * TODO: implement in ascendants
   */
  virtual void set_concentration_vector(Vec &vec){};
  
  //@}
  
  /// @name Inherited and not used.
  //@{
  virtual void choose_next_time(void);
  virtual void get_parallel_solution_vector(Vec &vc);
  virtual void get_solution_vector(double* &vector, unsigned int &size);
  //@}
                
  Element * get_element_for_dof_index(unsigned int idx);
                
protected:

  void initialize_substance_ids(const std::vector<string> &names, Input::Record in_rec);

  /**
   * Pointer to threedimensional array[mobile/immobile][species][elements] containing concentrations.
   */
  double **concentration_matrix;
  
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
  
  unsigned int n_substances_;   //< number of substances that take part in the reaction model
  unsigned int n_all_substances_;   //< number of all substances in the transport model
  std::map<unsigned int, unsigned int> substance_id;    //< mapping from local indexing of substances to global
};

#endif
