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
class OutputTime;


class Reaction: public EquationBase
{
public:
  
  /**
   * Static variable for definition of common input record in reactions.
   */
  static Input::Type::AbstractRecord input_type;
  
  /**
   * Specification of the output record. Need not to be used by all reaction models, but they should
   * allow output of similar fields.
   */
  static Input::Type::Record input_type_output_record;
    
  /**
   *  Constructor with parameter for initialization of a new declared class member
   *  TODO: parameter description
   */
  Reaction(Mesh &init_mesh, Input::Record in_rec, const std::vector<string> &names);
  /**
   * Destructor.
   */
  ~Reaction(void);
  
  /** Some of the ascendants need to do some job after setting time, mesh, distribution etc.
   * This method can be overriden just for that purpose.
   */
  virtual void initialize(void) {};
  
  /**
   * For simulation of chemical reaction in one element only. 
   * Inputs should be loc_el and local copies of concentrations of the element, which is then returned.
   */
  virtual double **compute_reaction(double **concentrations, int loc_el);
  /**
   * Virtual method that is reimplemented in ascendants. Computes new solution of the reaction model.
   */
  virtual void update_solution(void) = 0;
  
  /** Output method.
   * Some reaction models have their own data to output (sorption, dual porosity) - this is where it must be solved.
   * On the other hand, some do not have (linear reaction, pade approximant) - that is why it is not pure virtual.
   */
  virtual void output_data(void){};
  /**
   * Communicate parallel concentration vectors into sequential output vector.
   */
  virtual void output_vector_gather(void){};
  
  ///@name Setters
  //@{
  /**
   * Sets the concentration matrix for the mobile zone, all substances and on all elements.
   */
  virtual void set_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc, int *row_4_el);
  
  /// Set mesh used by the model.
  inline void set_mesh(Mesh &mesh) 
    { mesh_ = &mesh; };
  
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
                
protected:
  /** Initialize data from record in input file.
   * It is intended to use in ascendants.
   */
  virtual void init_from_input(Input::Record in_rec) {};
  
  void initialize_substance_ids(const std::vector<string> &names, Input::Record in_rec);

  /**
   * Pointer to two-dimensional array[species][elements] containing mobile concentrations.
   */
  double **concentration_matrix;
  
  /**
   * Indices of elements belonging to local dofs.
   */
  int *el_4_loc;
  /**
   * Indices of rows belonging to elements.
   */
  int *row_4_el;
  
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
  
  /// Record with output specification.
  Input::Record output_rec;

  OutputTime *output_stream;

};

#endif
