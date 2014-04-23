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


class ReactionTerm: public EquationBase
{
public:
  
  /**
   * Static variable for definition of common input record in reaction term.
   */
  static Input::Type::AbstractRecord input_type;
  
  /// Specification of the output record. 
  /**
   * Need not to be used by all reaction models, but they should
   * allow output of similar fields.
   */
  static Input::Type::Record input_type_output_record;

  /// Constructor.
  /** @param init_mesh is the reference to the computational mesh
   * @param in_rec is the input record
   */
  ReactionTerm(Mesh &init_mesh, Input::Record in_rec);

  /// Destructor.
  ~ReactionTerm(void);
  

  ///@name Setters
  //@{
  ///Sets the names of substances considered in transport.
  ReactionTerm &names(const std::vector<string> &names)
  {names_=names; return *this;}

  ///Sets the output stream which is given from transport class.
  ReactionTerm &output_stream(OutputTime &ostream)
  {output_stream_=&ostream; return *this;}

  /**
   * Sets the pointer to concentration matrix for the mobile zone, 
   * all substances and on all elements (given by transport).
   */
  ReactionTerm &concentration_matrix(double **concentration, Distribution *conc_distr, 
                                     int *el_4_loc, int *row_4_el)
  {
    concentration_matrix_ = concentration;
    distribution = conc_distr;
    this->el_4_loc = el_4_loc;
    this->row_4_el = row_4_el;
    return *this;
  }
  //@}

  /** @brief Output method.
   * 
   * Some reaction models have their own data to output (sorption, dual porosity) - this is where it must be solved.
   * On the other hand, some do not have (linear reaction, pade approximant) - that is why it is not pure virtual.
   */
  virtual void output_data(void){};

  /// Disable changes in TimeGovernor by empty method.
  void choose_next_time(void) override;

protected:
  /**
   * Communicate parallel concentration vectors into sequential output vector.
   */
  virtual void output_vector_gather(void){};

  /**
   * Computation of reaction term on a single element.
   * Inputs should be loc_el and local copies of concentrations of the element, which is then returned.
   */
  virtual double **compute_reaction(double **concentrations, int loc_el);

  /** Initialize data from record in input file.
   * It is intended to use in ascendants.
   */
  virtual void init_from_input(Input::Record in_rec) {};

  /**
   * Pointer to two-dimensional array[species][elements] containing concentrations.
   */
  double **concentration_matrix_;
  
  /// Indices of elements belonging to local dofs.
  int *el_4_loc;
  /// Indices of rows belonging to elements.
  int *row_4_el;
  
  /// Pointer to reference to distribution of elements between processors.
  Distribution *distribution;
  
  /// Names belonging to substances.
  /**  
   * Must be same as in the transport.
   */
  vector<string> names_;

  /// Pointer to a transport output stream.
  OutputTime *output_stream_;

};

#endif  // REACT
