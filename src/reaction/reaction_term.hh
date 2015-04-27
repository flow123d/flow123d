/** @brief Class ReactionTerm is an abstract class representing reaction term in transport.
 *
 * Descending classes implements different physical models - dual porosity, sorption, decays, 
 * chemical reactions.
 */
#ifndef REACTION_TERM_H
#define REACTION_TERM_H

#include "coupling/equation.hh"
#include "transport/substance.hh"

class Mesh;
class Distribution;
class OutputTime;


class ReactionTerm: public EquationBase
{
public:
    TYPEDEF_ERR_INFO( EI_Substance, std::string);
    TYPEDEF_ERR_INFO( EI_Model, std::string);
    DECLARE_INPUT_EXCEPTION( ExcUnknownSubstance, << "Unknown substance name: " << EI_Substance::qval);
    DECLARE_INPUT_EXCEPTION( ExcWrongDescendantModel, << "Impossible descendant model: " << EI_Model::qval);
    
  /**
   * Static variable for definition of common input record in reaction term.
   */
  static Input::Type::AbstractRecord & get_input_type();
  
  /// Specification of the output record. 
  /**
   * Need not to be used by all reaction models, but they should
   * allow output of similar fields.
   */
  static Input::Type::Record input_type_output_record;

  /// Constructor.
  /** @param init_mesh is the reference to the computational mesh
   *  @param in_rec is the input record
   */
  ReactionTerm(Mesh &init_mesh, Input::Record in_rec);

  /// Destructor.
  ~ReactionTerm(void);
  

  ///@name Setters
  //@{
  ///Sets the names of substances considered in transport.
  ReactionTerm &substances(SubstanceList &substances)
  {substances_.initialize(substances); return *this;}

  ///Sets the output stream which is given from transport class.
  ReactionTerm &output_stream(std::shared_ptr<OutputTime> ostream)
  {output_stream_=ostream; return *this;}

  /**
   * Sets the pointer to concentration matrix for the mobile zone, 
   * all substances and on all elements (given by transport).
   */
  ReactionTerm &concentration_matrix(double **concentration, Distribution *conc_distr, 
                                     int *el_4_loc, int *row_4_el)
  {
    concentration_matrix_ = concentration;
    distribution_ = conc_distr;
    el_4_loc_ = el_4_loc;
    row_4_el_ = row_4_el;
    return *this;
  }
  //@}

  /** @brief Output method.
   * 
   * Some reaction models have their own data to output (sorption, dual porosity) 
   * - this is where it must be reimplemented.
   * On the other hand, some do not have (linear reaction, pade approximant) 
   * - that is why it is not pure virtual.
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
  virtual double **compute_reaction(double **concentrations, int loc_el) =0;

  /**
   * Pointer to two-dimensional array[species][elements] containing concentrations.
   */
  double **concentration_matrix_;
  
  /// Indices of elements belonging to local dofs.
  int *el_4_loc_;
  /// Indices of rows belonging to elements.
  int *row_4_el_;
  
  /// Pointer to reference to distribution of elements between processors.
  Distribution *distribution_;
  
  /// Names belonging to substances.
  /**  
   * Must be same as in the transport.
   */
  SubstanceList substances_;

  /// Pointer to a transport output stream.
  std::shared_ptr<OutputTime> output_stream_;

};

#endif  // REACTION_TERM_H
