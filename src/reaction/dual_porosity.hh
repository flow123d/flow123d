/** @brief Class Dual_por_exchange implements the model of dual porosity.
 *
 * It can be part of the transport model and it computes the concentrations of substances both in 
 * mobile and immobile zone. This model can also work above the sorption model - the sorbed concentration
 * is then computed both from mobile and immobile concentrations. Linear reactions can be define 
 * also in both zones.
 *
 */
#ifndef DUAL_POROSITY_H_
#define DUAL_POROSITY_H_

#include <vector>
#include "input/accessors.hh"

#include "fields/field_algo_base.hh"
#include "fields/field_set.hh"
#include "fields/multi_field.hh"
#include "fields/vec_seq_double.hh"
#include "./reaction/reaction_term.hh"

class Mesh;

/// Class representing dual porosity model in transport.
class DualPorosity:  public ReactionTerm
{
public:
  typedef ReactionTerm FactoryBaseType;

  /**
   * Static variable for new input data types input
   */
  static const Input::Type::Record & get_input_type();

  /// DualPorosity data
  class EqData : public FieldSet
  {
  public:

    /// Collect all fields
    EqData();

    Field<3, FieldValue<3>::Vector > diffusion_rate_immobile;   ///< Mass transfer coefficients between mobile and immobile pores.
    Field<3, FieldValue<3>::Scalar > porosity_immobile;    ///< Immobile porosity field.
    
    Field<3, FieldValue<3>::Vector> init_conc_immobile; ///< Initial concentrations in the immobile zone. 

    Field<3, FieldValue<3>::Scalar > porosity; ///< Porosity field.
    
    MultiField<3, FieldValue<3>::Scalar>  conc_immobile;    ///< Calculated concentrations in the immobile zone.

    /// Fields indended for output, i.e. all input fields plus those representing solution.
    FieldSet output_fields;

    static const Input::Type::Selection & get_output_selection();
  };

  /// Constructor.
  DualPorosity(Mesh &init_mesh, Input::Record in_rec);
   
  ///Destructor.
  ~DualPorosity(void);

  /// Prepares the object to usage.
  /**
   * Allocating memory, reading input, initialization of fields.
   */
  void initialize() override;
  
  /**
   * Does first computation after initialization process.
   * The time is set and initial condition is set and output.
   */
  void zero_time_step() override;
  
  /**
   * Updates the solution according to the dual porosity model.
   */
  void update_solution(void) override;
  
  /// Main output routine.
  void output_data(void) override;
  
protected:
  /**
   * This method disables to use constructor without parameters.
   */
  DualPorosity();

  /// Resolves construction of following reactions.
  void make_reactions();
  
  /// Sets initial condition from input.
  void set_initial_condition();
  /// Initializes field sets.
  void initialize_fields();
  /// Allocates petsc vectors, prepares them for output and creates output vector scatter.
  void allocate_output_mpi(void);
  
  double **compute_reaction(double **concentrations, int loc_el) override;
  
  /// Gathers all the parallel vectors to enable them to be output.
  void output_vector_gather(void) override;
  
  /**
   * Pointer to twodimensional array[substance][elements] containing concentrations either in immobile.
   */
  double **conc_immobile;

  /**
   * Equation data - all data fields are in this set.
   */
  EqData data_;

  Input::Array output_array;

  /**
   * Input data set - fields in this set are read from the input file.
   */
  FieldSet input_data_set_;
  
  std::shared_ptr<ReactionTerm> reaction_mobile;       ///< Reaction running in mobile zone
  std::shared_ptr<ReactionTerm> reaction_immobile;     ///< Reaction running in immobile zone
  
  /// Dual porosity computational scheme tolerance. 
  /** According to this tolerance the analytical solution of dual porosity concentrations or
   * simple forward difference approximation of concentrations is chosen for computation.
   */
  double scheme_tolerance_;
  
  ///@name members used in output routines
  //@{
  VecScatter vconc_out_scatter; ///< Output vector scatter.
  Vec *vconc_immobile; ///< PETSC concentration vector for immobile phase (parallel).
  std::vector<VectorSeqDouble> conc_immobile_out; ///< concentration array output for immobile phase (gathered - sequential)
  //@}
  
private:
  /// Registrar of class to factory
  static const int registrar;

};

#endif  //DUAL_POROSITY_H_
