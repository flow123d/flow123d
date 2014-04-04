/** @brief Class Dual_por_exchange implements the model of dual porosity.
 *
 * It can be part of the transport model and it computes the concentrations of substances both in 
 * mobile and immobile zone. This model can also work above the sorption model - the sorbed concentration
 * is then computed both from mobile and immobile concentrations. Linear reactions can be define 
 * also in both zones.
 *
 */
#ifndef DUAL_POROSITY
#define DUAL_POROSITY

#include <vector>
#include <input/input_type.hh>

#include "fields/field_base.hh"
#include "fields/field_set.hh"
#include "./reaction/reaction.hh"

class Mesh;
class Distribution;
class SorptionBase;

class DualPorosity:  public ReactionTerm
{
public:
  /**
   * Static variable for new input data types input
   */
  static Input::Type::Record input_type;

  class EqData : public FieldSet // should be written in class Sorption
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

    static Input::Type::Selection output_selection;
  };

  DualPorosity(Mesh &init_mesh, Input::Record in_rec, vector<string> &names);
  /**
   * Destructor.
   */
  ~DualPorosity(void);
                
  /**
   * Updates the solution according to the dual porosity model.
   */
  void update_solution(void) override;
  
  /**
   * Initialization routines after all necessary members have been set.
   * It also sets and initializes possible following reaction models.
   */
  void initialize(OutputTime *stream) override;
  
  void output_data(void) override;
  void output_vector_gather(void) override;
  
  /**
   * Set the porosity field which is passed from transport.
   */
  void set_porosity(Field<3, FieldValue<3>::Scalar > &por_m);
  
  double **compute_reaction(double **concentrations, int loc_el) override;
  
protected:
  /**
   * This method disables to use constructor without parameters.
   */
  DualPorosity();

  void allocate_output_mpi(void);
  
  /**
   * Pointer to thwodimensional array[species][elements] containing concentrations either in immobile.
   */
  double **conc_immobile;

  /**
   * Equation data - all data field are in this set.
   */
  EqData data_;

  Input::Array output_array;

  /**
   * Input data set - fields in this set are read from the input file.
   */
  FieldSet input_data_set_;
  
  ReactionTerm *reaction_mobile;       ///< Reaction running in mobile zone
  ReactionTerm *reaction_immobile;     ///< Reaction running in immobile zone
  
  /** Minimal time for which the analytical solution of dual porosity concentrations are evaluated.
   * Else it is replaced with simple forward difference approximation.
   */
  static const double min_dt;
  
  ///@name members used in output routines
  //@{
  Vec *vconc_immobile; ///< PETSC concentration vector for immobile phase (parallel).
  Vec *vconc_immobile_out; ///< PETSC concentration vector output for immobile phase (gathered - sequential)
  double **conc_immobile_out; ///< concentration array output for immobile phase (gathered - sequential)  
  //@}
  
};

#endif  //DUAL_POROSITY
