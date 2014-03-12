/** @brief class Dual_por_exchange is used to enable simulation of sorption described by either linear or Langmuir isotherm in combination with limited solubility under consideration.
 *
 * Class in this file makes it possible to handle the dataset describing solid phase as either precipitated or sorbed species.
 *
 */
#ifndef DUAL_POROSITY
#define DUAL_POROSITY

#include <vector>
#include <input/input_type.hh>

#include "fields/field_base.hh"
#include "./reaction/reaction.hh"

/// TODO: incorporate index mapping for substances indices

class Mesh;
class Distribution;
class Reaction;
class SorptionBase;

typedef Field<3, FieldValue<3>::Scalar > * pScalar;

class Dual_por_exchange:  public Reaction
{
public:
  /**
   * Static variable for new input data types input
   */
  static Input::Type::Record input_type;

  class EqData : public EqDataBase // should be written in class Sorption
  {
  public:

    /// Collect all fields
    EqData();

    Field<3, FieldValue<3>::Vector > alpha;            ///< Mass transfer coefficients between mobile and immobile pores.
    Field<3, FieldValue<3>::Scalar > immob_porosity;    ///< Immobile porosity
    
    Field<3, FieldValue<3>::Vector> init_conc_immobile; ///< Initial concentrations in the immobile zone. 

    pScalar porosity; ///< Pointer to mobile porosity
  };

  Dual_por_exchange(Mesh &init_mesh, Input::Record in_rec, vector<string> &names);
  /**
   * Destructor.
   */
  ~Dual_por_exchange(void);
                
  /**
   * Updates the solution according to the dual porosity model.
   */
  void update_solution(void) override;
  
  void set_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc_) override;
  /**
   *
   */
  void set_porosity(pScalar porosity);
  
  /// Initialize from input interface.
  void init_from_input(Input::Record in_rec) override;
  
  double **compute_reaction(double **concentrations, int loc_el) override;
  
protected:
  /**
   * This method disables to use constructor without parameters.
   */
  Dual_por_exchange();

  /**
   * Pointer to thwodimensional array[species][elements] containing concentrations either in immobile.
   */
  double **immob_concentration_matrix;

  /**
   *
   */
  EqData data_;
  
  
  Reaction *reaction_mob;       //< Reaction running in mobile zone
  Reaction *reaction_immob;     //< Reaction running in immobile zone
  
  /** Minimal time for which the analytical solution of dual porosity concentrations are evaluated.
   * Else it is replaced with simple forward difference approximation.
   */
  static const double min_dt;
  
};

#endif  //DUAL_POROSITY
