/** @brief class Sorption is used to enable simulation of sorption described by either linear or Langmuir isotherm in combination with limited solubility under consideration.
 *
 * Class in this file makes it possible to handle the dataset describing solid phase as either precipitated or sorbed species.
 *
 */
#ifndef SORPTION_IMMOB_H
#define SORPTION_IMMOB_H

#include "reaction/isotherm.hh"
#include "reaction/sorption_dual.hh"

class Mesh;
class Isotherm;

class SorptionImmob:  public SorptionDual
{
public:

	static Input::Type::Record input_type;

  /**
   *  Constructor with parameter for initialization of a new declared class member
   */
  SorptionImmob(Mesh &init_mesh, Input::Record in_rec);
  /**
   * Destructor.
   */
  ~SorptionImmob(void);

protected:
  /**
   *
   */
  void isotherm_reinit(std::vector<Isotherm> &isotherms_vec, const ElementAccessor<3> &elem) override;
  
  //double compute_sorbing_scale(double por_m, double por_imm) override;
};

#endif  //SORPTION_IMMOB_H
