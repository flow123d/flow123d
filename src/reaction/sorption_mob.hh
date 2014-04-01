/** @brief class Sorption is used to enable simulation of sorption described by either linear or Langmuir isotherm in combination with limited solubility under consideration.
 *
 * Class in this file makes it possible to handle the dataset describing solid phase as either precipitated or sorbed species.
 *
 */
#ifndef SORPTION_DPMOB
#define SORPTION_DPMOB

#include "reaction/isotherm.hh"
#include "reaction/sorption_dual.hh"

class Mesh;
class Isotherm;

class SorptionMob:  public SorptionDual
{
public:

	class EqData : public SorptionBase::EqData
	{
	public:
		EqData();

		static Input::Type::Selection output_selection;
	};

	static Input::Type::Record input_type;

  /**
   *  Constructor with parameter for initialization of a new declared class member
   */
  SorptionMob(Mesh &init_mesh, Input::Record in_rec, vector<string> &names); //, pScalar mob_porosity, pScalar immob_porosity);
  /**
   * Destructor.
   */
  ~SorptionMob(void);
  
protected:
  /**
   *
   */
  void isotherm_reinit(std::vector<Isotherm> &isotherms_vec, const ElementAccessor<3> &elem) override;
  
  ///Renaming output names of substances.
  void set_output_names(void) override;
  
  Input::Type::Selection get_output_selection() override
  { return data_.output_selection; }

  SorptionBase::EqData &data() override
  { return data_; }

  //double compute_sorbing_scale(double por_m, double por_imm) override;

  EqData data_;

};

#endif
