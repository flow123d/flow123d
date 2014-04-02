/** @brief class Sorption is used to enable simulation of sorption described by either linear or Langmuir isotherm in combination with limited solubility under consideration.
 *
 * Class in this file makes it possible to handle the dataset describing solid phase as either precipitated or sorbed species.
 *
 */
#ifndef SORPTION_SIMPLE
#define SORPTION_SIMPLE

#include "fields/field_base.hh"
#include "reaction/sorption_base.hh"

class Mesh;
class Distribution;
class Isotherm;

class SorptionSimple:  public SorptionBase
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
  SorptionSimple(Mesh &init_mesh, Input::Record in_rec, vector<string> &names); 
  /**
   * Destructor.
   */
  ~SorptionSimple(void);
  
protected:
  /**
   * This method disables to use constructor without parameters.
   */
  SorptionSimple();
  
  Input::Type::Selection get_output_selection() override
  { return data_.output_selection; }

  /**
   *
   */
  void isotherm_reinit(std::vector<Isotherm> &isotherms, const ElementAccessor<3> &elm) override;

  SorptionBase::EqData &data() override
  { return data_; }

  EqData data_;

};

#endif
