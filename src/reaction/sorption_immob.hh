/** @brief class Sorption is used to enable simulation of sorption described by either linear or Langmuir isotherm in combination with limited solubility under consideration.
 *
 * Class in this file makes it possible to handle the dataset describing solid phase as either precipitated or sorbed species.
 *
 */
#ifndef SORPTION_IMMOB
#define SORPTION_IMMOB

#include <vector>

#include "input/input_type.hh"
#include "reaction/isotherm.hh"
#include "reaction/sorption_dual.hh"

class Mesh;
class Isotherm;

class SorptionImmob:  public SorptionDual
{
public:
    /**
     *  Constructor with parameter for initialization of a new declared class member
     *  TODO: parameter description
     */
	SorptionImmob(Mesh &init_mesh, Input::Record in_rec, vector<string> &names); //, pScalar mob_porosity, pScalar immob_porosity);
	/**
	*	Destructor.
	*/
	~SorptionImmob(void);
	/**
	*
	*/
	void isotherm_reinit(std::vector<Isotherm> &isotherms_vec, const ElementAccessor<3> &elem) override;
	
	//double compute_sorbing_scale(double por_m, double por_imm) override;
};

#endif
