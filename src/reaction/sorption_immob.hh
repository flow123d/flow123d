/** @brief class Sorption is used to enable simulation of sorption described by either linear or Langmuir isotherm in combination with limited solubility under consideration.
 *
 * Class in this file makes it possible to handle the dataset describing solid phase as either precipitated or sorbed species.
 *
 *
 * TODO:
 * sorption.cc:
 * - Why tests from line 151- 162
 * - use just one switch according to isotherm type
 * - what about time dependent sorption parameters?
 * - line 260: existence of appropriate table should be tested, faster and simpler
 *   even better have method is_precomputed() of the Isotherm class.
 *
 * - consider make sorbed concentration internaly also kg (of substnance) /kg (of rock) and
 *   convert it to mol/kg only on output, should be faster and safer
 *
 * - move all code that computes only values of one isotherm (of any type) to some other class
 *   proposed IsothermFactory
 *
 * - Idea is to have one Isotherm factory object, all necessary fields availabel in Sorption::EqDAta object
 *   and have virtual isotherm_reinit method, called during table initialization and/or in iterative solution.
 *
 *   Then we can have various derived Sorption classes for various purposes. Changing just the isotherm_reinit method we can
 *   change the parameters actually used in isotherms.
 *
 *   Need prototype of dual-porosity reaction class to design precise placement of various involved fields.
 *
 */
#ifndef SORPTION_IMMOB
#define SORPTION_IMMOB

#include <vector>
#include <input/input_type.hh>

#include "fields/field_base.hh"
#include "reaction/sorption_dual.hh"


class Mesh;
class Distribution;
class Reaction;
class Isotherm;
//class SorptionDual;

typedef Field<3, FieldValue<3>::Scalar > * pScalar;

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
	void isotherm_reinit(std::vector<Isotherm> &isotherms_vec, const ElementAccessor<3> &elem);
	/**
	*
	*/
	void set_immob_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc_);
	/**
	*
	*/
	void make_tables(void);
	/**
	*
	*/
	void update_solution(void) override;
protected:
	/**
	*	Pointer to thwodimensional array[species][elements] containing concentrations either in immobile.
	*/
	double **immob_concentration_matrix;
};

#endif
