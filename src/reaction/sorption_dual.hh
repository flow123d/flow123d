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
#ifndef SORPTION_DUAL
#define SORPTION_DUAL

#include <vector>
#include <input/input_type.hh>

#include "fields/field_base.hh"
#include "reaction/isotherm.hh"
#include "reaction/sorption_base.hh"

class Mesh;
class Distribution;
class Reaction;
class Isotherm;


class SorptionDual:  public SorptionBase
{
public:
  
    Field<3, FieldValue<3>::Scalar > immob_porosity_; //< Immobile porosity field copied from transport
    /** Fraction of the mobile porosity and the whole porosity, 
     * it was meant to be fraction of the total sorption surface exposed to the mobile zone, in interval (0,1).
     */
    Field<3, FieldValue<3>::Scalar > phi_; 
                
                



    /**
     *  Constructor with parameter for initialization of a new declared class member
     *  TODO: parameter description
     */
    SorptionDual(Mesh &init_mesh, Input::Record in_rec, vector<string> &names); //, pScalar mob_porosity, pScalar immob_porosity);
    /**
     * Destructor.
     */
    ~SorptionDual(void);
		/**
		*	For simulation of sorption in just one element either inside of MOBILE or IMMOBILE pores.
		*/
		//double **compute_reaction(double **concentrations, int loc_el);
		/**
		*
		*/
		virtual void isotherm_reinit(std::vector<Isotherm> &isotherms, const ElementAccessor<3> &elm);
		/**
		*	Prepared to compute sorption inside all of considered elements. It calls compute_reaction(...) for all the elements controled by concrete processor, when the computation is paralelized.
		*/
		virtual void update_solution(void) override;



    /** Sets the phi field.
     */
    inline void set_phi(Field<3, FieldValue<3>::Scalar > phi)
      { phi_.copy_from(phi); }
    
    /** Sets the immobile porosity field.
     */
    inline void set_porosity_immobile(Field<3, FieldValue<3>::Scalar > por_imm)
                  { immob_porosity_.copy_from(por_imm); }

                  
		virtual void make_tables(void);

	protected:
		/**
		*	This method disables to use constructor without parameters.
		*/
		SorptionDual();
};

#endif
