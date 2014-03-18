/** @brief class Sorption is used to enable simulation of sorption described by either linear or Langmuir isotherm in combination with limited solubility under consideration.
 *
 * Class in this file makes it possible to handle the dataset describing solid phase as either precipitated or sorbed species.
 *
 *
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
     * This method will be implemented in descendants - it is different in each zone.
     */
    virtual void isotherm_reinit(std::vector<Isotherm> &isotherms, const ElementAccessor<3> &elm) = 0;

    
    /** Sets the immobile porosity field.
     */
    inline void set_porosity_immobile(Field<3, FieldValue<3>::Scalar > &por_imm)
      { 
        immob_porosity_.copy_from(por_imm); 
        data_.operator+=(immob_porosity_);
      }

protected:
    /**
     * This method disables to use constructor without parameters.
     */
    SorptionDual();
    
    Field<3, FieldValue<3>::Scalar > immob_porosity_; //< Immobile porosity field copied from transport
    
    //virtual double compute_sorbing_scale(double por_m, double por_imm) = 0;
};

#endif
