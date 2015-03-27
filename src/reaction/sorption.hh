/** @brief This file contains classes representing sorption model. 
 * 
 * Sorption model can be computed both in case the dual porosity is considered or not.
 * The difference is only in the isotherm_reinit method. 
 * Passing immobile porosity from dual porosity model is solved in abstract class SorptionDual.
 * 
 * TODO:
 * It seems that the methods isotherm_reinit() are different only at computation of scale_aqua and scale_sorbed.
 * So it could be moved to SorptionDual and the only method which would be virtual would be 
 * compute_sorbing_scale(). It is prepared in comment code.
 */
#ifndef SORPTION_H
#define SORPTION_H

#include "fields/field_algo_base.hh"
#include "reaction/sorption_base.hh"
#include "input/factory.hh"

class Mesh;
class Isotherm;

/** @brief Simple sorption model without dual porosity.
 * 
 */
class SorptionSimple:  public SorptionBase
{
public:
	typedef ReactionTerm FactoryBaseType;

    static Input::Type::Record input_type;

    /// Constructor.
    SorptionSimple(Mesh &init_mesh, Input::Record in_rec);
    
    /// Destructor.
    ~SorptionSimple(void);
  
protected:
    /// Reinitializes the isotherm.
    void isotherm_reinit(std::vector<Isotherm> &isotherms, const ElementAccessor<3> &elm) override;

private:
    /// Registrar of class to factory
    static const int registrar;
};


/** @brief Abstract class of sorption model in case dual porosity is considered.
 * 
 */
class SorptionDual:  public SorptionBase
{
public:
    /// Constructor.
    SorptionDual(Mesh &init_mesh, Input::Record in_rec,
                const string &output_conc_name,
                const string &output_selection_name);

    /// Destructor.
    ~SorptionDual(void);
    
    /// Sets the immobile porosity field.
    inline void set_porosity_immobile(Field<3, FieldValue<3>::Scalar > &por_imm)
      { 
        immob_porosity_.copy_from(por_imm); 
      }

protected:
    /// Reinitializes the isotherm.
    virtual void isotherm_reinit(std::vector<Isotherm> &isotherms, const ElementAccessor<3> &elm) = 0;
    
    Field<3, FieldValue<3>::Scalar > immob_porosity_; //< Immobile porosity field copied from transport

    //virtual double compute_sorbing_scale(double por_m, double por_imm) = 0;
};


/** @brief Sorption model in mobile zone in case dual porosity is considered.
 * 
 */
class SorptionMob:  public SorptionDual
{
public:
	typedef ReactionTerm FactoryBaseType;

    static Input::Type::Record input_type;

    /// Constructor.
    SorptionMob(Mesh &init_mesh, Input::Record in_rec);
    
    /// Destructor.
    ~SorptionMob(void);
  
protected:
    /// Reinitializes the isotherm.
    void isotherm_reinit(std::vector<Isotherm> &isotherms_vec, const ElementAccessor<3> &elem) override;

    //double compute_sorbing_scale(double por_m, double por_imm) override;

private:
    /// Registrar of class to factory
    static const int registrar;
};


/** @brief Sorption model in immobile zone in case dual porosity is considered.
 * 
 */
class SorptionImmob:  public SorptionDual
{
public:
	typedef ReactionTerm FactoryBaseType;

    static Input::Type::Record input_type;

    /// Constructor.
    SorptionImmob(Mesh &init_mesh, Input::Record in_rec);
    
    /// Destructor.
    ~SorptionImmob(void);

protected:
    /// Reinitializes the isotherm.
    void isotherm_reinit(std::vector<Isotherm> &isotherms_vec, const ElementAccessor<3> &elem) override;

    //double compute_sorbing_scale(double por_m, double por_imm) override;

private:
    /// Registrar of class to factory
    static const int registrar;
};


#endif  // SORPTION_H
