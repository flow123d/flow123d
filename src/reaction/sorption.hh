/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    sorption.hh
 * @brief   This file contains classes representing sorption model.
 *          Sorption model can be computed both in case the dual porosity is considered or not.
 * 
 * The difference is only in the computation of scale_aqua and scale_sorbed.
 * Passing immobile porosity from dual porosity model is solved in abstract class SorptionDual.
 */

#ifndef SORPTION_H
#define SORPTION_H


#include <string>                     // for string
#include <vector>                     // for vector
#include "fields/field.hh"            // for Field
#include "fields/field_values.hh"     // for FieldValue<>::Scalar, FieldValue
#include "input/type_base.hh"         // for Array
#include "input/type_generic.hh"      // for Instance
#include "reaction/reaction_term.hh"  // for ReactionTerm
#include "reaction/sorption_base.hh"

class Mesh;
class Isotherm;
namespace Input {
	class Record;
	namespace Type { class Record; }
}
template <int spacedim> class ElementAccessor;


/** @brief Simple sorption model without dual porosity.
 * 
 */
class SorptionSimple:  public SorptionBase
{
public:
	typedef ReactionTerm FactoryBaseType;

    static const Input::Type::Record & get_input_type();

    /// Constructor.
    SorptionSimple(Mesh &init_mesh, Input::Record in_rec);
    
    /// Destructor.
    ~SorptionSimple(void);
  
protected:
    /// Implements @p SorptionBase::init_field_models
    void init_field_models() override;

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
	class EqFields : public SorptionBase::EqFields {
	public:
		EqFields(const string &output_field_name, const string &output_field_desc);

		Field<3, FieldValue<3>::Scalar > immob_porosity_; //< Immobile porosity field copied from transport
	};

	/// Constructor.
    SorptionDual(Mesh &init_mesh, Input::Record in_rec,
                const string &output_conc_name,
                const string &output_conc_desc);

    /// Destructor.
    ~SorptionDual(void);
    
    /// Sets the immobile porosity field.
    inline void set_porosity_immobile(Field<3, FieldValue<3>::Scalar > &por_imm)
    {
        eq_fields_dual_->immob_porosity_.copy_from(por_imm);
    }

protected:
    std::shared_ptr<EqFields> eq_fields_dual_;  ///< Overwrites SorptionBase::eq_fields_.
};


/** @brief Sorption model in mobile zone in case dual porosity is considered.
 * 
 */
class SorptionMob:  public SorptionDual
{
public:
	typedef ReactionTerm FactoryBaseType;

    static const Input::Type::Record & get_input_type();

    /// Constructor.
    SorptionMob(Mesh &init_mesh, Input::Record in_rec);
    
    /// Destructor.
    ~SorptionMob(void);
  
protected:
    /// Implements @p SorptionBase::init_field_models
    void init_field_models() override;

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

    static const Input::Type::Record & get_input_type();

    /// Constructor.
    SorptionImmob(Mesh &init_mesh, Input::Record in_rec);
    
    /// Destructor.
    ~SorptionImmob(void);

protected:
    /// Implements @p SorptionBase::init_field_models
    void init_field_models() override;

private:
    /// Registrar of class to factory
    static const int registrar;
};


#endif  // SORPTION_H
