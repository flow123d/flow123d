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
 * @file    reaction_term.hh
 * @brief   Class ReactionTerm is an abstract class representing reaction term in transport.
 *
 * Descending classes implements different physical models - dual porosity, sorption, decays, 
 * chemical reactions.
 */

#ifndef REACTION_TERM_H
#define REACTION_TERM_H

#include <memory>                    // for shared_ptr
#include <string>                    // for string
#include "coupling/equation.hh"      // for EquationBase
#include "system/index_types.hh"     // for LongInt
#include "input/input_exception.hh"  // for DECLARE_INPUT_EXCEPTION, Exception
#include "system/exceptions.hh"      // for ExcStream, operator<<, EI, TYPED...
#include "transport/substance.hh"    // for SubstanceList
#include "fem/dofhandler.hh"         // for DOFHandlerMultiDim
#include "fields/field_fe.hh"        // for FieldFE
#include "fields/field_set.hh"        // for FieldSet

class Distribution;
class Mesh;
class OutputTime;
namespace Input {
	class Record;
	namespace Type {
		class Abstract;
	}
}


class ReactionTerm: public EquationBase
{
public:
    TYPEDEF_ERR_INFO( EI_Substance, std::string);
    TYPEDEF_ERR_INFO( EI_Model, std::string);
    DECLARE_INPUT_EXCEPTION( ExcUnknownSubstance, << "Unknown substance name: " << EI_Substance::qval);
    DECLARE_INPUT_EXCEPTION( ExcWrongDescendantModel, << "Impossible descendant model: " << EI_Model::qval);
    
    typedef std::vector<std::shared_ptr<FieldFE< 3, FieldValue<3>::Scalar>>> FieldFEScalarVec;

    /**
     * Static variable for definition of common input record in reaction term.
     */
    static Input::Type::Abstract & it_abstract_term();
    static Input::Type::Abstract & it_abstract_mobile_term();
    static Input::Type::Abstract & it_abstract_immobile_term();
    static Input::Type::Abstract & it_abstract_reaction();
  
    /// Reaction data
    class EqFields : public FieldSet
    {
    public:
        /// Constructor
    	EqFields()
    	{}

        /// FieldFEs representing P0 interpolation of mobile concentration (passed from transport).
        FieldFEScalarVec conc_mobile_fe;
    };

    /// Reaction data
    class EqData
    {
    public:
        /// Constructor
        EqData();

        /**
         * Names belonging to substances.
         *
         * Must be same as in the transport.
         */
        SubstanceList substances_;

        /// Pointer to DOF handler used through the reaction tree
        std::shared_ptr<DOFHandlerMultiDim> dof_handler_;
    };

    /**
     * Constructor.
     *
     * @param init_mesh is the reference to the computational mesh
     * @param in_rec is the input record
     */
    ReactionTerm(Mesh &init_mesh, Input::Record in_rec);

    /// Destructor.
    ~ReactionTerm(void);
  

    ///@name Setters
    //@{
    ///Sets the names of substances considered in transport.
    ReactionTerm &substances(SubstanceList &substances)
    {eq_data_base_->substances_.initialize(substances); return *this;}

    ///Sets the output stream which is given from transport class.
    ReactionTerm &output_stream(std::shared_ptr<OutputTime> ostream)
    {output_stream_=ostream; return *this;}
  
    /**
     * Sets the pointer to concentration matrix for the mobile zone,
     * all substances and on all elements (given by transport).
     */
    ReactionTerm &concentration_fields(FieldFEScalarVec& conc_mobile)
    {
        ASSERT_PTR(eq_fields_base_);
        eq_fields_base_->conc_mobile_fe = conc_mobile;
        eq_data_base_->dof_handler_ = eq_fields_base_->conc_mobile_fe[0]->get_dofhandler();
        return *this;
    }
    //@}

    /** @brief Output method.
     *
     * Some reaction models have their own data to output (sorption, dual porosity)
     * - this is where it must be reimplemented.
     * On the other hand, some do not have (linear reaction, pade approximant)
     * - that is why it is not pure virtual.
     */
    virtual void output_data(void) override {};

    /// Disable changes in TimeGovernor by empty method.
    void choose_next_time(void) override;

protected:
    /// Compute reaction on a single element.
    virtual void compute_reaction(const DHCellAccessor& dh_cell) = 0;

    /// Pointer to a transport output stream.
    std::shared_ptr<OutputTime> output_stream_;

    /// Equation data - all fields needs in assembly class.
    std::shared_ptr<EqFields> eq_fields_base_;

    /// Equation data - all data needs in assembly class.
    std::shared_ptr<EqData> eq_data_base_;
};

#endif  // REACTION_TERM_H
