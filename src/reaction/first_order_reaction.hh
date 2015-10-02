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
 * @file    first_order_reaction.hh
 * @brief   
 */

#ifndef FIRST_ORDER_REACTION_H_
#define FIRST_ORDER_REACTION_H_

#include <vector>

#include "reaction/first_order_reaction_base.hh"
#include "input/accessors.hh"

class Mesh;

/** @brief Class implements the linear reactions.
 * 
 * This class implements the user interface for linear reactions and prepares the reaction matrix.
 * Common features are inherited from the @p FirstOrderReactionBase class.
 */
class FirstOrderReaction: public FirstOrderReactionBase
{
public:
	typedef ReactionTerm FactoryBaseType;

    static const Input::Type::Record & get_input_type();                 ///< Input record for class FirstOrderReaction.
    static const Input::Type::Record & get_input_type_single_reaction(); ///< Input record which defines particular reaction.
    static const Input::Type::Record & get_input_type_reactant();        ///< Input record for a reactant of a reaction.
    static const Input::Type::Record & get_input_type_product();         ///< Input record for a product of a reaction.

    /// Constructor.
    FirstOrderReaction(Mesh &init_mesh, Input::Record in_rec);

    /// Destructor.
    ~FirstOrderReaction(void);
    
protected:

    /// Implements the assembly of the system matrix of the ODEs.
    void assemble_ode_matrix(void) override;
    
    /// Initializes private members of sorption from the input record.
    void initialize_from_input() override;
    
    std::vector<double> reaction_rates_;    ///< Vector of reaction rates of the transported substances.

private:
    /// Registrar of class to factory
    static const int registrar;
};

#endif  // FIRST_ORDER_REACTION_H_
