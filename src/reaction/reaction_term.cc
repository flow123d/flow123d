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
 * @file    reaction_term.cc
 * @brief   
 */

#include "reaction/reaction_term.hh"
#include "system/global_defs.h"

#include "io/output_time.hh"

using namespace Input::Type;
        
Abstract & ReactionTerm::it_abstract_term() {
	return Abstract("ReactionTerm",
			"Abstract equation for a reaction term (dual porosity, sorption, reactions). Can be part of coupling with a transport equation via. operator splitting.")
			.close();
}

Abstract & ReactionTerm::it_abstract_mobile_term() {
    return Abstract("ReactionTermMobile",
            "Abstract equation for a reaction term of the MOBILE pores (sorption, reactions). Is part of dual porosity model.")
            .close();
}

Abstract & ReactionTerm::it_abstract_immobile_term() {
    return Abstract("ReactionTermImmobile",
            "Abstract equation for a reaction term of the IMMOBILE pores (sorption, reactions). Is part of dual porosity model.")
            .close();
}

Abstract & ReactionTerm::it_abstract_reaction() {
    return Abstract("GenericReaction",
            "Abstract equation for a reaction of species in single compartment (e.g. mobile solid)."
            "It can be part of: direct operator splitting coupling, dual porosity model, any sorption.")
            .close();
}


ReactionTerm::ReactionTerm(Mesh &init_mesh, Input::Record in_rec)
    : EquationBase(init_mesh, in_rec)
{
}

ReactionTerm::~ReactionTerm()
{
}




void ReactionTerm::choose_next_time(void)
{
  ASSERT(0).error("ReactionTerm does not change TimeGovernor.\n");
}
