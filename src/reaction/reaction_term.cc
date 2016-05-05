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
#include "mesh/mesh.h"

using namespace Input::Type;
        
Abstract & ReactionTerm::get_input_type() {
	return Abstract("ReactionTerm",
			"Equation for reading information about simple chemical reactions.")
			.close();
}

ReactionTerm::ReactionTerm(Mesh &init_mesh, Input::Record in_rec)
    : EquationBase(init_mesh, in_rec),
      concentration_matrix_(nullptr),
      el_4_loc_(nullptr),
      row_4_el_(nullptr),
      distribution_(nullptr)

{
}

ReactionTerm::~ReactionTerm()
{
}




void ReactionTerm::choose_next_time(void)
{
  OLD_ASSERT(0,"ReactionTerm does not change TimeGovernor.\n");
}
