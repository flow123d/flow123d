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
 * @file    isotherm.cc
 * @brief   
 */

#include <utility>

#include "reaction/isotherm.hh"
#include "system/sys_profiler.hh"

Isotherm::Isotherm()
: table_limit_(0.0)
{
}

void Isotherm::clear_table()
{
    table_limit_ = 0.0;
    interpolation_table.clear();
}

void Isotherm::make_table(unsigned int n_points, double table_limit)
{
    START_TIMER("Isotherm::make_table");
        table_limit_ = table_limit;
	if(table_limit_ > 0.0) switch(adsorption_type_)
	{
		case 0: // none
		 {
			 None obj_isotherm;
			 make_table(obj_isotherm, 1);
		 }
		break;
		case 1: //  linear:
	 	 {
		 	Linear obj_isotherm(mult_coef_);
			make_table(obj_isotherm, n_points);
	 	 }
	 	 break;
	 	 case 2: // freundlich:
	 	 {
		 	Freundlich obj_isotherm(mult_coef_, second_coef_);
			make_table(obj_isotherm, n_points);
	 	 }
	 	 break;
	 	 case 3: // langmuir:
	 	 {
		 	Langmuir obj_isotherm(mult_coef_, second_coef_);
			make_table(obj_isotherm, n_points);
	 	 }
	 	 break;
	 	 default:
	 	 {
		 	 ;
	 	 }
	 	 break;
	}
	else
        {
           clear_table();
        }
}


template<> Isotherm::ConcPair Isotherm::solve_conc(Isotherm::ConcPair c_pair, const None &isotherm)
{
    return c_pair;
}

template<> void Isotherm::make_table(const None &isotherm, int n_steps)
{
    // Solve_conc returns the same, so we need to do that also in compute_projection.
    // We set size of the table to 1, so it follow the conditions into solve_conc again.
    
    limited_solubility_on_ = false; // so it cannot go in precipitate function
    
    total_mass_step_ = 1;            // set just one step in the table, so we void zero division
    interpolation_table.resize(1,0); // set one value in the table so the condition in compute_projection fails
    return;
}
