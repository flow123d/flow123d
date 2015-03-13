/*
 * isotherm.cc
 *
 *  Created on: Mar 15, 2013
 *      Author: lukas
 */
#include <utility>

#include "reaction/isotherm.hh"

void Isotherm::make_table(int nr_of_points)
{
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
			make_table(obj_isotherm, nr_of_points);
	 	 }
	 	 break;
	 	 case 2: // freundlich:
	 	 {
		 	Freundlich obj_isotherm(mult_coef_, second_coef_);
			make_table(obj_isotherm, nr_of_points);
	 	 }
	 	 break;
	 	 case 3: // langmuir:
	 	 {
		 	Langmuir obj_isotherm(mult_coef_, second_coef_);
			make_table(obj_isotherm, nr_of_points);
	 	 }
	 	 break;
	 	 default:
	 	 {
		 	 ;
	 	 }
	 	 break;
	}
	return;
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