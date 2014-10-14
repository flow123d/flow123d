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
			 Linear obj_isotherm(0.0);
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

