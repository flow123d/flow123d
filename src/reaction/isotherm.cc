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
	xprintf(MsgDbg,"adsorption_type %d\n",adsorption_type_);
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
	xprintf(MsgDbg,"interpolation_table.size() is %d\n", interpolation_table.size());
	return;
}



/*
SorptionType Isotherm::get_sorption_type(void)
{
	return adsorption_type_;
}

void Isotherm::set_iso_params(SorptionType sorp_type, double mult_coef, double second_coef)
{
	adsorption_type_ = sorp_type;
	mult_coef_ = mult_coef;
	second_coef_ = second_coef;
	return;
}

void Isotherm::set_kind_of_pores(int kind_of_pores)
{
	kind_of_pores_ = kind_of_pores;
	return;
}
*/
