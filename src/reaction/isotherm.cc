/*
 * isotherm.cc
 *
 *  Created on: Mar 15, 2013
 *      Author: lukas
 */
#include <utility>

//#include "system/sys_profiler.hh"
//#include "transport/transport.h"
#include "reaction/isotherm.hh"









void Isotherm::make_table(int nr_of_points)
{
	xprintf(Msg,"adsorption_type %d\n",adsorption_type_);
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
	xprintf(Msg,"interpolation_table.size() is %d\n", interpolation_table.size());
	return;
}



template<class Func>
void Isotherm::make_table(const Func &isotherm, int n_steps)
{
    double mass_limit;
    mass_limit = scale_aqua_ * table_limit_ + scale_sorbed_ * const_cast<Func &>(isotherm)(table_limit_ / this->rho_aqua_);
    if(mass_limit < 0.0)
    {
    	xprintf(UsrErr,"Isotherm mass_limit has negative value.\n");
    	//cout << "isotherm mass_limit has negative value " << mass_limit << ", scale_aqua "  << scale_aqua_ << ", c_aq_limit " << table_limit_ << ", scale_sorbed " << scale_sorbed_ << endl;
    }
    total_mass_step_ = mass_limit / n_steps;
    double mass = 0.0;
    for(int i=0; i<= n_steps; i++) {
    	 // aqueous concentration (original coordinates c_a) corresponding to i-th total_mass_step_
        ConcPair c_pair( mass/scale_aqua_, 0.0 );

        ConcPair result = solve_conc( c_pair, isotherm);
    	double c_sorbed_rot = ( result.solid * scale_aqua_ - result.fluid * scale_sorbed_);
        interpolation_table.push_back(c_sorbed_rot);
        mass = mass+total_mass_step_;
    }

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
