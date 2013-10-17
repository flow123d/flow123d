/*
 * isotherm.cc
 *
 *  Created on: Mar 15, 2013
 *      Author: lukas
 */
#include <utility>

#include "transport/transport.h"
#include "reaction/isotherm.hh"

void Linear::reinit(double mult_coef)
{
	mult_coef_ = mult_coef;
	return;
}

void Freundlich::reinit(double mult_coef, double exponent)
{
	mult_coef_ = mult_coef;
	exponent_ = exponent;
	return;
}

void Langmuir::reinit(double mult_coef, double alpha)
{
	mult_coef_ = mult_coef;
	alpha_ = alpha;
	return;
}

void Isotherm::reinit(enum SorptionType adsorption_type, double rho_aqua, double scale_aqua, double scale_sorbed, double c_aqua_limit, double mult_coef, double second_coef)
{
	adsorption_type_ = adsorption_type;
	rho_aqua_ = rho_aqua;
	scale_aqua_ = scale_aqua;
	scale_sorbed_ = scale_sorbed;
    inv_scale_aqua_ = scale_aqua_/(scale_aqua_*scale_aqua_ + scale_sorbed_*scale_sorbed_);
    inv_scale_sorbed_ = scale_sorbed_/(scale_aqua_*scale_aqua_ + scale_sorbed_*scale_sorbed_);
    c_aqua_limit_ = c_aqua_limit;
    mult_coef_ = mult_coef;
    second_coef_ = second_coef;
}

bool Isotherm::compute_projection(double &c_aqua, double &c_sorbed)
{
    double total_mass = (scale_aqua_* c_aqua + scale_sorbed_ * c_sorbed);
    double total_mass_steps = total_mass / total_mass_step_;
    int total_mass_idx = static_cast <int>(std::floor(total_mass_steps));
    xprintf(Msg,"total_mass %f, total_mass_idx %d, total_mass_step_ %f, scale_aqua_ %f, scale_sorbed_ %f, c_aqua %f, c_sorbed %f\n", total_mass, total_mass_idx, total_mass_step_, scale_aqua_, scale_sorbed_, c_aqua, c_sorbed);
    if ( total_mass_idx < 0 ) {xprintf(UsrErr,"total_mass %f\n", total_mass); }
    if ( (unsigned int)(total_mass_idx) < (interpolation_table.size() - 1) ) {
    	double rot_sorbed = interpolation_table[total_mass_idx] + (total_mass_steps - total_mass_idx)*(interpolation_table[total_mass_idx+1] - interpolation_table[total_mass_idx]);
        c_aqua = (total_mass * inv_scale_aqua_ - rot_sorbed * inv_scale_sorbed_);
        c_sorbed = (total_mass * inv_scale_sorbed_ + rot_sorbed * inv_scale_aqua_);
        return true;
    } else {
    	if (c_aqua_limit_ > 0.0) {
    		precipitate(c_aqua, c_sorbed);
    	} else {
    		ConcPair conc(c_aqua, c_sorbed);
    		conc = solve_conc(conc);
    		c_aqua = conc.first;
    		c_sorbed = conc.second;
    	}
    }

    return true;
}

template<class Func>
void Isotherm::solve_conc(double &c_aqua, double &c_sorbed, const Func &isotherm)
{
    boost::uintmax_t max_iter = 20;
    tolerance<double> toler(30);
	double total_mass = (scale_aqua_*c_aqua + scale_sorbed_ * c_sorbed);
	double critic_total_mass = c_aqua_limit_*scale_aqua_ + const_cast<Func &>(isotherm)(c_aqua_limit_ / this->rho_aqua_)*scale_sorbed_;

	const double upper_solution_bound = critic_total_mass / scale_aqua_ + 0.00001;

	if(total_mass < critic_total_mass)
	{
		// equation describing one point on the isotherm
		CrossFunction<Func> eq_func(isotherm, total_mass, scale_aqua_, scale_sorbed_, this->rho_aqua_);
		pair<double,double> solution = boost::math::tools::toms748_solve(eq_func, 0.0, upper_solution_bound, toler, max_iter);
		c_aqua = (solution.first + solution.second)/2;
		c_sorbed = (total_mass - scale_aqua_ * c_aqua)/scale_sorbed_;
	}else{
		precipitate(c_aqua, c_sorbed);
	}

    return;
}

template void Isotherm::solve_conc<Linear>(double &c_aqua, double &c_sorbed, const Linear &isotherm);

template void Isotherm::solve_conc<Langmuir>(double &c_aqua, double &c_sorbed, const Langmuir &isotherm);

template void Isotherm::solve_conc<Freundlich>(double &c_aqua, double &c_sorbed, const Freundlich &isotherm);

ConcPair Isotherm::solve_conc(ConcPair conc)
{
	double c_aqua = conc.first;
	double c_sorbed = conc.second;

	switch(adsorption_type_)
	{
		case 1: //  linear:
		{
			Linear obj_isotherm(mult_coef_);
			solve_conc(c_aqua, c_sorbed, obj_isotherm);
		}
		break;
		case 2: // freundlich
		{
			Freundlich obj_isotherm(mult_coef_, second_coef_);
			solve_conc(c_aqua, c_sorbed, obj_isotherm);
		}
		break;
		case 3:  // langmuir:
		{
			Langmuir obj_isotherm(mult_coef_, second_coef_);
			solve_conc(c_aqua, c_sorbed, obj_isotherm);
		}
		break;
		default:
		{
			;
		}
		break;
	}
	conc.first = c_aqua;
	conc.second = c_sorbed;

	return conc;
}

void Isotherm::make_table(int nr_of_points)
{
	switch(adsorption_type_)
	{
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

void Isotherm::precipitate(double &c_aqua, double &c_sorbed)
{
	double total_mass = (scale_aqua_*c_aqua + scale_sorbed_ * c_sorbed);

	c_aqua = c_aqua_limit_;
	c_sorbed = (total_mass - scale_aqua_ * c_aqua_limit_)/scale_sorbed_;

	return;
}

template<class Func>
void Isotherm::make_table(const Func &isotherm, int n_steps)
{
    double mass_limit;
    if (c_aqua_limit_ > 0.0) {
        mass_limit = scale_aqua_ * c_aqua_limit_ + scale_sorbed_ * const_cast<Func &>(isotherm)(c_aqua_limit_ / this->rho_aqua_);
        if(mass_limit < 0.0)
        {
        	cout << "isotherm mass_limit has negative value " << mass_limit << ", scale_aqua "  << scale_aqua_ << ", c_aq_limit " << c_aqua_limit_ << ", scale_sorbed " << scale_sorbed_ << endl;
        }
    } else {
        cout << "Solubility limit has to be higher than 0.0" << endl;
        return;
    }
    total_mass_step_ = mass_limit / n_steps;
    double mass = 0.0;
    for(int i=0; i<= n_steps; i++) {
        double c_aqua = mass/(scale_aqua_); // aqueous concentration (original coordinates c_a) corresponding to i-th total_mass_step_
        double c_sorbed = 0.0;
        solve_conc(c_aqua, c_sorbed, isotherm);
    	double c_sorbed_rot = ( c_sorbed * scale_aqua_ - c_aqua * scale_sorbed_);
        interpolation_table.push_back(c_sorbed_rot);
        mass = mass+total_mass_step_;
    }

    return;
}

template void Isotherm::make_table<Linear>(const Linear &isotherm, int n_steps);

template void Isotherm::make_table<Langmuir>(const Langmuir &isotherm, int n_steps);

template void Isotherm::make_table<Freundlich>(const Freundlich &isotherm, int n_steps);

int Isotherm::is_precomputed()
{
	return  interpolation_table.size();
}

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
