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



void Isotherm::reinit(enum SorptionType sorp_type, double rock_density, double rho_aqua, double por_m, double por_imm, double phi, double molar_mass, double c_aqua_limit)
{
	this->set_sorption_type(sorp_type);
	this->set_rho_aqua(rho_aqua);
	this->set_scales(por_m, por_imm, phi, rock_density, molar_mass);
    this->set_inv_scale_aqua(scale_aqua_/((scale_aqua_*scale_aqua_ + scale_sorbed_*scale_sorbed_)));
    this->set_inv_scale_sorbed(scale_sorbed_/((scale_aqua_*scale_aqua_ + scale_sorbed_*scale_sorbed_)));
    this->set_caq_limmit(c_aqua_limit);
}

bool Isotherm::compute_projection(double &c_aqua, double &c_sorbed)
{
    double total_mass = (scale_aqua_* c_aqua + scale_sorbed_ * c_sorbed);
    double total_mass_steps = total_mass / total_mass_step_;
    int total_mass_idx = static_cast <int>(std::floor(total_mass_steps));
    if ( total_mass_idx < 0 ) return false;
    if ( (unsigned int)(total_mass_idx) < interpolation_table.size() ) {
    	double rot_sorbed = interpolation_table[total_mass_idx] + (total_mass_steps - total_mass_idx)*(interpolation_table[total_mass_idx+1] - interpolation_table[total_mass_idx]);
        c_aqua = (total_mass * inv_scale_aqua_ - rot_sorbed * inv_scale_sorbed_); // (scale_aqua_ * total_mass * inv_scale_aqua_ - scale_sorbed_ * rot_sorbed * inv_scale_sorbed_);
        c_sorbed = (total_mass * inv_scale_sorbed_ + rot_sorbed * inv_scale_aqua_); // (scale_sorbed_ * total_mass * inv_scale_aqua_ + scale_aqua_ * rot_sorbed * inv_scale_sorbed_);
        return true;
    } else {
    	if (c_aqua_limit_ > 0.0) {
    		precipitate(c_aqua, c_sorbed); //, scale_aqua_, scale_sorbed_);
    	} else
        {
        	cout << "c_aqua_limit_ has the value " << c_aqua_limit_ << endl;
        	return false;
        }
    }

    return true; //false;
}

template<class Func>
void Isotherm::solve_conc(double &c_aqua, double &c_sorbed, const Func &isotherm)
{
    boost::uintmax_t max_iter = 20;
    dekl_tolerance<double> toler(30);
	double total_mass = (scale_aqua_*c_aqua + scale_sorbed_ * c_sorbed);
	double critic_total_mass = c_aqua_limit_*scale_aqua_ + const_cast<Func &>(isotherm)(c_aqua_limit_ / this->rho_aqua_)*scale_sorbed_;

	const double upper_solution_bound = critic_total_mass / scale_aqua_ + 0.00001; // corresponds to c_a^max, where substance stsarts to precipitate + 1.0

	if(total_mass < critic_total_mass)
	{
		// equation describing one point on the isotherm
		CrossFunction<Func> eq_func(isotherm, total_mass, scale_aqua_, scale_sorbed_, this->rho_aqua_);
		pair<double,double> solution = boost::math::tools::toms748_solve(eq_func, 0.0, upper_solution_bound, toler, max_iter);
		c_aqua = (solution.first + solution.second)/2; // = average of the pair solution defined above, midpoint
		c_sorbed = (total_mass - scale_aqua_ * c_aqua)/scale_sorbed_; //const_cast<Func &>(isotherm)(c_aqua);
	}else{
		precipitate(c_aqua, c_sorbed); //, scale_aqua_, scale_sorbed_);
	}

    return;
}

template void Isotherm::solve_conc<Linear>(double &c_aqua, double &c_sorbed, const Linear &isotherm);

template void Isotherm::solve_conc<Langmuir>(double &c_aqua, double &c_sorbed, const Langmuir &isotherm);

template void Isotherm::solve_conc<Freundlich>(double &c_aqua, double &c_sorbed, const Freundlich &isotherm);

void Isotherm::precipitate(double &c_aqua, double &c_sorbed) //, double scale_aqua, double scale_sorbed)
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

void Isotherm::make_one_point_table(void)
{
	interpolation_table.resize(1);
	interpolation_table[0] = 1.0;
	return;
}

void Isotherm::set_sorption_type(SorptionType sorp_type)
{
	sorption_type_ = sorp_type;
	return;
}

SorptionType Isotherm::get_sorption_type(void)
{
	return sorption_type_;
}

void Isotherm::set_mult_coef_(double mult_coef)
{
	mult_coef_ = mult_coef;
	return;
}

double Isotherm::get_mult_coef_(void)
{
	return mult_coef_;
}

void Isotherm::set_second_coef_(double second_coef)
{
	second_coef_ = second_coef;
	return;
}

double Isotherm::get_second_coef_(void)
{
	return second_coef_;
}

void Isotherm::set_rho_aqua(double rho_aqua)
{
	rho_aqua_ = rho_aqua;
	return;
}

void Isotherm::set_scale_aqua(double scale_aqua)
{
	scale_aqua_ = scale_aqua;
	return;
}

void Isotherm::set_inv_scale_aqua(double inv_scale_aqua)
{
	inv_scale_aqua_ = inv_scale_aqua;
	return;
}

void Isotherm::set_scale_sorbed(double scale_sorbed)
{
	scale_sorbed_ = scale_sorbed;
	return;
}

void Isotherm::set_inv_scale_sorbed(double inv_scale_sorbed)
{
	inv_scale_sorbed_ = inv_scale_sorbed;
	return;
}

void Isotherm::set_caq_limmit(double caq_limmit)
{
	c_aqua_limit_ = caq_limmit;
	return;
}

void Isotherm::set_scales(double por_m, double por_imm, double phi, double rock_density, double molar_mass)
{
	switch (kind_of_pores_)
	{
		case IMMOBILE :
		 scale_aqua_ = por_imm;
	 	 scale_sorbed_ = (1 - phi) * (1 - por_m - por_imm) * rock_density * molar_mass; // (1 - phi) * (1 - por_m - por_imm)
	 	 break;
		case MOBILE :
		 scale_aqua_ = por_m;
	 	 scale_sorbed_ = phi * (1 - por_m - por_imm) * rock_density * molar_mass;
	 	break;
		default :
			xprintf(UsrErr,"Unknown type of pores inside a rock matrix.\n");
	}
	return;
}

void Isotherm::set_kind_of_pores(int kind_of_pores)
{
	kind_of_pores_ = kind_of_pores;
	return;
}
