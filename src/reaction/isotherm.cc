/*
 * isotherm.cc
 *
 *  Created on: Mar 15, 2013
 *      Author: lukas
 */
#include <utility>

#include "reaction/isotherm.hh"

void Linear::reinit(double mult_coef)
{
	mult_coef_ = mult_coef;
	return;
}

/*Isotherm::Isotherm(void)
{
}*/

void Isotherm::reinit(enum SorptionType sorp_type, double rock_density, double rho_aqua, double porosity, double molar_mass, double c_aqua_limit)
{
    // set class variables
	sorption_type_ = sorp_type;
    scale_aqua_ = porosity * rho_aqua;
    scale_sorbed_ = (1-porosity) * rock_density * molar_mass;
    inv_scale_aqua_ = 1/((scale_aqua_*scale_aqua_ + scale_sorbed_*scale_sorbed_)); // scale_aqua/(scale_aqua*scale_aqua + scale_sorbed*scale_sorbed);
    inv_scale_sorbed_ = 1/((scale_aqua_*scale_aqua_ + scale_sorbed_*scale_sorbed_)); // scale_sorbed/(scale_aqua*scale_aqua + scale_sorbed*scale_sorbed);
    c_aqua_limit_ = c_aqua_limit;
    /*cout << "sorp_type " << sorption_type << endl;
    cout << "scale_aqua " << scale_aqua << endl;
    cout << "scale_sorbed " << scale_sorbed << endl;
    cout << "inv_scale_aqua " << inv_scale_aqua << endl;
    cout << "inv_scale_sorbed " << inv_scale_sorbed << endl;
    cout << "c_aqua_limit " << c_aqua_limit_ << endl;*/
}

//inline
bool Isotherm::compute_projection(double &c_aqua, double &c_sorbed) //clear as glass but the inline command makes troubles, probably
{
    double total_mass = scale_aqua_* c_aqua + scale_sorbed_ * c_sorbed;
    double elem_volume = 1.0; // must be corrected later, this value is not correct, but the function is not used in this time
    //unsigned
    //if(total_mass < 0.0) total_mass = 0.0; // this must be solved somehow else, negative total mass is strange
    int i_total_mass = total_mass / total_mass_step;
    /*double min_max_total_mass = 0.0;
    if(total_mass < 0.0) //|| (i_total_mass < 0))
    {
    	if(total_mass < min_max_total_mass )
    	//cout << "i_total_mass is " << i_total_mass << " and total mass is " << total_mass << " and total_mass_step has the value " << total_mass_step << endl;
    }*/
    if (i_total_mass < 0) return false;
    //cout << "interpolation_table size is " << interpolation_table.size() << endl;
    if ( (unsigned int)(i_total_mass) < interpolation_table.size() ) {
    	int iso_ind_floor, iso_ind_ceil;
    	iso_ind_floor = (int)(total_mass/(total_mass_step)); iso_ind_ceil = iso_ind_floor + 1;
    	double rot_sorbed = interpolation_table[iso_ind_floor] + (total_mass - iso_ind_floor*total_mass_step)*(interpolation_table[iso_ind_ceil] - interpolation_table[iso_ind_floor])/total_mass_step;
        c_aqua = (total_mass * inv_scale_aqua_ - rot_sorbed*inv_scale_sorbed_);
        c_sorbed = (total_mass * inv_scale_sorbed_ + rot_sorbed * inv_scale_aqua_);
        return true;
    } else {
    	//cout << "c_aqua_limit_ has the value " << c_aqua_limit_ << endl;
        if (c_aqua_limit_ > 0.0) {
    		precipitate(c_aqua, c_sorbed, scale_aqua_, scale_sorbed_); // , elem_volume); // last two parameters are probably not correct, they are given for the function to make it compilable
            /*c_sorbed = (total_mass - scale_aqua_* c_aqua_limit_) * inv_scale_sorbed_;
            c_aqua = c_aqua_limit_;*/
        } else
        {
        	//cout << "c_aqua_limit_ has the value " << c_aqua_limit_ << endl;
        	return false;
        }
    }

    return true; //false;
}

template<class Func>
void Isotherm::solve_conc(double &c_aqua, double &c_sorbed, const Func &isotherm) // , double elem_volume) // Probably not used at this time. CrossFunction needs to be redefined.
{
    //double mass_limit;
    boost::uintmax_t max_iter=100;
    boost::math::tools::eps_tolerance<double> toler(60);
    // the condition written bellow seems to be obsolete or should be placed
    /*double f_max = const_cast<Func &>(isotherm)(c_aqua_limit_); //iso_hlp(c_aqua_limit_);
    if (c_aqua_limit_ >0) {
        mass_limit = scale_aqua*c_aqua_limit_ + scale_sorbed*f_max; // isotherm(c_aqua_limit_);
    } else {
        cout << "Solubility limit has to positive " << c_aqua_limit_ << endl;
    }*/
    // Following 5 lines are modificated temporarily by multiplication by elem_volume.
	double scale_aqua = scale_aqua_; // * elem_volume;
	double scale_sorbed = scale_sorbed_; // * elem_volume;
	double total_mass = (scale_aqua*c_aqua + scale_sorbed * c_sorbed);
	double inv_scale_aqua = inv_scale_aqua_; // /(elem_volume * elem_volume);
	double inv_scale_sorbed = inv_scale_sorbed_; // /(elem_volume * elem_volume);
	double critic_total_mass = c_aqua_limit_*scale_aqua + const_cast<Func &>(isotherm)(c_aqua_limit_)*scale_sorbed;

	const double upper_solution_bound = critic_total_mass * inv_scale_aqua + 1.0;

	if(total_mass < critic_total_mass)
	{
		CrossFunction<Func> eq_func(isotherm, total_mass, scale_aqua, scale_sorbed); // equation describing one point on the isotherm
		pair<double,double> solution = boost::math::tools::toms748_solve(eq_func, 0.0, upper_solution_bound, toler, max_iter); // (eq_func, 0.0, 10.0, toler, max_iter);
		//SOLUTION IS AN INTERVAL CONTAINING SOLUTION, because of that following two lines are commented
		//toms748_solve returns interval bounds
		c_aqua = (solution.first + solution.second)/2; // = average of the pair solution defined above, midpoint
		//cout << "aqueous concentration is " << c_aqua << endl;
		c_sorbed = const_cast<Func &>(isotherm)(c_aqua); // = f(midpoint)
	}else{
		precipitate(c_aqua, c_sorbed, scale_aqua, scale_sorbed); // , elem_volume);
	}

    return;
}

template void Isotherm::solve_conc<Linear>(double &c_aqua, double &c_sorbed, const Linear &isotherm); // , double elem_volume);

template void Isotherm::solve_conc<Langmuir>(double &c_aqua, double &c_sorbed, const Langmuir &isotherm); // , double elem_volume);

template void Isotherm::solve_conc<Freundlich>(double &c_aqua, double &c_sorbed, const Freundlich &isotherm); //, double elem_volume);

void Isotherm::precipitate(double &c_aqua, double &c_sorbed, double scale_aqua, double scale_sorbed) //, double elem_volume)
{
	double total_mass = (scale_aqua*c_aqua + scale_sorbed * c_sorbed);
	double inv_scale_sorbed = 1/((scale_aqua_*scale_aqua + scale_sorbed*scale_sorbed_)); // *(elem_volume * elem_volume));

	c_aqua = c_aqua_limit_;
	c_sorbed = (total_mass - scale_aqua * c_aqua_limit_)*inv_scale_sorbed;

	return;
}

template<class Func>
void Isotherm::make_table(const Func &isotherm, int n_steps) { //const Func &isotherm, int n_steps
    double mass_limit; // c_aqua, c_sorbed;
    //Func &iso_hlp = const_cast<Func &>(isotherm);
    //double f_max = isotherm(c_aqua_limit_);
    //interpolation_table.resize(n_steps); // obsolete
    double f_max = const_cast<Func &>(isotherm)(c_aqua_limit_);
    //if(f_max < 0.0) cout << "Functor for isotherm returns negative value " << f_max << endl; else cout << "Functor for isotherm returns positive value " << f_max << endl; // correct
    SorptionType sorpt_type = this->get_sorption_type();
    //cout << "Sorption type is " << sorpt_type << endl; // correct
    if (c_aqua_limit_ > 0.0) {
        mass_limit = scale_aqua_*c_aqua_limit_ + scale_sorbed_*f_max; //isotherm(c_aqua_limit_);
        if(mass_limit < 0.0)
        {
        	cout << "isotherm type " << sorpt_type << ", mass_limit has negative value " << mass_limit << ", scale_aqua "  << scale_aqua_ << ", c_aq_limit " << c_aqua_limit_ << ", scale_sorbed " << scale_sorbed_ << ", f_max " << f_max << endl;
        }
    } else {
        cout << "Solubility limit has to be higher than 0.0" << endl;
        return;
    	//mass_limit = scale_aqua + scale_sorbed;// set mass_limit from max conc = 1, needs to be computed somehow else
    }
    total_mass_step = mass_limit / n_steps;
    double mass = 0.0; // total_mass_step;
    for(int i=0; i<= n_steps;i++, mass+=total_mass_step) {
        double c_aqua = mass * inv_scale_aqua_; // aqueous concentration (original coordinates c_a) corresponding to total mass
        double c_sorbed = const_cast<Func &>(isotherm)(c_aqua); // functional value appropriate to f(c_a)
        //solve_conc(c_aqua, c_sorbed, isotherm); // this line seems to be obsolete
    	double c_sorbed_rot = (c_sorbed * scale_aqua_ - c_aqua * scale_sorbed_); // const_cast<Func &>(isotherm)(mass);
        interpolation_table.push_back(c_sorbed_rot);
    }

    return;
}

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

template void Isotherm::make_table<Linear>(const Linear &isotherm, int n_steps);

template void Isotherm::make_table<Langmuir>(const Langmuir &isotherm, int n_steps);

template void Isotherm::make_table<Freundlich>(const Freundlich &isotherm, int n_steps);


double Isotherm::get_scale_aqua(void)
{
	return scale_aqua_;
}

double Isotherm::get_scale_sorbed(void)
{
	return scale_sorbed_;
}

int Isotherm::get_interpolation_table_size(void)
{
	return interpolation_table.size();
}

/*template <class Func>
CrossFunction::CrossFunction(const Func &func_,  double total_mass, double scale_aqua, double scale_sorbed)
: func(func_), total_mass_(total_mass), scale_aqua(scale_aqua), scale_sorbed(scale_sorbed)
{
	return;
};*/

