/*
 * isotherm.cc
 *
 *  Created on: Mar 15, 2013
 *      Author: lukas
 */

#include "reaction/isotherm.hh"

void Linear::reinit(double mult_coef)
{
	mult_coef_ = mult_coef;
	return;
};

void Isotherm::reinit(enum SorptionType sorp_type, double rock_density, double rho_aqua, double porosity, double molar_mass, double c_aqua_limit)
{
    // set class variables
	sorption_type = sorp_type;
    scale_aqua = porosity * rho_aqua;
    scale_sorbed = (1-porosity) * rock_density * molar_mass;
    inv_scale_aqua = scale_aqua/(scale_aqua*scale_aqua + scale_sorbed*scale_sorbed);
    inv_scale_sorbed = scale_sorbed/(scale_aqua*scale_aqua + scale_sorbed*scale_sorbed);
    c_aqua_limit_=c_aqua_limit;
};

//inline
bool Isotherm::compute_projection(double &c_aqua, double &c_sorbed) //clear as glass but the inline command makes troubles, probably
{
    double total_mass = scale_aqua* c_aqua + scale_sorbed * c_sorbed;
    unsigned int i_total_mass = total_mass / total_mass_step;
    if (i_total_mass < 0) return false;
    if (i_total_mass < interpolation_table.size()) {
    	int iso_ind_floor, iso_ind_ceil;
    	iso_ind_floor = (int)(total_mass/(total_mass_step)); iso_ind_ceil = iso_ind_floor + 1;
    	double rot_sorbed = interpolation_table[iso_ind_floor] + (total_mass - iso_ind_floor*total_mass_step)*(interpolation_table[iso_ind_ceil] - interpolation_table[iso_ind_floor])/total_mass_step;
        c_aqua = (total_mass * inv_scale_aqua - rot_sorbed*inv_scale_sorbed);
        c_sorbed = (total_mass * inv_scale_sorbed + rot_sorbed * inv_scale_aqua);
        return true;
    } else {
        if (c_aqua_limit_ > 0.0) {
            c_sorbed = (total_mass - scale_aqua* c_aqua_limit_)*inv_scale_sorbed;
            c_aqua = c_aqua_limit_;
        } else return false;
    }

    return false;
};

template<class Func>
void Isotherm::solve_conc(double &c_aqua, double &c_sorbed, const Func &isotherm) // Probably not used at this time. CrossFunction needs to be redefined.
{
    double mass_limit;
    boost::uintmax_t max_iter=100;
    boost::math::tools::eps_tolerance<double> toler(60);
    //Func &iso_hlp = const_cast<Func &>(isotherm);
    double f_max = const_cast<Func &>(isotherm)(c_aqua_limit_); //iso_hlp(c_aqua_limit_);
    if (c_aqua_limit_ >0) {
        mass_limit = scale_aqua*c_aqua_limit_ + scale_sorbed*f_max; // isotherm(c_aqua_limit_);
    } else {
        mass_limit = scale_aqua + scale_sorbed;// set mass_limit from max conc = 1, needs to be computed somehow else
    }
	double total_mass = scale_aqua*c_aqua + scale_sorbed * c_sorbed;
    CrossFunction<Func> eq_func(isotherm, total_mass, scale_aqua, scale_sorbed); // equation describing one point on the isotherm
    pair<double,double> solution = boost::math::tools::toms748_solve(eq_func, 0.0, 10.0, toler, max_iter);
    //PROBABLY COMPLETELY WRONG, SOLUTION IS AN INTERVAL CONTAINING SOLUTION, because of that following two lines are commented
    //toms748_solve returns interval bounds
    //c_aqua = (total_mass - scale_sorbed * solution.first) / scale_aqua;
    //c_sorbed = (total_mass - scale_aqua * solution.second) / scale_sorbed;
    //MUST BE REPARED, LATER.
    c_aqua = 1.0;
    c_sorbed = 1.0;

    return;
};

template void Isotherm::solve_conc<Linear>(double &c_aqua, double &c_sorbed, const Linear &isotherm);

template void Isotherm::solve_conc<Langmuir>(double &c_aqua, double &c_sorbed, const Langmuir &isotherm);

//template void Isotherm::solve_conc<Freundlich>(double &c_aqua, double &c_sorbed, const Freundlich &isotherm);

template<class Func>
void Isotherm::make_table(const Func &isotherm, int n_steps) { //const Func &isotherm, int n_steps
    double mass_limit; // c_aqua, c_sorbed;
    //Func &iso_hlp = const_cast<Func &>(isotherm);
    //double f_max = isotherm(c_aqua_limit_);
    interpolation_table.resize(n_steps);
    double f_max = const_cast<Func &>(isotherm)(c_aqua_limit_);
    if (c_aqua_limit_ >0) {
        mass_limit = scale_aqua*c_aqua_limit_ + scale_sorbed*f_max; //isotherm(c_aqua_limit_);
    } else {
        mass_limit = scale_aqua + scale_sorbed;// set mass_limit from max conc = 1, needs to be computed somehow else
    }
    total_mass_step = mass_limit / n_steps;
    double mass = total_mass_step; // we need not to save value for zero mass, it is zero
    for(int i=0; i< n_steps;i++, mass+=total_mass_step) {
        /*double c_aqua = mass * inv_scale_aqua; // aqueous concentration (original coordinates c_a) corresponding to total mass
        double c_sorbed = const_cast<Func &>(isotherm)(c_aqua); // mass * inv_scale_sorbed;
        solve_conc(c_aqua, c_sorbed, isotherm);
        interpolation_table.push_back( c_sorbed * scale_sorbed - c_aqua * scale_aqua);*/
    	double c_sorbed_rot = const_cast<Func &>(isotherm)(mass);
        interpolation_table.push_back(c_sorbed_rot);
    }

    return;
};

template void Isotherm::make_table<Linear>(const Linear &isotherm, int n_steps);

template void Isotherm::make_table<Langmuir>(const Langmuir &isotherm, int n_steps);

//template void Isotherm::make_table<Freundlich>(const Freundlich &isotherm, int n_steps);


double Isotherm::get_scale_aqua(void)
{
	return scale_aqua;
};

double Isotherm::get_scale_sorbed(void)
{
	return scale_sorbed;
};


/*template <class Func>
CrossFunction::CrossFunction(const Func &func_,  double total_mass, double scale_aqua, double scale_sorbed)
: func(func_), total_mass_(total_mass), scale_aqua(scale_aqua), scale_sorbed(scale_sorbed)
{
	return;
};*/

