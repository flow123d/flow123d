/*
 * isotherms.hh
 *
 *  Created on: Jan 15, 2013
 *      Author: jb
 */

#ifndef ISOTHERMS_HH_
#define ISOTHERMS_HH_

#include <boost/math/tools/roots.hpp>
//#include <boost/math/tools/toms748_solve.hpp>

// Functors for isotherms.


class Langmuir {
public:
    Langmuir( double mult_coef_, double alpha_)
    : alpha(alpha_) , mult_coef(mult_coef_) {}

    double operator()( double x) {
        return mult_coef*(alpha * x)/(1+alpha*x);
    }
private:
    double mult_coef;
    double alpha;
};

class Linear {
public:
	Linear(double direct_)
	: direct(direct_) {}

	double operator()(double x){
		return direct*x;
	}
private:
	double direct;
};

// Equation for intersection points for any isotherm
template <class Func>
class CrossFunction {
public:
    CrossFunction(const Func &func_,  double total_mass_, double slope_)
    : func(func_), total_mass(total_mass_), slope(slope_) {}


    double operator()( double conc_aqua) {
        return total_mass  - conc_aqua * slope - func(conc_aqua); // not shure about slope sign
    }
private:
    Func func;
    double total_mass, slope;
};

// Equation determines total_mass parameter for the line going through the point [c_k,f(c_k)], it should probably rather produce array or vector output instead of scalar value
template <class Func>
class TotalMass {
public:
	TotalMass(const Func &func_, double c_aq_max_, double slope_, int n_points_)
	: func(func_), c_aq_max(c_aq_max_), n_points(n_points_), slope(slope_) {}

	std::vector<double> total_mass; //size needs to be specified
	//double *operator ()(void){
	std::vector<double> *operator ()(void){
		double MaxTotalMass;
		MaxTotalMass = slope*c_aq_max + func(c_aq_max);

		for(int point = 0; point < n_points; point++){
			total_mass[point] = point * MaxTotalMass/n_points;
		}
		return total_mass;
	}

private:
	Func func;
	double c_aq_max, slope;
	int n_points;
};

template <class Func>
class CrossSolve { // where to place this function, probably inside of determine_crossection instead in sorption.cc
public:
	CrossSolve(const Func &func_, double porosity_, double rho_rock_, double rho_water_, double adsorbent_molar_mass_, int n_points_, double slope_)
	: func(func_), porosity(porosity_), rho_rock(rho_rock_), rho_water(rho_water_), adsorbent_molar_mass(adsorbent_molar_mass_), n_points(n_points_), slope(slope_) {}

	std::vector<double> conc_aqua; //size need to be specified
	double *operator()(double c_aq_max){
		 //Langmuir isotherm( mult_coef, alpha );
		 slope = (rho_water * porosity)/(adsorbent_molar_mass * rho_rock *(1 - porosity));
		 // total mass is not a mass from physical point of wiev but it is mass balance line crossection with y-axis
		 TotalMass<Func> total_mass(func, c_aq_max, slope, n_points);
		 for(int i = 0; i < n_points; i++ ){
			 CrossFunction<Func> funct(func, total_mass[i], slope);
			 conc_aqua = boost::math::tools::toms748_solve( funct, 0, c_aq_max, boost::math::tools::eps_tolerance<double>(60), 100) . first;
		 }
		 return conc_aqua;
	}

private:
	Func func;
	//double *conc_aqua;
	double porosity, rho_rock, rho_water, adsorbent_molar_mass, slope;
	int n_points;
};


#endif /* ISOTHERMS_HH_ */
