/*
 * isotherm.hh
 *
 *  Created on: Mar 7, 2013
 *      Author: jb
 *
 *

 */

#ifndef SORPTION_IMPL_HH_
#define SORPTION_IMPL_HH_

#include <vector>
#include <input/input_type.hh>
#include <boost/math/tools/roots.hpp>
#include "fields/field.hh"



/**
 * Convergence criteria for interval based nonlinear solver. It is functor, that
 * returns true if bounds a,b of the solution are close enough.
 * We use relative criteria.
 */
template <class T>
class tolerance
{
public:
   tolerance(unsigned bits)
   {
      BOOST_MATH_STD_USING
      eps = T(ldexp(1.0F, 1-bits));
   }
   bool operator()(const T& a, const T& b)
   {
      BOOST_MATH_STD_USING
      return fabs(a - b) <= (eps * (std::max)(fabs(a), fabs(b)));
   }
private:
   T eps;
};



/**
 * Functor for linear isotherm
 */
class Linear {
public:
	/// Constructor to set parameters
    Linear(double mult_coef) : mult_coef_(mult_coef) {}
    /// Isotherm definition.
    inline double operator()(double x) {
    	return (mult_coef_*x);
    }
private:
    /// Parameters of the isotherm.
    double mult_coef_;
};



/**
 * Functor for Langmuir isotherm.
 */
class Langmuir {
public:
	/// Constructor to set parameters
	Langmuir( double mult_coef, double alpha) : mult_coef_(mult_coef), alpha_(alpha) {}
    /// Isotherm definition.
    inline double operator()( double x) {
    	return (mult_coef_*(alpha_ * x)/(alpha_ *x + 1));
    }

private:
    /// Parameters of the isotherm.
    double mult_coef_;
    double alpha_;
};



/**
 * Functor for Freundlich isotherm.
 */
class Freundlich {
public:
	/// Constructor to set parameters
	Freundlich(double mult_coef, double exponent) : mult_coef_(mult_coef), exponent_(exponent){}
    /// Isotherm definition.
	inline double operator()(double x) {
		return (mult_coef_*pow(x, exponent_));
	}

private:
    /// Parameters of the isotherm.
	double mult_coef_;
	double exponent_;
};


/**
* Class describing one isotherm with possibly precalculated interpolation table.
*/
class Isotherm {
public:

	/// Type of adsorption isotherm.
	enum SorptionType {
		none = 0,
		linear = 1,
		freundlich = 2,
		langmuir = 3
	};

	/// Pair of soluted and adsorbed concentration.
	struct ConcPair {
		ConcPair(double x, double y) : fluid(x), solid(y) {}
		double fluid;
		double solid;
	};

    /**
     * Setting adsorption parameters for general isotherm. These parameters are then used either
     * for creation of the interpolation table via @p make_table method or just one adsorption is computed
     * through @p compute method.  Provided parameters are:
     * @param sorption_type - type of isotherm
     * @param limited_solubility_on - true if @p c_aqua_limit is solubility limit
     * @param aqua_density - density of the liquid phase
     * @param scale_aqua - generalized porosity, fraction of the space with liquid phase
     * @param scale_sorbed  - fraction of the space with the solid to which we adsorp
     * @param c_aqua_limit - limit for interpolation table, possibly solubility limit
     * @param mult_coef - multiplicative coefficient of the isotherm (all isotherms have one)
     * @param secodn_coef - possibly second parameter of the isotherm
     */
	inline void reinit(enum SorptionType sorption_type, bool limited_solubility_on,
			double aqua_density, double scale_aqua, double scale_sorbed,
			double c_aqua_limit, double mult_coef, double second_coef);
    /**
     * Create interpolation table for isotherm in rotated coordinate system with X axes given by total mass in
     * both phases. Size of the table is the only parameter. Currently we support only linear interpolation.
     * @p reinit has to be called just before this method.
     */
    void make_table(int n_points);

    /**
    * Direct calculation of the equilibrium adsorption using a non-linear solver.
    * @p reinit has to be called just before this method.
    */
    inline void compute(double &c_aqua, double &c_sorbed);

    /**
     * Use interpolation to determine equilibrium state.
     * Assumes previous call to @p make_table. If total mass is larger then table limit we either
     * call @p precipitate (limit_solubility_on) or use direct computation.
     */
     inline void interpolate(double &c_aqua, double &c_sorbed);

    /**
     * Returns true if interpolation table is created.
     */
    inline bool is_precomputed(void) {
        return interpolation_table.size() != 0;
    }  
    
protected:
	/**
     * Implementation of interpolation construction for particular isotherm functor.
     */
    template<class Func>
    void make_table(const Func &isotherm, int n_points);
    /**
     * Find new values for concentrations in @p c_pair that has same total mass and lies on the
     * @p isotherm (functor object).
     */
    template<class Func>
    inline ConcPair solve_conc(ConcPair c_pair, const Func &isotherm);
    /**
     * Dispatch isotherm type and use appropriate template.
     */
    inline ConcPair solve_conc(ConcPair conc);
    /**
     * Update concentrations using interopolation.
     */
    inline ConcPair compute_projection( ConcPair conc );
    /**
     * Modify concentrations after adsorption for limited solubility.
     */
    inline ConcPair precipitate( ConcPair conc );

    /****************************************
     * Data
     */

    /// Type of isotherm
    enum SorptionType adsorption_type_;

    /// Multiplication parameter of the isotherm
    double mult_coef_;

    /// Optional secod parameter of the isotherm
    double second_coef_;

    /*  Concentration in liquid phase for limit of the interpolation table, or
     *  solubility limit.
     */
    double table_limit_;

    /// Solubility limit flag
    bool limited_solubility_on_;

    /// density of the solvent
    double rho_aqua_;
    /// coefficient that convert soluted concentration to mass; porosity = k_W, originally rho_aqua*porosity = k_W
    double scale_aqua_;
    /// coefficient that convert adsorbed molar concentration to mass; molar_weight * rho_rock * (1 - porosity) = k_H
    double scale_sorbed_;
    /// reciprocal values
    double inv_scale_aqua_, inv_scale_sorbed_;
    /**
     * Interpolation table of isotherm in the rotated coordinates.
     * The X axes of rotated system is total mass, the Y axes is perpendicular.
     */
    vector<double> interpolation_table;
    /**
     * Step on the rotated X axes (total mass).
     */
    double total_mass_step_;

};


/**
 *  Functor for solved equation in form F(x) ==0.
 *  Function @p func is an isotherm functor object in concentration based coordinated system.
 *  We solve the equation in modified system (scaled and rotated) for better numerical stability.
 *  The solved equation reads:
 *  F(X) -Y =0, where
 *  X is total mass , Y
 */
template <class Func>
class CrossFunction
{
public:
    CrossFunction(const Func &func_,  double total_mass, double scale_aqua, double scale_sorbed, double rho_aqua)
    : func(func_), total_mass_(total_mass),
      scale_sorbed_(scale_sorbed), scale_aqua_(scale_aqua), rho_aqua_(rho_aqua)
    {}

    double operator()( double conc_aqua)
    {
    	// that is the  selected isotherm // scale_sorbed_ * func( conc_aqua ) + scale_aqua_ * conc_aqua - total_mass_
        return scale_sorbed_*func( conc_aqua/rho_aqua_) + (scale_aqua_) * conc_aqua - total_mass_;
    }
private:
    Func func;
    double total_mass_, scale_sorbed_, scale_aqua_, rho_aqua_;
};


/*****************************************************************************************************************
 * IMPLEMENTATION
 */


inline void Isotherm::reinit(enum SorptionType adsorption_type, bool limited_solubility_on,
		              double rho_aqua, double scale_aqua, double scale_sorbed,
		              double c_aqua_limit, double mult_coef, double second_coef)
{
	adsorption_type_ = adsorption_type;
	rho_aqua_ = rho_aqua;
	scale_aqua_ = scale_aqua;
	scale_sorbed_ = scale_sorbed;
    inv_scale_aqua_ = scale_aqua_/(scale_aqua_*scale_aqua_ + scale_sorbed_*scale_sorbed_);
    inv_scale_sorbed_ = scale_sorbed_/(scale_aqua_*scale_aqua_ + scale_sorbed_*scale_sorbed_);
    table_limit_ = c_aqua_limit;
    limited_solubility_on_ = limited_solubility_on;
    mult_coef_ = mult_coef;
    second_coef_ = second_coef;

}



inline void Isotherm::compute( double &c_aqua, double &c_sorbed ) {
	ConcPair c_pair(c_aqua, c_sorbed);
	ConcPair result(0,0);

    if (limited_solubility_on_ && (c_pair.fluid > table_limit_)) {
        result = precipitate( c_pair );
    } else {
       	result = solve_conc( c_pair );
    }
    c_aqua=result.fluid;
    c_sorbed=result.solid;
}


inline void Isotherm::interpolate( double &c_aqua, double &c_sorbed ) {
	ConcPair c_pair(c_aqua, c_sorbed);
	ConcPair result(0,0);

	result = compute_projection( c_pair );

	c_aqua=result.fluid;
    c_sorbed=result.solid;
}



inline Isotherm::ConcPair Isotherm::compute_projection( Isotherm::ConcPair c_pair ) {
  double total_mass = (scale_aqua_* c_pair.fluid + scale_sorbed_ * c_pair.solid);
  double total_mass_steps = total_mass / total_mass_step_;
  int total_mass_idx = static_cast <int>(std::floor(total_mass_steps));

  if ( total_mass_idx < 0 ) {xprintf(UsrErr,"total_mass %f seems to have negative value.\n", total_mass); }
  if ((unsigned int)(total_mass_idx) < (interpolation_table.size() - 1) ) {
      double rot_sorbed = interpolation_table[total_mass_idx] + (total_mass_steps - total_mass_idx)*(interpolation_table[total_mass_idx+1] - interpolation_table[total_mass_idx]);
      return ConcPair( (total_mass * inv_scale_aqua_ - rot_sorbed * inv_scale_sorbed_),
                       (total_mass * inv_scale_sorbed_ + rot_sorbed * inv_scale_aqua_) );
  } else {
      if (limited_solubility_on_) {
              return precipitate( c_pair );
      } else {
              return solve_conc( c_pair );
      }
  }
}



inline Isotherm::ConcPair Isotherm::precipitate( Isotherm::ConcPair c_pair) {
	double total_mass = (scale_aqua_*c_pair.fluid + scale_sorbed_ * c_pair.solid);
	return ConcPair(	table_limit_,
						(total_mass - scale_aqua_ * table_limit_)/scale_sorbed_  );
}



template<class Func>
inline Isotherm::ConcPair Isotherm::solve_conc(Isotherm::ConcPair c_pair, const Func &isotherm)
{
	//START_TIMER("new-sorption, solve_conc, toms748_solve");
	boost::uintmax_t max_iter = 20;
	tolerance<double> toler(30);
	double total_mass = (scale_aqua_*c_pair.fluid + scale_sorbed_ * c_pair.solid);
	double mass_limit = table_limit_*scale_aqua_ + const_cast<Func &>(isotherm)(table_limit_ / this->rho_aqua_)*scale_sorbed_;
	double upper_solution_bound;

	if(total_mass >= mass_limit)
	{
		mass_limit = total_mass;
	}
	upper_solution_bound = mass_limit / scale_aqua_;
	//xprintf(Msg,"%s upper_solution_bound %f, total_mass %e \n",typeid(Func).name(), upper_solution_bound, total_mass);
	CrossFunction<Func> eq_func(isotherm, total_mass, scale_aqua_, scale_sorbed_, this->rho_aqua_);
	//xprintf(Msg,"CrossFunction returns %e, %f, scale_aqua_ %f, scale_sorbed_ %f, c_aqua %e, c_sorbed %e\n", eq_func(0), eq_func(upper_solution_bound), scale_aqua_, scale_sorbed_, c_aqua, c_sorbed);
	pair<double,double> solution;
	if (total_mass > 0) // here should be probably some kind of tolerance instead of "0"
		solution = boost::math::tools::toms748_solve(eq_func, 0.0, upper_solution_bound, toler, max_iter);
	double difference;
	difference = (solution.second - solution.first)/2;
	double c_aqua = solution.first + difference;
	return ConcPair( c_aqua,
					 (total_mass - scale_aqua_ * c_aqua)/scale_sorbed_);
}


Isotherm::ConcPair Isotherm::solve_conc(Isotherm::ConcPair conc)
{
        switch(adsorption_type_)
        {
                case 0: // none
                {
                        Linear obj_isotherm(0.0);
                        return solve_conc( conc, obj_isotherm);
                }
                break;
                case 1: //  linear:
                {
                        Linear obj_isotherm(mult_coef_);
                        return solve_conc( conc, obj_isotherm);
                }
                break;
                case 2: // freundlich
                {
                        Freundlich obj_isotherm(mult_coef_, second_coef_);
                        return solve_conc( conc, obj_isotherm);
                }
                break;
                case 3:  // langmuir:
                {
                        Langmuir obj_isotherm(mult_coef_, second_coef_);
                        return solve_conc( conc, obj_isotherm);
                }
                break;
        }
        return conc;
}


template<class Func>
void Isotherm::make_table(const Func &isotherm, int n_steps)
{
    double mass_limit = scale_aqua_ * table_limit_ + scale_sorbed_ * const_cast<Func &>(isotherm)(table_limit_ / this->rho_aqua_);
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


#endif /* SORPTION_IMPL_HH_ */
