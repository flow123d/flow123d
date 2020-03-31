/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    isotherm.hh
 * @brief   
 *
 * Other possible transformation of coordinates:
 *
 * c_l - conc. liquid
 * c_s - conc. solid
 * c_t = h_l * c_l + h_s * c_s = X' + Y'
 * X = c_t
 * Y = X' - Y' = h_l * c_l - h_s * c_s
 *
 * A) make table for function c_t -> Y
 * 1) for given c_t solve nonlinear eq.
 *    h_l * c_l + h_s * f(c_l) = c_t
 *
 * 2) from c_l compute
 *    Y = c_t - 2* h_l * c_l
 *
 *
 * B) calculation of new c_l, c_s from c_t using table:
 *
 * 1) use table to get value of Y for given c_t
 * 2) compute:
 *    c_l = (c_t - Y) / (2 * h_l)
 *    c_s = (c_t + Y) / (2 * h_s)
 *
 * ========================
 * The transformation currently in use transforms
 * pair (c_l, c_s) directly to (c_t, W) in ortogonal way.
 * Proposed transformation first scale (c_l, c_s) to (X',Y')
 * and then transform scaled coordinates in ortogonal way.
 */

#ifndef ISOTHERM_IMPL_HH_
#define ISOTHERM_IMPL_HH_

#include <stdint.h>                                           // for uintmax_t
#include <algorithm>                                          // for max
#include <vector>
#include <cmath>                                              // for floor, pow
#include <complex>                                            // for fabs
#include <ostream>                                            // for operator<<
#include <string>                                             // for string
#include "input/input_exception.hh"                           // for EI_Address
#include "system/exceptions.hh"                               // for operator<<
#include <boost/math/special_functions/detail/round_fwd.hpp>  // for BOOST_M...

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
 * Complementary functor for isotherm of type "none"
 */
class None {
public:
    /// Constructor to set parameters
    None(){}
    /// Isotherm definition.
    inline double operator()(double) {
        return (0.0);
    }
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
    TYPEDEF_ERR_INFO( EI_BoostMessage, std::string);
    DECLARE_EXCEPTION( ExcBoostSolver, << 
    "The Boost solver of nonlinear equation did not converge in sorption model.\n" 
    << "Check input at the following address for possible user error: \t" << Input::EI_Address::val << "\n"
    "Solver message: \n" << EI_BoostMessage::val);
    
    TYPEDEF_ERR_INFO( EI_TotalMass, double);
    DECLARE_EXCEPTION( ExcNegativeTotalMass, << 
    "The total mass in sorption model became negative during the computation (value: " 
    << EI_TotalMass::val << ").\n"
    << "Check input at the following address for possible user error: \t" << Input::EI_Address::val << "\n");
    
    
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

    /// Default constructor.
    Isotherm();

    /**
     * Setting adsorption parameters for general isotherm. These parameters are then used either
     * for creation of the interpolation table via @p make_table method or just one adsorption is computed
     * through @p compute method.  Provided parameters are:
     * @param sorption_type - type of isotherm
     * @param limited_solubility_on - true if @p c_aqua_limit is solubility limit
     * @param aqua_density - density of the liquid phase
     * @param scale_aqua - generalized porosity, fraction of the space with liquid phase
     * @param scale_sorbed  - fraction of the space with the solid to which we adsorp
     * @param c_aqua_limit - solubility limit (limit aqueous concentration)
     * @param mult_coef - multiplicative coefficient of the isotherm (all isotherms have one)
     * @param second_coef - possibly second parameter of the isotherm
     */
    void reinit(enum SorptionType sorption_type, bool limited_solubility_on,
                double aqua_density, double scale_aqua, double scale_sorbed,
                double c_aqua_limit, double mult_coef, double second_coef);

    /**
     * Create interpolation table for isotherm in rotated coordinate system with X axes given by total mass in
     * both phases. Currently we support only linear interpolation.
     * @p reinit has to be called just before this method.
     * @param n_points is the size of the table
     * @param table_limit is the limit value of aqueous concentration
     */
    void make_table(unsigned int n_points, double table_limit);

    /**
     * Clears the interpolation table and resets the table limit.
     * (That means, no interpolation takes place in the computation.)
     */
    void clear_table();

    /**
    * Direct calculation of the equilibrium adsorption using a non-linear solver.
    * @p reinit has to be called just before this method.
    */
    void compute(double &c_aqua, double &c_sorbed);

    /**
     * Use interpolation to determine equilibrium state.
     * Assumes previous call to @p make_table. If total mass is larger then table limit we either
     * call @p precipitate (limit_solubility_on) or use direct computation.
     */
    void interpolate(double &c_aqua, double &c_sorbed);

    /**
     * Returns true if interpolation table is created.
     */
    inline bool is_precomputed(void) {
        return interpolation_table.size() != 0;
    }

    /// Getter for table limit (limit aqueous concentration).
    inline double table_limit(void) const
    { return table_limit_;}
    
protected:
	/**
     * Implementation of interpolation construction for particular isotherm functor.
     * Specialized for sorption 'none' - creates empty table.
     */
    template<class Func>
    void make_table(const Func &isotherm, int n_points);
    /**
     * Find new values for concentrations in @p c_pair that has same total mass and lies on the
     * @p isotherm (functor object).
     * Specialized for sorption 'none' - returns unchanged @p c_pair.
     */
    template<class Func>
    ConcPair solve_conc(ConcPair c_pair, const Func &isotherm);
    /**
     * Dispatch isotherm type and use appropriate template.
     */
    ConcPair solve_conc(ConcPair conc);
    /**
     * Update concentrations using interopolation.
     */
    ConcPair compute_projection( ConcPair conc );
    /**
     * Modify concentrations after adsorption for limited solubility.
     */
    ConcPair precipitate( ConcPair conc );

    double get_total_mass( ConcPair conc );

    /****************************************
     * Data
     */

    /// Type of isotherm
    enum SorptionType adsorption_type_;

    /// Multiplication parameter of the isotherm
    double mult_coef_;

    /// Optional secod parameter of the isotherm
    double second_coef_;

    /// Concentration in liquid phase for limit of the interpolation table.
    double table_limit_;

    /// Solubility limit flag
    bool limited_solubility_on_;
    /// Concentration limit in liquid phase (solubility limit).
    double solubility_limit_;

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
    std::vector<double> interpolation_table;
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

#endif /* ISOTHERM_IMPL_HH_ */
