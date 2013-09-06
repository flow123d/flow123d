/*
 * isotherm.hh
 *
 *  Created on: Mar 7, 2013
 *      Author: jb
 */

#ifndef SORPTION_IMPL_HH_
#define SORPTION_IMPL_HH_

#include <vector>
#include <input/input_type.hh>
#include <boost/math/tools/roots.hpp>

enum SorptionType {
	none = 0,
	linear = 1,
	freundlich = 2,
	langmuir = 3
};

template <class T>
class dekl_tolerance
{
public:
   dekl_tolerance(unsigned bits)
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
 * Functor for Langmuir isotherm.
 */
class Langmuir {
public:
	/**
	* 	Original constructor, Langmuir( double mult_coef, double alpha) : alpha(alpha), mult_coef_(mult_coef) {}
	*/
    Langmuir( double mult_coef, double alpha) : mult_coef_(mult_coef), alpha(alpha) {}
    /**
    * 	Destructor.
    */
    ~Langmuir(void){}
    /**
    * 	Operator.
    */
    double operator()( double x) { return (mult_coef_*(alpha * x)/(alpha*x + 1)); }

private:
    double mult_coef_;
    double alpha;
};


/**
 * Functor for linear isotherm
 */
class Linear {
public:
	/**
	* 	Original constructor, Linear(double mult_coef) : mult_coef_(mult_coef) {}
	*/
    Linear(double mult_coef) : mult_coef_(mult_coef) {}
    /**
    * Destructor.
    */
    ~Linear(void) {}
	/**
	* 	Just the test to define multiplication coefficient other way.
	*/
	void reinit(double mult_coef);
    /**
    * 	Operator.
    */
    double operator()(double x) { return (mult_coef_*x); }

private:
    double mult_coef_;
};

class Freundlich {
public:
	/**
	* 	Constructor.
	*/
	Freundlich(double mult_coef, double exponent) : mult_coef_(mult_coef), exponent_(exponent){}
	/**
	* 	Destructor.
	*/
	~Freundlich(void){}
	/**
	* 	Operator.
	*/
	double operator()(double x){ return (mult_coef_*pow(x, exponent_)); }
private:
	double mult_coef_;
	double exponent_;
};


/**
* Class describing one isotherm with possibly precalculated interpolation table.
*/
class Isotherm {
public:
    /**
     * Initialization of the isotherm approximation.
     * @p isotherm is a functor object representing the isotherm. @p rock_density and @p porosity are
     * material parameters and final parameter is the @p molar_density of the adsorbed substance.
     */
    void reinit(enum SorptionType sorption_type, double rock_density, double aqua_density, double porosity, double molar_mass, double c_aqua_limit, bool dual_porosity_on, double phi);
    /**
     *
     */
    template<class Func>
    void make_table(const Func &isotherm, int n_points); // const Func &isotherm
    /**
     * Find new values for concentrations @p c_aqua, @p c_sorbed that has same total mass and lies on the
     * @p isotherm (functor object).
     */
    template<class Func>
    void solve_conc(double &c_aqua, double &c_sorbed, const Func &isotherm); // , double elem_volume); // const Func &isotherm
    /**
     * Update concentrations.
     */
    bool compute_projection(double &c_aqua, double &c_sorbed);
    /**
    *	Enables to set private parameter.
    */
    void set_inv_scale_aqua(double inv_scale_aqua);
    /**
    *	Enables to set private parameter.
    */
    void set_inv_scale_sorbed(double inv_scale_sorbed);
    /**
    *
    */
    void set_rho_aqua(double rho_aqua);
    /**
    *
    */
    double get_rho_aqua(void);
    /**
    *	Enables to set private parameter.
    */
    void set_scale_aqua(double scale_aqua);
    /**
    *	Enables to get private parameter.
    */
    double get_scale_aqua(void);
    /**
    *	Enables to set private parameter.
    */
    void set_scale_sorbed(double scale_sorbed);
    /**
    *	Enables to get private parameter.
    */
    double get_scale_sorbed(void);
    /**
    *	Enables to set private parameter.
    */
    void set_caq_limmit(double caq_limmit);
    /**
    * 	Verifies how big interpolation table is defined
    */
    int get_interpolation_table_size(void);
    /**
    * 	Creates interpolation table containing just one point
    */
    void make_one_point_table(void);
    /**
    *  Returns sorption type
    */
    SorptionType get_sorption_type(void);
    /**
    *  Sets sorption type
    */
    void set_sorption_type(SorptionType sorp_type);
    /**
    *
    */
    void set_mult_coef_(double mult_coef);
    /**
    *
    */
    double get_mult_coef_(void);
    /**
    *
    */
    void set_second_coef_(double second_coef);
    /**
    *
    */
    double get_second_coef_(void);
    /**
    *
    */
    void precipitate(double &c_aqua, double &c_sorbed); //, double scale_aqua, double scale_sorbed);
private:
    /// density of the solvent
    double rho_aqua_;
    /// coefficient that convert soluted concentration to mass; porosity = k_W, originally rho_aqua*porosity = k_W
    double scale_aqua_;
    /// coefficient that convert adsorbed molar concentration to mass; molar_weight * rho_rock * (1 - porosity) = k_H
    double scale_sorbed_;
    /// reciprocal values
    double inv_scale_aqua_, inv_scale_sorbed_;
    /// Limit concentration in solution, we model coagulation as adsorption
    double c_aqua_limit_;
    /// Type of isotherm
    SorptionType sorption_type_;
    /**
    * 	Multiplication parameter of the isotherm
    */
    double mult_coef_;
    /**
    * 	Second potential parameter of the isotherm
    */
    double second_coef_;
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
    CrossFunction(const Func &func_,  double total_mass, double scale_aqua, double scale_sorbed)
    : func(func_), total_mass_(total_mass), scale_sorbed_(scale_sorbed), scale_aqua_(scale_aqua) {}

    double operator()( double conc_aqua)
    {
        return scale_sorbed_*func( conc_aqua ) + (scale_aqua_) * conc_aqua - total_mass_; // that is the  selected isotherm // scale_sorbed_ * func( conc_aqua ) + scale_aqua_ * conc_aqua - total_mass_
    }
private:
    Func func;
    double total_mass_, scale_sorbed_, scale_aqua_;
};

#endif /* SORPTION_IMPL_HH_ */
