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
    * 	Mysterious operator.
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
    * 	Mysterious operator.
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

/*void Linear::reinit(double mult_coef)
{
	mult_coef_ = mult_coef;
	return;
};*/


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
    void reinit(enum SorptionType sorption_type, double rock_density, double aqua_density, double porosity, double molar_mass, double c_aqua_limit);
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
    //inline
    bool compute_projection(double &c_aqua, double &c_sorbed);
    //bool compute_projection(double &c_aqua);
    /**
    *	Enables to get private parameter.
    */
    double get_scale_aqua(void);
    /**
    *	Enables to get private parameter.
    */
    double get_scale_sorbed(void);
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
    void precipitate(double &c_aqua, double &c_sorbed, double scale_aqua, double scale_sorbed); // , double elem_volume);
private:
    /**
    * 	Suppresses the use of implicit constructor.
    */
    //Isotherm();
    /// coefficient that convert soluted concentration to mass; rho_aqua*porosity = k_W
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

/*void Isotherm::reinit(enum SorptionType sorp_type, double rock_density, double rho_aqua, double porosity, double molar_mass, double c_aqua_limit)
{
    // set class variables
	sorption_type = sorp_type;
    scale_aqua = porosity * rho_aqua;
    scale_sorbed = (1-porosity) * rock_density * molar_mass;
    inv_scale_aqua = 1/scale_aqua/2;
    inv_scale_sorbed = 1/scale_sorbed/2;
    c_aqua_limit_=c_aqua_limit;
};



inline bool Isotherm::compute_projection(double &c_aqua, double &c_sorbed) //clear as glass
{
    double total_mass = scale_aqua* c_aqua + scale_sorbed * c_sorbed;
    unsigned int i_total_mass = total_mass / total_mass_step;
    if (i_total_mass < 0) return false;
    if (i_total_mass < interpolation_table.size()) {
    	int iso_ind_floor, iso_ind_ceil;
    	iso_ind_floor = (int)(total_mass/(total_mass_step)); iso_ind_ceil = iso_ind_floor + 1;
    	double rot_sorbed = interpolation_table[iso_ind_floor] + (total_mass - iso_ind_floor*total_mass_step)*(interpolation_table[iso_ind_ceil] - interpolation_table[iso_ind_floor])/total_mass_step;
        c_aqua = (total_mass + rot_sorbed) * inv_scale_aqua;
        c_sorbed = (total_mass - rot_sorbed) * inv_scale_sorbed;
        return true;
    } else {
        if (c_aqua_limit_ > 0.0) {
            c_sorbed = (total_mass - scale_aqua* c_aqua_limit_)*inv_scale_sorbed;
            c_aqua = c_aqua_limit_;
        } else return false;
    }
    return false;
};*/


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


/*template<class Func>
void Isotherm::solve_conc(double &c_aqua, double &c_sorbed, const Func &isotherm) // const Func &isotherm
{
    double mass_limit;
    boost::uintmax_t max_iter=100;
    boost::math::tools::eps_tolerance<double> toler(60);
    Func &iso_hlp = const_cast<Func &>(isotherm);
    double f_max = iso_hlp(c_aqua_limit_);
    if (c_aqua_limit_ >0) {
        mass_limit = scale_aqua*c_aqua_limit_ + scale_sorbed*f_max; // isotherm(c_aqua_limit_);
    } else {
        mass_limit = scale_aqua + scale_sorbed;// set mass_limit from max conc = 1, needs to be computed somehow else
    }
	double total_mass = scale_aqua*c_aqua + scale_sorbed * c_sorbed;
    CrossFunction<Func> eq_func(isotherm, total_mass, scale_aqua, scale_sorbed); // equation describing one point on the isotherm
    pair<double,double> solution = boost::math::tools::toms748_solve(eq_func, 0.0, 10.0, toler, max_iter);
    //c_sorbed = (total_mass - scale_aqua * solution.first) / scale_sorbed;
    //PROBABLY COMPLETELY WRONG, SOLUTION IS AN INTERVAL CONTAINING SOLUTION, because of that following two lines are commented
    //c_aqua = (total_mass - scale_sorbed * solution.first) / scale_aqua;
    //c_sorbed = (total_mass - scale_aqua * solution.second) / scale_sorbed;
    //MUST BE REPARED, LATER.
    c_aqua = 1.0;
    c_sorbed = 1.0;
};

template<class Func>
void Isotherm::make_table(const Func &isotherm, int n_steps) { //const Func &isotherm, int n_steps
    double mass_limit, c_aqua, c_sorbed;
    Func &iso_hlp = const_cast<Func &>(isotherm);
    //double f_max = isotherm(c_aqua_limit_);
    double f_max = iso_hlp(c_aqua_limit_);
    if (c_aqua_limit_ >0) {
        mass_limit = scale_aqua*c_aqua_limit_ + scale_sorbed*f_max; //isotherm(c_aqua_limit_);
    } else {
        mass_limit = scale_aqua + scale_sorbed;// set mass_limit from max conc = 1, needs to be computed somehow else
    }
    total_mass_step = mass_limit / n_steps;
    double mass = total_mass_step; // we need not to save value for zero mass, it is zero
    for(int i=0; i< n_steps;i++, mass+=total_mass_step) {
        double c_aqua = mass * inv_scale_aqua;
        double c_sorbed = mass * inv_scale_sorbed;
        solve_conc(c_aqua, c_sorbed, isotherm); //isotherm);
        interpolation_table.push_back( c_sorbed * scale_sorbed - c_aqua * scale_aqua);
    }
};

double Isotherm::get_scale_aqua(void)
{
	return scale_aqua;
};

double Isotherm::get_scale_sorbed(void)
{
	return scale_sorbed;
};*/

#endif /* SORPTION_IMPL_HH_ */
