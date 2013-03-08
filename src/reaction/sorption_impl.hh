/*
 * sorption_impl.hh
 *
 *  Created on: Mar 7, 2013
 *      Author: jb
 */

#ifndef SORPTION_IMPL_HH_
#define SORPTION_IMPL_HH_


/**
 * Functor for Langmuir isotherm.
 */
class Langmuir {
public:
    Langmuir( double mult_coef, double alpha) : alpha(alpha), mult_coef(mult_coef) {}

    double operator()( double x) { return mult_coef*(alpha * x)/(1+alpha*x); }

private:
    double mult_coef;
    double alpha;
};


/**
 * Functor for linear isotherm
 */
class Linear {
public:
    Linear(double mult_coef) : direct(mult_coef) {}

    double operator()(double x) { return direct*x; }

private:
    double direct;
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
    void reinit(double rock_density, double porosity, double molar_mass);

    /**
     *
     */
    template<class Func>
    void make_table(const Func &isotherm);

    /**
     * Find new values for concentrations @p c_aqua, @p c_sorbed that has same total mass and lies on the
     * @p isotherm (functor object).
     */
    template<class Func>
    void solve_conc(double &c_aqua, double &c_sorbed, const Func &isotherm);

    /**
     * Update concentrations.
     */
    inline void update_conc(double &c_aqua, double &c_sorbed);
private:

    /// coefficient that convert soluted concentration to mass; rho_aqua*porosity
    double scale_aqua;
    /// coefficient that convert adsorbed molar concentration to mass; molar_weight * rho_rock * (1 - porosity)
    double scale_sorbed;
    /// reciprocal values divided by 2
    double inv_scale_aqua, inv_scale_sorbed;
    /// Limit concentration in solution, we model coagulation as adsorption
    double c_aqua_limit_;
    /**
     * Interpolation table of isotherm in the rotated coordinates.
     * The X axes of rotated system is total mass, the Y axes is perpendicular.
     */
    vector<double> interpolation_table;
    /**
     * Step on the rotated X axes (total mass).
     */
    double total_mass_step;

};



void Isotherm::reinit(double rock_density, double porosity, double molar_mass, double c_aqua_limit) {
    // set class variables
    scale_aqua = porosity * rho_aqua;
    scale_sorbed = (1-porosity) * rock_density * molar_mass;
    inv_scale_aqua = 1/scale_aqua/2;
    inv_scale_sorbed = 1/scale_sorbed/2;
    c_aqua_limit_=c_aqua_limit;
}



inline bool Isotherm::update_conc(double &c_aqua, double &c_sorbed) {

    double total_mass = scale_aqua* c_aqua + scale_sorbed * c_sorbed;
    unsigned int i_total_mass = total_mass / total_mass_step;
    if (i_total_mass < 0) return false;
    if (i_total_mass < interpolation_table.size()) {
        double rot_sorbed = interpolation_table[i_total_mass]; // TODO: linear interpolaton
        c_aqua = (total_mass - rot_sorbed) * inv_scale_aqua;
        c_sorbed = (total_mass + rot_sorbed) * inv_scale_sorbed;
        return true;
    } else {
        if (c_aqua_limit > 0.0) {
            c_sorbed = (total_mass - scale_c_aqua* c_aqua_limit)/scale_sorbed;
            c_aqua = c_aqua_limit;
        } else return false;
    }
    return false;
}


/**
 *  Functor for solved equation in form F(x) ==0.
 *  Function @p func is an isotherm functor object in concentration based coordinated system.
 *  We solve the equation in modified system (scaled and rotated) for better numerical stability.
 *  The solved equation reads:
 *  F(X) -Y =0, where
 *  X is total mass , Y
 */
template <class Func>
class CrossFunction {
public:
    CrossFunction(const Func &func,  double total_mass, double scale_aqua, double scale_sorbed)
    : func(func_), total_mass(total_mass_), scale_aqua(scale_aqua), scale_sorbed(scale_sorbed) {}


    double operator()( double conc_aqua) {
        return scale_sorbed * func( conc_aqua ) - total_mass + scale_aqua * conc_aqua;
    }
private:
    Func func;
    double total_mass, scale_sorbed, scale_aqua;
};


template<class Func>
void solve_conc(double &c_aqua, double &c_sorbed, const Func &isotherm) {

    double total_mass = scale_aqua*c_aqua + scale_sorbed * c_sorbed;
    CrossFunction<Func> eq_func(isotherm, total_mass, scale_aqua, scale_sorbed);
    pair<double,double> solution = boost::math::tools::toms748_solve(
                                                eq_func, 0, total_mass / scale_aqua,
                                                boost::math::tools::eps_tolerance<double>(60), 100);
    c_sorbed = (total_mass - scale_aqua * solution.first) / c_sorbed;
}




template<class Func>
void Isotherm::make_table(const Func &isotherm) {
    double mass_limit;
    if (c_aqua_limit >0) {
        // set mass_limit from c_aqua_limit
    } else {
        // set mass_limit from max conc = 1
    }
    total_mass_step = mass_limit / n_steps;
    double mass = total_mass_step; // we need not to save value for zero mass, it is zero
    for(int i=0; i< n_steps;i++, mass+=total_mass_step) {
        double c_aqua = mass * inv_scale_aqua;
        double c_sorbed = mass * inv_scale_sorbed;
        solve_conc(c_aqua, c_sorbed, isotherm);
        interpolation_table.push_back( c_sorbed * scale_sorbed - c_aqua * scale_aqua);
    }
}


#endif /* SORPTION_IMPL_HH_ */
