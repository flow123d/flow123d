/*
 * isotherms.hh
 *
 *  Created on: Jan 15, 2013
 *      Author: jb
 */

#ifndef ISOTHERMS_HH_
#define ISOTHERMS_HH_


// Functors for isotherms.


class Langmuir {
public:
    Langmuir( double alpha)
    : alpha_(alpha) {}

    double operator()( double x) {
        return (alpha * x)/(1+alpha*x);
    }
private:
    double alpha_;
};


// Equation for intersection points for any isotherm
template <class Func>
class CrossFunction {
public:
    CrossFunction(const Func &func,  double total_mass, double V_rho)
    : func_(func), total_mass(total_mass), V_rho(V_rho) {}


    double operator()( double conc_aqua) {
        return total_mass -con_aqua * V_rho - func_(conc_aqua);
    }
private:
    Func func_;
    double total_mass, V_rho;
};


void solve() {
     Langmuir isotherm( alpha );
     V_rho=... // zavisi na latce a diskretizaci V_elem * rho_latka
     for( total_mass ) {
         CrossFunction<Langmuir> func(isotherm, total_mass, V_rho);
         conc_aqua = boost::toms748_solve( func, 0, conc_max, boost::eps_tolerance(60), 100) . first;
     }

}


#endif /* ISOTHERMS_HH_ */
