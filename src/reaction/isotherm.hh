/*
 * isotherm.hh
 *
 *  Created on: Mar 7, 2013
 *      Author: jb
 *
 *
 *  TODO:
 *  - fix identificators to something meaningfull (in english):
 *    dekl_tolerance (changed to tolerance), precipitate is correct because These are 1) precipitation of an element due to supersaturation of that element in the soil solution, 2) surface precipitation, and 3) co-precipitation of elements. These chemical reactions can potentially control the mobility and toxicity of trace metals in soil systems.,
 *    isotherm.cc - iso_ind_floor (someone changed it), etc.
 *  - isotherm.cc, compute_projection:
 *    - (done) index must be less then size -1, in order to get correct linear interpolation
 *    - make this function never fail - use solve_conc when outside of the interval
 *      This is not so easy  since solve_cons is function template. Maybe we can
 *      store appropriate function pointer together with the table
 *    - proposed minor optimization
 *  - setters/getters are useless
 *  - reinit - do not pass phi and dual_porosity, but aqua_fracture (porosity) and rock_fracture, aqua_fraction, rock_fraction
 *    better: pass all parameters of one isotherm through reinit, make swith according type there
 *
 *
 *
 *
 *

 /// Pair of soluted concentration and adsorbed concentration
 typedef pair<double, double> ConcPair;
 class IsothermFactory {
     void reint(int adsorption_type, ... ) {
         tmp_isotherm.reinit(...);
         // save functor parameters
     }
     Conc Pair solve_conc( ConcPair conc ) {
          switch (adsoption_type_) {
            case Isotherm::linear:
                Linear functor(...);
                return tmp_isotherm.solve_conc(conc, functor);
                break;
            case ....

          }
     void make_table(Isotherm &table) {
          switch (adsoption_type_) {
            case Isotherm::linear:
                Linear functor(...);
                table = tmp_isotherm; // copy reinit
                table.make_table(functor);
                //table.set_factory(this); // allows compute values outside of the table
                break;
            case ....
     }

     Isotherm tmp_isotherm;

 }

 */

#ifndef SORPTION_IMPL_HH_
#define SORPTION_IMPL_HH_

#include <vector>
#include <input/input_type.hh>
#include <boost/math/tools/roots.hpp>
#include "fields/field_base.hh"

typedef Field<3, FieldValue<3>::Scalar > * pScalar;
typedef pair<double, double> ConcPair;

enum SorptionType {
	none = 0,
	linear = 1,
	freundlich = 2,
	langmuir = 3
};

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
 * Functor for Langmuir isotherm.
 */
class Langmuir {
public:
	/**
	* 	Original constructor, Langmuir( double mult_coef, double alpha) : alpha(alpha), mult_coef_(mult_coef) {}
	*/
    Langmuir( double mult_coef, double alpha) : mult_coef_(mult_coef), alpha_(alpha) {}
    /**
    * 	Destructor.
    */
    ~Langmuir(void){}
	/**
	* 	Just the test to define coefficients other way.
	*/
	void reinit(double mult_coef, double alpha);
    /**
    * 	Operator.
    */
    double operator()( double x) { return (mult_coef_*(alpha_ * x)/(alpha_ *x + 1)); }

private:
    double mult_coef_;
    double alpha_;
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
	* 	Just the test to define coefficients other way.
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
	* 	Just the test to define multiplication coefficient other way.
	*/
	void reinit(double mult_coef, double exponent);
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
	void reinit(enum SorptionType sorption_type, double aqua_density, double scale_aqua, double scale_sorbed, double c_aqua_limit, double mult_coef, double second_coef);
	/**
     *
     */
    template<class Func>
    void make_table(const Func &isotherm, int n_points);
    /**
    *
    */
    void make_table(int n_points);
    /**
     * Find new values for concentrations @p c_aqua, @p c_sorbed that has same total mass and lies on the
     * @p isotherm (functor object).
     */
    template<class Func>
    void solve_conc(double &c_aqua, double &c_sorbed, const Func &isotherm);
    /**
    *
    */
    ConcPair solve_conc(ConcPair conc);
    /**
     * Update concentrations.
     */
    bool compute_projection(double &c_aqua, double &c_sorbed);
    /**
    *  Returns sorption type
    */
    SorptionType get_sorption_type(void);
    /**
    *
    */
    void set_iso_params(SorptionType sorp_type, double mult_coef, double second_coef);
    /**
    *
    */
    void set_kind_of_pores(int kind_of_pores);
    /**
    *
    */
    void precipitate(double &c_aqua, double &c_sorbed); //double &c_aqua, double &c_sorbed); //
    /**
    * Informs ifever the interpolation table is precomputed and has more than one cell
    */
    int is_precomputed(void);
    /// Type of isotherm
    enum SorptionType adsorption_type_;
    /**
    * 	Multiplication parameter of the isotherm
    */
    double mult_coef_;
    /**
    * 	Second potential parameter of the isotherm
    */
    double second_coef_;
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
    /**
     * Interpolation table of isotherm in the rotated coordinates.
     * The X axes of rotated system is total mass, the Y axes is perpendicular.
     */
    vector<double> interpolation_table;
    /**
     * Step on the rotated X axes (total mass).
     */
    double total_mass_step_;
    /**
    *
    */
    int kind_of_pores_;
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
    : func(func_), total_mass_(total_mass), scale_sorbed_(scale_sorbed), scale_aqua_(scale_aqua), rho_aqua_(rho_aqua) {}

    double operator()( double conc_aqua)
    {
        return scale_sorbed_*func( conc_aqua/rho_aqua_) + (scale_aqua_) * conc_aqua - total_mass_; // that is the  selected isotherm // scale_sorbed_ * func( conc_aqua ) + scale_aqua_ * conc_aqua - total_mass_
    }
private:
    Func func;
    double total_mass_, scale_sorbed_, scale_aqua_, rho_aqua_;
};

/*class IsothermFactory
{
public:
	void reinit(Isotherm &isotherm, enum SorptionType adsorption_type, double rho_aqua, double scale_aqua, double scale_sorbed, double c_aqua_limit, double mult_coef, double second_coef)
	{
        //tmp_isotherm = isotherm;
        if((mult_coef_ = mult_coef) < 0.0 ) xprintf(UsrErr,"Multiplication coefficient for an isotherm must equal to or higher then 0.0");
        second_coef_ = second_coef;
        // save functor parameters
		isotherm.reinit(adsorption_type, rho_aqua, scale_aqua, scale_sorbed, c_aqua_limit, mult_coef, second_coef);
    }
    ConcPair solve_conc( Isotherm &isotherm, ConcPair conc ) {
    	double c_aqua = conc.first;
    	double c_sorbed = conc.second;
    	switch (isotherm.adsorption_type_) {
    	   case none:
    	   {
    		   xprintf(Msg,"IsothermFactory solve_conc not computed. No sorption is considered.\n");
    		   break;
    	   }
    	   case linear:
    	   {
    		   Linear iso_functor(mult_coef_);
    		   isotherm.solve_conc(c_aqua, c_sorbed, iso_functor);
    		   conc.first = c_aqua;
    		   conc.second = c_sorbed;
               return conc;
               break;
    	   }
           case freundlich:
           {
        	   Freundlich iso_functor(mult_coef_, second_coef_);
        	   isotherm.solve_conc(c_aqua, c_sorbed, iso_functor);
    		   conc.first = c_aqua;
    		   conc.second = c_sorbed;
               return conc;
               break;
           }
           case langmuir:
           {
        	   Langmuir iso_functor(mult_coef_, second_coef_);
        	   isotherm.solve_conc(c_aqua, c_sorbed, iso_functor);
    		   conc.first = c_aqua;
    		   conc.second = c_sorbed;
               return conc;
               break;
           }
           default:
           {
        	   xprintf(UsrErr,"IsothermFactory solve_conc failed."); // Sorption of %d-th specie in %d-th region has unknown type %d.", i_subst, reg_idx, hlp_iso_type)
        	   break;
           }
         }
    }

    void make_table(Isotherm &isotherm)
    {
    	int n_points=100; //gives the number of samples in the interpolation table

         switch (isotherm.adsorption_type_)
         {
           case none:
           {
        	   isotherm.make_one_point_table();
        	   break;
           }
           case linear:
           {
        	   Linear iso_functor(mult_coef_);
               isotherm.make_table(iso_functor, n_points);
               break;
           }
           case freundlich:
           {
        	   Freundlich iso_functor(mult_coef_, second_coef_);
               isotherm.make_table(iso_functor, n_points);
        	   break;
           }
           case langmuir:
           {
        	   Langmuir iso_functor(mult_coef_, second_coef_);
               isotherm.make_table(iso_functor, n_points);
        	   break;
           }
           default:
           {
        	   xprintf(UsrErr,"IsothermFactory make_table failed. Unknown type of an isotherm."); // Sorption of %d-th specie in %d-th region has unknown type %d.", i_subst, reg_idx, hlp_iso_type)
        	   break;
           }
         }
    }

    //Isotherm tmp_isotherm; // o jeden parametr min pri reinitu
    double mult_coef_, second_coef_;
};*/

#endif /* SORPTION_IMPL_HH_ */
