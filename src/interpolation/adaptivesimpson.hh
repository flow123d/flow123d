#ifndef ADAPTIVESIMPSON_H
#define ADAPTIVESIMPSON_H

#include "functors.hh"

#define MAX_RECURSION 1e+7


///Static class implementing integration using Simpson's rule.
/** Uses 3-point Simpson's rule to evaluate an intergral on interval a,b.
  * Divides interval in halves in recusion until the difference
  * between values of Simpson's rule before and after division is
  * smaller then $15\epsilon$.
  */
class AdaptiveSimpson
{
  private:
    ///Evaluates the Simpson's rule.
    static double Simpson ( const double& h, const double &fa, const double &fc, const double &fb );
  
    ///the recursive method
    static double SimpsonAd( FunctorBase<double> &func, 
		      const double& h, const double &a, const double &c, const double &b,
		      const double &fa, const double &fc, const double &fb, 
		      const double &sx, const double &tol, long &recursion );
      
  public:
	       
    ///main method that starts the evaluation and calls the recursion
    static double AdaptSimpson( FunctorBase<double> &func,
		     const double& a, const double& b, 
		     const double& tol );	    		
};
  


#endif	//ADAPTIVESIMPSON_H