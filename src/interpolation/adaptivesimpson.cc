
#include "adaptivesimpson.hh"
#include "functors_impl.hh"
#include "system/xio.h"

#include <math.h>

double AdaptiveSimpson::Simpson(const double& h, const double &fa, 
				const double &fc, const double &fb)
{
  return (fa + fb + 4.0*fc)*(h/6.0);
}

double AdaptiveSimpson::SimpsonAd(Functor<double> &func, 
				  const double& h, const double &a, 
				  const double &c, const double &b,
				  const double &fa, const double &fc, 
				  const double &fb, const double &sx, 
				  const double &tol, long &recursion)
{
  recursion++;
  double ca = 0.5*(a+c);
  double cb = 0.5*(c+b);
  double fca = func(ca);
  double fcb = func(cb);
   
  double h2 = 0.5*h;
  double sa = Simpson(h2,fa,fca,fc);
  double sb = Simpson(h2,fc,fcb,fb);
  
  double err_est = (sa+sb-sx)/15.0;
  if (fabs(err_est) <= tol || recursion >= MAX_RECURSION)
  {
    return sa+sb+err_est;
  }
  else
  {
    //std::cout << "simpsonad -else " << recursion << std::endl;
    //std::cout << h2 << "  " << a << "  " << ca << "  " << c << "  " << fa  << "  " << fca  << "  " << fc  << "  " << sa << std::endl;
    return SimpsonAd(func,h2,a,ca,c,fa,fca,fc,sa,tol,recursion) 
	   + SimpsonAd(func,h2,c,cb,b,fc,fcb,fb,sb,tol,recursion);
  }
}

double AdaptiveSimpson::AdaptSimpson( Functor<double> &func,
				      const double& a, const double& b, 
				      const double& tol )
{
  //DBGMSG("AdaptiveSimpson: a(%f) b(%f)\n",a,b);
  double c = 0.5*(b+a);
  double fa,fb,fc;
  double sx,res;
  fa = func(a);
  fb = func(b);
  fc = func(c);
  long recursion = 0;
  //DBGMSG("fa=%f \tfb=%f \tfc=%f\n",fa,fb,fc);
  sx = Simpson(b-a, fa, fc, fb);
  res = SimpsonAd(func,b-a,a,c,b,fa,fc,fb,sx,tol,recursion);
  //DBGMSG("\trecursions: %f\n",recursion);
  return res;
}




