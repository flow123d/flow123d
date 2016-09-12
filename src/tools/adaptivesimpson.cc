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
 * @file    adaptivesimpson.cc
 * @brief   
 */

#include "adaptivesimpson.hh"
#include "functors_impl.hh"

#include <math.h>

double AdaptiveSimpson::Simpson(const double& h, const double &fa, 
				const double &fc, const double &fb)
{
  return (fa + fb + 4.0*fc)*(h/6.0);
}

double AdaptiveSimpson::SimpsonAd(FunctorBase<double> &func, 
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
    //DebugOut().fmt("simpsonad, recursion:  {} \t {}\n", recursion, err_est);
    //std::cout << "simpsonad -else " << recursion << std::endl;
    //std::cout << h2 << "  " << a << "  " << ca << "  " << c << "  " << fa  << "  " << fca  << "  " << fc  << "  " << sa << std::endl;
    return SimpsonAd(func,h2,a,ca,c,fa,fca,fc,sa,tol,recursion) 
	   + SimpsonAd(func,h2,c,cb,b,fc,fcb,fb,sb,tol,recursion);
  }
}

double AdaptiveSimpson::AdaptSimpson( FunctorBase<double> &func,
				      const double& a, const double& b, 
				      const double& tol )
{
  //DebugOut().fmt("AdaptiveSimpson: a({}) b({})\n",a,b);
  double c = 0.5*(b+a);
  double fa,fb,fc;
  double sx,res;
  fa = func(a);
  fb = func(b);
  fc = func(c);
  long recursion = 0;
  //DebugOut().fmt("fa={} \tfb={} \tfc={}\n",fa,fb,fc);
  sx = Simpson(b-a, fa, fc, fb);
  res = SimpsonAd(func,b-a,a,c,b,fa,fc,fb,sx,tol,recursion);
  //DebugOut().fmt("\trecursions: {}\n",recursion);
  return res;
}




