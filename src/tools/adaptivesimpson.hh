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
 * @file    adaptivesimpson.hh
 * @brief   
 */

#ifndef ADAPTIVESIMPSON_H
#define ADAPTIVESIMPSON_H

#include "functors.hh"

#define MAX_RECURSION 1e+7


///Static class implementing integration using Simpson's rule.
/** Uses 3-point Simpson's rule to evaluate an intergral on interval a,b.
  * Divides interval in halves in recusion until the difference
  * between values of Simpson's rule before and after division is
  * smaller then \f$ 15\epsilon \f$.
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