/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief Definitions of particular quadrature rules on simplices.
 *  @author Jan Stebel
 */

#include "quadrature_lib.hh"

template <>
QGauss<1>::QGauss (const unsigned int n) :
        Quadrature<1> ((const unsigned int)n)
{
  if (n == 0)
    return;

  const unsigned int m = (n+1)/2;

                                   // tolerance for the Newton
                                   // iteration below. we need to make
                                   // it adaptive since on some
                                   // machines (for example PowerPC)
                                   // long double is the same as
                                   // double -- in that case we can
                                   // only get to a certain multiple
                                   // of the accuracy of double there,
                                   // while on other machines we'd
                                   // like to go further down
                   //
                   // the situation is complicated by
                   // the fact that even if long
                   // double exists and is described
                   // by std::numeric_limits, we may
                   // not actually get the additional
                   // precision. One case where this
                   // happens is on x86, where one can
                   // set hardware flags that disable
                   // long double precision even for
                   // long double variables. these
                   // flags are not usually set, but
                   // for example matlab sets them and
                   // this then breaks deal.II code
                   // that is run as a subroutine to
                   // matlab...
                                   //
                                   // a similar situation exists, btw,
                                   // when running programs under
                                   // valgrind up to and including at
                                   // least version 3.3: valgrind's
                                   // emulator only supports 64 bit
                                   // arithmetic, even for 80 bit long
                                   // doubles.
#ifdef HAVE_STD_NUMERIC_LIMITS
  const long double
    long_double_eps = static_cast<long double>(std::numeric_limits<long double>::epsilon()),
    double_eps      = static_cast<long double>(std::numeric_limits<double>::epsilon());
#else
  const long double
    long_double_eps = 1.09e-19L,
    double_eps      = 2.23e-16L;
#endif

                   // now check whether long double is more
                   // accurate than double, and set
                   // tolerances accordingly. generate a one
                   // that really is generated at run-time
                   // and is not optimized away by the
                   // compiler. that makes sure that the
                   // tolerance is set at run-time with the
                   // current behavior, not at compile-time
                   // (not doing so leads to trouble with
                   // valgrind for example).
  volatile long double runtime_one = 1.0;
  const long double tolerance
    = (runtime_one + long_double_eps != runtime_one
       ?
       max (double_eps / 100, long_double_eps * 5)
       :
       double_eps * 5
       );


  for (unsigned int i=1; i<=m; ++i)
    {
      long double z = cos(math::pi() * (i-.25)/(n+.5));

      long double pp;
      long double p1, p2, p3;

                                       // Newton iteration
      do
    {
                       // compute L_n (z)
      p1 = 1.;
      p2 = 0.;
      for (unsigned int j=0;j<n;++j)
        {
          p3 = p2;
          p2 = p1;
          p1 = ((2.*j+1.)*z*p2-j*p3)/(j+1);
        }
      pp = n*(z*p1-p2)/(z*z-1);
      z = z-p1/pp;
    }
      while (abs(p1/pp) > tolerance);

      double x = .5*z;
      this->quadrature_points[i-1](0) = .5-x;
      this->quadrature_points[n-i](0) = .5+x;

      double w = 1./((1.-z*z)*pp*pp);
      this->weights[i-1] = w;
      this->weights[n-i] = w;
    }
}

