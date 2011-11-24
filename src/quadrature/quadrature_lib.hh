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

#ifndef QUADRATURE_LIB_HH_
#define QUADRATURE_LIB_HH_

#include "quadrature/quadrature.hh"


/**
 * Gauss-Legendre quadrature of arbitrary order.
 *
 * The coefficients are computed by a function from <tt>Numerical Recipes</tt>.
 *
 */
template<unsigned int dim>
class QGauss : public Quadrature<dim> {
public:
    /**
     * Generate a formula of given order (exact for polynomials of degree @p order).
     */
    QGauss(const unsigned int order);
};





#endif /* QUADRATURE_LIB_HH_ */
