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
 * @file    quadrature_lib.hh
 * @brief   Definitions of particular quadrature rules on simplices.
 * @author  Jan Stebel
 */

#ifndef QUADRATURE_LIB_HH_
#define QUADRATURE_LIB_HH_

#include <array>
#include "quadrature/quadrature.hh"



/**
 * @brief Symmetric Gauss-Legendre quadrature formulae on simplices.
 *
 * Symmetric Gauss-Legendre quadrature on a @p dim dimensional simplex.
 * The coefficients are taken from Parallel Hierarchical Grid project.
 *
 */
class QGauss : public Quadrature {
public:
    typedef std::array<QGauss, 4> array;
    /**
     * Create a
     */
    inline static std::array<QGauss, 4> make_array(unsigned int order)
    {
        return {QGauss(0, order), QGauss(1, order), QGauss(2, order), QGauss(3, order)};
    }

    /**
     * @brief Create a formula of given order.
     *
     * The formula is exact for polynomials of degree @p order.
     */
    QGauss(unsigned int dim, const unsigned int order);

    template<uint dimension>
    void init(uint order);
};





#endif /* QUADRATURE_LIB_HH_ */
