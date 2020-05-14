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
 * @file    quadrature_lib.cc
 * @brief   Definitions of particular quadrature rules on simplices.
 * @author  Jan Stebel
 */

#include "system/global_defs.h"
#include "system/system.hh"
#include "quadrature/quadrature_lib.hh"


#define FLOAT double
#define SHORT int
#define _F(n) n
#include "quad.c"
#undef FLOAT
#undef SHORT
#undef _F

typedef std::vector<QUAD *> DimQuadList;
std::vector<DimQuadList> __gauss_quadratures = {
        {},
        { QUAD_1D_P1,
                QUAD_1D_P1, QUAD_1D_P2, QUAD_1D_P3,
                QUAD_1D_P4, QUAD_1D_P5, QUAD_1D_P6,
                QUAD_1D_P7, QUAD_1D_P8, QUAD_1D_P9,
                QUAD_1D_P10,    QUAD_1D_P11,    QUAD_1D_P12,
                QUAD_1D_P13,    QUAD_1D_P14,    QUAD_1D_P15,
                QUAD_1D_P16,    QUAD_1D_P17,    QUAD_1D_P18,
                QUAD_1D_P19,    QUAD_1D_P20,    QUAD_1D_P21
              },
        { QUAD_2D_P1,
                QUAD_2D_P1, QUAD_2D_P2, QUAD_2D_P3, QUAD_2D_P4,
                QUAD_2D_P5, QUAD_2D_P6, QUAD_2D_P7, QUAD_2D_P8,
                QUAD_2D_P9, QUAD_2D_P10,    QUAD_2D_P11,    QUAD_2D_P12,
                QUAD_2D_P13,    QUAD_2D_P14,    QUAD_2D_P15,    QUAD_2D_P16,
                QUAD_2D_P17,    QUAD_2D_P18,    QUAD_2D_P19,    QUAD_2D_P20,
                QUAD_2D_P21
              },
        { QUAD_3D_P1,
                QUAD_3D_P1, QUAD_3D_P2, QUAD_3D_P3, QUAD_3D_P4,
                QUAD_3D_P5, QUAD_3D_P6, QUAD_3D_P7, QUAD_3D_P8,
                QUAD_3D_P9, QUAD_3D_P10,    QUAD_3D_P11,    QUAD_3D_P12,
                QUAD_3D_P13,    QUAD_3D_P14
              }
};

double __unit_cell_volume[] = { 1, 1, 0.5, 1./6 };

template<uint dimension>
void QGauss::init(uint order) {
    DimQuadList & quads = __gauss_quadratures[dimension];
    OLD_ASSERT(order < quads.size(), "Quadrature of given order is not implemented.");
    auto &point_list = quads[order];

    this->quadrature_points.reinit(point_list->npoints);
    this->weights.resize(0);
    for (int i=0; i<point_list->npoints; i++)
    {
        Armor::ArmaVec<double, dimension> p(& point_list->points[i*(dimension + 1)]);
        this->quadrature_points.append(p);
        this->weights.push_back(point_list->weights[i] * __unit_cell_volume[dimension]);
    }
}


QGauss::QGauss(unsigned int dim, unsigned int order)
: Quadrature(dim)
{

    switch (dim)
    {
    case 0:
        // Quadrature on 0-dim element have single quadrature point
        // with 0 local coordinates.
        // No way to append 0-dim arma vec, we just resize the Array.
        this->quadrature_points.reinit(1);
        this->quadrature_points.resize(1);
        this->weights.push_back(1);
        ASSERT_EQ_DBG(size(), 1);
        return;
    case 1:
        init<1>(order);
        break;
    case 2:
        init<2>(order);
        break;
    case 3:
        init<3>(order);
        break;
    }


}




