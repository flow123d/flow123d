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


using namespace arma;


template<unsigned int dim>
QGauss<dim>::QGauss(const unsigned int order)
{
    typedef QUAD* pQUAD;

    static const pQUAD
        q1d[] = { QUAD_1D_P1,
                QUAD_1D_P1, QUAD_1D_P2, QUAD_1D_P3,
                QUAD_1D_P4, QUAD_1D_P5, QUAD_1D_P6,
                QUAD_1D_P7, QUAD_1D_P8, QUAD_1D_P9,
                QUAD_1D_P10,    QUAD_1D_P11,    QUAD_1D_P12,
                QUAD_1D_P13,    QUAD_1D_P14,    QUAD_1D_P15,
                QUAD_1D_P16,    QUAD_1D_P17,    QUAD_1D_P18,
                QUAD_1D_P19,    QUAD_1D_P20,    QUAD_1D_P21
              },
        q2d[] = { QUAD_2D_P1,
                QUAD_2D_P1, QUAD_2D_P2, QUAD_2D_P3, QUAD_2D_P4,
                QUAD_2D_P5, QUAD_2D_P6, QUAD_2D_P7, QUAD_2D_P8,
                QUAD_2D_P9, QUAD_2D_P10,    QUAD_2D_P11,    QUAD_2D_P12,
                QUAD_2D_P13,    QUAD_2D_P14,    QUAD_2D_P15,    QUAD_2D_P16,
                QUAD_2D_P17,    QUAD_2D_P18,    QUAD_2D_P19,    QUAD_2D_P20,
                QUAD_2D_P21
              },
        q3d[] = { QUAD_3D_P1,
                QUAD_3D_P1, QUAD_3D_P2, QUAD_3D_P3, QUAD_3D_P4,
                QUAD_3D_P5, QUAD_3D_P6, QUAD_3D_P7, QUAD_3D_P8,
                QUAD_3D_P9, QUAD_3D_P10,    QUAD_3D_P11,    QUAD_3D_P12,
                QUAD_3D_P13,    QUAD_3D_P14
              };
    static const double unit_cell_volume[] = { 1, 1, 0.5, 1./6 };
    const pQUAD *q;
    int nquads;
    vec::fixed<dim> p;

    switch (dim)
    {
    case 0:
        this->quadrature_points.push_back(p);
        this->weights.push_back(1);
        return;
    case 1:
        q = q1d;
        nquads = sizeof(q1d) / sizeof(pQUAD);
        break;
    case 2:
        q = q2d;
        nquads = sizeof(q2d) / sizeof(pQUAD);
        break;
    case 3:
        q = q3d;
        nquads = sizeof(q3d) / sizeof(pQUAD);
        break;
    }

    ASSERT(order < nquads, "Quadrature of given order is not implemented.");

    for (int i=0; i<q[order]->npoints; i++)
    {
        for (int j=0; j<dim; j++)
            p(j) = q[order]->points[i*(dim+1)+j];

        this->quadrature_points.push_back(p);
        // The weights must be adjusted according to the volume of the unit cell:
        // 1D cell: volume 1
        // 2D cell: volume 1/2
        // 3D cell: volume 1/6
        this->weights.push_back(q[order]->weights[i]*unit_cell_volume[dim]);
    }
}


template class QGauss<0>;
template class QGauss<1>;
template class QGauss<2>;
template class QGauss<3>;
