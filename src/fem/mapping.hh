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
 * @brief Class Mapping calculates data related to the mapping
 *        of the reference cell to the actual cell, such as Jacobian
 *        and normal vectors.
 * @author Jan Stebel
 */

#ifndef MAPPING_HH_
#define MAPPING_HH_

#include <armadillo>
#include <vector>
#include "fem/dofhandler.hh"



template<unsigned int dim> class Quadrature;

using namespace arma;
using namespace std;


/**
 * Calculates determinant of a rectangular matrix.
 */
template<unsigned int nr, unsigned int nc>
double determinant(const mat::fixed<nr,nc> &M);

template<unsigned int n> inline double determinant(const mat::fixed<n,n> &M)
{
    return det(M);
}

template<> inline double determinant(const mat::fixed<1,2> &M)
{
    return M(0,0)*M(0,0)+M(1,1)*M(1,1);
}

template<> inline double determinant(const mat::fixed<1,3> &M)
{
    return M(0,0)*M(0,0)+M(0,1)*M(0,1)+M(0,2)*M(0,2);
}

template<> inline double determinant(const mat::fixed<2,3> &M)
{
    return (M(0,0)*M(0,0)+M(0,1)*M(0,1)+M(0,2)*M(0,2))*(M(1,0)*M(1,0)+M(1,1)*M(1,1)+M(1,2)*M(1,2))
          -(M(0,0)*M(1,0)+M(0,1)*M(1,1)+M(0,2)*M(1,2))*(M(0,0)*M(1,0)+M(0,1)*M(1,1)+M(0,2)*M(1,2));
}



/**
 * Class Mapping calculates data related to the mapping of the
 * reference cell to the actual cell, such as Jacobian and normal
 * vectors.
 */
template<unsigned int dim, unsigned int spacedim>
class Mapping
{
public:
    /**
     * Calculates the mapping data and stores them in the provided
     * structures.
     */
    virtual void fill_fe_values(const typename DOFHandler<dim>::CellIterator &cell,
                        const Quadrature<dim> &q,
                        vector< mat::fixed<dim,spacedim> > &jacobians,
                        vector<double> &JxW_values,
                        vector< mat > &inverse_jacobians,
                        vector< vec::fixed<spacedim> > &normal_vectors) = 0;

    virtual ~Mapping() {};
};










#endif /* MAPPING_HH_ */
