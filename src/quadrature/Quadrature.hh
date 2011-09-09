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
 * $Id: time_governor.hh 1323 2011-09-08 08:00:03Z jan.stebel $
 * $Revision: 1323 $
 * $LastChangedBy: jan.stebel $
 * $LastChangedDate: 2011-09-08 10:00:03 +0200 (Čt, 08 zář 2011) $
 *
 * @file
 * @brief Basic definitions of quadratures.
 *  @author Jan Stebel
 */

#ifndef QUADRATURE_HH_
#define QUADRATURE_HH_

#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

/**
 * Base class for quadrature formulae in arbitrary dimensions.
 * This class stores quadrature points and weights on the reference line,
 * triangle, tetrahedron etc.
 */
template<int dim>
class Quadrature {
public:
    /**
     * Constructor.
     */
    Quadrature(const unsigned int n_quadrature_points = 0);

    /**
     * Copy constructor.
     */
    Quadrature(const Quadrature<dim> &q);

    /**
     * Virtual destructor.
     */
    virtual ~Quadrature();

    /**
     * Number of quadrature points.
     */
    unsigned int size() const;

    /**
     * Return the <tt>i</tt>th quadrature point.
     */
    const vec::fixed<dim> & point(const unsigned int i) const;

    /**
     * Return a reference to the whole array of quadrature points.
     */
    const vector<vec::fixed<dim> > & get_points() const;

    /**
     * Return the <tt>i</tt>th weight.
     */
    double weight(const unsigned int i) const;

    /**
     * Return a reference to the whole array of weights.
     */
    const vector<double> & get_weights() const;

protected:
    /**
     * List of quadrature points.
     * To be filled by the constructors of the derived classes.
     */
    vector<vec::fixed<dim> > quadrature_points;

    /**
     * List of weights to the quadrature points.
     * To be filled by the constructors of the derived classes.
     */
    vector<double> weights;
};

template<int dim>
inline unsigned int Quadrature<dim>::size() const {
    return weights.size();
}

template<int dim>
inline const vec::fixed<dim> & Quadrature<dim>::point(
        const unsigned int i) const {
    return quadrature_points[i];
}

template<int dim>
inline const vector<vec::fixed<dim> > & Quadrature<dim>::get_points() const {
    return quadrature_points;
}

template<int dim>
inline double Quadrature<dim>::weight(const unsigned int i) const {
    return weights[i];
}

template<int dim>
inline const vector<double> & Quadrature<dim>::get_weights() const {
    return weights;
}

#endif /* QUADRATURE_HH_ */
