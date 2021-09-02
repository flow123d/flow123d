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
 * @file    quadrature.hh
 * @brief   Basic definitions of numerical quadrature rules.
 * @author  Jan Stebel
 */

#ifndef QUADRATURE_HH_
#define QUADRATURE_HH_

#include <armadillo>
#include <vector>

#include "system/armor.hh"
#include "mesh/ref_element.hh"



/**
 * @brief Base class for quadrature rules on simplices in arbitrary dimensions.
 *
 * This class stores quadrature points and weights on the reference line,
 * triangle or tetrahedron, respectively.
 * Quadrature rules are used for evaluation of integrals over elements.
 * In particular, for a reference element @f$E@f$ we have:
 * @f[
 * 		\int_E f(x)\,dx \approx \sum_{i=1}^{N} w_i f(p_i),
 * @f]
 * where @f$\{w_i\}@f$, @f$\{p_i\}@f$ are the quadrature weights and the quadrature points,
 * respectively.
 *
 * TODO:
 * - remove set_weight, set point; quadrature should be set by its descendants
 * - introduce Quadrature point which should store point coords together with the weight in raw double[1+dim] array
 *   and just return arma object on the flys
 */
class Quadrature {
public:

    /// Copy constructor.
    Quadrature(const Quadrature &q);
    
    /**
     * @brief Constructor.
     * @param n_quadrature_points Number of quadrature points to be allocated.
     */
    Quadrature(unsigned int dimension, unsigned int n_quadrature_points = 0);

    /** @brief Constructor from quadrature of lower dimension (e.g. for side integration).
     * @param sub_quadrature lower dimensional (dim-1) quadrature
     * @param sid local index of side
     * @param pid index of permutation of nodes on given side
     */
//     template<unsigned int quad_dim>
//     explicit Quadrature(const Quadrature &sub_quadrature, unsigned int sid, unsigned int pid);
    
    /// Virtual destructor.
    virtual ~Quadrature()
    {
    };
    
    inline unsigned int dim() const
    { return dim_; }

    /**
     * @brief Modify the number of quadrature points.
     * @param n_q_points New number of quadrature points.
     */
    inline void resize(unsigned int n_q_points)
    {
        quadrature_points.resize(n_q_points);
        weights.resize(n_q_points, 0);
    }

    /// Returns number of quadrature points.
    inline unsigned int size() const
    { return weights.size(); }

    /// Returns the <tt>i</tt>th quadrature point.
    template<unsigned int point_dim>
    inline Armor::ArmaVec<double, point_dim> point(unsigned int i) const
    {
        ASSERT_EQ_DBG(point_dim, dim_);
        return quadrature_points.vec<point_dim>(i);
    }

    inline Armor::Array<double>::ArrayMatSet set(uint i)
    {
        return quadrature_points.set(i);
    }

    /// Return a reference to the whole array of quadrature points.
    inline const Armor::Array<double> & get_points() const
    { return quadrature_points; }

    /// Returns the <tt>i</tt>th weight.
    inline double weight(unsigned int i) const
    { return weights[i]; }
    
    /// Returns the <tt>i</tt>th weight (non-const version).
    inline double &weight(unsigned int i)
    { return weights[i]; }

    /// Return a reference to the whole array of weights.
    inline const std::vector<double> & get_weights() const
    { return weights; }

    Quadrature &operator=(const Quadrature &q);
    
    /**
     * Create bulk quadrature from side quadrature.
     * 
     * Consider *this as quadrature on a side of an element and create
     * higher dimensional quadrature considering side and permutation index.
     */
    template<unsigned int bulk_dim>
    Quadrature make_from_side(unsigned int sid, unsigned int pid) const;
    

protected:
    
    /// Dimension of quadrature points.
    const unsigned int dim_;
    
    /**
     * @brief List of quadrature points.
     *
     * To be filled by the constructors of the derived classes.
     */
    Armor::Array<double> quadrature_points;

    /**
     * @brief List of weights to the quadrature points.
     *
     * To be filled by the constructors of the derived classes.
     */
    std::vector<double> weights;

};

#endif /* QUADRATURE_HH_ */
