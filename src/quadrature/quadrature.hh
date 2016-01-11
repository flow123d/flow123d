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
 */
template<unsigned int dim>
class Quadrature {
public:
    /**
     * @brief Constructor.
     * @param n_quadrature_points Number of quadrature points to be allocated.
     */
    Quadrature(const unsigned int n_quadrature_points = 0);

    /// Copy constructor.
    Quadrature(const Quadrature<dim> &q);

    /// Virtual destructor.
    virtual ~Quadrature();

    /**
     * @brief Modify the number of quadrature points.
     * @param n_q_points New number of quadrature points.
     */
    void resize(const unsigned int n_q_points);

    /// Returns number of quadrature points.
    const unsigned int size() const;

    /// Returns the <tt>i</tt>th quadrature point.
    const arma::vec::fixed<dim> & point(const unsigned int i) const;

    /// Return a reference to the whole array of quadrature points.
    const std::vector<arma::vec::fixed<dim> > & get_points() const;

    /**
     * @brief Sets individual quadrature point coordinates.
     * @param i Number of the quadrature point.
     * @param p New coordinates.
     */
    void set_point(const unsigned int i, const arma::vec::fixed<dim> &p);

    /// Returns the <tt>i</tt>th weight.
    double weight(const unsigned int i) const;

    /// Return a reference to the whole array of weights.
    const std::vector<double> & get_weights() const;

    /// Sets individual quadrature weight.
    void set_weight(const unsigned int i, const double w);

protected:
    /**
     * @brief List of quadrature points.
     *
     * To be filled by the constructors of the derived classes.
     */
    std::vector<arma::vec::fixed<dim> > quadrature_points;

    /**
     * @brief List of weights to the quadrature points.
     *
     * To be filled by the constructors of the derived classes.
     */
    std::vector<double> weights;

};



template<unsigned int dim>
Quadrature<dim>::Quadrature(const unsigned int n_q)
{
    resize(n_q);
}

template<unsigned int dim>
Quadrature<dim>::Quadrature(const Quadrature<dim> &q) :
        quadrature_points(q.quadrature_points),
        weights(q.weights)
{}

template<unsigned int dim>
void Quadrature<dim>::resize(const unsigned int n_q)
{
    arma::vec::fixed<dim> v;
    v.fill(0);
    quadrature_points.resize(n_q, v);
    weights.resize(n_q, 0);
}

template<unsigned int dim>
inline const unsigned int Quadrature<dim>::size() const {
    return weights.size();
}

template<unsigned int dim>
inline const arma::vec::fixed<dim> & Quadrature<dim>::point(
        const unsigned int i) const {
    return quadrature_points[i];
}

template<unsigned int dim>
inline const std::vector<arma::vec::fixed<dim> > & Quadrature<dim>::get_points() const {
    return quadrature_points;
}

template<unsigned int dim>
inline void Quadrature<dim>::set_point(const unsigned int i, const arma::vec::fixed<dim> &p)
{
    quadrature_points[i] = p;
}

template<unsigned int dim>
inline double Quadrature<dim>::weight(const unsigned int i) const {
    return weights[i];
}

template<unsigned int dim>
inline const std::vector<double> & Quadrature<dim>::get_weights() const {
    return weights;
}

template<unsigned int dim>
inline void Quadrature<dim>::set_weight(const unsigned int i, const double w)
{
    weights[i] = w;
}

template<unsigned int dim>
Quadrature<dim>::~Quadrature()
{}






#endif /* QUADRATURE_HH_ */
