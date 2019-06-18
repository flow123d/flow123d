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

    /** @brief Constructor from quadrature of lower dimension (e.g. for side integration).
     * @param sub_quadrature lower dimnesional (dim-1) quadrature
     * @param sid local index of side
     * @param pid index of permutation of nodes on given side
     */
    Quadrature(const Quadrature<dim-1> &sub_quadrature, unsigned int sid, unsigned int pid);
    
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

template<unsigned int dim> inline
Quadrature<dim>::Quadrature(const Quadrature<dim-1> &subq, unsigned int sid, unsigned int pid)
{
    resize(subq.size());

//     double lambda;
// 
//     // vectors of barycentric coordinates of quadrature points
//     arma::vec::fixed<dim+1> el_bar_coords;
//     arma::vec::fixed<dim> side_bar_coords;
// 
//     for (unsigned int k=0; k<subq.size(); k++)
//     {
//         const arma::vec::fixed<dim-1> &sub_point = subq.point(k);
//         // Calculate barycentric coordinates on the side of the k-th
//         // quadrature point.
//         el_bar_coords.zeros();
//         lambda = 0;
//         // Apply somewhere permutation of indices!
//         for (unsigned int j=0; j<dim-1; j++)
//         {
//             side_bar_coords(j) = sub_point(j);
//             lambda += sub_point(j);
//         }
//         side_bar_coords(dim-1) = 1.0 - lambda;
// 
//         // transform to element coordinates
//         auto side_nodes = RefElement<dim>::interact(Interaction<0, (dim - 1)>(sid));
//         for (unsigned int i=0; i<dim; i++) {
//             // TODO: use RefElement<>::interpolate to map coordinates from the subelement
//             unsigned int i_node = (side_nodes[RefElement<dim>::side_permutations[pid][i]]+dim)%(dim+1);
//             el_bar_coords(i_node) = side_bar_coords((i+dim-1)%dim);
//         }
//         quadrature_points[k] = el_bar_coords.subvec(0,dim-1);
//         weights[k] = subq.weight(k);
//     }
    
    arma::vec::fixed<dim+1> el_bar_coords, final_bar;
    
    for (unsigned int k=0; k<subq.size(); k++)
    {
        //compute barycentric coordinates on element
        arma::vec::fixed<dim> p = RefElement<dim-1>::local_to_bary(subq.point(k));
        arma::vec::fixed<dim> pp;
        
        //permute
        for (unsigned int i=0; i<RefElement<dim>::n_nodes_per_side; i++) {
            pp(RefElement<dim>::side_permutations[pid][i]) = p(i);
        }
        
        el_bar_coords = RefElement<dim>::template interpolate<dim-1>(pp,sid);
        
        //get local coordinates and set
        quadrature_points[k] = RefElement<dim>::bary_to_local(el_bar_coords);
        weights[k] = subq.weight(k);
    }
}

// Specialized subquadrature consructor for dim=1.
template<> inline
Quadrature<1>::Quadrature(const Quadrature<0> &subq, unsigned int sid, unsigned int pid)
{
    arma::vec::fixed<1> p({(double)sid});
    quadrature_points.push_back(p);
    weights.push_back(1);
}

#endif /* QUADRATURE_HH_ */
