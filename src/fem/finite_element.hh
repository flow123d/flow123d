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
 * @file    finite_element.hh
 * @brief   Abstract class for description of finite elements.
 * @author  Jan Stebel
 */

#ifndef FINITE_ELEMENT_HH_
#define FINITE_ELEMENT_HH_

#include <armadillo>
#include <map>
#include <vector>
#include <boost/assign/list_of.hpp>
#include "fem/update_flags.hh"



template<unsigned int dim, unsigned int spacedim> class FEValuesData;
template<unsigned int dim> class Quadrature;





/**
 * @brief Multiplicity of finite element dofs.
 *
 * Multiplicities describe groups of dofs whose order changes with
 * the configuration (e.g. the rotation or orientation) of
 * the geometrical entity relative to the actual cell.
 *
 * In each spatial dimension we accept the following dof multiplicities:
 *
 * 0) Point (1 configuration):
 *    - single dofs
 *
 * 1) Line (2 possible configurations=orientations):
 *    - single dofs
 *    - pairs
 *
 * 2) Triangle (2 orientations and 3 rotations=6 configurations):
 *    - single dofs
 *    - pairs
 *    - triples
 *    - sextuples
 *
 * 3) Tetrahedron (1 configuration, since it is always the cell):
 *    - single dofs
 */
enum DofMultiplicity {
    DOF_SINGLE = 1, DOF_PAIR = 2, DOF_TRIPLE = 3, DOF_SEXTUPLE = 6
};

const std::vector<DofMultiplicity> dof_multiplicities = boost::assign::list_of(
        DOF_SINGLE)(DOF_PAIR)(DOF_TRIPLE)(DOF_SEXTUPLE);

/**
 * @brief Structure for storing the precomputed finite element data.
 */
class FEInternalData
{
public:
    /**
     * @brief Precomputed values of basis functions at the quadrature points.
     */
    std::vector<arma::vec> basis_values;

    /**
     * @brief Precomputed gradients of basis functions at the quadrature points.
     */
    std::vector<arma::mat > basis_grads;


    /**
     * @brief Precomputed values of basis functions at the quadrature points.
     *
     * For vectorial finite elements.
     */
    std::vector<std::vector<arma::vec> > basis_vectors;

    /**
     * @brief Precomputed gradients of basis functions at the quadrature points.
     *
     * For vectorial finite elements:
     */
    std::vector<std::vector<arma::mat> > basis_grad_vectors;
};


/**
 * @brief Abstract class for the description of a general finite element on
 * a reference simplex in @p dim dimensions.
 *
 * Description of dofs:
 *
 * The reference cell consists of lower dimensional entities (points,
 * lines, triangles). Each dof is associated to one of these
 * entities. This means that if the entity is shared by 2 or more
 * neighbouring cells in the mesh then this dof is shared by the
 * finite elements on all of these cells. If a dof is associated
 * to the cell itself then it is not shared with neighbouring cells.
 * The ordering of nodes in the entity may not be appropriate for the
 * finite elements on the neighbouring cells, hence we need to
 * describe how the order of dofs changes with the relative
 * configuration of the entity with respect to the actual cell.
 * For this reason we define the dof multiplicity which allows to
 * group the dofs as described in \ref DofMultiplicity.
 *
 * Support points:
 *
 * Sometimes it is convenient to describe the function space using
 * a basis (called the raw basis) that is different from the set of
 * shape functions for the finite element (the actual basis). For
 * this reason we define the support points which play the role of
 * nodal functionals associated to the particular dofs. To convert
 * between the two bases one can use the @p node_matrix, which is
 * constructed by the method compute_node_matrix(). In the case of
 * non-Lagrangean finite elements the dofs are not associated to
 * nodal functionals but e.g. to derivatives or averages. For that
 * reason we distinguish between the unit support points which are
 * uniquely associated to the dofs and the generalized support
 * points that are auxiliary for the calculation of the dof
 * functionals.
 *
 *
 */
template<unsigned int dim, unsigned int spacedim>
class FiniteElement {
public:

    /**
     * @brief Constructor.
     */
    FiniteElement();

    /**
     * @brief Clears all internal structures.
     */
    void init();

    /**
     * @brief Returns the number of degrees of freedom needed by the finite
     * element.
     */
    const unsigned int n_dofs() const;

    /**
     * @brief Returns the number of single dofs/dof pairs/triples/sextuples
     * that lie on a single geometric entity of the dimension
     * @p object_dim.
     *
     * @param object_dim Dimension of the geometric entity.
     * @param multiplicity Multiplicity of dofs.
     */
    const unsigned int n_object_dofs(unsigned int object_dim,
            DofMultiplicity multiplicity);

    /**
     * @brief Calculates the value of the @p i-th raw basis function at the
     * point @p p on the reference element.
     *
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    virtual double basis_value(const unsigned int i,
            const arma::vec::fixed<dim> &p) const = 0;

    /**
     * @brief Variant of basis_value() for vectorial finite elements.
     *
     * Calculates the value @p i-th vector-valued raw basis function
     * at the point @p p on the reference element.
     *
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    virtual arma::vec::fixed<dim> basis_vector(const unsigned int i,
            const arma::vec::fixed<dim> &p) const = 0;

    /**
     * @brief Calculates the gradient of the @p i-th raw basis function at the
     * point @p p on the reference element.
     *
     * The gradient components
     * are relative to the reference cell coordinate system.
     *
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    virtual arma::vec::fixed<dim> basis_grad(const unsigned int i,
            const arma::vec::fixed<dim> &p) const = 0;

    /**
     * @brief Variant of basis_grad() for vectorial finite elements.
     *
     * Calculates the gradient of the @p i-th vector-valued raw basis
     * function at the point @p p on the reference element. The gradient
     * components are relative to the reference cell coordinate system.
     *
     * @param i Number of the basis function.
     * @param p Point of evaluation.
     */
    virtual arma::mat::fixed<dim,dim> basis_grad_vector(const unsigned int i,
            const arma::vec::fixed<dim> &p) const = 0;

    /**
     * @brief Initializes the @p node_matrix for computing the coefficients
     * of the raw basis functions from values at support points.
     *
     * The method is implemented for the case of Langrangean finite
     * element. In other cases it may be reimplemented.
     */
    virtual void compute_node_matrix();

    /**
     * @brief Calculates the data on the reference cell.
     *
     * @param q Quadrature rule.
     * @param flags Update flags.
     */
    virtual FEInternalData *initialize(const Quadrature<dim> &q, UpdateFlags flags);

    /**
     * @brief Decides which additional quantities have to be computed
     * for each cell.
     *
     * @param flags Computed update flags.
     */
    virtual UpdateFlags update_each(UpdateFlags flags);

    /**
     * @brief Computes the shape function values and gradients on the actual cell
     * and fills the FEValues structure.
     *
     * @param q Quadrature rule.
     * @param data Precomputed finite element data.
     * @param fv_data Data to be computed.
     */
    virtual void fill_fe_values(const Quadrature<dim> &q,
            FEInternalData &data,
            FEValuesData<dim,spacedim> &fv_data);

    /**
     * @brief Returns the maximum degree of
     * space of polynomials contained in the finite element space.
     *
     * For possible use in hp methods.
     */
    virtual const unsigned int polynomial_order() const {
        return order;
    };

    /**
     * @brief Indicates whether the finite element function space is scalar
     * or vectorial.
     */
    const bool is_scalar() const {
        return is_scalar_fe;
    };

    /**
     * @brief Returns either the generalized support points (if they are defined)
     * or the unit support points.
     */
    const std::vector<arma::vec::fixed<dim> > &get_generalized_support_points();

    /**
     * @brief Destructor.
     */
    virtual ~FiniteElement();

protected:

    /**
     * @brief Total number of degrees of freedom at one finite element.
     */
    unsigned int number_of_dofs;

    /**
     * @brief Number of single dofs at one geometrical entity of the given
     * dimension (point, line, triangle, tetrahedron).
     */
    unsigned int number_of_single_dofs[dim + 1];

    /**
     * @brief Number of pairs of dofs at one geometrical entity of the given
     * dimension (applicable to lines and triangles).
     */
    unsigned int number_of_pairs[dim + 1];

    /**
     * @brief Number of triples of dofs associated to one triangle.
     */
    unsigned int number_of_triples[dim + 1];

    /**
     * @brief Number of sextuples of dofs associated to one triangle.
     */
    unsigned int number_of_sextuples[dim + 1];

    /**
     * @brief Polynomial order - to be possibly used in hp methods.
     */
    unsigned int order;

    /**
     * @brief Indicator of scalar versus vectorial finite element.
     */
    bool is_scalar_fe;

    /**
     * @brief Matrix that determines the coefficients of the raw basis
     * functions from the values at the support points.
     */
    arma::mat node_matrix;

    /**
     * @brief Support points for Lagrangean finite elements.
     *
     * Support points are points in the reference element where
     * function values determine the dofs. In case of Lagrangean
     * finite elements the dof values are precisely the function
     * values at @p unit_support_points.
     *
     */
    std::vector<arma::vec::fixed<dim> > unit_support_points;

    /**
     * @brief Support points for non-Lagrangean finite elements.
     *
     * In case of non-Lagrangean finite elements the meaning of the
     * support points is different, hence we denote the structure
     * as @p generalized_support_points.
     */
    std::vector<arma::vec::fixed<dim> > generalized_support_points;
};




#endif /* FINITE_ELEMENT_HH_ */
