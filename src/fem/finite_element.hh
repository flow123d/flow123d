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
#include <stdio.h>                             // for sprintf
#include <string.h>                            // for memcpy
#include <algorithm>                           // for max, min
#include <boost/assign/list_of.hpp>            // for generic_list, list_of
#include <boost/exception/info.hpp>            // for error_info::error_info...
#include <new>                                 // for operator new[]
#include <ostream>                             // for operator<<
#include <string>                              // for operator<<
#include <vector>                              // for vector
#include "fem/update_flags.hh"                 // for operator&, operator|=
#include "system/exceptions.hh"                // for ExcAssertMsg::~ExcAsse...

template<unsigned int dim> class FESystem;
template<unsigned int spacedim> class FEValues;
template<unsigned int dim> class FE_P_disc;





// Possible types are: value, gradient, cell integral, ...
enum DofType { Value = 1 };

/**
 * Class Dof is a general description of functionals (degrees of freedom)
 * determining the finite element. We assume that the Dof is defined as
 * a linear combination of components of function value at a given point:
 * 
 *     Dof_value = a_1.f_1(x) + ... + a_n.f_n(x),
 * 
 * where (a_1, ... , a_n) are given by @p coefs, x is the support point
 * given by barycentric @p coords and (f_1, ..., f_n) is a generally
 * vector-valued function. For the simplest Dof, i.e. value of a scalar
 * function at point p, we set
 * 
 *    @p coords = p, @p coefs = { 1 }.
 * The member @p dim denotes the affiliation of Dof to n-face:
 *    Nodal dofs have      dim = 0,
 *    Dofs on lines:       dim = 1,
 *    Dofs on triangles:   dim = 2,
 *    Dofs in tetrahedron: dim = 3.
 * It means that when a node, line or triangle is shared by adjacent cells,
 * also the Dofs on this n-face are shared by the cells. Therefore
 * for DG finite elements we set for all dofs the highest possible dimension.
 * 
 * The class implements the method evaluate() which computes the Dof value
 * for a basis function from given FunctionSpace.
 */
class Dof {
public:
    
    Dof(unsigned int dim_,
        unsigned int n_face_idx_,
        arma::vec coords_,
        arma::vec coefs_,
        DofType type_)
    
        : dim(dim_),
          n_face_idx(n_face_idx_),
          coords(coords_),
          coefs(coefs_),
          type(type_)
    {}
    
    /// Evaulate dof for basis function of given function space.
    template<class FS>
    double evaluate(const FS &function_space, 
                    unsigned int basis_idx) const;
    
    /// Association to n-face of given dimension (point, line, triangle, tetrahedron.
    unsigned int dim;
    
    /// Index of n-face to which the dof is associated.
    unsigned int n_face_idx;
    
    /// Barycentric coordinates.
    arma::vec coords;
    
    /// Coefficients of linear combination of function value components.
    arma::vec coefs;
    
    /** 
     * Currently we consider only type=Value, possibly we could have Gradient,
     * CellIntegral, FaceIntegral or other types.
     */
    DofType type;
};


/**
 * FunctionSpace is an abstract class that is used to describe finite elements.
 * It is determined by the dimension of the field of definition (@p space_dim_),
 * by the dimension of the range (@p n_components_) and by the values and
 * gradients of a basis functions.
 * 
 * Combining FunctionSpace with Dof(s), the FiniteElement class constructs the
 * shape functions, i.e. basis of FunctionSpace for which the Dof(s) attain
 * the value 0 or 1.
 */
class FunctionSpace {
public:
    
    FunctionSpace(unsigned int space_dim,
                  unsigned int n_components)
    
        : space_dim_(space_dim),
          n_components_(n_components)
    {}
    
    
    /**
     * @brief Value of the @p i th basis function at point @p point.
     * @param basis_index  Index of the basis function.
     * @param point        Point coordinates.
     * @param comp_index   Index of component (>0 for vector-valued functions).
     */
    virtual double basis_value(unsigned int basis_index,
                               const arma::vec &point,
                               unsigned int comp_index = 0
                               ) const = 0;
    
    /**
     * @brief Gradient of the @p i th basis function at point @p point.
     * @param basis_index  Index of the basis function.
     * @param point        Point coordinates.
     * @param comp_index   Index of component (>0 for vector-valued functions).
     */
    virtual const arma::vec basis_grad(unsigned int basis_index,
                                       const arma::vec &point,
                                       unsigned int comp_index = 0
                                      ) const = 0;
    
    /// Dimension of function space (number of basis functions).
    virtual unsigned int dim() const = 0;
    
    /// Getter for space dimension.
    unsigned int space_dim() const { return space_dim_; }
    
    /// Getter for number of components.
    unsigned int n_components() const { return n_components_; }
    
    virtual ~FunctionSpace() {}
    
protected:
    
    /// Space dimension of function arguments (i.e. 1, 2 or 3).
    unsigned int space_dim_;
    
    /// Number of components of function values.
    unsigned int n_components_;
};


/**
 * Types of FiniteElement: scalar, vector-valued, tensor-valued or mixed system.
 * 
 * The type also indicates how the shape functions and their gradients are transformed
 * from reference element to arbitrary element. In particular:
 * 
 *     TYPE              OBJECT  EXPRESSION
 * -----------------------------------------------------------
 * FEScalar, FEVector,   value   ref_value
 * FETensor              grad    J^{-T} * ref_grad
 * -----------------------------------------------------------
 * FEVectorContravariant value   J * ref_value
 *                       grad    J^{-T} * ref_grad * J^T
 * -----------------------------------------------------------
 * FEVectorPiola         value   J * ref_value / |J|
 *                       grad    J^{-T} * ref_grad * J^T / |J|
 * -----------------------------------------------------------
 * FEMixedSystem         value   depends on sub-elements
 *                       grad    depends on sub-elements
 * 
 * The transformation itself is done in FEValuesBase::fill_..._data() methods.
 * 
 * Note that we use columnwise gradients, i.e. gradient of each component is a column vector.
 */
enum FEType {
  FEScalar = 0,
  FEVector = 1,
  FEVectorContravariant = 2,
  FEVectorPiola = 3,
  FETensor = 4,
  FEMixedSystem = 5
};



/**
 * @brief Abstract class for the description of a general finite element on
 * a reference simplex in @p dim dimensions.
 * 
 * The finite element is determined by a @p function_space_ and a set
 * of @p dofs_. Further it must be set whether the finite element
 * @p is_primitive_, which means that even if the functions in
 * @p function_space_ are vector-valued, the basis functions have
 * only one nonzero component (this is typical for tensor-product FE,
 * e.g. vector-valued polynomial FE, but not for Raviart-Thomas FE).
 *
 * Description of dofs:
 *
 * The reference cell consists of lower dimensional entities (points,
 * lines, triangles). Each dof is associated to one of these
 * entities. If the entity is shared by 2 or more neighbouring cells
 * in the mesh then this dof is shared by the finite elements on all
 * of these cells. If a dof is associated to the cell itself then it
 * is not shared with neighbouring cells.
 * 
 * 
 * Shape functions:
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
template<unsigned int dim>
class FiniteElement {
public:
  
    /**
     * @brief Constructor.
     */
    FiniteElement();
    
    /**
     * @brief Returns the number of degrees of freedom needed by the finite
     * element.
     */
    inline unsigned int n_dofs() const
    { return dofs_.size(); }

    /**
     * @brief Calculates the value of the @p comp-th component of
     * the @p i-th shape function at the
     * point @p p on the reference element.
     *
     * @param i    Number of the shape function.
     * @param p    Point of evaluation.
     * @param comp Number of vector component.
     */
    double shape_value(const unsigned int i,
                       const arma::vec::fixed<dim> &p,
                       const unsigned int comp = 0) const;

    /**
     * @brief Calculates the @p comp-th component of the gradient
     * of the @p i-th shape function at the point @p p on the
     * reference element.
     *
     * @param i    Number of the shape function.
     * @param p    Point of evaluation.
     * @param comp Number of vector component.
     */
    arma::vec::fixed<dim> shape_grad(const unsigned int i,
                                     const arma::vec::fixed<dim> &p,
                                     const unsigned int comp = 0) const;

    /// Returns numer of components of the basis function.    
    inline unsigned int n_components() const
    { return function_space_->n_components(); }
    
    /// Returns @p i -th degree of freedom.
    inline const Dof &dof(unsigned int i) const
    { return dofs_[i]; }
    
    /// Number of components of FE in a mapped space with dimension @p spacedim.
    unsigned int n_space_components(unsigned int spacedim);
    
    /// Get barycentric coordinates of the points on the reference element associated with the dofs.
    /// Used in BDDC for unknown reason.
    virtual std::vector< arma::vec::fixed<dim+1> > dof_points() const;

    /**
     * @brief Destructor.
     */
    virtual ~FiniteElement() {};

protected:
  
    /**
     * @brief Clears all internal structures.
     */
    void init(bool primitive = true,
              FEType type = FEScalar);
    
    /**
     * @brief Initialize vectors with information about components of basis functions.
     */
    void setup_components();
    
    /**
     * @brief Decides which additional quantities have to be computed
     * for each cell.
     *
     * @param flags Computed update flags.
     */
    virtual UpdateFlags update_each(UpdateFlags flags);

    /**
     * @brief Initializes the @p node_matrix for computing the coefficients
     * of the shape functions in the raw basis of @p functions_space_.
     * This is done by evaluating the @p dofs_ for the basis function
     * and by inverting the resulting matrix.
     */
    virtual void compute_node_matrix();
    
    /**
     * @brief Indicates whether the basis functions have one or more
     * nonzero components (scalar FE spaces are always primitive).
     */
    inline bool is_primitive() const
    { return is_primitive_; }
    
    /**
     * @brief Returns the component index for vector valued finite elements.
     * @param sys_idx Index of shape function.
     */
    inline unsigned int system_to_component_index(unsigned sys_idx) const
    { return component_indices_[sys_idx]; }
    
    /**
     * @brief Returns the mask of nonzero components for given basis function.
     * @param sys_idx Index of basis function.
     */
    inline const std::vector<bool> &get_nonzero_components(unsigned int sys_idx) const
    { return nonzero_components_[sys_idx]; }

    
    
    /// Type of FiniteElement.
    FEType type_;

    /**
     * @brief Primitive FE is using componentwise shape functions,
     * i.e. only one component is nonzero for each shape function.
     */
    bool is_primitive_;
    
    /// Indices of nonzero components of shape functions (for primitive FE).
    std::vector<unsigned int> component_indices_;
    
    /// Footprints of nonzero components of shape functions.
    std::vector<std::vector<bool> > nonzero_components_;

    /**
     * @brief Matrix that determines the coefficients of the raw basis
     * functions from the values at the support points.
     */
    arma::mat node_matrix;

    /// Function space defining the FE.
    std::shared_ptr<FunctionSpace> function_space_;
    
    /// Set of degrees of freedom (functionals) defining the FE.
    std::vector<Dof> dofs_;
    
    
    friend class FESystem<dim>;
    friend class FEValues<3>;
    friend class FE_P_disc<dim>;
    friend class SubDOFHandlerMultiDim;
};




#endif /* FINITE_ELEMENT_HH_ */
