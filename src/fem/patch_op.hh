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
 * @file    patch_op.hh
 * @brief   Base class of FE operations.
 * @author  David Flanderka
 */

#ifndef PATCH_OP_HH_
#define PATCH_OP_HH_

#include <Eigen/Dense>

#include "fem/eigen_tools.hh"
#include "fem/arena_resource.hh"
#include "fem/arena_vec.hh"
#include "mesh/ref_element.hh"


template<unsigned int spacedim> class PatchFEValues;


/// Distinguishes bulk and side domain
enum ElemDomain
{
    domain_bulk,
    domain_side
};


/**
 * @brief Class represents element or FE operations.
 */
template<unsigned int spacedim = 3>
class PatchOp {
public:
    /**
     * Constructor
     *
     * Set all data members.
     */
    PatchOp(uint dim, PatchFEValues<spacedim> &pfev, std::initializer_list<uint> shape, OpSizeType size_type, uint n_dofs = 1)
    : dim_(dim), shape_(set_shape_vec(shape)), size_type_(size_type), n_dofs_(n_dofs), patch_fe_(&pfev)
    {
        ASSERT_GT(n_dofs, 0);
        result_ = Eigen::Vector<ArenaVec<double>, Eigen::Dynamic>(shape_[0] * shape_[1] * n_dofs_);
    }

    /// Destructor
    virtual ~PatchOp()
    {}


    /// Aligns shape_vec to 2 items (equal to matrix number of dimensions)
    std::vector<uint> set_shape_vec(std::initializer_list<uint> shape) const {
        std::vector<uint> shape_vec(shape);
        if (shape_vec.size() == 1) shape_vec.push_back(1);
        ASSERT_EQ(shape_vec.size(), 2);
        return shape_vec;
    }

    /**
     * Return number of operation components
     *
     * Value is computed from shape_ vector
     */
    inline uint n_comp() const {
        return shape_[0] * shape_[1];
    }

//    /// Getter for dimension
//    inline uint dimension() const { // dim is in conflict with template of some descendants
//        return dim_;
//    }

    /// Getter for bulk_side flag
    inline ElemDomain domain() const {
        return domain_;
    }

    /// Getter for size_type_
    OpSizeType size_type() const {
        return size_type_;
    }

    /// Getter for n_dofs_
    inline uint n_dofs() const {
        return n_dofs_;
    }

    /// Return pointer to operation of i_op index in input operation vector.
    inline PatchOp<spacedim> *input_ops(uint i_op) const {
        return input_ops_[i_op];
    }

    /// Getter for shape_
    inline const std::vector<uint> &shape() const {
        return shape_;
    }

    /**
     * Format shape to string
     *
     * Method is used in output development method.
     */
    std::string format_shape() const {
        std::stringstream ss;
        ss << shape_[0] << "x" << shape_[1];
        return ss.str();
    }

    inline void allocate_result(size_t data_size, PatchArena &arena) {
        for (uint i=0; i<n_comp()*n_dofs_; ++i)
            result_(i) = ArenaVec<double>(data_size, arena);
    }

    inline void allocate_const_result(ArenaVec<double> &value_vec) {
        for (uint i=0; i<n_comp()*n_dofs_; ++i)
            result_(i) = value_vec;
    }

    /// Return map referenced result as Eigen::Vector
    inline Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> result_matrix() {
        return Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>>(result_.data(), shape_[0], shape_[1] * n_dofs_);
    }

    /// Return map referenced result of DOF values as Eigen::Matrix
    inline Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> result_sub_matrix(uint i_dof) {
        ASSERT_LT(i_dof, n_dofs_);
        uint n_dof_comps = shape_[0] * shape_[1];
        return Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>>(result_.data()+i_dof*n_dof_comps, shape_[0], shape_[1]);
    }

    /// Return map referenced result as Eigen::Vector
    inline Eigen::Vector<ArenaVec<double>, Eigen::Dynamic> &raw_result() {
        return result_;
    }

    /// Same as previous but return const reference
    inline const Eigen::Vector<ArenaVec<double>, Eigen::Dynamic> &raw_result() const {
        return result_;
    }

    /// Return reference of PatchPointValues
    inline PatchPointValues<spacedim> &ppv() {
        return patch_fe_->patch_point_vals_[domain_][this->dim_-1];
    }

    /// Reinit function of operation. Implementation in descendants.
    virtual void eval() =0;

    /**
     * Returns output value of data stored by elements.
     *
     * @param point_idx   Index of quadrature point in ElementCacheMap
     */
    template <class ValueType>
    ValueType elem_value(uint point_idx) const;

    /**
     * Returns output value on quadrature point.
     *
     * @param point_idx   Index of quadrature point in ElementCacheMap
     * @param i_dof       Index of DOF
     */
    template <class ValueType>
    ValueType point_value(uint point_idx, uint i_dof=0) const;


protected:
    /**
     * @brief Precomputed values of basis functions at the bulk quadrature points.
     *
     * Dimensions:   (no. of quadrature points)
     *             x (no. of dofs)
     *             x (no. of components in ref. cell)
     */
    template<unsigned int FE_dim>
    std::vector<std::vector<arma::vec> > ref_shape_values_bulk(Quadrature *q, std::shared_ptr<FiniteElement<FE_dim>> fe) {
        std::vector<std::vector<arma::vec> > ref_shape_vals( q->size(), vector<arma::vec>(fe->n_dofs()) );

        arma::mat shape_values(fe->n_dofs(), fe->n_components());
        for (unsigned int i=0; i<q->size(); i++)
        {
            for (unsigned int j=0; j<fe->n_dofs(); j++)
            {
                for (unsigned int c=0; c<fe->n_components(); c++)
                    shape_values(j,c) = fe->shape_value(j, q->point<FE_dim>(i), c);

                ref_shape_vals[i][j] = trans(shape_values.row(j));
            }
        }

        return ref_shape_vals;
    }

    /**
     * @brief Precomputed values of basis functions at the side quadrature points.
     *
     * Dimensions:   (sides)
     *             x (no. of quadrature points)
     *             x (no. of dofs)
     *             x (no. of components in ref. cell)
     */
    template<unsigned int FE_dim>
    std::vector< std::vector<std::vector<arma::vec> > > ref_shape_values_side(Quadrature *q, std::shared_ptr<FiniteElement<FE_dim>> fe) {
        std::vector< std::vector<std::vector<arma::vec> > > ref_shape_vals( FE_dim+1, std::vector<std::vector<arma::vec> >(q->size(), vector<arma::vec>(fe->n_dofs())) );

        arma::mat shape_values(fe->n_dofs(), fe->n_components());

        for (unsigned int sid=0; sid<FE_dim+1; sid++) {
            auto quad = q->make_from_side<FE_dim>(sid);
        	for (unsigned int i=0; i<quad.size(); i++)
            {
                for (unsigned int j=0; j<fe->n_dofs(); j++)
                {
                    for (unsigned int c=0; c<fe->n_components(); c++) {
                        shape_values(j,c) = fe->shape_value(j, quad.template point<FE_dim>(i), c);
                    }

                    ref_shape_vals[sid][i][j] = trans(shape_values.row(j));
                }
            }
        }

        return ref_shape_vals;
    }

    /**
     * @brief Precomputed gradients of basis functions at the bulk quadrature points.
     *
     * Dimensions:   (no. of quadrature points)
     *             x (no. of dofs)
     *             x ((dim_ of. ref. cell)x(no. of components in ref. cell))
     */
    template<unsigned int FE_dim>
    std::vector<std::vector<arma::mat> > ref_shape_gradients_bulk(Quadrature *q, std::shared_ptr<FiniteElement<FE_dim>> fe) {
    	std::vector<std::vector<arma::mat> > ref_shape_grads( q->size(), vector<arma::mat>(fe->n_dofs()) );

        arma::mat grad(FE_dim, fe->n_components());
        for (unsigned int i_pt=0; i_pt<q->size(); i_pt++)
        {
            for (unsigned int i_dof=0; i_dof<fe->n_dofs(); i_dof++)
            {
                grad.zeros();
                for (unsigned int c=0; c<fe->n_components(); c++)
                    grad.col(c) += fe->shape_grad(i_dof, q->point<FE_dim>(i_pt), c);

                ref_shape_grads[i_pt][i_dof] = grad;
            }
        }

        return ref_shape_grads;
    }

    /**
     * @brief Precomputed gradients of basis functions at the side quadrature points.
     *
     * Dimensions:   (sides)
     *             x (no. of quadrature points)
     *             x (no. of dofs)
     *             x ((dim_ of. ref. cell)x(no. of components in ref. cell))
     */
    template<unsigned int FE_dim>
    std::vector<std::vector<std::vector<arma::mat> > > ref_shape_gradients_side(Quadrature *q, std::shared_ptr<FiniteElement<FE_dim>> fe) {
        std::vector<std::vector<std::vector<arma::mat> > > ref_shape_grads( FE_dim+1, std::vector<std::vector<arma::mat> >(q->size(), vector<arma::mat>(fe->n_dofs())) );

        arma::mat grad(FE_dim, fe->n_components());
        for (unsigned int sid=0; sid<FE_dim+1; sid++) {
            auto quad = q->make_from_side<FE_dim>(sid);
            for (unsigned int i_pt=0; i_pt<quad.size(); i_pt++)
            {
                for (unsigned int i_dof=0; i_dof<fe->n_dofs(); i_dof++)
                {
                    grad.zeros();
                    for (unsigned int c=0; c<fe->n_components(); c++)
                        grad.col(c) += fe->shape_grad(i_dof, quad.template point<FE_dim>(i_pt), c);

                    ref_shape_grads[sid][i_pt][i_dof] = grad;
                }
            }
        }

        return ref_shape_grads;
    }

    uint dim_;                                    ///< Dimension
    ElemDomain domain_;                           ///< Flag: BulkOp = 0, SideOp = 1
    std::vector<uint> shape_;                     ///< Shape of stored data (size of vector or number of rows and cols of matrix)
    Eigen::Vector<ArenaVec<double>, Eigen::Dynamic> result_;    ///< Result matrix of operation
    OpSizeType size_type_;                        ///< Type of operation by size of vector (element, point or fixed size)
    std::vector<PatchOp<spacedim> *> input_ops_;  ///< Indices of operations in PatchPointValues::operations_ vector on which PatchOp is depended
    uint n_dofs_;                                 ///< Number of DOFs of FE operations (or 1 in case of element operations)
    PatchFEValues<spacedim> *patch_fe_;           ///< Pointer to PatchFEValues object

    friend class PatchFEValues<spacedim>;
};


#endif /* PATCH_OP_HH_ */
