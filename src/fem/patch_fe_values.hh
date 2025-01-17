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
 * @file    patch_fe_values.hh
 * @brief   Class FEValues calculates finite element data on the actual
 *          cells such as shape function values, gradients, Jacobian of
 *          the mapping from the reference cell etc.
 * @author  Jan Stebel, David Flanderka
 */

#ifndef PATCH_FE_VALUES_HH_
#define PATCH_FE_VALUES_HH_


#include <string.h>                           // for memcpy
#include <algorithm>                          // for swap
#include <new>                                // for operator new[]
#include <string>                             // for operator<<
#include <vector>                             // for vector
#include "fem/fe_system.hh"                   // for FESystem
#include "fem/eigen_tools.hh"
#include "fem/patch_point_values.hh"
#include "fem/mapping_p1.hh"
#include "mesh/ref_element.hh"                // for RefElement
#include "mesh/accessors.hh"
#include "fem/update_flags.hh"                // for UpdateFlags
#include "quadrature/quadrature_lib.hh"
#include "fields/eval_subset.hh"
#include "fem/arena_resource.hh"
#include "fem/arena_vec.hh"

template<unsigned int spacedim> class PatchFEValues;




template <class ValueType>
class ElQ {
public:
    /// Forbidden default constructor
    ElQ() = delete;

    /// Constructor
    ElQ(PatchPointValues<3> *patch_point_vals, unsigned int op_idx)
    : patch_point_vals_(patch_point_vals), op_idx_(op_idx) {}

    ValueType operator()(const BulkPoint &point) const;

    ValueType operator()(const SidePoint &point) const;

private:
    PatchPointValues<3> *patch_point_vals_; ///< Reference to PatchPointValues
    unsigned int op_idx_;                   ///< Index of operation in patch_point_vals_.operations vector
};


template <class ValueType>
class FeQ {
public:
    /// Forbidden default constructor
    FeQ() = delete;

    // Class similar to current FeView
    FeQ(PatchPointValues<3> *patch_point_vals, bool is_bulk, unsigned int op_idx)
    : patch_point_vals_bulk_(nullptr), patch_point_vals_side_(nullptr), op_idx_(op_idx), i_shape_fn_idx_(0) {
        if (is_bulk) patch_point_vals_bulk_ = patch_point_vals;
        else patch_point_vals_side_ = patch_point_vals;
    }

    /// Constructor used only in FeQArray::shape()
    FeQ(PatchPointValues<3> *patch_point_vals_bulk, PatchPointValues<3> *patch_point_vals_side,
            unsigned int op_idx, unsigned int i_shape_fn_idx)
    : patch_point_vals_bulk_(patch_point_vals_bulk), patch_point_vals_side_(patch_point_vals_side),
      op_idx_(op_idx), i_shape_fn_idx_(i_shape_fn_idx) {}


    ValueType operator()(const BulkPoint &point) const;

    ValueType operator()(const SidePoint &point) const;

    // Implementation for EdgePoint, SidePoint, and JoinPoint shoud have a common implementation
    // resolving to side values

private:
    PatchPointValues<3> *patch_point_vals_bulk_; ///< Pointer to bulk PatchPointValues
    PatchPointValues<3> *patch_point_vals_side_; ///< Pointer to side PatchPointValues
    unsigned int op_idx_;                        ///< Index of operation in patch_point_vals_.operations vector
    unsigned int i_shape_fn_idx_;                ///< Index of shape function
};


template <class ValueType>
class FeQArray {
public:
    /// Forbidden default constructor
    FeQArray() = delete;

    // Class similar to current FeView
    FeQArray(PatchPointValues<3> *patch_point_vals, bool is_bulk, unsigned int op_idx, unsigned int n_dofs)
    : patch_point_vals_bulk_(nullptr), patch_point_vals_side_(nullptr), op_idx_(op_idx), n_dofs_(n_dofs) {
        ASSERT_GT(n_dofs, 0).error("Invalid number of DOFs.\n");

        if (is_bulk) patch_point_vals_bulk_ = patch_point_vals;
        else patch_point_vals_side_ = patch_point_vals;
    }


    FeQ<ValueType> shape(unsigned int i_shape_fn_idx) const {
        ASSERT_LT(i_shape_fn_idx, n_dofs_);
        return FeQ<ValueType>(patch_point_vals_bulk_, patch_point_vals_side_, op_idx_, i_shape_fn_idx);
    }

    /// Return number of DOFs
    inline unsigned int n_dofs() const {
        return n_dofs_;
    }

private:
    PatchPointValues<3> *patch_point_vals_bulk_; ///< Reference to bulk PatchPointValues
    PatchPointValues<3> *patch_point_vals_side_; ///< Reference to side PatchPointValues
    unsigned int op_idx_;                        ///< Index of operation in patch_point_vals_.operations vector
    unsigned int n_dofs_;                        ///< Number of DOFs
};


template <class ValueType>
class FeQJoin {
public:
    /// Default constructor
    FeQJoin()
    : patch_point_vals_bulk_(nullptr), patch_point_vals_side_(nullptr) {}

    /**
     * Constructor
     *
     * @param patch_point_vals_bulk  Pointer to PatchPointValues bulk object.
     * @param patch_point_vals_side  Pointer to PatchPointValues side object.
     * @param begin                  Index of the first component of the bulk Quantity.
     * @param begin_side             Index of the first component of the side Quantity.
     * @param n_dofs_bulk            Number of DOFs of bulk (lower-dim) element.
     * @param n_dofs_side            Number of DOFs of side (higher-dim) element.
     */
    FeQJoin(PatchPointValues<3> *patch_point_vals_bulk, PatchPointValues<3> *patch_point_vals_side, unsigned int n_dofs_bulk,
            unsigned int n_dofs_side, unsigned int op_idx_bulk, unsigned int op_idx_side)
    : patch_point_vals_bulk_(patch_point_vals_bulk), patch_point_vals_side_(patch_point_vals_side),
      n_dofs_high_(n_dofs_side), n_dofs_low_(n_dofs_bulk), op_idx_bulk_(op_idx_bulk), op_idx_side_(op_idx_side) {}


    inline unsigned int n_dofs_low() const {
        return n_dofs_low_;
    }

    inline unsigned int n_dofs_high() const {
        return n_dofs_high_;
    }

    inline unsigned int n_dofs_both() const {
        return n_dofs_high_ + n_dofs_low_;
    }

    /// Return local index of DOF (on low / high-dim) - should be private method
    inline unsigned int local_idx(unsigned int i_join_idx) const {
        if (this->is_high_dim(i_join_idx)) return (i_join_idx - n_dofs_low());
        else return i_join_idx;
    }

    inline bool is_high_dim(unsigned int i_join_idx) const {
        return (i_join_idx >= n_dofs_low());
    }

    FeQ<ValueType> shape(unsigned int i_join_idx) const {
        ASSERT_LT(i_join_idx, n_dofs_both());

        if (this->is_high_dim(i_join_idx))
            return FeQ<ValueType>(patch_point_vals_side_->zero_values(), patch_point_vals_side_, op_idx_side_, i_join_idx - n_dofs_low());
        else
            return FeQ<ValueType>(patch_point_vals_bulk_, patch_point_vals_bulk_->zero_values(), op_idx_bulk_, i_join_idx);
    }


private:
    // attributes:
    PatchPointValues<3> *patch_point_vals_bulk_;  ///< Pointer to bulk PatchPointValues
    PatchPointValues<3> *patch_point_vals_side_;  ///< Pointer to side PatchPointValues
    unsigned int n_dofs_high_;                    ///< Number of DOFs on high-dim element
    unsigned int n_dofs_low_;                     ///< Number of DOFs on low-dim element
    unsigned int op_idx_bulk_;                    ///< Index of operation in patch_point_vals_bulk_.operations vector
    unsigned int op_idx_side_;                    ///< Index of operation in patch_point_vals_side_.operations vector
};



template<unsigned int dim>
class BaseValues
{
protected:
	// Default constructor
	BaseValues()
	{}

    /// Return FiniteElement of \p component_idx for FESystem or \p fe for other types
    template<unsigned int FE_dim>
    std::shared_ptr<FiniteElement<FE_dim>> fe_comp(std::shared_ptr< FiniteElement<FE_dim> > fe, uint component_idx) {
        if (fe->fe_type() == FEMixedSystem) {
            FESystem<FE_dim> *fe_sys = dynamic_cast<FESystem<FE_dim>*>( fe.get() );
            return fe_sys->fe()[component_idx];
        } else {
            ASSERT_EQ(component_idx, 0).warning("Non-zero component_idx can only be used for FESystem.");
            return fe;
        }
    }

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

        arma::mat grad(dim, fe->n_components());
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

};

template<unsigned int dim>
class BulkValues : public BaseValues<dim>
{
public:
	/// Constructor
	BulkValues(PatchPointValues<3> *patch_point_vals, MixedPtr<FiniteElement> fe)
	: BaseValues<dim>(), patch_point_vals_(patch_point_vals) {
	    ASSERT_EQ(patch_point_vals->dim(), dim);
	    fe_ = fe[Dim<dim>{}];
	}

    /**
     * @brief Register the product of Jacobian determinant and the quadrature
     * weight at bulk quadrature points.
     *
     * @param quad Quadrature.
     */
    inline FeQ<Scalar> JxW()
    {
        return FeQ<Scalar>(patch_point_vals_, true, FeBulk::BulkOps::opJxW);
    }

	/// Create bulk accessor of coords entity
    inline ElQ<Vector> coords()
    {
        return ElQ<Vector>(patch_point_vals_, FeBulk::BulkOps::opCoords);
    }

//    inline ElQ<Tensor> jacobian(std::initializer_list<Quadrature *> quad_list)
//    {}

    /// Create bulk accessor of jac determinant entity
    inline ElQ<Scalar> determinant()
    {
        return ElQ<Scalar>(patch_point_vals_, FeBulk::BulkOps::opJacDet);
    }

    inline FeQArray<Scalar> ref_scalar(uint component_idx = 0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");

        uint n_points = patch_point_vals_->get_quadrature()->size();
        uint n_dofs = fe_component->n_dofs();

        auto *ref_scalar_op = patch_point_vals_->make_fixed_fe_op(FeBulk::BulkOps::opRefScalar, {n_dofs}, &common_reinit::op_base, n_dofs);

        auto ref_shape_vals = this->ref_shape_values_bulk(patch_point_vals_->get_quadrature(), fe_component);
        ref_scalar_op->allocate_result(n_points, patch_point_vals_->asm_arena());
        auto ref_scalar_value = ref_scalar_op->result_matrix();
        for (unsigned int i_p = 0; i_p < n_points; i_p++)
            for (unsigned int i_dof = 0; i_dof < n_dofs; i_dof++) {
                ref_scalar_value(i_dof)(i_p) = ref_shape_vals[i_p][i_dof][0];
            }

        return FeQArray<Scalar>(patch_point_vals_, true, FeBulk::BulkOps::opRefScalar, fe_component->n_dofs());
    }

    inline FeQArray<Vector> ref_vector(uint component_idx = 0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT((fe_component->fe_type() == FEType::FEVector) | (fe_component->fe_type() == FEType::FEVectorPiola) | (fe_component->fe_type() == FEType::FEVectorContravariant))
                .error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");

        uint n_points = patch_point_vals_->get_quadrature()->size();
        uint n_dofs = fe_component->n_dofs();

        auto *vector_ref_op = patch_point_vals_->make_fixed_fe_op(FeBulk::BulkOps::opRefVector, {3, n_dofs}, &common_reinit::op_base, n_dofs);

        vector_ref_op->allocate_result(n_points, patch_point_vals_->asm_arena());
        auto vector_ref_result = vector_ref_op->result_matrix();
        auto ref_shape_vals = this->ref_shape_values_bulk(patch_point_vals_->get_quadrature(), fe_component);

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
            for (uint i_p=0; i_p<n_points; ++i_p)
                for (uint c=0; c<3; ++c)
                    vector_ref_result(c, i_dof)(i_p) = ref_shape_vals[i_p][i_dof](c);

        return FeQArray<Vector>(patch_point_vals_, true, FeBulk::BulkOps::opRefVector, n_dofs);
    }

    inline FeQArray<Vector> ref_scalar_grad(uint component_idx = 0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");

        uint n_points = patch_point_vals_->get_quadrature()->size();
        uint n_dofs = fe_component->n_dofs();

        auto *ref_scalar_op = patch_point_vals_->make_fixed_fe_op(FeBulk::BulkOps::opRefScalarGrad, {dim, n_dofs}, &common_reinit::op_base, n_dofs);

        std::vector<std::vector<arma::mat> > ref_shape_grads = this->ref_shape_gradients_bulk(patch_point_vals_->get_quadrature(), fe_component);
        ref_scalar_op->allocate_result(n_points, patch_point_vals_->asm_arena());
        auto ref_scalar_value = ref_scalar_op->result_matrix();
        for (uint i=0; i<ref_scalar_value.rows(); ++i) {
            for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
                for (uint i_p=0; i_p<n_points; ++i_p)
                    ref_scalar_value(i, i_dof)(i_p) = ref_shape_grads[i_p][i_dof](i);
        }

        return FeQArray<Vector>(patch_point_vals_, true, FeBulk::BulkOps::opRefScalarGrad, fe_component->n_dofs());
    }

    inline FeQArray<Tensor> ref_vector_grad(uint component_idx = 0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT((fe_component->fe_type() == FEType::FEVector) | (fe_component->fe_type() == FEType::FEVectorPiola) | (fe_component->fe_type() == FEType::FEVectorContravariant))
                .error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");

        uint n_points = patch_point_vals_->get_quadrature()->size();
        uint n_dofs = fe_component->n_dofs();

        auto *ref_vector_op = patch_point_vals_->make_fixed_fe_op(FeBulk::BulkOps::opRefVectorGrad, {dim, 3*n_dofs}, &common_reinit::op_base, n_dofs);

        std::vector<std::vector<arma::mat> > ref_shape_grads = this->ref_shape_gradients_bulk(patch_point_vals_->get_quadrature(), fe_component);
        ref_vector_op->allocate_result(n_points, patch_point_vals_->asm_arena());
        auto ref_vector_value = ref_vector_op->result_matrix();
        for (uint i_c=0; i_c<3; ++i_c) {
            for (uint i_dim=0; i_dim<dim; ++i_dim)
                for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
                    for (uint i_p=0; i_p<n_points; ++i_p)
                        ref_vector_value(i_dim,3*i_dof+i_c)(i_p) = ref_shape_grads[i_p][i_dof](i_dim, i_c);
        }

        return FeQArray<Tensor>(patch_point_vals_, true, FeBulk::BulkOps::opRefVectorGrad, n_dofs);
    }

    /**
     * @brief Return the value of the @p function_no-th shape function at
     * the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Scalar> scalar_shape(uint component_idx = 0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");

        uint n_dofs = fe_component->n_dofs();
        patch_point_vals_->make_fe_op(FeBulk::BulkOps::opScalarShape, {n_dofs}, bulk_reinit::ptop_scalar_shape, n_dofs);

        return FeQArray<Scalar>(patch_point_vals_, true, FeBulk::BulkOps::opScalarShape, n_dofs);
    }

    inline FeQArray<Vector> vector_shape(uint component_idx = 0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);

        uint n_dofs = fe_component->n_dofs();
        uint vector_shape_op_idx = FeBulk::BulkOps::opVectorShape;

        switch (fe_component->fe_type()) {
            case FEVector:
            {
                patch_point_vals_->make_fe_op(vector_shape_op_idx, {3, n_dofs}, bulk_reinit::ptop_vector_shape, n_dofs);
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                patch_point_vals_->make_fe_op(vector_shape_op_idx, {3, n_dofs}, bulk_reinit::ptop_vector_contravariant_shape, n_dofs);
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
                patch_point_vals_->make_fe_op(vector_shape_op_idx, {3, n_dofs}, bulk_reinit::ptop_vector_piola_shape, n_dofs);
                break;
            }
            default:
                ASSERT(false).error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");
        }

        return FeQArray<Vector>(patch_point_vals_, true, vector_shape_op_idx, fe_component->n_dofs());
    }

//    inline FeQArray<Tensor> tensor_shape(uint component_idx = 0)
//    {}

    /**
     * @brief Return the value of the @p function_no-th gradient shape function at
     * the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Vector> grad_scalar_shape(uint component_idx=0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape accessor must be FEScalar!\n");

        uint n_dofs = fe_component->n_dofs();
        patch_point_vals_->make_fe_op(FeBulk::BulkOps::opGradScalarShape, {3, n_dofs}, bulk_reinit::ptop_scalar_shape_grads<dim>, n_dofs);

        return FeQArray<Vector>(patch_point_vals_, true, FeBulk::BulkOps::opGradScalarShape, n_dofs);
    }

    /**
     * @brief Return the value of the @p function_no-th gradient vector shape function
     * at the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> grad_vector_shape(uint component_idx=0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);

        uint n_dofs = fe_component->n_dofs();
        uint vector_shape_grads_op_idx = FeBulk::BulkOps::opGradVectorShape;

        switch (fe_component->fe_type()) {
            case FEVector:
            {
                patch_point_vals_->make_fe_op(vector_shape_grads_op_idx,
                                             {3, 3*n_dofs},
                                             bulk_reinit::ptop_vector_shape_grads<dim>,
                                             n_dofs);
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Grad vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                patch_point_vals_->make_fe_op(vector_shape_grads_op_idx,
                                             {3, 3*n_dofs},
                                             bulk_reinit::ptop_vector_contravariant_shape_grads<dim>,
                                             n_dofs);
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Grad vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
                patch_point_vals_->make_fe_op(vector_shape_grads_op_idx,
                                              {3, 3*n_dofs},
                                              bulk_reinit::ptop_vector_piola_shape_grads<dim>,
                                              n_dofs);
                break;
            }
            default:
                ASSERT(false).error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");
        }

        return FeQArray<Tensor>(patch_point_vals_, true, vector_shape_grads_op_idx, n_dofs);
    }

    /**
     * @brief Return the value of the @p function_no-th vector symmetric gradient
     * at the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> vector_sym_grad(uint component_idx=0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        uint n_dofs = fe_component->n_dofs();
        //ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape accessor must be FEScalar!\n");

        patch_point_vals_->make_fe_op(FeBulk::BulkOps::opVectorSymGrad, {3,3*n_dofs}, common_reinit::ptop_vector_sym_grad, n_dofs);

        return FeQArray<Tensor>(patch_point_vals_, true, FeBulk::BulkOps::opVectorSymGrad, n_dofs);
    }

    /**
     * @brief Return the value of the @p function_no-th vector divergence at
     * the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Scalar> vector_divergence(uint component_idx=0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        uint n_dofs = fe_component->n_dofs();
        //ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape accessor must be FEScalar!\n");

        patch_point_vals_->make_fe_op(FeBulk::BulkOps::opVectorDivergence, {n_dofs}, common_reinit::ptop_vector_divergence, n_dofs);

        return FeQArray<Scalar>(patch_point_vals_, true, FeBulk::BulkOps::opVectorDivergence, n_dofs);
    }

private:
    PatchPointValues<3> *patch_point_vals_;
    std::shared_ptr< FiniteElement<dim> > fe_;
};


template<unsigned int dim>
class SideValues : public BaseValues<dim>
{
public:
	/// Constructor
	SideValues(PatchPointValues<3> *patch_point_vals, MixedPtr<FiniteElement> fe)
	: BaseValues<dim>(), patch_point_vals_(patch_point_vals) {
	    ASSERT_EQ(patch_point_vals->dim(), dim);
	    fe_ = fe[Dim<dim>{}];
	}

    /// Same as BulkValues::JxW but register at side quadrature points.
    inline FeQ<Scalar> JxW()
    {
        return FeQ<Scalar>(patch_point_vals_, false, FeSide::SideOps::opJxW);
    }

    /**
     * @brief Register the normal vector to a side at side quadrature points.
     *
     * @param quad Quadrature.
     */
	inline ElQ<Vector> normal_vector()
	{
        return ElQ<Vector>(patch_point_vals_, FeSide::SideOps::opNormalVec);
	}

	/// Create side accessor of coords entity
    inline ElQ<Vector> coords()
    {
        return ElQ<Vector>(patch_point_vals_, FeSide::SideOps::opCoords);
    }

    /// Create bulk accessor of jac determinant entity
    inline ElQ<Scalar> determinant()
    {
        return ElQ<Scalar>(patch_point_vals_, FeSide::SideOps::opSideJacDet);
    }

    inline FeQArray<Scalar> ref_scalar(uint component_idx = 0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");

        uint n_points = patch_point_vals_->get_quadrature()->size();
        uint n_dofs = fe_component->n_dofs();

        auto *ref_scalar_op = patch_point_vals_->make_fixed_fe_op(FeSide::SideOps::opRefScalar, {dim+1, n_dofs}, &common_reinit::op_base, n_dofs);

        auto ref_shape_vals = this->ref_shape_values_side(patch_point_vals_->get_quadrature(), fe_component);
        ref_scalar_op->allocate_result(n_points, patch_point_vals_->asm_arena());
        auto ref_scalar_value = ref_scalar_op->result_matrix();
        for (unsigned int s=0; s<dim+1; ++s) {
            for (unsigned int i_p = 0; i_p < n_points; i_p++)
                for (unsigned int i_dof = 0; i_dof < n_dofs; i_dof++) {
                    ref_scalar_value(s, i_dof)(i_p) = ref_shape_vals[s][i_p][i_dof][0];
                }
        }

        return FeQArray<Scalar>(patch_point_vals_, true, FeSide::SideOps::opRefScalar, fe_component->n_dofs());
    }

    inline FeQArray<Vector> ref_vector(uint component_idx = 0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT((fe_component->fe_type() == FEType::FEVector) | (fe_component->fe_type() == FEType::FEVectorPiola) | (fe_component->fe_type() == FEType::FEVectorContravariant))
                .error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");

        uint n_points = patch_point_vals_->get_quadrature()->size();
        uint n_dofs = fe_component->n_dofs();

        auto *ref_vector_op = patch_point_vals_->make_fixed_fe_op(FeSide::SideOps::opRefVector, {dim+1,3*n_dofs}, &common_reinit::op_base, n_dofs);
        ref_vector_op->allocate_result(n_points, patch_point_vals_->asm_arena());
        auto ref_vector_value = ref_vector_op->result_matrix();
        auto ref_shape_vals = this->ref_shape_values_side(patch_point_vals_->get_quadrature(), fe_component);
        for (unsigned int s=0; s<dim+1; ++s)
            for (unsigned int i_p = 0; i_p < n_points; i_p++)
                for (unsigned int i_dof = 0; i_dof < n_dofs; i_dof++)
                    for (uint c=0; c<3; ++c)
                        ref_vector_value(s,3*i_dof+c)(i_p) = ref_shape_vals[s][i_p][i_dof][c];

        return FeQArray<Vector>(patch_point_vals_, true, FeSide::SideOps::opRefVector, fe_component->n_dofs());
    }

    inline FeQArray<Vector> ref_scalar_grad(uint component_idx = 0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");

        uint n_points = patch_point_vals_->get_quadrature()->size();
        uint n_dofs = fe_component->n_dofs();

        auto *ref_scalar_op = patch_point_vals_->make_fixed_fe_op(FeSide::SideOps::opRefScalarGrad, {dim+1,dim*n_dofs}, &common_reinit::op_base, n_dofs);

        std::vector<std::vector<std::vector<arma::mat> > > ref_shape_grads = this->ref_shape_gradients_side(patch_point_vals_->get_quadrature(), fe_component);
        ref_scalar_op->allocate_result(n_points, patch_point_vals_->asm_arena());
        auto ref_scalar_value = ref_scalar_op->result_matrix();
        for (unsigned int s=0; s<dim+1; ++s)
            for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
                for (uint i_p=0; i_p<n_points; ++i_p)
                    for (uint c=0; c<dim; ++c)
                        ref_scalar_value(s,dim*i_dof+c)(i_p) = ref_shape_grads[s][i_p][i_dof](c);

        return FeQArray<Vector>(patch_point_vals_, true, FeSide::SideOps::opRefScalarGrad, fe_component->n_dofs());
    }

    inline FeQArray<Tensor> ref_vector_grad(uint component_idx = 0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT((fe_component->fe_type() == FEType::FEVector) | (fe_component->fe_type() == FEType::FEVectorPiola) | (fe_component->fe_type() == FEType::FEVectorContravariant))
                .error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");

        uint n_points = patch_point_vals_->get_quadrature()->size();
        uint n_dofs = fe_component->n_dofs();
        uint n_sides = dim+1;

        auto *ref_vector_op = patch_point_vals_->make_fixed_fe_op(FeSide::SideOps::opRefVectorGrad, {n_sides*dim, 3*n_dofs}, &common_reinit::op_base, n_dofs);

        std::vector<std::vector<std::vector<arma::mat> > > ref_shape_grads = this->ref_shape_gradients_side(patch_point_vals_->get_quadrature(), fe_component);
        ref_vector_op->allocate_result(n_points, patch_point_vals_->asm_arena());
        auto ref_vector_value = ref_vector_op->result_matrix();
        for (uint i_sd=0; i_sd<n_sides; ++i_sd)
            for (uint i_c=0; i_c<3; ++i_c)
                for (uint i_dim=0; i_dim<dim; ++i_dim)
                    for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
                        for (uint i_p=0; i_p<n_points; ++i_p) {
                            ref_vector_value(i_sd*dim+i_dim, 3*i_dof+i_c)(i_p) = ref_shape_grads[i_sd][i_p][i_dof](i_dim, i_c);
                        }

        return FeQArray<Tensor>(patch_point_vals_, true, FeSide::SideOps::opRefVectorGrad, fe_component->n_dofs());
    }

    /// Same as BulkValues::scalar_shape but register at side quadrature points.
    inline FeQArray<Scalar> scalar_shape(uint component_idx = 0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");

        uint n_dofs = fe_component->n_dofs();
        patch_point_vals_->make_fe_op(FeSide::SideOps::opScalarShape, {n_dofs}, side_reinit::ptop_scalar_shape, n_dofs);

        return FeQArray<Scalar>(patch_point_vals_, false, FeSide::SideOps::opScalarShape, n_dofs);
    }

    /// Same as BulkValues::vector_shape but register at side quadrature points.
    inline FeQArray<Vector> vector_shape(uint component_idx = 0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        //ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");

        uint n_dofs = fe_component->n_dofs();
        uint vector_shape_op_idx = FeSide::SideOps::opVectorShape;

        switch (fe_component->fe_type()) {
            case FEVector:
            {
                patch_point_vals_->make_fe_op(vector_shape_op_idx, {3, n_dofs}, side_reinit::ptop_vector_shape, n_dofs);
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                patch_point_vals_->make_fe_op(vector_shape_op_idx, {3, n_dofs}, side_reinit::ptop_vector_contravariant_shape, n_dofs);
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
                patch_point_vals_->make_fe_op(vector_shape_op_idx, {3, n_dofs}, side_reinit::ptop_vector_piola_shape, n_dofs);
                break;
            }
            default:
                ASSERT(false).error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");
        }
        return FeQArray<Vector>(patch_point_vals_, false, vector_shape_op_idx, n_dofs);
    }

    /// Same as BulkValues::grad_scalar_shape but register at side quadrature points.
    inline FeQArray<Vector> grad_scalar_shape(uint component_idx=0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape accessor must be FEScalar!\n");

        uint n_dofs = fe_component->n_dofs();
        patch_point_vals_->make_fe_op(FeSide::SideOps::opGradScalarShape, {3, n_dofs}, side_reinit::ptop_scalar_shape_grads<dim>, n_dofs);

        return FeQArray<Vector>(patch_point_vals_, false, FeSide::SideOps::opGradScalarShape, n_dofs);
    }

    /**
     * @brief Return the value of the @p function_no-th gradient vector shape function
     * at the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> grad_vector_shape(uint component_idx=0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);

        uint n_dofs = fe_component->n_dofs();
        uint vector_shape_grads_op_idx = FeSide::SideOps::opGradVectorShape;

        switch (fe_component->fe_type()) {
            case FEVector:
            {
                patch_point_vals_->make_fe_op(vector_shape_grads_op_idx,
                                             {3, 3*n_dofs},
                                             side_reinit::ptop_vector_shape_grads<dim>,
                                             n_dofs);
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Grad vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                patch_point_vals_->make_fe_op(vector_shape_grads_op_idx,
                                             {3, 3*n_dofs},
                                             side_reinit::ptop_vector_contravariant_shape_grads<dim>,
                                             n_dofs);
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Grad vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
                patch_point_vals_->make_fe_op(vector_shape_grads_op_idx,
                                              {3, 3*n_dofs},
                                              side_reinit::ptop_vector_piola_shape_grads<dim>,
                                              n_dofs);
                break;
            }
            default:
                ASSERT(false).error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");
        }

        return FeQArray<Tensor>(patch_point_vals_, false, vector_shape_grads_op_idx, n_dofs);
    }

    /**
     * @brief Return the value of the @p function_no-th vector symmetric gradient
     * at the @p p side quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Tensor> vector_sym_grad(uint component_idx=0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        uint n_dofs = fe_component->n_dofs();
        //ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape accessor must be FEScalar!\n");

        patch_point_vals_->make_fe_op(FeSide::SideOps::opVectorSymGrad, {3,3*n_dofs}, common_reinit::ptop_vector_sym_grad, n_dofs);

        return FeQArray<Tensor>(patch_point_vals_, false, FeSide::SideOps::opVectorSymGrad, n_dofs);
    }

    /**
     * @brief Return the value of the @p function_no-th vector divergence at
     * the @p p side quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQArray<Scalar> vector_divergence(uint component_idx=0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        uint n_dofs = fe_component->n_dofs();
        //ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape accessor must be FEScalar!\n");

        patch_point_vals_->make_fe_op(FeSide::SideOps::opVectorDivergence, {n_dofs}, common_reinit::ptop_vector_divergence, n_dofs);

        return FeQArray<Scalar>(patch_point_vals_, false, FeSide::SideOps::opVectorDivergence, n_dofs);
    }

private:
    PatchPointValues<3> *patch_point_vals_;
    std::shared_ptr< FiniteElement<dim> > fe_;
};


template<unsigned int dim>
class JoinValues : public BaseValues<dim>
{
public:
	/// Constructor
	JoinValues(PatchPointValues<3> *patch_point_vals_bulk, PatchPointValues<3> *patch_point_vals_side, MixedPtr<FiniteElement> fe)
	: BaseValues<dim>(), patch_point_vals_bulk_(patch_point_vals_bulk), patch_point_vals_side_(patch_point_vals_side) {
	    ASSERT_EQ(patch_point_vals_bulk->dim(), dim-1);
	    ASSERT_EQ(patch_point_vals_side->dim(), dim);
	    fe_high_dim_ = fe[Dim<dim>{}];
	    fe_low_dim_ = fe[Dim<dim-1>{}];
	}

    inline FeQJoin<Scalar> scalar_join_shape(uint component_idx = 0)
    {
    	// element of lower dim (bulk points)
        auto fe_component_low = this->fe_comp(fe_low_dim_, component_idx);
        ASSERT_EQ(fe_component_low->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");
        uint n_dofs_low = fe_component_low->n_dofs();
        patch_point_vals_bulk_->make_fe_op(FeBulk::BulkOps::opScalarShape, {1}, bulk_reinit::ptop_scalar_shape, n_dofs_low);
        patch_point_vals_bulk_->zero_values_needed();

    	// element of higher dim (side points)
        auto fe_component_high = this->fe_comp(fe_high_dim_, component_idx);
        ASSERT_EQ(fe_component_high->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");
        uint n_dofs_high = fe_component_high->n_dofs();
        patch_point_vals_side_->make_fe_op(FeSide::SideOps::opScalarShape, {1}, side_reinit::ptop_scalar_shape, n_dofs_high);
        patch_point_vals_side_->zero_values_needed();

        return FeQJoin<Scalar>(patch_point_vals_bulk_, patch_point_vals_side_, n_dofs_low, n_dofs_high,
                               FeBulk::BulkOps::opScalarShape, FeSide::SideOps::opScalarShape);
    }

    inline FeQJoin<Vector> vector_join_shape(uint component_idx = 0)
    {
    	// element of lower dim (bulk points)
        auto fe_component_low = this->fe_comp(fe_low_dim_, component_idx);
        uint op_idx_bulk = FeBulk::BulkOps::opVectorShape;

        // element of higher dim (side points)
        auto fe_component_high = this->fe_comp(fe_high_dim_, component_idx);
        uint op_idx_side = FeSide::SideOps::opVectorShape;

        ASSERT_EQ(fe_component_high->fe_type(), fe_component_low->fe_type()).error("Type of FiniteElement of low and high element must be same!\n");
        switch (fe_component_low->fe_type()) {
            case FEVector:
            {
                patch_point_vals_bulk_->make_fe_op(op_idx_bulk, {3}, bulk_reinit::ptop_vector_shape, fe_component_low->n_dofs());
                patch_point_vals_side_->make_fe_op(op_idx_side, {3}, side_reinit::ptop_vector_shape, fe_component_high->n_dofs());
                patch_point_vals_bulk_->zero_values_needed();
                patch_point_vals_side_->zero_values_needed();
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                //patch_point_vals_.make_fe_op({3}, bulk_reinit::ptop_vector_contravariant_shape, fe_component->n_dofs());
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
                //patch_point_vals_.make_fe_op({3}, bulk_reinit::ptop_vector_piola_shape, fe_component->n_dofs());
                break;
            }
            default:
                ASSERT(false).error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");
        }

        return FeQJoin<Vector>(patch_point_vals_bulk_, patch_point_vals_side_, fe_component_low->n_dofs(),
                fe_component_high->n_dofs(), op_idx_bulk, op_idx_side);
    }

    inline FeQJoin<Tensor> gradient_vector_join_shape(uint component_idx = 0)
    {
    	// element of lower dim (bulk points)
        auto fe_component_low = this->fe_comp(fe_low_dim_, component_idx);
        uint op_idx_bulk = FeBulk::BulkOps::opGradVectorShape;

        // element of higher dim (side points)
        auto fe_component_high = this->fe_comp(fe_high_dim_, component_idx);
        uint op_idx_side = FeSide::SideOps::opGradVectorShape;

        ASSERT_EQ(fe_component_high->fe_type(), fe_component_low->fe_type()).error("Type of FiniteElement of low and high element must be same!\n");
        switch (fe_component_low->fe_type()) {
            case FEVector:
            {
                patch_point_vals_bulk_->make_fe_op(op_idx_bulk,
                                                  {3, 3},
                                                  bulk_reinit::ptop_vector_shape_grads<dim-1>,
                                                  fe_component_low->n_dofs());

                patch_point_vals_side_->make_fe_op(op_idx_side,
                                                  {3, 3},
                                                  side_reinit::ptop_vector_shape_grads<dim>,
                                                  fe_component_high->n_dofs());

                patch_point_vals_bulk_->zero_values_needed();
                patch_point_vals_side_->zero_values_needed();
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
//                patch_point_vals_bulk_->make_fe_op(op_idx_bulk,
//                                                  {3, 3},
//                                                  bulk_reinit::ptop_vector_contravariant_shape_grads<dim-1>,
//                                                  fe_component_low->n_dofs());
//
//                patch_point_vals_side_->make_fe_op(op_idx_side,
//                                                  {3, 3},
//                                                  side_reinit::ptop_vector_contravariant_shape_grads<dim>,
//                                                  fe_component_high->n_dofs());
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
//                patch_point_vals_bulk_->make_fe_op(op_idx_bulk,
//                                                  {3, 3},
//                                                  bulk_reinit::ptop_vector_piola_shape_grads<dim-1>,
//                                                  fe_component_low->n_dofs());
//
//                patch_point_vals_side_->make_fe_op(op_idx_side,
//                                                  {3, 3},
//                                                  side_reinit::ptop_vector_piola_shape_grads<dim>,
//                                                  fe_component_high->n_dofs());
                break;
            }
            default:
                ASSERT(false).error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");
        }

        return FeQJoin<Tensor>(patch_point_vals_bulk_, patch_point_vals_side_, fe_component_low->n_dofs(),
                fe_component_high->n_dofs(), op_idx_bulk, op_idx_side);
    }

private:
    PatchPointValues<3> *patch_point_vals_bulk_;
    PatchPointValues<3> *patch_point_vals_side_;
    std::shared_ptr< FiniteElement<dim> > fe_high_dim_;
    std::shared_ptr< FiniteElement<dim-1> > fe_low_dim_;
};

/// Template specialization of dim = 1
template <>
class JoinValues<1> : public BaseValues<1>
{
public:
	/// Constructor
	JoinValues(FMT_UNUSED PatchPointValues<3> *patch_point_vals_bulk, FMT_UNUSED PatchPointValues<3> *patch_point_vals_side, FMT_UNUSED MixedPtr<FiniteElement> fe)
	: BaseValues<1>() {}

    inline FeQJoin<Scalar> scalar_join_shape(FMT_UNUSED uint component_idx = 0)
    {
        return FeQJoin<Scalar>();
    }

    inline FeQJoin<Vector> vector_join_shape(FMT_UNUSED uint component_idx = 0)
    {
        return FeQJoin<Vector>();
    }

    inline FeQJoin<Tensor> gradient_vector_join_shape(FMT_UNUSED uint component_idx = 0)
    {
        return FeQJoin<Tensor>();
    }
};


template<unsigned int spacedim = 3>
class PatchFEValues {
public:
    typedef typename PatchPointValues<spacedim>::PatchFeData PatchFeData;

    /// Struct for pre-computing number of elements, sides, bulk points and side points on each dimension.
    struct TableSizes {
    public:
        /// Constructor
        TableSizes() {
            elem_sizes_ = std::vector<std::vector<uint> >(2, std::vector<uint>(spacedim,0));
            point_sizes_ = std::vector<std::vector<uint> >(2, std::vector<uint>(spacedim,0));
        }

        /// Set all values to zero
        void reset() {
            std::fill(elem_sizes_[0].begin(), elem_sizes_[0].end(), 0);
            std::fill(elem_sizes_[1].begin(), elem_sizes_[1].end(), 0);
            std::fill(point_sizes_[0].begin(), point_sizes_[0].end(), 0);
            std::fill(point_sizes_[1].begin(), point_sizes_[1].end(), 0);
        }

        /// Copy values of other TableSizes instance
        void copy(const TableSizes &other) {
            elem_sizes_[0] = other.elem_sizes_[0];
            elem_sizes_[1] = other.elem_sizes_[1];
            point_sizes_[0] = other.point_sizes_[0];
            point_sizes_[1] = other.point_sizes_[1];
        }

        /**
         * Holds number of elements and sides on each dimension
         * Format:
         *  { {n_elements_1D, n_elements_2D, n_elements_3D },
         *    {n_sides_1D, n_sides_2D, n_sides_3D } }
         */
        std::vector<std::vector<uint> > elem_sizes_;

        /**
         * Holds number of bulk and side points on each dimension
         * Format:
         *  { {n_bulk_points_1D, n_bulk_points_2D, n_bulk_points_3D },
         *    {n_side_points_1D, n_side_points_2D, n_side_points_3D } }
         */
        std::vector<std::vector<uint> > point_sizes_;
    };

    PatchFEValues()
    : patch_fe_data_(1024 * 1024, 256),
      patch_point_vals_bulk_{ {FeBulk::PatchPointValues(1, 0, patch_fe_data_),
                               FeBulk::PatchPointValues(2, 0, patch_fe_data_),
                               FeBulk::PatchPointValues(3, 0, patch_fe_data_)} },
      patch_point_vals_side_{ {FeSide::PatchPointValues(1, 0, patch_fe_data_),
                               FeSide::PatchPointValues(2, 0, patch_fe_data_),
                               FeSide::PatchPointValues(3, 0, patch_fe_data_)} }
    {
        used_quads_[0] = false; used_quads_[1] = false;
    }

    PatchFEValues(unsigned int quad_order, MixedPtr<FiniteElement> fe)
    : patch_fe_data_(1024 * 1024, 256),
      patch_point_vals_bulk_{ {FeBulk::PatchPointValues(1, quad_order, patch_fe_data_),
    	                       FeBulk::PatchPointValues(2, quad_order, patch_fe_data_),
                               FeBulk::PatchPointValues(3, quad_order, patch_fe_data_)} },
      patch_point_vals_side_{ {FeSide::PatchPointValues(1, quad_order, patch_fe_data_),
                               FeSide::PatchPointValues(2, quad_order, patch_fe_data_),
                               FeSide::PatchPointValues(3, quad_order, patch_fe_data_)} },
      fe_(fe)
    {
        used_quads_[0] = false; used_quads_[1] = false;

        // TODO move initialization zero_vec_ to patch_fe_data_ constructor when we will create separate ArenaVec of DOshape functions
        uint max_n_dofs = std::max(fe_[Dim<1>{}]->n_dofs(), std::max(fe_[Dim<2>{}]->n_dofs(), fe_[Dim<3>{}]->n_dofs()) );
        uint zero_vec_size = 300 * max_n_dofs;
        patch_fe_data_.zero_vec_ = ArenaVec<double>(zero_vec_size, patch_fe_data_.asm_arena_);
        for (uint i=0; i<zero_vec_size; ++i) patch_fe_data_.zero_vec_(i) = 0.0;
    }


    /// Destructor
    ~PatchFEValues()
    {}

    /// Return bulk or side quadrature of given dimension
    Quadrature *get_quadrature(uint dim, bool is_bulk) const {
        if (is_bulk) return patch_point_vals_bulk_[dim-1].get_quadrature();
        else return patch_point_vals_side_[dim-1].get_quadrature();
    }

    /**
	 * @brief Initialize structures and calculates cell-independent data.
	 *
	 * @param _quadrature The quadrature rule for the cell associated
     *                    to given finite element or for the cell side.
	 * @param _flags The update flags.
	 */
    template<unsigned int DIM>
    void initialize(Quadrature &_quadrature)
    {
        if ( _quadrature.dim() == DIM ) {
            used_quads_[0] = true;
            patch_point_vals_bulk_[DIM-1].initialize(); // bulk
        } else {
            used_quads_[1] = true;
            patch_point_vals_side_[DIM-1].initialize(); // side
        }
    }

    /// Finalize initialization, creates child (patch) arena and passes it to PatchPointValue objects
    void init_finalize() {
        patch_fe_data_.patch_arena_ = patch_fe_data_.asm_arena_.get_child_arena();
    }

    /// Reset PatchpointValues structures
    void reset()
    {
        for (unsigned int i=0; i<3; ++i) {
            if (used_quads_[0]) patch_point_vals_bulk_[i].reset();
            if (used_quads_[1]) patch_point_vals_side_[i].reset();
        }
        patch_fe_data_.patch_arena_->reset();
    }

    /// Reinit data.
    void reinit_patch()
    {
        for (unsigned int i=0; i<3; ++i) {
            if (used_quads_[0]) patch_point_vals_bulk_[i].reinit_patch();
            if (used_quads_[1]) patch_point_vals_side_[i].reinit_patch();
        }
    }

    /**
     * @brief Returns the number of shape functions.
     */
    template<unsigned int dim>
    inline unsigned int n_dofs() const {
        ASSERT((dim>=0) && (dim<=3))(dim).error("Dimension must be 0, 1, 2 or 3.");
        return fe_[Dim<dim>{}]->n_dofs();
    }

    /// Return BulkValue object of dimension given by template parameter
    template<unsigned int dim>
    BulkValues<dim> bulk_values() {
    	ASSERT((dim>0) && (dim<=3))(dim).error("Dimension must be 1, 2 or 3.");
        return BulkValues<dim>(&patch_point_vals_bulk_[dim-1], fe_);
    }

    /// Return SideValue object of dimension given by template parameter
    template<unsigned int dim>
    SideValues<dim> side_values() {
    	ASSERT((dim>0) && (dim<=3))(dim).error("Dimension must be 1, 2 or 3.");
        return SideValues<dim>(&patch_point_vals_side_[dim-1], fe_);
    }

    /// Return JoinValue object of dimension given by template parameter
    template<unsigned int dim>
    JoinValues<dim> join_values() {
    	//ASSERT((dim>1) && (dim<=3))(dim).error("Dimension must be 2 or 3.");
        return JoinValues<dim>(&patch_point_vals_bulk_[dim-2], &patch_point_vals_side_[dim-1], fe_);
    }

    /** Following methods are used during update of patch. **/

    /// Resize tables of patch_point_vals_
    void resize_tables(TableSizes table_sizes) {
        for (uint i=0; i<3; ++i) {
            if (used_quads_[0]) patch_point_vals_bulk_[i].resize_tables(table_sizes.elem_sizes_[0][i], table_sizes.point_sizes_[0][i]);
            if (used_quads_[1]) patch_point_vals_side_[i].resize_tables(table_sizes.elem_sizes_[1][i], table_sizes.point_sizes_[1][i]);
        }
    }

    /// Register element to patch_point_vals_ table by dimension of element
    uint register_element(DHCellAccessor cell, uint element_patch_idx) {
        arma::mat coords;
        switch (cell.dim()) {
        case 1:
            coords = MappingP1<1,spacedim>::element_map(cell.elm());
            return patch_point_vals_bulk_[0].register_element(coords, element_patch_idx);
            break;
        case 2:
        	coords = MappingP1<2,spacedim>::element_map(cell.elm());
            return patch_point_vals_bulk_[1].register_element(coords, element_patch_idx);
            break;
        case 3:
        	coords = MappingP1<3,spacedim>::element_map(cell.elm());
            return patch_point_vals_bulk_[2].register_element(coords, element_patch_idx);
            break;
        default:
        	ASSERT(false);
        	return 0;
            break;
        }
    }

    /// Register side to patch_point_vals_ table by dimension of side
    uint register_side(DHCellSide cell_side) {
        arma::mat side_coords(spacedim, cell_side.dim());
        for (unsigned int n=0; n<cell_side.dim(); n++)
            for (unsigned int c=0; c<spacedim; c++)
                side_coords(c,n) = (*cell_side.side().node(n))[c];

        arma::mat elm_coords;
        DHCellAccessor cell = cell_side.cell();
        switch (cell.dim()) {
        case 1:
            elm_coords = MappingP1<1,spacedim>::element_map(cell.elm());
            return patch_point_vals_side_[0].register_side(elm_coords, side_coords, cell_side.side_idx());
            break;
        case 2:
            elm_coords = MappingP1<2,spacedim>::element_map(cell.elm());
            return patch_point_vals_side_[1].register_side(elm_coords, side_coords, cell_side.side_idx());
            break;
        case 3:
            elm_coords = MappingP1<3,spacedim>::element_map(cell.elm());
            return patch_point_vals_side_[2].register_side(elm_coords, side_coords, cell_side.side_idx());
            break;
        default:
            ASSERT(false);
            return 0;
            break;
        }
    }

    /// Register bulk point to patch_point_vals_ table by dimension of element
    uint register_bulk_point(DHCellAccessor cell, uint elem_table_row, uint value_patch_idx, uint i_point_on_elem) {
        return patch_point_vals_bulk_[cell.dim()-1].register_bulk_point(elem_table_row, value_patch_idx, cell.elm_idx(), i_point_on_elem);
    }

    /// Register side point to patch_point_vals_ table by dimension of side
    uint register_side_point(DHCellSide cell_side, uint elem_table_row, uint value_patch_idx, uint i_point_on_side) {
        return patch_point_vals_side_[cell_side.dim()-1].register_side_point(elem_table_row, value_patch_idx, cell_side.elem_idx(),
                cell_side.side_idx(), i_point_on_side);
    }

    /// Temporary development method
    void print_data_tables(ostream& stream, bool points, bool ints, bool only_bulk=true) const {
        stream << endl << "Table of patch FE data:" << endl;
        for (uint i=0; i<3; ++i) {
            stream << std::setfill('-') << setw(100) << "" << endl;
            stream << "Bulk, dimension " << (i+1) << endl;
            patch_point_vals_bulk_[i].print_data_tables(stream, points, ints);
        }
        if (!only_bulk)
            for (uint i=0; i<3; ++i) {
                stream << std::setfill('-') << setw(100) << "" << endl;
                stream << "Side, dimension " << (i+1) << endl;
                patch_point_vals_side_[i].print_data_tables(stream, points, ints);
            }
        stream << std::setfill('=') << setw(100) << "" << endl;
    }

    /// Temporary development method
    void print_operations(ostream& stream) const {
        stream << endl << "Table of patch FE operations:" << endl;
        for (uint i=0; i<3; ++i) {
            stream << std::setfill('-') << setw(100) << "" << endl;
            stream << "Bulk, dimension " << (i+1) << endl;
            patch_point_vals_bulk_[i].print_operations(stream, 0);
        }
        for (uint i=0; i<3; ++i) {
            stream << std::setfill('-') << setw(100) << "" << endl;
            stream << "Side, dimension " << (i+1) << endl;
            patch_point_vals_side_[i].print_operations(stream, 1);
        }
        stream << std::setfill('=') << setw(100) << "" << endl;
    }

private:
    PatchFeData patch_fe_data_;
    std::array<FeBulk::PatchPointValues<spacedim>, 3> patch_point_vals_bulk_;  ///< Sub objects of bulk data of dimensions 1,2,3
    std::array<FeSide::PatchPointValues<spacedim>, 3> patch_point_vals_side_;  ///< Sub objects of side data of dimensions 1,2,3

    MixedPtr<FiniteElement> fe_;   ///< Mixed of shared pointers of FiniteElement object
    bool used_quads_[2];           ///< Pair of flags signs holds info if bulk and side quadratures are used

    template <class ValueType>
    friend class ElQ;
    template <class ValueType>
    friend class FeQ;
};


template <class ValueType>
ValueType ElQ<ValueType>::operator()(const BulkPoint &point) const {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_->scalar_elem_value(op_idx_, value_cache_idx);
}

template <>
inline Vector ElQ<Vector>::operator()(const BulkPoint &point) const {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_->vector_elem_value(op_idx_, value_cache_idx);
}

template <>
inline Tensor ElQ<Tensor>::operator()(const BulkPoint &point) const {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_->tensor_elem_value(op_idx_, value_cache_idx);
}

template <class ValueType>
ValueType ElQ<ValueType>::operator()(const SidePoint &point) const {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_->scalar_elem_value(op_idx_, value_cache_idx);
}

template <>
inline Vector ElQ<Vector>::operator()(const SidePoint &point) const {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_->vector_elem_value(op_idx_, value_cache_idx);
}

template <>
inline Tensor ElQ<Tensor>::operator()(const SidePoint &point) const {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_->tensor_elem_value(op_idx_, value_cache_idx);
}

template <class ValueType>
ValueType FeQ<ValueType>::operator()(const BulkPoint &point) const {
    ASSERT_PTR(patch_point_vals_bulk_);

    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_bulk_->scalar_value(op_idx_, value_cache_idx, i_shape_fn_idx_);
}

template <>
inline Vector FeQ<Vector>::operator()(const BulkPoint &point) const {
	ASSERT_PTR(patch_point_vals_bulk_);

    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_bulk_->vector_value(op_idx_, value_cache_idx, i_shape_fn_idx_);
}

template <>
inline Tensor FeQ<Tensor>::operator()(const BulkPoint &point) const {
	ASSERT_PTR(patch_point_vals_bulk_);

    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_bulk_->tensor_value(op_idx_, value_cache_idx, i_shape_fn_idx_);
}

template <class ValueType>
ValueType FeQ<ValueType>::operator()(const SidePoint &point) const {
	ASSERT_PTR(patch_point_vals_side_);

    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_side_->scalar_value(op_idx_, value_cache_idx, i_shape_fn_idx_);
}

template <>
inline Vector FeQ<Vector>::operator()(const SidePoint &point) const {
	ASSERT_PTR(patch_point_vals_side_);

    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_side_->vector_value(op_idx_, value_cache_idx, i_shape_fn_idx_);
}

template <>
inline Tensor FeQ<Tensor>::operator()(const SidePoint &point) const {
	ASSERT_PTR(patch_point_vals_side_);

    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_side_->tensor_value(op_idx_, value_cache_idx, i_shape_fn_idx_);
}


#endif /* PATCH_FE_VALUES_HH_ */
