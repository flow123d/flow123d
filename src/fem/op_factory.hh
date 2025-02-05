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
 * @file    op_factory.hh
 * @brief   Declares top level factory classes of FE operations.
 * @author  David Flanderka
 */

#ifndef OP_FACTORY_HH_
#define OP_FACTORY_HH_

#include "fem/eigen_tools.hh"
#include "fem/patch_point_values.hh"
#include "fem/op_function.hh"
#include "fem/op_accessors_impl.hh"

template<unsigned int spacedim> class PatchFEValues;


template<unsigned int dim>
class BaseValues
{
protected:
	// Default constructor
	BaseValues(PatchFEValues<3> &pfev) : patch_fe_values_(pfev)
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

    PatchFEValues<3> &patch_fe_values_;
};

template<unsigned int dim>
class BulkValues : public BaseValues<dim>
{
public:
	/// Constructor
	BulkValues(PatchPointValues<3> *patch_point_vals, PatchFEValues<3> &pfev, MixedPtr<FiniteElement> fe)
	: BaseValues<dim>(pfev), patch_point_vals_(patch_point_vals) {
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
        PatchOp<3> *op = this->patch_fe_values_.template get< Op::Bulk::Pt::OpJxW<dim, 3> >(dim);
        return FeQ<Scalar>(op, true);
    }

	/// Create bulk accessor of coords entity
    inline FeQ<Vector> coords()
    {
        PatchOp<3> *op = this->patch_fe_values_.template get< Op::Bulk::Pt::OpCoords >(dim);
        return FeQ<Vector>(op, true);
    }

//    inline ElQ<Tensor> jacobian(std::initializer_list<Quadrature *> quad_list)
//    {}

    /// Create bulk accessor of jac determinant entity
    inline ElQ<Scalar> determinant()
    {
        PatchOp<3> *op = this->patch_fe_values_.template get< Op::Bulk::El::OpJacDet<dim, 3> >(dim);
        return ElQ<Scalar>(op);
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
//    inline FeQArray<Vector> grad_scalar_shape(uint component_idx=0)
//    {
//        auto fe_component = this->fe_comp(fe_, component_idx);
//        ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape accessor must be FEScalar!\n");
//
//        uint n_dofs = fe_component->n_dofs();
//        patch_point_vals_->make_fe_op(FeBulk::BulkOps::opGradScalarShape, {3, n_dofs}, bulk_reinit::ptop_scalar_shape_grads<dim>, n_dofs);
//
//        return FeQArray<Vector>(patch_point_vals_, true, FeBulk::BulkOps::opGradScalarShape, n_dofs);
//    }
    inline FeQArray<Vector> grad_scalar_shape(uint component_idx=0)
    {
        PatchOp<3> *op = this->patch_fe_values_.template get< Op::Bulk::Pt::OpGradScalarShape<dim, 3> >(dim, component_idx);
        return FeQArray<Vector>(op, true);
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
	SideValues(PatchPointValues<3> *patch_point_vals, PatchFEValues<3> &pfev, MixedPtr<FiniteElement> fe)
	: BaseValues<dim>(pfev), patch_point_vals_(patch_point_vals) {
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
	JoinValues(PatchPointValues<3> *patch_point_vals_bulk, PatchPointValues<3> *patch_point_vals_side, PatchFEValues<3> &pfev, MixedPtr<FiniteElement> fe)
	: BaseValues<dim>(pfev), patch_point_vals_bulk_(patch_point_vals_bulk), patch_point_vals_side_(patch_point_vals_side) {
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
        patch_point_vals_bulk_->make_fe_op(FeBulk::BulkOps::opScalarShape, {n_dofs_low}, bulk_reinit::ptop_scalar_shape, n_dofs_low);
        patch_point_vals_bulk_->zero_values_needed();

    	// element of higher dim (side points)
        auto fe_component_high = this->fe_comp(fe_high_dim_, component_idx);
        ASSERT_EQ(fe_component_high->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");
        uint n_dofs_high = fe_component_high->n_dofs();
        patch_point_vals_side_->make_fe_op(FeSide::SideOps::opScalarShape, {n_dofs_high}, side_reinit::ptop_scalar_shape, n_dofs_high);
        patch_point_vals_side_->zero_values_needed();

        return FeQJoin<Scalar>(patch_point_vals_bulk_, patch_point_vals_side_, n_dofs_low, n_dofs_high,
                               FeBulk::BulkOps::opScalarShape, FeSide::SideOps::opScalarShape);
    }

    inline FeQJoin<Vector> vector_join_shape(uint component_idx = 0)
    {
    	// element of lower dim (bulk points)
        auto fe_component_low = this->fe_comp(fe_low_dim_, component_idx);
        uint op_idx_bulk = FeBulk::BulkOps::opVectorShape;
        uint n_dofs_low = fe_component_low->n_dofs();

        // element of higher dim (side points)
        auto fe_component_high = this->fe_comp(fe_high_dim_, component_idx);
        uint op_idx_side = FeSide::SideOps::opVectorShape;
        uint n_dofs_high = fe_component_high->n_dofs();

        ASSERT_EQ(fe_component_high->fe_type(), fe_component_low->fe_type()).error("Type of FiniteElement of low and high element must be same!\n");
        switch (fe_component_low->fe_type()) {
            case FEVector:
            {
                patch_point_vals_bulk_->make_fe_op(op_idx_bulk, {3, n_dofs_low}, bulk_reinit::ptop_vector_shape, n_dofs_low);
                patch_point_vals_side_->make_fe_op(op_idx_side, {3, n_dofs_high}, side_reinit::ptop_vector_shape, n_dofs_high);
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

        return FeQJoin<Vector>(patch_point_vals_bulk_, patch_point_vals_side_, n_dofs_low, n_dofs_high, op_idx_bulk, op_idx_side);
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
                                                  {3, 3-fe_component_low->n_dofs()},
                                                  bulk_reinit::ptop_vector_shape_grads<dim-1>,
                                                  fe_component_low->n_dofs());

                patch_point_vals_side_->make_fe_op(op_idx_side,
                                                  {3, 3*fe_component_high->n_dofs()},
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
	JoinValues(FMT_UNUSED PatchPointValues<3> *patch_point_vals_bulk, FMT_UNUSED PatchPointValues<3> *patch_point_vals_side,
	        PatchFEValues<3> &pfev, FMT_UNUSED MixedPtr<FiniteElement> fe)
	: BaseValues<1>(pfev) {}

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


#endif /* OP_FACTORY_HH_ */
