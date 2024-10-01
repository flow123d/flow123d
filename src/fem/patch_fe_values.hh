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
#include "fem/element_values.hh"              // for ElementValues
#include "fem/fe_values.hh"                   // for FEValuesBase
#include "fem/fe_values_views.hh"             // for FEValuesViews
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
    ElQ(PatchPointValues<3> &patch_point_vals, unsigned int op_idx)
    : patch_point_vals_(patch_point_vals), op_idx_(op_idx) {}

    ValueType operator()(FMT_UNUSED const BulkPoint &point);

    ValueType operator()(FMT_UNUSED const SidePoint &point);

private:
    PatchPointValues<3> &patch_point_vals_; ///< Reference to PatchPointValues
    unsigned int op_idx_;                   ///< Index of operation in patch_point_vals_.operations vector
};


template <class ValueType>
class FeQ {
public:
    /// Forbidden default constructor
    FeQ() = delete;

    // Class similar to current FeView
    FeQ(PatchPointValues<3> &patch_point_vals, unsigned int op_idx, unsigned int n_dofs)
    : patch_point_vals_(patch_point_vals), op_idx_(op_idx), n_dofs_(n_dofs) {}


    ValueType operator()(FMT_UNUSED unsigned int shape_idx, FMT_UNUSED const BulkPoint &point);

    ValueType operator()(FMT_UNUSED unsigned int shape_idx, FMT_UNUSED const SidePoint &point);

    // Implementation for EdgePoint, SidePoint, and JoinPoint shoud have a common implementation
    // resolving to side values

private:
    PatchPointValues<3> &patch_point_vals_; ///< Reference to PatchPointValues
    unsigned int op_idx_;                   ///< Index of operation in patch_point_vals_.operations vector
    unsigned int n_dofs_;                   ///< Number of DOFs
};


template <class ValueType>
class JoinShapeAccessor {
public:
    /// Default constructor
    JoinShapeAccessor()
    : patch_point_vals_bulk_(nullptr), patch_point_vals_side_(nullptr), join_idx_(-1) {}

    /**
     * Constructor
     *
     * @param patch_point_vals_bulk  Pointer to PatchPointValues bulk object.
     * @param patch_point_vals_side  Pointer to PatchPointValues side object.
     * @param begin                  Index of the first component of the bulk Quantity.
     * @param begin_side             Index of the first component of the side Quantity.
     * @param n_dofs_bulk            Number of DOFs of bulk (lower-dim) element.
     * @param n_dofs_side            Number of DOFs of side (higher-dim) element.
     * @param join_idx               Index function.
     */
    JoinShapeAccessor(PatchPointValues<3> *patch_point_vals_bulk, PatchPointValues<3> *patch_point_vals_side, unsigned int n_dofs_bulk,
            unsigned int n_dofs_side, unsigned int op_idx_bulk, unsigned int op_idx_side, unsigned int join_idx)
    : patch_point_vals_bulk_(patch_point_vals_bulk), patch_point_vals_side_(patch_point_vals_side),
	  n_dofs_high_(n_dofs_side), n_dofs_low_(n_dofs_bulk), op_idx_bulk_(op_idx_bulk), op_idx_side_(op_idx_side), join_idx_(join_idx) {
        //ASSERT( (patch_point_vals_bulk->dim()==2) || (patch_point_vals_bulk->dim()==3) )(patch_point_vals_bulk->dim() ).error("Invalid dimension, must be 2 or 3!");
    }

    /// Return global index of DOF
    inline unsigned int join_idx() const {
        return join_idx_;
    }

    /// Return local index of DOF (on low / high-dim) - should be private method
    inline unsigned int local_idx() const {
        if (this->is_high_dim()) return (join_idx_ - n_dofs_low_);
        else return join_idx_;
    }

    inline unsigned int n_dofs_low() const {
        return n_dofs_low_;
    }

    inline unsigned int n_dofs_high() const {
        return n_dofs_high_;
    }

    inline unsigned int n_dofs_both() const {
        return n_dofs_high_ + n_dofs_low_;
    }

    inline bool is_high_dim() const {
        return (join_idx_ >= n_dofs_low_);
    }

    /// Iterates to next item.
    inline void inc() {
        join_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const JoinShapeAccessor<ValueType>& other) const {
    	return (join_idx_ == other.join_idx_);
    }


    ValueType operator()(const BulkPoint &point);

    ValueType operator()(const SidePoint &point);

private:
    // attributes:
    PatchPointValues<3> *patch_point_vals_bulk_;  ///< Pointer to bulk PatchPointValues
    PatchPointValues<3> *patch_point_vals_side_;  ///< Pointer to side PatchPointValues
    unsigned int n_dofs_high_;                    ///< Number of DOFs on high-dim element
    unsigned int n_dofs_low_;                     ///< Number of DOFs on low-dim element
    unsigned int op_idx_bulk_;                    ///< Index of operation in patch_point_vals_bulk_.operations vector
    unsigned int op_idx_side_;                    ///< Index of operation in patch_point_vals_side_.operations vector
    unsigned int join_idx_;                       ///< Index of processed DOF
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
            auto quad = q->make_from_side<dim>(sid);
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

};

template<unsigned int dim>
class BulkValues : public BaseValues<dim>
{
public:
	/// Constructor
	BulkValues(PatchPointValues<3> &patch_point_vals, MixedPtr<FiniteElement> fe)
	: BaseValues<dim>(), patch_point_vals_(patch_point_vals) {
	    ASSERT_EQ(patch_point_vals.dim(), dim);
	    fe_ = fe[Dim<dim>{}];
	}

    /**
     * @brief Register the product of Jacobian determinant and the quadrature
     * weight at bulk quadrature points.
     *
     * @param quad Quadrature.
     */
    inline ElQ<Scalar> JxW()
    {
        return ElQ<Scalar>(patch_point_vals_, FeBulk::BulkOps::opJxW);
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

    /**
     * @brief Return the value of the @p function_no-th shape function at
     * the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQ<Scalar> scalar_shape(uint component_idx = 0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");

        // use lambda reinit function
        std::vector< std::vector<double> > shape_values( fe_component->n_dofs(), vector<double>(patch_point_vals_.get_quadrature()->size()) );
        auto ref_shape_vals = this->ref_shape_values_bulk(patch_point_vals_.get_quadrature(), fe_component);
        for (unsigned int i = 0; i < patch_point_vals_.get_quadrature()->size(); i++)
            for (unsigned int j = 0; j < fe_component->n_dofs(); j++) {
            	shape_values[j][i] = ref_shape_vals[i][j][0];
            }
        uint scalar_shape_op_idx = patch_point_vals_.operations_.size(); // index in operations_ vector
        auto lambda_scalar_shape = [shape_values, scalar_shape_op_idx](std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
                bulk_reinit::ptop_scalar_shape(operations, shape_values, scalar_shape_op_idx);
            };
        patch_point_vals_.make_fe_op({1}, lambda_scalar_shape, {}, fe_component->n_dofs());
        patch_point_vals_.set_fe_op(FEOps::opScalarShape, scalar_shape_op_idx);

        return FeQ<Scalar>(patch_point_vals_, scalar_shape_op_idx, fe_component->n_dofs());
    }

    inline FeQ<Vector> vector_shape(uint component_idx = 0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);

        // use lambda reinit function
        std::vector< std::vector<arma::vec3> > shape_values( fe_component->n_dofs(), vector<arma::vec3>(patch_point_vals_.get_quadrature()->size()) );
        auto ref_shape_vals = this->ref_shape_values_bulk(patch_point_vals_.get_quadrature(), fe_component);
        for (unsigned int i = 0; i < patch_point_vals_.get_quadrature()->size(); i++)
            for (unsigned int j = 0; j < fe_component->n_dofs(); j++) {
            	shape_values[j][i] = ref_shape_vals[i][j];
            }
        uint vector_shape_op_idx = patch_point_vals_.operations_.size(); // index in operations_ vector

        switch (fe_component->fe_type()) {
            case FEVector:
            {
                auto lambda_vector_shape = [shape_values, vector_shape_op_idx](std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
                        bulk_reinit::ptop_vector_shape(operations, shape_values, vector_shape_op_idx);
                    };
                patch_point_vals_.make_fe_op({3}, lambda_vector_shape, {}, fe_component->n_dofs());
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                auto lambda_contravariant_shape = [shape_values, vector_shape_op_idx](std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
                        bulk_reinit::ptop_vector_contravariant_shape(operations, shape_values, vector_shape_op_idx);
                    };
                patch_point_vals_.make_fe_op({3}, lambda_contravariant_shape, {FeBulk::BulkOps::opJac}, fe_component->n_dofs());
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
                auto lambda_piola_shape = [shape_values, vector_shape_op_idx](std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
                        bulk_reinit::ptop_vector_piola_shape(operations, shape_values, vector_shape_op_idx);
                    };
                patch_point_vals_.make_fe_op({3}, lambda_piola_shape, {FeBulk::BulkOps::opJac, FeBulk::BulkOps::opJacDet}, fe_component->n_dofs());
                break;
            }
            default:
                ASSERT(false).error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");
        }
        patch_point_vals_.set_fe_op(FEOps::opVectorShape, vector_shape_op_idx);

        return FeQ<Vector>(patch_point_vals_, vector_shape_op_idx, fe_component->n_dofs());
    }

//    inline FeQ<Tensor> tensor_shape(uint component_idx = 0)
//    {}

    /**
     * @brief Return the value of the @p function_no-th gradient shape function at
     * the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQ<Vector> grad_scalar_shape(uint component_idx=0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape accessor must be FEScalar!\n");

        // use lambda reinit function
        auto ref_shape_grads = this->ref_shape_gradients(fe_component);
        uint scalar_shape_grads_op_idx = patch_point_vals_.operations_.size(); // index in operations_ vector
        auto lambda_scalar_shape_grad = [ref_shape_grads, scalar_shape_grads_op_idx](std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
                bulk_reinit::ptop_scalar_shape_grads<dim>(operations, ref_shape_grads, scalar_shape_grads_op_idx);
            };
        patch_point_vals_.make_fe_op({3}, lambda_scalar_shape_grad, {FeBulk::BulkOps::opInvJac}, fe_component->n_dofs());
        patch_point_vals_.set_fe_op(FEOps::opGradScalarShape, scalar_shape_grads_op_idx);

        return FeQ<Vector>(patch_point_vals_, scalar_shape_grads_op_idx, fe_component->n_dofs());
    }

    /**
     * @brief Return the value of the @p function_no-th gradient vector shape function
     * at the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQ<Tensor> grad_vector_shape(uint component_idx=0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);

        // use lambda reinit function
        auto ref_shape_grads = this->ref_shape_gradients(fe_component);
        uint vector_shape_grads_op_idx = patch_point_vals_.operations_.size(); // index in operations_ vector

        switch (fe_component->fe_type()) {
            case FEVector:
            {
                auto lambda_vector_shape_grad = [ref_shape_grads, vector_shape_grads_op_idx](std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
                        bulk_reinit::ptop_vector_shape_grads<dim>(operations, ref_shape_grads, vector_shape_grads_op_idx);
                    };
                patch_point_vals_.make_fe_op({3, 3},
                                             lambda_vector_shape_grad,
                                             {FeBulk::BulkOps::opInvJac},
                                             fe_component->n_dofs());
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Grad vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                auto lambda_contravariant_shape_grad = [ref_shape_grads, vector_shape_grads_op_idx](std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
                        bulk_reinit::ptop_vector_contravariant_shape_grads<dim>(operations, ref_shape_grads, vector_shape_grads_op_idx);
                    };
                patch_point_vals_.make_fe_op({3, 3},
                                             lambda_contravariant_shape_grad,
                                             {FeBulk::BulkOps::opInvJac, FeBulk::BulkOps::opJac},
                                             fe_component->n_dofs());
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Grad vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
                auto lambda_piola_shape_grad = [ref_shape_grads, vector_shape_grads_op_idx](std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
                        bulk_reinit::ptop_vector_piola_shape_grads<dim>(operations, ref_shape_grads, vector_shape_grads_op_idx);
                    };
                patch_point_vals_.make_fe_op({3, 3},
                                              lambda_piola_shape_grad,
                                              {FeBulk::BulkOps::opInvJac, FeBulk::BulkOps::opJac, FeBulk::BulkOps::opJacDet},
                                              fe_component->n_dofs());
                break;
            }
            default:
                ASSERT(false).error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");
        }
        patch_point_vals_.set_fe_op(FEOps::opGradVectorShape, vector_shape_grads_op_idx);

        return FeQ<Tensor>(patch_point_vals_, vector_shape_grads_op_idx, fe_component->n_dofs());
    }

    /**
     * @brief Return the value of the @p function_no-th vector symmetric gradient
     * at the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQ<Tensor> vector_sym_grad(uint component_idx=0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        //ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape accessor must be FEScalar!\n");

        // use lambda reinit function
        uint vector_sym_grad_op_idx = patch_point_vals_.operations_.size(); // index in operations_ vector
        uint grad_vector_op_idx = patch_point_vals_.get_fe_op(FEOps::opGradVectorShape);
        auto lambda_vector_sym_grad = [vector_sym_grad_op_idx](std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
                bulk_reinit::ptop_vector_sym_grad(operations, vector_sym_grad_op_idx);
            };
        patch_point_vals_.make_fe_op({3,3}, lambda_vector_sym_grad, {grad_vector_op_idx}, fe_component->n_dofs());
        patch_point_vals_.set_fe_op(FEOps::opVectorSymGrad, vector_sym_grad_op_idx);

        return FeQ<Tensor>(patch_point_vals_, vector_sym_grad_op_idx, fe_component->n_dofs());
    }

    /**
     * @brief Return the value of the @p function_no-th vector divergence at
     * the @p p bulk quadrature point.
     *
     * @param component_idx Number of the shape function.
     */
    inline FeQ<Scalar> vector_divergence(uint component_idx=0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        //ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape accessor must be FEScalar!\n");

        // use lambda reinit function
        uint vector_divergence_op_idx = patch_point_vals_.operations_.size(); // index in operations_ vector
        uint grad_vector_op_idx = patch_point_vals_.get_fe_op(FEOps::opGradVectorShape);
        auto lambda_vector_divergence = [vector_divergence_op_idx](std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
                bulk_reinit::ptop_vector_divergence(operations, vector_divergence_op_idx);
            };
        patch_point_vals_.make_fe_op({1}, lambda_vector_divergence, {grad_vector_op_idx}, fe_component->n_dofs());
        patch_point_vals_.set_fe_op(FEOps::opVectorDivergence, vector_divergence_op_idx);

        return FeQ<Scalar>(patch_point_vals_, vector_divergence_op_idx, fe_component->n_dofs());
    }

private:
    /**
     * @brief Precomputed gradients of basis functions at the quadrature points.
     *
     * Dimensions:   (no. of quadrature points)
     *             x (no. of dofs)
     *             x ((dim_ of. ref. cell)x(no. of components in ref. cell))
     */
    std::vector<std::vector<arma::mat> > ref_shape_gradients(std::shared_ptr<FiniteElement<dim>> fe) {
        Quadrature *q = patch_point_vals_.get_quadrature();
    	std::vector<std::vector<arma::mat> > ref_shape_grads( q->size(), vector<arma::mat>(fe->n_dofs()) );

        arma::mat grad(dim, fe->n_components());
        for (unsigned int i_pt=0; i_pt<q->size(); i_pt++)
        {
            for (unsigned int i_dof=0; i_dof<fe->n_dofs(); i_dof++)
            {
                grad.zeros();
                for (unsigned int c=0; c<fe->n_components(); c++)
                    grad.col(c) += fe->shape_grad(i_dof, q->point<dim>(i_pt), c);

                ref_shape_grads[i_pt][i_dof] = grad;
            }
        }

        return ref_shape_grads;
    }

    PatchPointValues<3> &patch_point_vals_;
    std::shared_ptr< FiniteElement<dim> > fe_;
};


template<unsigned int dim>
class SideValues : public BaseValues<dim>
{
public:
	/// Constructor
	SideValues(PatchPointValues<3> &patch_point_vals, MixedPtr<FiniteElement> fe)
	: BaseValues<dim>(), patch_point_vals_(patch_point_vals) {
	    ASSERT_EQ(patch_point_vals.dim(), dim);
	    fe_ = fe[Dim<dim>{}];
	}

    /// Same as BulkValues::JxW but register at side quadrature points.
    inline ElQ<Scalar> JxW()
    {
        return ElQ<Scalar>(patch_point_vals_, FeSide::SideOps::opJxW);
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

    /// Same as BulkValues::scalar_shape but register at side quadrature points.
    inline FeQ<Scalar> scalar_shape(uint component_idx = 0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");

        // use lambda reinit function
        std::vector< std::vector< std::vector<double> > > shape_values(
                dim+1,
                std::vector< std::vector<double> >(patch_point_vals_.get_quadrature()->size(), vector<double>(fe_component->n_dofs()) )
				);
        auto ref_shape_vals = this->ref_shape_values_side(patch_point_vals_.get_quadrature(), fe_component);
        for (unsigned int s=0; s<dim+1; ++s)
            for (unsigned int i = 0; i < patch_point_vals_.get_quadrature()->size(); i++)
                for (unsigned int j = 0; j < fe_component->n_dofs(); j++) {
            	    shape_values[s][i][j] = ref_shape_vals[s][i][j][0];
                }
        uint scalar_shape_op_idx = patch_point_vals_.operations_.size(); // index in operations_ vector
        auto lambda_scalar_shape = [shape_values, scalar_shape_op_idx](std::vector<ElOp<3>> &operations, IntTableArena &el_table) {
                side_reinit::ptop_scalar_shape(operations, el_table, shape_values, scalar_shape_op_idx);
            };
        patch_point_vals_.make_fe_op({1}, lambda_scalar_shape, {}, fe_component->n_dofs());
        patch_point_vals_.set_fe_op(FEOps::opScalarShape, scalar_shape_op_idx);

        return FeQ<Scalar>(patch_point_vals_, scalar_shape_op_idx, fe_component->n_dofs());
    }

    /// Same as BulkValues::grad_scalar_shape but register at side quadrature points.
    inline FeQ<Vector> grad_scalar_shape(uint component_idx=0)
    {
        auto fe_component = this->fe_comp(fe_, component_idx);
        ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape accessor must be FEScalar!\n");

        // use lambda reinit function
        auto ref_shape_grads = this->ref_shape_gradients(fe_component);
        uint scalar_shape_grads_op_idx = patch_point_vals_.operations_.size(); // index in operations_ vector
        auto lambda_scalar_shape_grad = [ref_shape_grads, scalar_shape_grads_op_idx](std::vector<ElOp<3>> &operations, IntTableArena &el_table) {
                side_reinit::ptop_scalar_shape_grads<dim>(operations, el_table, ref_shape_grads, scalar_shape_grads_op_idx);
            };
        patch_point_vals_.make_fe_op({3}, lambda_scalar_shape_grad, {FeSide::SideOps::opElInvJac}, fe_component->n_dofs());
        patch_point_vals_.set_fe_op(FEOps::opGradScalarShape, scalar_shape_grads_op_idx);

        return FeQ<Vector>(patch_point_vals_, scalar_shape_grads_op_idx, fe_component->n_dofs());
    }

private:
    /**
     * @brief Precomputed gradients of basis functions at the quadrature points.
     *
     * Dimensions:   (sides)
     *             x (no. of quadrature points)
     *             x (no. of dofs)
     *             x ((dim_ of. ref. cell)x(no. of components in ref. cell))
     */
    std::vector<std::vector<std::vector<arma::mat> > > ref_shape_gradients(std::shared_ptr<FiniteElement<dim>> fe) {
        Quadrature *q = patch_point_vals_.get_quadrature();
        std::vector<std::vector<std::vector<arma::mat> > > ref_shape_grads( dim+1, std::vector<std::vector<arma::mat> >(q->size(), vector<arma::mat>(fe->n_dofs())) );

        arma::mat grad(dim, fe->n_components());
        for (unsigned int sid=0; sid<dim+1; sid++) {
            auto quad = q->make_from_side<dim>(sid);
            for (unsigned int i_pt=0; i_pt<quad.size(); i_pt++)
            {
                for (unsigned int i_dof=0; i_dof<fe->n_dofs(); i_dof++)
                {
                    grad.zeros();
                    for (unsigned int c=0; c<fe->n_components(); c++)
                        grad.col(c) += fe->shape_grad(i_dof, quad.template point<dim>(i_pt), c);

                    ref_shape_grads[sid][i_pt][i_dof] = grad;
                }
            }
        }

        return ref_shape_grads;
    }

    PatchPointValues<3> &patch_point_vals_;
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

    inline Range< JoinShapeAccessor<Scalar> > scalar_join_shape(uint component_idx = 0)
    {
    	// element of lower dim (bulk points)
        auto fe_component_low = this->fe_comp(fe_low_dim_, component_idx);
        ASSERT_EQ(fe_component_low->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");
        // use lambda reinit function
        std::vector< std::vector<double> > shape_values_bulk( fe_component_low->n_dofs(), vector<double>(patch_point_vals_bulk_->get_quadrature()->size()) );
        auto ref_shape_vals_bulk = this->ref_shape_values_bulk(patch_point_vals_bulk_->get_quadrature(), fe_component_low);
        for (unsigned int i = 0; i < patch_point_vals_bulk_->get_quadrature()->size(); i++)
            for (unsigned int j = 0; j < fe_component_low->n_dofs(); j++) {
            	shape_values_bulk[j][i] = ref_shape_vals_bulk[i][j][0];
            }
        uint scalar_shape_op_idx_bulk = patch_point_vals_bulk_->operations_.size(); // index in operations_ vector
        auto lambda_scalar_shape_bulk = [shape_values_bulk, scalar_shape_op_idx_bulk](std::vector<ElOp<3>> &operations, FMT_UNUSED IntTableArena &el_table) {
                bulk_reinit::ptop_scalar_shape(operations, shape_values_bulk, scalar_shape_op_idx_bulk);
            };
        patch_point_vals_bulk_->make_fe_op({1}, lambda_scalar_shape_bulk, {}, fe_component_low->n_dofs());
        uint op_idx_bulk = patch_point_vals_bulk_->operations_.size()-1;

    	// element of higher dim (side points)
        auto fe_component_high = this->fe_comp(fe_high_dim_, component_idx);
        ASSERT_EQ(fe_component_high->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");
        // use lambda reinit function
        std::vector< std::vector< std::vector<double> > > shape_values_side(
                dim+1,
                std::vector< std::vector<double> >(patch_point_vals_side_->get_quadrature()->size(), vector<double>(fe_component_high->n_dofs()) )
				);
        auto ref_shape_vals_side = this->ref_shape_values_side(patch_point_vals_side_->get_quadrature(), fe_component_high);
        for (unsigned int s=0; s<dim+1; ++s)
            for (unsigned int i = 0; i < patch_point_vals_side_->get_quadrature()->size(); i++)
                for (unsigned int j = 0; j < fe_component_high->n_dofs(); j++) {
            	    shape_values_side[s][i][j] = ref_shape_vals_side[s][i][j][0];
                }
        uint scalar_shape_op_idx_side = patch_point_vals_side_->operations_.size(); // index in operations_ vector
        auto lambda_scalar_shape_side = [shape_values_side, scalar_shape_op_idx_side](std::vector<ElOp<3>> &operations, IntTableArena &el_table) {
                side_reinit::ptop_scalar_shape(operations, el_table, shape_values_side, scalar_shape_op_idx_side);
            };
        patch_point_vals_side_->make_fe_op({1}, lambda_scalar_shape_side, {}, fe_component_high->n_dofs());
        uint op_idx_side = patch_point_vals_side_->operations_.size()-1;

        auto bgn_it = make_iter<JoinShapeAccessor<Scalar>>( JoinShapeAccessor<Scalar>(patch_point_vals_bulk_, patch_point_vals_side_,
                fe_component_low->n_dofs(), fe_component_high->n_dofs(), op_idx_bulk, op_idx_side, 0) );
        unsigned int end_idx = fe_component_low->n_dofs() + fe_component_high->n_dofs();
        auto end_it = make_iter<JoinShapeAccessor<Scalar>>( JoinShapeAccessor<Scalar>(patch_point_vals_bulk_, patch_point_vals_side_,
                fe_component_low->n_dofs(), fe_component_high->n_dofs(), op_idx_bulk, op_idx_side, end_idx) );
        return Range<JoinShapeAccessor<Scalar>>(bgn_it, end_it);
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

    inline Range< JoinShapeAccessor<Scalar> > scalar_join_shape(FMT_UNUSED uint component_idx = 0)
    {
        return Range<JoinShapeAccessor<Scalar>>(
                make_iter<JoinShapeAccessor<Scalar>>(JoinShapeAccessor<Scalar>()),
                make_iter<JoinShapeAccessor<Scalar>>(JoinShapeAccessor<Scalar>()) );
    }
};


template<unsigned int spacedim = 3>
class PatchFEValues {
public:
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
    : asm_arena_(1024 * 1024, 256),
      patch_arena_(nullptr),
      patch_point_vals_bulk_{ {FeBulk::PatchPointValues(1, 0, asm_arena_),
                               FeBulk::PatchPointValues(2, 0, asm_arena_),
                               FeBulk::PatchPointValues(3, 0, asm_arena_)} },
      patch_point_vals_side_{ {FeSide::PatchPointValues(1, 0, asm_arena_),
                               FeSide::PatchPointValues(2, 0, asm_arena_),
                               FeSide::PatchPointValues(3, 0, asm_arena_)} }
    {
        used_quads_[0] = false; used_quads_[1] = false;
    }

    PatchFEValues(unsigned int quad_order, MixedPtr<FiniteElement> fe)
    : asm_arena_(1024 * 1024, 256),
      patch_arena_(nullptr),
      patch_point_vals_bulk_{ {FeBulk::PatchPointValues(1, quad_order, asm_arena_),
    	                       FeBulk::PatchPointValues(2, quad_order, asm_arena_),
                               FeBulk::PatchPointValues(3, quad_order, asm_arena_)} },
      patch_point_vals_side_{ {FeSide::PatchPointValues(1, quad_order, asm_arena_),
                               FeSide::PatchPointValues(2, quad_order, asm_arena_),
                               FeSide::PatchPointValues(3, quad_order, asm_arena_)} },
      fe_(fe)
    {
        used_quads_[0] = false; used_quads_[1] = false;
    }


    /// Destructor
    ~PatchFEValues()
    {
        if (patch_arena_!=nullptr)
            delete patch_arena_;
    }

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
        patch_arena_ = asm_arena_.get_child_arena();
        for (unsigned int i=0; i<3; ++i) {
            if (used_quads_[0]) patch_point_vals_bulk_[i].init_finalize(patch_arena_);
            if (used_quads_[1]) patch_point_vals_side_[i].init_finalize(patch_arena_);
        }
    }

    /// Reset PatchpointValues structures
    void reset()
    {
        for (unsigned int i=0; i<3; ++i) {
            if (used_quads_[0]) patch_point_vals_bulk_[i].reset();
            if (used_quads_[1]) patch_point_vals_side_[i].reset();
        }
        patch_arena_->reset();
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
        return BulkValues<dim>(patch_point_vals_bulk_[dim-1], fe_);
    }

    /// Return SideValue object of dimension given by template parameter
    template<unsigned int dim>
    SideValues<dim> side_values() {
    	ASSERT((dim>0) && (dim<=3))(dim).error("Dimension must be 1, 2 or 3.");
        return SideValues<dim>(patch_point_vals_side_[dim-1], fe_);
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
    AssemblyArena asm_arena_;
    PatchArena *patch_arena_;
    std::array<FeBulk::PatchPointValues<spacedim>, 3> patch_point_vals_bulk_;  ///< Sub objects of bulk data of dimensions 1,2,3
    std::array<FeSide::PatchPointValues<spacedim>, 3> patch_point_vals_side_;  ///< Sub objects of side data of dimensions 1,2,3

    MixedPtr<FiniteElement> fe_;   ///< Mixed of shared pointers of FiniteElement object
    bool used_quads_[2];           ///< Pair of flags signs holds info if bulk and side quadratures are used

    template <class ValueType>
    friend class ElQ;
    template <class ValueType>
    friend class FeQ;
    template <class ValueType>
    friend class JoinShapeAccessor;
};


template <class ValueType>
ValueType ElQ<ValueType>::operator()(const BulkPoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.scalar_value(op_idx_, value_cache_idx);
}

template <>
inline Vector ElQ<Vector>::operator()(const BulkPoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.vector_value(op_idx_, value_cache_idx);
}

template <>
inline Tensor ElQ<Tensor>::operator()(const BulkPoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.tensor_value(op_idx_, value_cache_idx);
}

template <class ValueType>
ValueType ElQ<ValueType>::operator()(const SidePoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.scalar_value(op_idx_, value_cache_idx);
}

template <>
inline Vector ElQ<Vector>::operator()(const SidePoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.vector_value(op_idx_, value_cache_idx);
}

template <>
inline Tensor ElQ<Tensor>::operator()(const SidePoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.tensor_value(op_idx_, value_cache_idx);
}

template <class ValueType>
ValueType FeQ<ValueType>::operator()(unsigned int shape_idx, const BulkPoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.scalar_value(op_idx_, value_cache_idx, shape_idx);
}

template <>
inline Vector FeQ<Vector>::operator()(unsigned int shape_idx, const BulkPoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.vector_value(op_idx_, value_cache_idx, shape_idx);
}

template <>
inline Tensor FeQ<Tensor>::operator()(unsigned int shape_idx, const BulkPoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.tensor_value(op_idx_, value_cache_idx, shape_idx);
}

template <class ValueType>
ValueType FeQ<ValueType>::operator()(unsigned int shape_idx, const SidePoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.scalar_value(op_idx_, value_cache_idx, shape_idx);
}

template <>
inline Vector FeQ<Vector>::operator()(unsigned int shape_idx, const SidePoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.vector_value(op_idx_, value_cache_idx, shape_idx);
}

template <>
inline Tensor FeQ<Tensor>::operator()(FMT_UNUSED unsigned int shape_idx, FMT_UNUSED const SidePoint &point) {
    unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
    return patch_point_vals_.tensor_value(op_idx_, value_cache_idx, shape_idx);
}


template <class ValueType>
ValueType JoinShapeAccessor<ValueType>::operator()(const BulkPoint &point) {
    if (this->is_high_dim()) {
        return 0.0;
    } else {
        unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
        return patch_point_vals_bulk_->scalar_value(op_idx_bulk_, value_cache_idx, this->local_idx());
    }
}

template <>
inline Vector JoinShapeAccessor<Vector>::operator()(FMT_UNUSED const BulkPoint &point) {
    Vector vect; vect.zeros();
    return vect;
}

template <>
inline Tensor JoinShapeAccessor<Tensor>::operator()(FMT_UNUSED const BulkPoint &point) {
	Tensor tens; tens.zeros();
    return tens;
}

template <class ValueType>
ValueType JoinShapeAccessor<ValueType>::operator()(const SidePoint &point) {
    if (this->is_high_dim()) {
        unsigned int value_cache_idx = point.elm_cache_map()->element_eval_point(point.elem_patch_idx(), point.eval_point_idx());
        return patch_point_vals_side_->scalar_value(op_idx_side_, value_cache_idx, this->local_idx());
    } else {
        return 0.0;
    }
}

template <>
inline Vector JoinShapeAccessor<Vector>::operator()(FMT_UNUSED const SidePoint &point) {
    Vector vect; vect.zeros();
    return vect;
}

template <>
inline Tensor JoinShapeAccessor<Tensor>::operator()(FMT_UNUSED const SidePoint &point) {
	Tensor tens; tens.zeros();
    return tens;
}


#endif /* PATCH_FE_VALUES_HH_ */
