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
 * @file    op_function.hh
 * @brief   Store finite element reinit functions.
 * @author  David Flanderka
 */

#ifndef OP_FUNCTION_HH_
#define OP_FUNCTION_HH_

#include "fem/patch_op.hh"
#include "fem/patch_fe_values.hh"



namespace Op {

/// Base class of bulk and side vector symmetric gradients
template<unsigned int dim, unsigned int spacedim>
class VectorSymGradBase : public PatchOp<spacedim> {
public:
    /// Constructor
    VectorSymGradBase(uint _dim, PatchFEValues<spacedim> &pfev, std::initializer_list<uint> shape, OpSizeType size_type)
    : PatchOp<spacedim>(_dim, pfev, shape, size_type)
    {
        ASSERT_EQ(this->dim_, dim);
    }

    void eval() override {
        for (uint i_dof=0; i_dof<this->n_dofs(); ++i_dof) {
            Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > grad_vector_dof = this->input_ops(0)->result_sub_matrix(i_dof);
            Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > sym_grad_dof = this->result_sub_matrix(i_dof);
            sym_grad_dof = 0.5 * (grad_vector_dof.transpose() + grad_vector_dof);
        }
    }
};

/// Base class of bulk and side vector divergence
template<unsigned int dim, unsigned int spacedim>
class VectorDivergenceBase : public PatchOp<spacedim> {
public:
    /// Constructor
    VectorDivergenceBase(uint _dim, PatchFEValues<spacedim> &pfev, std::initializer_list<uint> shape, OpSizeType size_type)
    : PatchOp<spacedim>(_dim, pfev, shape, size_type)
    {
        ASSERT_EQ(this->dim_, dim);
    }

    void eval() override {
        auto divergence_value = this->result_matrix();

        for (uint i_dof=0; i_dof<this->n_dofs(); ++i_dof) {
            Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > grad_vector_dof = this->input_ops(0)->result_sub_matrix(i_dof);
            divergence_value(i_dof) = grad_vector_dof(0,0) + grad_vector_dof(1,1) + grad_vector_dof(2,2);
        }
    }
};


namespace Bulk {

/// Common ancestor of all bulk operations.
template<unsigned int spacedim>
class Base : public PatchOp<spacedim> {
public:
    /// Constructor
	Base(uint dim, PatchFEValues<spacedim> &pfev, std::initializer_list<uint> shape, OpSizeType size_type)
    : PatchOp<spacedim>(dim, pfev, shape, size_type)
    {
        this->bulk_side_ = 0;
    }
};

namespace El {

template<unsigned int spacedim>
class OpCoords : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
    OpCoords(uint dim, PatchFEValues<spacedim> &pfev)
    : Op::Bulk::Base<spacedim>(dim, pfev, {spacedim, dim+1}, OpSizeType::elemOp) {}

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        this->allocate_result( ppv.n_elems_, *ppv.patch_fe_data_.patch_arena_ );
        auto result = this->result_matrix();

        for (uint i_elm=0; i_elm<ppv.elem_list_.size(); ++i_elm)
            for (uint i_col=0; i_col<this->dim_+1; ++i_col)
                for (uint i_row=0; i_row<spacedim; ++i_row) {
                    result(i_row, i_col)(i_elm) = ( *ppv.elem_list_[i_elm].node(i_col) )(i_row);
                }
    }
};

template<unsigned int spacedim>
class OpJac : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
    OpJac(uint dim, PatchFEValues<spacedim> &pfev)
    : Op::Bulk::Base<spacedim>(dim, pfev, {spacedim, dim}, OpSizeType::elemOp)
    {
        this->input_ops_.push_back( pfev.template get< OpCoords<spacedim> >(dim) );
    }

    void eval() override {
        auto jac_value = this->result_matrix();
        auto coords_value = this->input_ops(0)->result_matrix();
        for (unsigned int i=0; i<spacedim; i++)
            for (unsigned int j=0; j<this->dim_; j++)
                jac_value(i,j) = coords_value(i,j+1) - coords_value(i,0);
    }
};

template<unsigned int dim, unsigned int spacedim>
class OpInvJac : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
    OpInvJac(uint _dim, PatchFEValues<spacedim> &pfev)
    : Op::Bulk::Base<spacedim>(_dim, pfev, {dim, spacedim}, OpSizeType::elemOp)
    {
        ASSERT_EQ(this->dim_, dim);
        this->input_ops_.push_back( pfev.template get< OpJac<spacedim> >(dim) );
    }

    void eval() override {
        auto inv_jac_value = this->result_matrix();
        auto jac_value = this->input_ops(0)->result_matrix();
        inv_jac_value = eigen_arena_tools::inverse<spacedim, dim>(jac_value);
    }
};

template<unsigned int dim, unsigned int spacedim>
class OpJacDet : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
	OpJacDet(uint _dim, PatchFEValues<spacedim> &pfev)
	: Op::Bulk::Base<spacedim>(_dim, pfev, {1}, OpSizeType::elemOp)
	{
	    ASSERT_EQ(this->dim_, dim);
	    this->input_ops_.push_back( pfev.template get< OpJac<spacedim> >(dim) );
	}

    void eval() override {
        auto jac_det_value = this->result_matrix();
        auto jac_value = this->input_ops(0)->result_matrix();
        jac_det_value(0) = eigen_arena_tools::determinant<spacedim, dim>(jac_value).abs();
    }
};

} // end of namespace Op::Bulk::El

namespace Pt {

/// Fixed operation points weights
template<unsigned int spacedim>
class OpWeights : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
    OpWeights(uint dim, PatchFEValues<spacedim> &pfev)
    : Op::Bulk::Base<spacedim>(dim, pfev, {1}, OpSizeType::fixedSizeOp)
    {
        // create result vector of weights operation in assembly arena
        const std::vector<double> &point_weights_vec = pfev.get_bulk_quadrature(dim)->get_weights();
        this->create_result();
        this->allocate_result(point_weights_vec.size(), pfev.asm_arena());
        for (uint i=0; i<point_weights_vec.size(); ++i)
            this->result_(0)(i) = point_weights_vec[i];
    }

    void eval() override {}
};

/// Evaluates coordinates of quadrature points
template<unsigned int spacedim>
class OpCoords : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
    OpCoords(uint dim, PatchFEValues<spacedim> &pfev)
    : Op::Bulk::Base<spacedim>(dim, pfev, {spacedim}, OpSizeType::pointOp){}

    void eval() override {}
};

/// Evaluates JxW on quadrature points
template<unsigned int dim, unsigned int spacedim>
class OpJxW : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
    OpJxW(uint _dim, PatchFEValues<spacedim> &pfev)
    : Op::Bulk::Base<spacedim>(_dim, pfev, {1}, OpSizeType::pointOp)
    {
        this->input_ops_.push_back( pfev.template get< OpWeights<spacedim> >(dim) );
        this->input_ops_.push_back( pfev.template get< Op::Bulk::El::OpJacDet<dim, spacedim> >(dim) );
    }

    void eval() override {
        auto weights_value = this->input_ops(0)->result_matrix();
        auto jac_det_value = this->input_ops(1)->result_matrix();
        ArenaOVec<double> weights_ovec( weights_value(0,0) );
        ArenaOVec<double> jac_det_ovec( jac_det_value(0,0) );
        ArenaOVec<double> jxw_ovec = jac_det_ovec * weights_ovec;
        this->result_(0) = jxw_ovec.get_vec();
    }
};

/// Fixed operation of  scalar shape reference values
template<unsigned int dim, unsigned int spacedim>
class OpRefScalar : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
    OpRefScalar(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Bulk::Base<spacedim>(_dim, pfev, {1}, OpSizeType::fixedSizeOp)
    {
        ASSERT_EQ(this->dim_, dim);

        uint n_points = pfev.get_bulk_quadrature(dim)->size();
        uint n_dofs = fe->n_dofs();
        this->n_dofs_ = n_dofs;

        auto ref_shape_vals = this->template ref_shape_values_bulk<dim>(pfev.get_bulk_quadrature(dim), fe);
        this->create_result();
        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_scalar_value = this->result_matrix();
        for (unsigned int i_p = 0; i_p < n_points; i_p++)
            for (unsigned int i_dof = 0; i_dof < n_dofs; i_dof++) {
                ref_scalar_value(i_dof)(i_p) = ref_shape_vals[i_p][i_dof][0];
            }
    }

    void eval() override {}
};

/// Fixed operation of  scalar shape reference values
template<unsigned int dim, unsigned int spacedim>
class OpRefVector : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
	OpRefVector(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
	: Op::Bulk::Base<spacedim>(_dim, pfev, {spacedim}, OpSizeType::fixedSizeOp)
	{
        ASSERT_EQ(this->dim_, dim);

        uint n_points = pfev.get_bulk_quadrature(dim)->size();
        uint n_dofs = fe->n_dofs();
        this->n_dofs_ = n_dofs;

        auto ref_shape_vals = this->template ref_shape_values_bulk<dim>(pfev.get_bulk_quadrature(dim), fe);
        this->create_result();
        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_vector_value = this->result_matrix();

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
            for (uint i_p=0; i_p<n_points; ++i_p)
                for (uint c=0; c<spacedim; ++c)
                    ref_vector_value(c, i_dof)(i_p) = ref_shape_vals[i_p][i_dof](c);
	}

    void eval() override {}
};

/// Fixed operation of gradient scalar reference values
template<unsigned int dim, unsigned int spacedim>
class OpRefGradScalar : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
    OpRefGradScalar(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Bulk::Base<spacedim>(_dim, pfev, {dim, 1}, OpSizeType::fixedSizeOp)
    {
        ASSERT_EQ(this->dim_, dim);

        uint n_points = pfev.get_bulk_quadrature(dim)->size();
        uint n_dofs = fe->n_dofs();
        this->n_dofs_ = n_dofs;

        std::vector<std::vector<arma::mat> > ref_shape_grads = this->template ref_shape_gradients_bulk<dim>(pfev.get_bulk_quadrature(dim), fe);
        this->create_result();
        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_scalar_value = this->result_matrix();
        for (uint i_row=0; i_row<ref_scalar_value.rows(); ++i_row) {
            for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
                for (uint i_p=0; i_p<n_points; ++i_p)
                    ref_scalar_value(i_row, i_dof)(i_p) = ref_shape_grads[i_p][i_dof](i_row);
        }
    }

    void eval() override {}
};

/// Fixed operation of gradient scalar reference values
template<unsigned int dim, unsigned int spacedim>
class OpRefGradVector : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
	OpRefGradVector(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Bulk::Base<spacedim>(_dim, pfev, {dim, spacedim}, OpSizeType::fixedSizeOp)
    {
        ASSERT_EQ(this->dim_, dim);

        uint n_points = pfev.get_bulk_quadrature(dim)->size();
        uint n_dofs = fe->n_dofs();
        this->n_dofs_ = n_dofs;

        std::vector<std::vector<arma::mat> > ref_shape_grads = this->template ref_shape_gradients_bulk<dim>(pfev.get_bulk_quadrature(dim), fe);
        this->create_result();
        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_vector_value = this->result_matrix();
        for (uint i_c=0; i_c<spacedim; ++i_c) {
            for (uint i_dim=0; i_dim<dim; ++i_dim)
                for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
                    for (uint i_p=0; i_p<n_points; ++i_p)
                        ref_vector_value(i_dim,3*i_dof+i_c)(i_p) = ref_shape_grads[i_p][i_dof](i_dim, i_c);
        }
    }

    void eval() override {}
};

/// Evaluates scalar values
template<unsigned int dim, unsigned int spacedim>
class OpScalarShape : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
	OpScalarShape(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
	: Op::Bulk::Base<spacedim>(_dim, pfev, {1}, OpSizeType::pointOp)
	{
	    ASSERT_EQ(this->dim_, dim);

	    this->n_dofs_ = fe->n_dofs();
	    this->input_ops_.push_back( pfev.template get< OpRefScalar<dim, spacedim>, dim >(fe) );
	}

    void eval() override {
        auto ref_vec = this->input_ops(0)->result_matrix();
        auto result_vec = this->result_matrix();

        uint n_dofs = this->n_dofs();
        uint n_elem = this->ppv().n_elems_;

        ArenaVec<double> elem_vec(n_elem, this->patch_fe_->patch_arena());
        for (uint i=0; i<n_elem; ++i) {
            elem_vec(i) = 1.0;
        }
        ArenaOVec<double> elem_ovec(elem_vec);

        Eigen::Vector<ArenaOVec<double>, Eigen::Dynamic> ref_ovec(n_dofs);
        for (uint i=0; i<n_dofs; ++i) {
            ref_ovec(i) = ArenaOVec<double>( ref_vec(i) );
        }

        Eigen::Vector<ArenaOVec<double>, Eigen::Dynamic> result_ovec = elem_ovec * ref_ovec;
        for (uint i=0; i<n_dofs; ++i) {
            result_vec(i) = result_ovec(i).get_vec();
        }
    }
};

/// Evaluates vector values (FEType == FEVector)
template<unsigned int dim, unsigned int spacedim>
class OpVectorShape : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
	OpVectorShape(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : Op::Bulk::Base<spacedim>(_dim, pfev, {spacedim}, OpSizeType::pointOp), dispatch_op_(dispatch_op)
    {
        ASSERT_EQ(this->dim_, dim);

        this->n_dofs_ = fe->n_dofs();
        this->input_ops_.push_back( pfev.template get< OpRefVector<dim, spacedim>, dim >(fe) );
	}

    void eval() override {
        auto ref_shape_vec = this->input_ops(0)->result_matrix();
        auto result_vec = dispatch_op_.result_matrix();

        uint n_dofs = this->n_dofs();
        uint n_elem = this->ppv().n_elems_;

        ArenaVec<double> elem_vec(n_elem, this->patch_fe_->patch_arena());
        for (uint i=0; i<n_elem; ++i) {
            elem_vec(i) = 1.0;
        }
        ArenaOVec<double> elem_ovec(elem_vec);

        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> ref_shape_ovec(3, n_dofs);
        for (uint c=0; c<3*n_dofs; ++c) {
            ref_shape_ovec(c) = ArenaOVec(ref_shape_vec(c));
        }

        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> result_ovec = elem_ovec * ref_shape_ovec;
        for (uint c=0; c<3*n_dofs; ++c)
            result_vec(c) = result_ovec(c).get_vec();
    }

private:
    PatchOp<spacedim> &dispatch_op_;
};

// class OpVectorCovariantShape
// class OpVectorPiolaShape

/// Dispatch class of vector values
template<unsigned int dim, unsigned int spacedim>
class DispatchVectorShape : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
    DispatchVectorShape(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Bulk::Base<spacedim>(_dim, pfev, {spacedim}, OpSizeType::pointOp), in_op_(nullptr)
    {
        ASSERT_EQ(this->dim_, dim);
        this->n_dofs_ = fe->n_dofs();

        switch (fe->fe_type()) {
            case FEVector:
            {
                in_op_ = new OpVectorShape<dim, spacedim>(_dim, pfev, fe, *this);
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpVectorCovariantShape<dim, spacedim>(_dim, pfev, fe, *this);
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpVectorPiolaShape<dim, spacedim>(_dim, pfev, fe, *this);
                break;
            }
            default:
                ASSERT(false).error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");
        }

	}

    void eval() override {
        in_op_->eval();
    }

private:
    PatchOp<spacedim> *in_op_;
};

/// Evaluates gradient scalar values
template<unsigned int dim, unsigned int spacedim>
class OpGradScalarShape : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
	OpGradScalarShape(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
	: Op::Bulk::Base<spacedim>(_dim, pfev, {spacedim, 1}, OpSizeType::pointOp)
	{
	    ASSERT_EQ(this->dim_, dim);

	    this->n_dofs_ = fe->n_dofs();
	    this->input_ops_.push_back( pfev.template get< Op::Bulk::El::OpInvJac<dim, spacedim> >(dim) );
	    this->input_ops_.push_back( pfev.template get< OpRefGradScalar<dim, spacedim>, dim >(fe) );
	}

    void eval() override {
        auto inv_jac_vec = this->input_ops(0)->result_matrix();    // dim x spacedim=3
        auto ref_grads_vec = this->input_ops(1)->result_matrix();  // dim x n_dofs

        uint n_dofs = this->n_dofs();

        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> ref_grads_ovec(this->dim_, n_dofs);
        for (uint i=0; i<this->dim_*n_dofs; ++i) {
            ref_grads_ovec(i) = ArenaOVec(ref_grads_vec(i));
        }

        Eigen::Matrix<ArenaOVec<double>, dim, 3> inv_jac_ovec;
        for (uint i=0; i<this->dim_*spacedim; ++i) {
            inv_jac_ovec(i) = ArenaOVec(inv_jac_vec(i));
        }

        auto result_vec = this->result_matrix();
        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> result_ovec = inv_jac_ovec.transpose() * ref_grads_ovec;
        for (uint i=0; i<spacedim*n_dofs; ++i) {
            result_vec(i) = result_ovec(i).get_vec();
        }
    }
};

/// Evaluates vector values (FEType == FEVector)
template<unsigned int dim, unsigned int spacedim>
class OpGradVectorShape : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
	OpGradVectorShape(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : Op::Bulk::Base<spacedim>(_dim, pfev, {spacedim, spacedim}, OpSizeType::pointOp), dispatch_op_(dispatch_op)
    {
        ASSERT_EQ(this->dim_, dim);

        this->n_dofs_ = fe->n_dofs();
        this->input_ops_.push_back( pfev.template get< Op::Bulk::El::OpInvJac<dim, spacedim> >(dim) );
        this->input_ops_.push_back( pfev.template get< OpRefGradVector<dim, spacedim>, dim >(fe) );
	}

    void eval() override {
        auto inv_jac_vec = this->input_ops(0)->result_matrix();    // dim x spacedim
        auto ref_grads_vec = this->input_ops(1)->result_matrix();  // dim x spacedim
        auto result_vec = dispatch_op_.result_matrix();            // spacedim x spacedim

        uint n_dofs = this->n_dofs();

        Eigen::Matrix<ArenaOVec<double>, dim, 3> inv_jac_ovec;
        for (uint i=0; i<dim*3; ++i) {
            inv_jac_ovec(i) = ArenaOVec(inv_jac_vec(i));
        }

        Eigen::Matrix<ArenaOVec<double>, dim, 3> ref_grads_ovec;
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i=0; i<3*dim; ++i) {
                ref_grads_ovec(i) = ArenaOVec(ref_grads_vec(i_dof*3*dim + i));
            }

            Eigen::Matrix<ArenaOVec<double>, 3, 3> result_ovec = inv_jac_ovec.transpose() * ref_grads_ovec;
            for (uint i=0; i<3; ++i) {
                for (uint j=0; j<3; ++j) {
                   result_vec(j,i+3*i_dof) = result_ovec(i,j).get_vec();
                }
            }
        }
    }

private:
    PatchOp<spacedim> &dispatch_op_;
};

// class OpGradVectorCovariantShape
// class OpGradVectorPiolaShape

/// Dispatch class of vector values
template<unsigned int dim, unsigned int spacedim>
class DispatchGradVectorShape : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
	DispatchGradVectorShape(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Bulk::Base<spacedim>(_dim, pfev, {spacedim, spacedim}, OpSizeType::pointOp), in_op_(nullptr)
    {
        ASSERT_EQ(this->dim_, dim);
        this->n_dofs_ = fe->n_dofs();

        switch (fe->fe_type()) {
            case FEVector:
            {
                in_op_ = new OpGradVectorShape<dim, spacedim>(_dim, pfev, fe, *this);
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpGradVectorCovariantShape<dim, spacedim>(_dim, pfev, fe, *this);
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpGradVectorPiolaShape<dim, spacedim>(_dim, pfev, fe, *this);
                break;
            }
            default:
                ASSERT(false).error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");
        }

	}

    void eval() override {
        in_op_->eval();
    }

private:
    PatchOp<spacedim> *in_op_;
};

/// Evaluates vector symmetric gradients
template<unsigned int dim, unsigned int spacedim>
class OpVectorSymGrad : public Op::VectorSymGradBase<dim, spacedim> {
public:
    /// Constructor
	OpVectorSymGrad(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::VectorSymGradBase<dim, spacedim>(_dim, pfev, {spacedim, spacedim}, OpSizeType::pointOp)
    {
        this->bulk_side_ = 0;
        this->n_dofs_ = fe->n_dofs();
        this->input_ops_.push_back( pfev.template get< DispatchGradVectorShape<dim, spacedim>, dim >(fe) );
    }
};

/// Evaluates vector divergence
template<unsigned int dim, unsigned int spacedim>
class OpVectorDivergence : public Op::VectorDivergenceBase<dim, spacedim> {
public:
    /// Constructor
	OpVectorDivergence(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::VectorDivergenceBase<dim, spacedim>(_dim, pfev, {1}, OpSizeType::pointOp)
    {
        this->bulk_side_ = 0;
        this->n_dofs_ = fe->n_dofs();
        this->input_ops_.push_back( pfev.template get< DispatchGradVectorShape<dim, spacedim>, dim >(fe) );
    }
};

} // end of namespace Op::Bulk::Pt

} // end of namespace Op::Bulk

namespace Side {

/// Common ancestor of all side operations.
template<unsigned int spacedim>
class Base : public PatchOp<spacedim> {
public:
    /// Constructor
	Base(uint dim, PatchFEValues<spacedim> &pfev, std::initializer_list<uint> shape, OpSizeType size_type)
    : PatchOp<spacedim>(dim, pfev, shape, size_type)
    {
        this->bulk_side_ = 1;
    }
};

namespace El {

template<unsigned int spacedim>
class OpElCoords : public Op::Side::Base<spacedim> {
public:
    /// Constructor
    OpElCoords(uint dim, PatchFEValues<spacedim> &pfev)
    : Op::Side::Base<spacedim>(dim, pfev, {spacedim, dim+1}, OpSizeType::elemOp) {}

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        this->allocate_result( ppv.n_elems_, *ppv.patch_fe_data_.patch_arena_ );
        auto result = this->result_matrix();

        for (uint i_elm=0; i_elm<ppv.elem_list_.size(); ++i_elm)
            for (uint i_col=0; i_col<this->dim_+1; ++i_col)
                for (uint i_row=0; i_row<spacedim; ++i_row) {
                    result(i_row, i_col)(i_elm) = ( *ppv.elem_list_[i_elm].node(i_col) )(i_row);
                }
    }

};

template<unsigned int spacedim>
class OpElJac : public Op::Side::Base<spacedim> {
public:
    /// Constructor
    OpElJac(uint dim, PatchFEValues<spacedim> &pfev)
    : Op::Side::Base<spacedim>(dim, pfev, {spacedim, dim}, OpSizeType::elemOp)
    {
        this->input_ops_.push_back( pfev.template get< OpElCoords<spacedim> >(dim) );
    }

    void eval() override {
        auto jac_value = this->result_matrix();
        auto coords_value = this->input_ops(0)->result_matrix();
        for (unsigned int i=0; i<spacedim; i++)
            for (unsigned int j=0; j<this->dim_; j++)
                jac_value(i,j) = coords_value(i,j+1) - coords_value(i,0);
    }
};

template<unsigned int dim, unsigned int spacedim>
class OpElInvJac : public Op::Side::Base<spacedim> {
public:
    /// Constructor
	OpElInvJac(uint _dim, PatchFEValues<spacedim> &pfev)
    : Op::Side::Base<spacedim>(_dim, pfev, {dim, spacedim}, OpSizeType::elemOp)
    {
        ASSERT_EQ(this->dim_, dim);
        this->input_ops_.push_back( pfev.template get< OpElJac<spacedim> >(dim) );
    }

    void eval() override {
        auto inv_jac_value = this->result_matrix();
        auto jac_value = this->input_ops(0)->result_matrix();
        inv_jac_value = eigen_arena_tools::inverse<spacedim, dim>(jac_value);
    }
};

template<unsigned int spacedim>
class OpSideCoords : public Op::Side::Base<spacedim> {
public:
    /// Constructor
	OpSideCoords(uint dim, PatchFEValues<spacedim> &pfev)
    : Op::Side::Base<spacedim>(dim, pfev, {spacedim, dim}, OpSizeType::elemOp) {}

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        this->allocate_result( ppv.n_elems_, *ppv.patch_fe_data_.patch_arena_ );
        auto result = this->result_matrix();

        for (uint i_side=0; i_side<ppv.side_list_.size(); ++i_side)
            for (uint i_col=0; i_col<this->dim_; ++i_col)
                for (uint i_row=0; i_row<spacedim; ++i_row) {
                    result(i_row, i_col)(i_side) = ( *ppv.side_list_[i_side].node(i_col) )(i_row);
                }
    }

};

template<unsigned int spacedim>
class OpSideJac : public Op::Side::Base<spacedim> {
public:
    /// Constructor
	OpSideJac(uint dim, PatchFEValues<spacedim> &pfev)
    : Op::Side::Base<spacedim>(dim, pfev, {spacedim, dim-1}, OpSizeType::elemOp)
    {
        this->input_ops_.push_back( pfev.template get< OpSideCoords<spacedim> >(dim) );
    }

    void eval() override {
        auto jac_value = this->result_matrix();
        auto coords_value = this->input_ops(0)->result_matrix();
        for (unsigned int i=0; i<spacedim; i++)
            for (unsigned int j=0; j<this->dim_-1; j++)
                jac_value(i,j) = coords_value(i,j+1) - coords_value(i,0);
    }
};

template<unsigned int dim, unsigned int spacedim>
class OpSideJacDet : public Op::Side::Base<spacedim> {
public:
    /// Constructor
	OpSideJacDet(uint _dim, PatchFEValues<spacedim> &pfev)
	: Op::Side::Base<spacedim>(_dim, pfev, {1}, OpSizeType::elemOp)
	{
	    ASSERT_EQ(this->dim_, dim);
	    this->input_ops_.push_back( pfev.template get< OpSideJac<spacedim> >(dim) );
	}

    void eval() override {
        auto jac_det_value = this->result_matrix();
        auto jac_value = this->input_ops(0)->result_matrix();
        jac_det_value(0) = eigen_arena_tools::determinant<spacedim, dim-1>(jac_value).abs();
    }
};

/// Template specialization of previous: dim=1
template<unsigned int spacedim>
class OpSideJacDet<1, spacedim> : public Op::Side::Base<spacedim> {
public:
    /// Constructor
	OpSideJacDet(uint _dim, PatchFEValues<spacedim> &pfev)
	: Op::Side::Base<spacedim>(_dim, pfev, {1}, OpSizeType::elemOp)
	{
	    ASSERT_EQ(this->dim_, 1);
	}

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        this->allocate_result( ppv.n_elems_, *ppv.patch_fe_data_.patch_arena_ );
        auto jac_det_value = this->result_matrix();
        for (uint i=0;i<ppv.n_elems_; ++i) {
            jac_det_value(0,0)(i) = 1.0;
        }
    }
};

} // end of namespace Op::Side::El

namespace Pt {

/// Fixed operation points weights
template<unsigned int spacedim>
class OpWeights : public Op::Side::Base<spacedim> {
public:
    /// Constructor
    OpWeights(uint dim, PatchFEValues<spacedim> &pfev)
    : Op::Side::Base<spacedim>(dim, pfev, {1}, OpSizeType::fixedSizeOp)
    {
        // create result vector of weights operation in assembly arena
        const std::vector<double> &point_weights_vec = pfev.get_side_quadrature(dim)->get_weights();
        this->create_result();
        this->allocate_result(point_weights_vec.size(), pfev.asm_arena());
        for (uint i=0; i<point_weights_vec.size(); ++i)
            this->result_(0)(i) = point_weights_vec[i];
    }

    void eval() override {}
};

/// Evaluates coordinates of quadrature points
template<unsigned int spacedim>
class OpCoords : public Op::Side::Base<spacedim> {
public:
    /// Constructor
    OpCoords(uint dim, PatchFEValues<spacedim> &pfev)
    : Op::Side::Base<spacedim>(dim, pfev, {spacedim}, OpSizeType::pointOp){}

    void eval() override {}
};

/// Evaluates JxW on quadrature points
template<unsigned int dim, unsigned int spacedim>
class OpJxW : public Op::Side::Base<spacedim> {
public:
    /// Constructor
    OpJxW(uint _dim, PatchFEValues<spacedim> &pfev)
    : Op::Side::Base<spacedim>(_dim, pfev, {1}, OpSizeType::pointOp)
    {
        this->input_ops_.push_back( pfev.template get< OpWeights<spacedim> >(dim) );
        this->input_ops_.push_back( pfev.template get< Op::Side::El::OpSideJacDet<dim, spacedim> >(dim) );
    }

    void eval() override {
        auto weights_value = this->input_ops(0)->result_matrix();
        auto jac_det_value = this->input_ops(1)->result_matrix();
        ArenaOVec<double> weights_ovec( weights_value(0,0) );
        ArenaOVec<double> jac_det_ovec( jac_det_value(0,0) );
        ArenaOVec<double> jxw_ovec = jac_det_ovec * weights_ovec;
        this->result_(0) = jxw_ovec.get_vec();
    }
};

/// Evaluates normal vector on quadrature points
template<unsigned int dim, unsigned int spacedim>
class OpNormalVec : public Op::Side::Base<spacedim> {
public:
    /// Constructor
    OpNormalVec(uint _dim, PatchFEValues<spacedim> &pfev)
    : Op::Side::Base<spacedim>(_dim, pfev, {spacedim}, OpSizeType::elemOp)
    {
        this->input_ops_.push_back( pfev.template get< Op::Side::El::OpElInvJac<dim, spacedim> >(dim) );
    }

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        auto normal_value = this->result_matrix();
        auto inv_jac_value = this->input_ops(0)->result_matrix();
        normal_value = inv_jac_value.transpose() * RefElement<dim>::normal_vector_array( ppv.int_table_(3) );

        ArenaVec<double> norm_vec( ppv.n_elems_, *ppv.patch_fe_data_.patch_arena_ );
        Eigen::VectorXd A(3);
        for (uint i=0; i<normal_value(0).data_size(); ++i) {
            A(0) = normal_value(0)(i);
            A(1) = normal_value(1)(i);
            A(2) = normal_value(2)(i);
            norm_vec(i) = A.norm();
        }

        for (uint i=0; i<3; ++i) {
            normal_value(i) = normal_value(i) / norm_vec;
        }
    }
};

/// Fixed operation of  scalar shape reference values
template<unsigned int dim, unsigned int spacedim>
class OpRefScalar : public Op::Side::Base<spacedim> {
public:
    /// Constructor
    OpRefScalar(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Side::Base<spacedim>(_dim, pfev, {dim+1}, OpSizeType::fixedSizeOp)
    {
        ASSERT_EQ(this->dim_, dim);

        uint n_points = pfev.get_side_quadrature(dim)->size();
        uint n_dofs = fe->n_dofs();
        this->n_dofs_ = n_dofs;

        auto ref_shape_vals = this->template ref_shape_values_side<dim>(pfev.get_side_quadrature(dim), fe);
        this->create_result();
        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_scalar_value = this->result_matrix();
        for (unsigned int s=0; s<dim+1; ++s) {
            for (unsigned int i_p = 0; i_p < n_points; i_p++)
                for (unsigned int i_dof = 0; i_dof < n_dofs; i_dof++) {
                    ref_scalar_value(s, i_dof)(i_p) = ref_shape_vals[s][i_p][i_dof][0];
                }
        }
	}

    void eval() override {}
};

/// Fixed operation of vector shape reference values
template<unsigned int dim, unsigned int spacedim>
class OpRefVector : public Op::Side::Base<spacedim> {
public:
    /// Constructor
    OpRefVector(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Side::Base<spacedim>(_dim, pfev, {dim+1, spacedim}, OpSizeType::fixedSizeOp)
    {
        ASSERT_EQ(this->dim_, dim);

        uint n_points = pfev.get_side_quadrature(dim)->size();
        uint n_dofs = fe->n_dofs();
        this->n_dofs_ = n_dofs;

        this->create_result();
        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_vector_value = this->result_matrix();
        auto ref_shape_vals = this->template ref_shape_values_side<dim>(pfev.get_side_quadrature(dim), fe);

        for (unsigned int s=0; s<dim+1; ++s)
            for (unsigned int i_p = 0; i_p < n_points; i_p++)
                for (unsigned int i_dof = 0; i_dof < n_dofs; i_dof++)
                    for (uint c=0; c<spacedim; ++c)
                        ref_vector_value(s,3*i_dof+c)(i_p) = ref_shape_vals[s][i_p][i_dof][c];
	}

    void eval() override {}
};

/// Fixed operation of gradient scalar reference values
template<unsigned int dim, unsigned int spacedim>
class OpRefGradScalar : public Op::Side::Base<spacedim> {
public:
    /// Constructor
    OpRefGradScalar(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Side::Base<spacedim>(_dim, pfev, {dim+1, dim}, OpSizeType::fixedSizeOp)
    {
        ASSERT_EQ(this->dim_, dim);

        uint n_points = pfev.get_side_quadrature(dim)->size();
        uint n_dofs = fe->n_dofs();
        this->n_dofs_ = n_dofs;

        std::vector<std::vector<std::vector<arma::mat> > > ref_shape_grads = this->template ref_shape_gradients_side<dim>(pfev.get_side_quadrature(dim), fe);
        this->create_result();
        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_scalar_value = this->result_matrix();
        for (unsigned int s=0; s<dim+1; ++s)
            for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
                for (uint i_p=0; i_p<n_points; ++i_p)
                    for (uint c=0; c<dim; ++c)
                        ref_scalar_value(s,dim*i_dof+c)(i_p) = ref_shape_grads[s][i_p][i_dof](c);
    }

    void eval() override {}
};

/// Fixed operation of gradient scalar reference values
template<unsigned int dim, unsigned int spacedim>
class OpRefGradVector : public Op::Side::Base<spacedim> {
public:
    /// Constructor
	OpRefGradVector(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Side::Base<spacedim>(_dim, pfev, {(dim+1)*dim, spacedim}, OpSizeType::fixedSizeOp)
    {
        ASSERT_EQ(this->dim_, dim);

        uint n_points = pfev.get_side_quadrature(dim)->size();
        uint n_dofs = fe->n_dofs();
        this->n_dofs_ = n_dofs;
        uint n_sides = dim+1;

        std::vector<std::vector<std::vector<arma::mat> > > ref_shape_grads = this->template ref_shape_gradients_side<dim>(pfev.get_side_quadrature(dim), fe);
        this->create_result();
        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_vector_value = this->result_matrix();
        for (uint i_sd=0; i_sd<n_sides; ++i_sd)
            for (uint i_c=0; i_c<spacedim; ++i_c)
                for (uint i_dim=0; i_dim<dim; ++i_dim)
                    for (uint i_dof=0; i_dof<n_dofs; ++i_dof)
                        for (uint i_p=0; i_p<n_points; ++i_p) {
                            ref_vector_value(i_sd*dim+i_dim, 3*i_dof+i_c)(i_p) = ref_shape_grads[i_sd][i_p][i_dof](i_dim, i_c);
                        }
    }

    void eval() override {}
};

/// Evaluates scalar values
template<unsigned int dim, unsigned int spacedim>
class OpScalarShape : public Op::Side::Base<spacedim> {
public:
    /// Constructor
	OpScalarShape(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
	: Op::Side::Base<spacedim>(_dim, pfev, {1}, OpSizeType::pointOp)
	{
	    ASSERT_EQ(this->dim_, dim);

	    this->n_dofs_ = fe->n_dofs();
	    this->input_ops_.push_back( pfev.template get< OpRefScalar<dim, spacedim>, dim >(fe) );
	}

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        this->allocate_result(ppv.n_points_, ppv.patch_arena());

        auto ref_vec = this->input_ops(0)->result_matrix();
        auto result_vec = this->result_matrix();

        uint n_dofs = this->n_dofs();
        uint n_sides = ppv.int_table_(3).data_size();        // number of sides on patch
        uint n_patch_points = ppv.int_table_(4).data_size(); // number of points on patch

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i_pt=0; i_pt<n_patch_points; ++i_pt) {
                result_vec(i_dof)(i_pt) = ref_vec(ppv.int_table_(4)(i_pt), i_dof)(i_pt / n_sides);
            }
        }
    }
};

/// Evaluates vector values (FEType == FEVector)
template<unsigned int dim, unsigned int spacedim>
class OpVectorShape : public Op::Side::Base<spacedim> {
public:
    /// Constructor
	OpVectorShape(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : Op::Side::Base<spacedim>(_dim, pfev, {spacedim}, OpSizeType::pointOp), dispatch_op_(dispatch_op)
    {
        ASSERT_EQ(this->dim_, dim);

        this->n_dofs_ = fe->n_dofs();
        this->input_ops_.push_back( pfev.template get< OpRefVector<dim, spacedim>, dim >(fe) );
	}

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        dispatch_op_.allocate_result(ppv.n_points_, ppv.patch_arena());

        auto ref_shape_vec = this->input_ops(0)->result_matrix();  // dim+1 x spacedim
        auto result_vec = dispatch_op_.result_matrix();            // spacdim x 1

        uint n_dofs = this->n_dofs();
        uint n_sides = ppv.int_table_(3).data_size();
        uint n_patch_points = ppv.int_table_(4).data_size();

        for (uint c=0; c<3*n_dofs; c++)
        	result_vec(c) = ArenaVec<double>(n_patch_points, ppv.patch_arena());

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i_pt=0; i_pt<n_patch_points; ++i_pt)
                for (uint c=0; c<3; c++)
                    result_vec(c,i_dof)(i_pt) = ref_shape_vec(ppv.int_table_(4)(i_pt),3*i_dof+c)(i_pt / n_sides);
        }
    }

private:
    PatchOp<spacedim> &dispatch_op_;
};

// class OpVectorCovariantShape
// class OpVectorPiolaShape

/// Dispatch class of vector values
template<unsigned int dim, unsigned int spacedim>
class DispatchVectorShape : public Op::Side::Base<spacedim> {
public:
    /// Constructor
    DispatchVectorShape(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Side::Base<spacedim>(_dim, pfev, {spacedim}, OpSizeType::pointOp), in_op_(nullptr)
    {
        ASSERT_EQ(this->dim_, dim);
        this->n_dofs_ = fe->n_dofs();

        switch (fe->fe_type()) {
            case FEVector:
            {
                in_op_ = new OpVectorShape<dim, spacedim>(_dim, pfev, fe, *this);
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpVectorCovariantShape<dim, spacedim>(_dim, pfev, fe, *this);
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpVectorPiolaShape<dim, spacedim>(_dim, pfev, fe, *this);
                break;
            }
            default:
                ASSERT(false).error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");
        }

	}

    void eval() override {
        in_op_->eval();
    }

private:
    PatchOp<spacedim> *in_op_;
};

/// Evaluates gradient scalar values
template<unsigned int dim, unsigned int spacedim>
class OpGradScalarShape : public Op::Side::Base<spacedim> {
public:
    /// Constructor
	OpGradScalarShape(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
	: Op::Side::Base<spacedim>(_dim, pfev, {spacedim, 1}, OpSizeType::pointOp)
	{
	    ASSERT_EQ(this->dim_, dim);

	    this->n_dofs_ = fe->n_dofs();
	    this->input_ops_.push_back( pfev.template get< Op::Side::El::OpElInvJac<dim, spacedim> >(dim) );
	    this->input_ops_.push_back( pfev.template get< OpRefGradScalar<dim, spacedim>, dim >(fe) );
	}

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        this->allocate_result(ppv.n_points_, ppv.patch_arena());

        auto ref_shape_grads = this->input_ops(1)->result_matrix();
        auto grad_scalar_shape_value = this->result_matrix();

        uint n_dofs = this->n_dofs();
        uint n_points = ref_shape_grads(0).data_size();
        uint n_sides = ppv.int_table_(3).data_size();
        uint n_patch_points = ppv.int_table_(4).data_size();

        // Expands inverse jacobian to inv_jac_expd_value
        auto inv_jac_value = this->input_ops(0)->result_matrix();
        Eigen::Matrix<ArenaVec<double>, dim, 3> inv_jac_expd_value;
        for (uint i=0; i<dim*3; ++i) {
        	inv_jac_expd_value(i) = ArenaVec<double>( n_patch_points, ppv.patch_arena() );
        	for (uint j=0; j<n_patch_points; ++j)
        	    inv_jac_expd_value(i)(j) = inv_jac_value(i)(j%n_sides);
        }

        // Fill ref shape gradients by q_point. DOF and side_idx
        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> ref_shape_grads_expd(dim, n_dofs);
        for (uint i=0; i<dim*n_dofs; ++i) {
            ref_shape_grads_expd(i) = ArenaVec<double>( n_patch_points, ppv.patch_arena() );
        }
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i_pt=0; i_pt<n_points; ++i_pt) {
                uint i_begin = i_pt * n_sides;
                for (uint i_sd=0; i_sd<n_sides; ++i_sd) {
                    for (uint i_c=0; i_c<dim; ++i_c) {
                        ref_shape_grads_expd(i_c, i_dof)(i_begin + i_sd) = ref_shape_grads(ppv.int_table_(3)(i_sd), i_dof*dim+i_c)(i_pt);
                    }
                }
            }
        }

        // computes operation result
        grad_scalar_shape_value = inv_jac_expd_value.transpose() * ref_shape_grads_expd;
    }
};

/// Evaluates vector values (FEType == FEVector)
template<unsigned int dim, unsigned int spacedim>
class OpGradVectorShape : public Op::Side::Base<spacedim> {
public:
    /// Constructor
	OpGradVectorShape(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : Op::Side::Base<spacedim>(_dim, pfev, {spacedim, spacedim}, OpSizeType::pointOp), dispatch_op_(dispatch_op)
    {
        ASSERT_EQ(this->dim_, dim);

        this->n_dofs_ = fe->n_dofs();
        this->input_ops_.push_back( pfev.template get< Op::Side::El::OpElInvJac<dim, spacedim> >(dim) );
        this->input_ops_.push_back( pfev.template get< OpRefGradVector<dim, spacedim>, dim >(fe) );
	}

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        auto inv_jac_value = this->input_ops(0)->result_matrix();    // dim x spacedim
        auto ref_vector_grad = this->input_ops(1)->result_matrix();  // n_sides*dim x spacedim

        uint n_dofs = this->n_dofs();
        uint n_points = ref_vector_grad(0).data_size();
        uint n_patch_sides = ppv.int_table_(3).data_size();
        uint n_patch_points = ppv.int_table_(4).data_size();

        // Expands inverse jacobian to inv_jac_expd_value
        Eigen::Matrix<ArenaVec<double>, dim, 3> inv_jac_expd_value;
        for (uint i=0; i<dim*3; ++i) {
        	inv_jac_expd_value(i) = ArenaVec<double>( n_patch_points, ppv.patch_arena() );
        	for (uint j=0; j<n_patch_points; ++j)
        	    inv_jac_expd_value(i)(j) = inv_jac_value(i)(j%n_patch_sides);
        }

        // Fill ref shape gradients by q_point. DOF and side_idx
        Eigen::Matrix<ArenaVec<double>, dim, 3> ref_shape_grads_expd;
        for (uint i=0; i<3*dim; ++i) {
            ref_shape_grads_expd(i) = ArenaVec<double>( n_patch_points, ppv.patch_arena() );
        }
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {

            for (uint i_pt=0; i_pt<n_points; ++i_pt) {
                uint i_begin = i_pt * n_patch_sides;
                for (uint i_sd=0; i_sd<n_patch_sides; ++i_sd) {
                    for (uint i_dim=0; i_dim<dim; ++i_dim) {
                        for (uint i_c=0; i_c<3; ++i_c) {
                            ref_shape_grads_expd(i_dim, i_c)(i_begin + i_sd) = ref_vector_grad(ppv.int_table_(3)(i_sd)*dim+i_dim, 3*i_dof+i_c)(i_pt);
                        }
                    }
                }
            }

            // computes operation result
            Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > res_submat = dispatch_op_.result_sub_matrix(i_dof);
            res_submat = (inv_jac_expd_value.transpose() * ref_shape_grads_expd).transpose();
        }
    }

private:
    PatchOp<spacedim> &dispatch_op_;
};

// class OpGradVectorCovariantShape
// class OpGradVectorPiolaShape

/// Dispatch class of vector values
template<unsigned int dim, unsigned int spacedim>
class DispatchGradVectorShape : public Op::Side::Base<spacedim> {
public:
    /// Constructor
	DispatchGradVectorShape(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Side::Base<spacedim>(_dim, pfev, {spacedim, spacedim}, OpSizeType::pointOp), in_op_(nullptr)
    {
        ASSERT_EQ(this->dim_, dim);
        this->n_dofs_ = fe->n_dofs();

        switch (fe->fe_type()) {
            case FEVector:
            {
                in_op_ = new OpGradVectorShape<dim, spacedim>(_dim, pfev, fe, *this);
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpGradVectorCovariantShape<dim, spacedim>(_dim, pfev, fe, *this);
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpGradVectorPiolaShape<dim, spacedim>(_dim, pfev, fe, *this);
                break;
            }
            default:
                ASSERT(false).error("Type of FiniteElement of grad_vector_shape accessor must be FEVector, FEVectorPiola or FEVectorContravariant!\n");
        }

	}

    void eval() override {
        in_op_->eval();
    }

private:
    PatchOp<spacedim> *in_op_;
};

/// Evaluates vector symmetric gradients
template<unsigned int dim, unsigned int spacedim>
class OpVectorSymGrad : public Op::VectorSymGradBase<dim, spacedim> {
public:
    /// Constructor
	OpVectorSymGrad(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::VectorSymGradBase<dim, spacedim>(_dim, pfev, {spacedim, spacedim}, OpSizeType::pointOp)
    {
        this->bulk_side_ = 1;
        this->n_dofs_ = fe->n_dofs();
        this->input_ops_.push_back( pfev.template get< DispatchGradVectorShape<dim, spacedim>, dim >(fe) );
    }
};

/// Evaluates vector divergence
template<unsigned int dim, unsigned int spacedim>
class OpVectorDivergence : public Op::VectorDivergenceBase<dim, spacedim> {
public:
    /// Constructor
	OpVectorDivergence(uint _dim, PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::VectorDivergenceBase<dim, spacedim>(_dim, pfev, {1}, OpSizeType::pointOp)
    {
        this->bulk_side_ = 1;
        this->n_dofs_ = fe->n_dofs();
        this->input_ops_.push_back( pfev.template get< DispatchGradVectorShape<dim, spacedim>, dim >(fe) );
    }
};

} // end of namespace Op::Side::Pt

} // end of namespace Side

} // end of namespace Op


#endif /* OP_FUNCTION_HH_ */
