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
	OpRefScalar(uint _dim, PatchFEValues<spacedim> &pfev, uint component_idx)
	: Op::Bulk::Base<spacedim>(_dim, pfev, {1}, OpSizeType::fixedSizeOp)
	{
	    ASSERT_EQ(this->dim_, dim);
	    auto fe_component = pfev.template fe_comp<dim>(component_idx);
	    ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");

	    uint n_points = pfev.get_bulk_quadrature(dim)->size();
	    uint n_dofs = fe_component->n_dofs();
	    this->n_dofs_ = n_dofs;

	    auto ref_shape_vals = this->template ref_shape_values_bulk<dim>(pfev.get_bulk_quadrature(dim), fe_component);
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

/// Fixed operation of gradient scalar reference values
template<unsigned int dim, unsigned int spacedim>
class OpRefGradScalar : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
    OpRefGradScalar(uint _dim, PatchFEValues<spacedim> &pfev, uint component_idx)
    : Op::Bulk::Base<spacedim>(_dim, pfev, {dim, 1}, OpSizeType::fixedSizeOp)
    {
        ASSERT_EQ(this->dim_, dim);
        auto fe_component = pfev.template fe_comp<dim>(component_idx);
        ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape accessor must be FEScalar!\n");

        uint n_points = pfev.get_bulk_quadrature(dim)->size();
        uint n_dofs = fe_component->n_dofs();
        this->n_dofs_ = n_dofs;

        std::vector<std::vector<arma::mat> > ref_shape_grads = this->template ref_shape_gradients_bulk<dim>(pfev.get_bulk_quadrature(dim), fe_component);
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

/// Evaluates scalar values
template<unsigned int dim, unsigned int spacedim>
class OpScalarShape : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
	OpScalarShape(uint _dim, PatchFEValues<spacedim> &pfev, uint component_idx)
	: Op::Bulk::Base<spacedim>(_dim, pfev, {1}, OpSizeType::pointOp)
	{
	    ASSERT_EQ(this->dim_, dim);
	    auto fe_component = pfev.template fe_comp<dim>(component_idx);
	    ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape accessor must be FEScalar!\n");

	    uint n_dofs = fe_component->n_dofs();
	    this->n_dofs_ = n_dofs;

	    this->input_ops_.push_back( pfev.template get< OpRefScalar<dim, spacedim> >(dim, component_idx) );
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

/// Evaluates gradient scalar values
template<unsigned int dim, unsigned int spacedim>
class OpGradScalarShape : public Op::Bulk::Base<spacedim> {
public:
    /// Constructor
	OpGradScalarShape(uint _dim, PatchFEValues<spacedim> &pfev, uint component_idx)
	: Op::Bulk::Base<spacedim>(_dim, pfev, {spacedim, 1}, OpSizeType::pointOp)
	{
	    ASSERT_EQ(this->dim_, dim);
	    auto fe_component = pfev.template fe_comp<dim>(component_idx);
	    ASSERT_EQ(fe_component->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape accessor must be FEScalar!\n");

	    uint n_dofs = fe_component->n_dofs();
	    this->n_dofs_ = n_dofs;

	    this->input_ops_.push_back( pfev.template get< Op::Bulk::El::OpInvJac<dim, spacedim> >(dim) );
	    this->input_ops_.push_back( pfev.template get< OpRefGradScalar<dim, spacedim> >(dim, component_idx) );
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
        this->allocate_result( ppv.n_elems_, *this->patch_fe_->patch_fe_data_.patch_arena_ );
        auto result = this->result_matrix();

        for (uint i_elm=0; i_elm<ppv.elem_list_.size(); ++i_elm)
            for (uint i_col=0; i_col<this->dim_+1; ++i_col)
                for (uint i_row=0; i_row<spacedim; ++i_row) {
                    result(i_row, i_col)(i_elm) = ( *ppv.elem_list_[i_elm].node(i_col) )(i_row);
                }
    }

};

template<unsigned int spacedim>
class OpSdCoords : public Op::Side::Base<spacedim> {
public:
    /// Constructor
	OpSdCoords(uint dim, PatchFEValues<spacedim> &pfev)
    : Op::Side::Base<spacedim>(dim, pfev, {spacedim, dim}, OpSizeType::elemOp) {}

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        this->allocate_result( ppv.n_elems_, *this->patch_fe_->patch_fe_data_.patch_arena_ );
        auto result = this->result_matrix();

        for (uint i_side=0; i_side<ppv.side_list_.size(); ++i_side)
            for (uint i_col=0; i_col<this->dim_; ++i_col)
                for (uint i_row=0; i_row<spacedim; ++i_row) {
                    result(i_row, i_col)(i_side) = ( *ppv.side_list_[i_side].node(i_col) )(i_row);
                }
    }

};

} // end of namespace Op::Side::El

namespace Pt {

} // end of namespace Op::Side::Pt

} // end of namespace Side

} // end of namespace Op


#endif /* OP_FUNCTION_HH_ */
