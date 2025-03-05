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

/// Class used as template type for type resolution Bulk / Side
class BulkDomain {
public:
    static uint domain() {
        return 0;
    }

    static inline constexpr uint n_nodes(uint dim) {
        return dim+1;
    }
};

/// Class used as template type for type resolution Bulk / Side
class SideDomain {
public:
    static uint domain() {
        return 1;
    }

    static inline constexpr uint n_nodes(uint dim) {
        return dim;
    }
};


/**
 * Evaluates node coordinates on Bulk (Element) / Side
 *
 * Template parameters:
 *   dim       Dimension of operation
 *   ElDomain  Target domain - data are evaluated on Bulk / Side domain
 *   Domain    Source domain - operation is called from Bulk / Side domain
 *   spacedim  Dimension of the solved task
 */
template<unsigned int dim, class ElDomain, class Domain, unsigned int spacedim = 3>
class Coords : public PatchOp<spacedim> {
public:
    /// Constructor
    Coords(PatchFEValues<spacedim> &pfev)
    : PatchOp<spacedim>(dim, pfev, {spacedim, ElDomain::n_nodes(dim)}, OpSizeType::elemOp) {
        this->bulk_side_ = Domain::domain();
    }

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        this->allocate_result( ppv.n_elems_, *ppv.patch_fe_data_.patch_arena_ );
        auto result = this->result_matrix();

        for (uint i_elm=0; i_elm<ppv.elem_list_.size(); ++i_elm)
            for (uint i_col=0; i_col<ElDomain::n_nodes(dim); ++i_col)
                for (uint i_row=0; i_row<spacedim; ++i_row) {
                    result(i_row, i_col)(i_elm) = ( *ppv.template node<ElDomain>(i_elm, i_col) )(i_row);
                }
    }

};

/// Evaluates Jacobians on Bulk (Element) / Side
template<unsigned int dim, class ElDomain, class Domain, unsigned int spacedim = 3>
class Jac : public PatchOp<spacedim> {
public:
    /// Constructor
    Jac(PatchFEValues<spacedim> &pfev)
    : PatchOp<spacedim>(dim, pfev, {spacedim, ElDomain::n_nodes(dim)-1}, OpSizeType::elemOp)
    {
        this->bulk_side_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< Op::Coords<dim, ElDomain, Domain, spacedim>, dim >() );
    }

    void eval() override {
        auto jac_value = this->result_matrix();
        auto coords_value = this->input_ops(0)->result_matrix();
        for (unsigned int i=0; i<spacedim; i++)
            for (unsigned int j=0; j<ElDomain::n_nodes(dim)-1; j++)
                jac_value(i,j) = coords_value(i,j+1) - coords_value(i,0);
    }
};

/// Evaluates Jacobian determinants on Bulk (Element) / Side
template<unsigned int dim, class ElDomain, class Domain, unsigned int spacedim = 3>
class JacDet : public PatchOp<spacedim> {
public:
    /// Constructor
	JacDet(PatchFEValues<spacedim> &pfev)
	: PatchOp<spacedim>(dim, pfev, {1}, OpSizeType::elemOp)
	{
        this->bulk_side_ = Domain::domain();
	    this->input_ops_.push_back( pfev.template get< Op::Jac<dim, ElDomain, Domain, spacedim>, dim >() );
	}

    void eval() override {
        auto jac_det_value = this->result_matrix();
        auto jac_value = this->input_ops(0)->result_matrix();
        jac_det_value(0) = eigen_arena_tools::determinant<spacedim, ElDomain::n_nodes(dim)-1>(jac_value).abs();
    }
};

/// Template specialization of previous: dim=1, domain=Side
template<>
class JacDet<1, Op::SideDomain, Op::SideDomain, 3> : public PatchOp<3> {
public:
    /// Constructor
    JacDet(PatchFEValues<3> &pfev)
    : PatchOp<3>(1, pfev, {1}, OpSizeType::elemOp)
    {
        this->bulk_side_ = 1;
    }

    void eval() override {
        PatchPointValues<3> &ppv = this->ppv();
        this->allocate_result( ppv.n_elems_, *ppv.patch_fe_data_.patch_arena_ );
        auto jac_det_value = this->result_matrix();
        for (uint i=0;i<ppv.n_elems_; ++i) {
            jac_det_value(0,0)(i) = 1.0;
        }
    }
};

/**
 * Evaluates Inverse Jacobians on Bulk (Element) / Side
 * ElDomain (target) is always Bulk
 */
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class InvJac : public PatchOp<spacedim> {
public:
    /// Constructor
	InvJac(PatchFEValues<spacedim> &pfev)
    : PatchOp<spacedim>(dim, pfev, {dim, spacedim}, OpSizeType::elemOp)
    {
        this->bulk_side_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< Op::Jac<dim, BulkDomain, Domain, spacedim>, dim >() );
    }

    void eval() override {
        auto inv_jac_value = this->result_matrix();
        auto jac_value = this->input_ops(0)->result_matrix();
        inv_jac_value = eigen_arena_tools::inverse<spacedim, dim>(jac_value);
    }
};

/// Evaluates coordinates of quadrature points - not implemented yet
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class PtCoords : public PatchOp<spacedim> {
public:
    /// Constructor
    PtCoords(PatchFEValues<spacedim> &pfev)
    : PatchOp<spacedim>(dim, pfev, {spacedim}, OpSizeType::pointOp)
    {
        this->bulk_side_ = Domain::domain();
    }

    void eval() override {}
};

/**
 * Fixed operation points weights
 * ElDomain (target) is equivalent with Domain (source)
 */
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class Weights : public PatchOp<spacedim> {
public:
    /// Constructor
    Weights(PatchFEValues<spacedim> &pfev)
    : PatchOp<spacedim>(dim, pfev, {1}, OpSizeType::fixedSizeOp)
    {
        this->bulk_side_ = Domain::domain();
        // create result vector of weights operation in assembly arena
        const std::vector<double> &point_weights_vec = this->ppv().get_quadrature()->get_weights();
        this->allocate_result(point_weights_vec.size(), pfev.asm_arena());
        for (uint i=0; i<point_weights_vec.size(); ++i)
            this->result_(0)(i) = point_weights_vec[i];
    }

    void eval() override {}
};

/**
 * Evaluates JxW on quadrature points
 * ElDomain (target) is equivalent with Domain (source)
 */
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class JxW : public PatchOp<spacedim> {
public:
    /// Constructor
    JxW(PatchFEValues<spacedim> &pfev)
    : PatchOp<spacedim>(dim, pfev, {1}, OpSizeType::pointOp)
    {
        this->bulk_side_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< Op::Weights<dim, Domain, spacedim>, dim >() );
        this->input_ops_.push_back( pfev.template get< Op::JacDet<dim, Domain, Domain, spacedim>, dim >() );
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

/// Fixed operation of scalar shape reference values
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class RefScalar : public PatchOp<spacedim> {
public:
    /// Constructor
    RefScalar(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, {1}, OpSizeType::fixedSizeOp, fe->n_dofs())
    {
        this->bulk_side_ = Domain::domain();
        uint n_points = pfev.get_bulk_quadrature(dim)->size();

        auto ref_shape_vals = this->template ref_shape_values_bulk<dim>(pfev.get_bulk_quadrature(dim), fe);
        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_scalar_value = this->result_matrix();
        for (unsigned int i_p = 0; i_p < n_points; i_p++)
            for (unsigned int i_dof = 0; i_dof < this->n_dofs_; i_dof++) {
                ref_scalar_value(i_dof)(i_p) = ref_shape_vals[i_p][i_dof][0];
            }
    }

    void eval() override {}
};

/// Template specialization of previous: Domain=SideDomain
template<unsigned int dim, unsigned int spacedim>
class RefScalar<dim, Op::SideDomain, spacedim> : public PatchOp<spacedim> {
public:
    /// Constructor
    RefScalar(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, {dim+1}, OpSizeType::fixedSizeOp, fe->n_dofs())
    {
        this->bulk_side_ = Op::SideDomain::domain();
        uint n_points = pfev.get_side_quadrature(dim)->size();

        auto ref_shape_vals = this->template ref_shape_values_side<dim>(pfev.get_side_quadrature(dim), fe);
        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_scalar_value = this->result_matrix();
        for (unsigned int s=0; s<dim+1; ++s) {
            for (unsigned int i_p = 0; i_p < n_points; i_p++)
                for (unsigned int i_dof = 0; i_dof < this->n_dofs_; i_dof++) {
                    ref_scalar_value(s, i_dof)(i_p) = ref_shape_vals[s][i_p][i_dof][0];
                }
        }
	}

    void eval() override {}
};

/// Fixed operation of vector shape reference values
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class RefVector : public PatchOp<spacedim> {
public:
    /// Constructor
	RefVector(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
	: PatchOp<spacedim>(dim, pfev, {spacedim}, OpSizeType::fixedSizeOp, fe->n_dofs())
	{
        this->bulk_side_ = Domain::domain();
        uint n_points = pfev.get_bulk_quadrature(dim)->size();

        auto ref_shape_vals = this->template ref_shape_values_bulk<dim>(pfev.get_bulk_quadrature(dim), fe);
        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_vector_value = this->result_matrix();

        for (uint i_dof=0; i_dof<this->n_dofs_; ++i_dof)
            for (uint i_p=0; i_p<n_points; ++i_p)
                for (uint c=0; c<spacedim; ++c)
                    ref_vector_value(c, i_dof)(i_p) = ref_shape_vals[i_p][i_dof](c);
	}

    void eval() override {}
};

/// Template specialization of previous: Domain=SideDomain
template<unsigned int dim, unsigned int spacedim>
class RefVector<dim, Op::SideDomain, spacedim> : public PatchOp<spacedim> {
public:
    /// Constructor
    RefVector(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, {dim+1, spacedim}, OpSizeType::fixedSizeOp, fe->n_dofs())
    {
        this->bulk_side_ = Op::SideDomain::domain();
        uint n_points = pfev.get_side_quadrature(dim)->size();

        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_vector_value = this->result_matrix();
        auto ref_shape_vals = this->template ref_shape_values_side<dim>(pfev.get_side_quadrature(dim), fe);

        for (unsigned int s=0; s<dim+1; ++s)
            for (unsigned int i_p = 0; i_p < n_points; i_p++)
                for (unsigned int i_dof = 0; i_dof < this->n_dofs_; i_dof++)
                    for (uint c=0; c<spacedim; ++c)
                        ref_vector_value(s,3*i_dof+c)(i_p) = ref_shape_vals[s][i_p][i_dof][c];
	}

    void eval() override {}
};

/// Fixed operation of gradient scalar reference values
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class RefGradScalar : public PatchOp<spacedim> {
public:
    /// Constructor
    RefGradScalar(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, {dim, 1}, OpSizeType::fixedSizeOp, fe->n_dofs())
    {
        this->bulk_side_ = Domain::domain();
        uint n_points = pfev.get_bulk_quadrature(dim)->size();

        std::vector<std::vector<arma::mat> > ref_shape_grads = this->template ref_shape_gradients_bulk<dim>(pfev.get_bulk_quadrature(dim), fe);
        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_scalar_value = this->result_matrix();
        for (uint i_row=0; i_row<ref_scalar_value.rows(); ++i_row) {
            for (uint i_dof=0; i_dof<this->n_dofs_; ++i_dof)
                for (uint i_p=0; i_p<n_points; ++i_p)
                    ref_scalar_value(i_row, i_dof)(i_p) = ref_shape_grads[i_p][i_dof](i_row);
        }
    }

    void eval() override {}
};

/// Template specialization of previous: Domain=SideDomain
template<unsigned int dim, unsigned int spacedim>
class RefGradScalar<dim, Op::SideDomain, spacedim> : public PatchOp<spacedim> {
public:
    /// Constructor
    RefGradScalar(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, {dim+1, dim}, OpSizeType::fixedSizeOp, fe->n_dofs())
    {
        this->bulk_side_ = Op::SideDomain::domain();
        uint n_points = pfev.get_side_quadrature(dim)->size();

        std::vector<std::vector<std::vector<arma::mat> > > ref_shape_grads = this->template ref_shape_gradients_side<dim>(pfev.get_side_quadrature(dim), fe);
        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_scalar_value = this->result_matrix();
        for (unsigned int s=0; s<dim+1; ++s)
            for (uint i_dof=0; i_dof<this->n_dofs_; ++i_dof)
                for (uint i_p=0; i_p<n_points; ++i_p)
                    for (uint c=0; c<dim; ++c)
                        ref_scalar_value(s,dim*i_dof+c)(i_p) = ref_shape_grads[s][i_p][i_dof](c);
    }

    void eval() override {}
};

/// Fixed operation of gradient scalar reference values
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class RefGradVector : public PatchOp<spacedim> {
public:
    /// Constructor
	RefGradVector(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, {dim, spacedim}, OpSizeType::fixedSizeOp, fe->n_dofs())
    {
        this->bulk_side_ = Domain::domain();
        uint n_points = pfev.get_bulk_quadrature(dim)->size();

        std::vector<std::vector<arma::mat> > ref_shape_grads = this->template ref_shape_gradients_bulk<dim>(pfev.get_bulk_quadrature(dim), fe);
        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_vector_value = this->result_matrix();
        for (uint i_c=0; i_c<spacedim; ++i_c) {
            for (uint i_dim=0; i_dim<dim; ++i_dim)
                for (uint i_dof=0; i_dof<this->n_dofs_; ++i_dof)
                    for (uint i_p=0; i_p<n_points; ++i_p)
                        ref_vector_value(i_dim,3*i_dof+i_c)(i_p) = ref_shape_grads[i_p][i_dof](i_dim, i_c);
        }
    }

    void eval() override {}
};

/// Template specialization of previous: Domain=SideDomain
template<unsigned int dim, unsigned int spacedim>
class RefGradVector<dim, Op::SideDomain, spacedim> : public PatchOp<spacedim> {
public:
    /// Constructor
	RefGradVector(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, {(dim+1)*dim, spacedim}, OpSizeType::fixedSizeOp, fe->n_dofs())
    {
        this->bulk_side_ = Op::SideDomain::domain();
        uint n_points = pfev.get_side_quadrature(dim)->size();
        uint n_sides = dim+1;

        std::vector<std::vector<std::vector<arma::mat> > > ref_shape_grads = this->template ref_shape_gradients_side<dim>(pfev.get_side_quadrature(dim), fe);
        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_vector_value = this->result_matrix();
        for (uint i_sd=0; i_sd<n_sides; ++i_sd)
            for (uint i_c=0; i_c<spacedim; ++i_c)
                for (uint i_dim=0; i_dim<dim; ++i_dim)
                    for (uint i_dof=0; i_dof<this->n_dofs_; ++i_dof)
                        for (uint i_p=0; i_p<n_points; ++i_p) {
                            ref_vector_value(i_sd*dim+i_dim, 3*i_dof+i_c)(i_p) = ref_shape_grads[i_sd][i_p][i_dof](i_dim, i_c);
                        }
    }

    void eval() override {}
};


/* ************************************************* */
/// Base class of bulk and side vector symmetric gradients
template<unsigned int dim, unsigned int spacedim = 3>
class VectorSymGradBase : public PatchOp<spacedim> {
public:
    /// Constructor
    VectorSymGradBase(PatchFEValues<spacedim> &pfev, std::initializer_list<uint> shape, OpSizeType size_type, uint n_dofs)
    : PatchOp<spacedim>(dim, pfev, shape, size_type, n_dofs)
    {}

    void eval() override {
        for (uint i_dof=0; i_dof<this->n_dofs(); ++i_dof) {
            Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > grad_vector_dof = this->input_ops(0)->result_sub_matrix(i_dof);
            Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > sym_grad_dof = this->result_sub_matrix(i_dof);
            sym_grad_dof = 0.5 * (grad_vector_dof.transpose() + grad_vector_dof);
        }
    }
};

/// Base class of bulk and side vector divergence
template<unsigned int dim, unsigned int spacedim = 3>
class VectorDivergenceBase : public PatchOp<spacedim> {
public:
    /// Constructor
    VectorDivergenceBase(PatchFEValues<spacedim> &pfev, std::initializer_list<uint> shape, OpSizeType size_type, uint n_dofs)
    : PatchOp<spacedim>(dim, pfev, shape, size_type, n_dofs)
    {}

    void eval() override {
        auto divergence_value = this->result_matrix();

        for (uint i_dof=0; i_dof<this->n_dofs(); ++i_dof) {
            Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > grad_vector_dof = this->input_ops(0)->result_sub_matrix(i_dof);
            divergence_value(i_dof) = grad_vector_dof(0,0) + grad_vector_dof(1,1) + grad_vector_dof(2,2);
        }
    }
};

/// Class represents zero operation of Join quantities.
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class OpZero : public PatchOp<spacedim> {
public:
    /// Constructor
	OpZero(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, {spacedim, spacedim}, OpSizeType::fixedSizeOp, fe->n_dofs())
    {
        this->bulk_side_ = Domain::domain();
        this->allocate_const_result( this->ppv().patch_fe_data_.zero_vec_ );
    }

    void eval() override {}
};


namespace Bulk {

/// Common ancestor of all bulk operations.
template<unsigned int dim, unsigned int spacedim = 3>
class Base : public PatchOp<spacedim> {
public:
    /// Constructor
	Base(PatchFEValues<spacedim> &pfev, std::initializer_list<uint> shape, OpSizeType size_type, uint n_dofs = 1)
    : PatchOp<spacedim>(dim, pfev, shape, size_type, n_dofs)
    {
        this->bulk_side_ = 0;
    }
};

namespace Pt {

/// Evaluates scalar values
template<unsigned int dim, unsigned int spacedim = 3>
class OpScalarShape : public Op::Bulk::Base<dim, spacedim> {
public:
    /// Constructor
    OpScalarShape(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Bulk::Base<dim, spacedim>(pfev, {1}, OpSizeType::pointOp, fe->n_dofs())
    {
        ASSERT_EQ(fe->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape must be FEScalar!\n");
        this->input_ops_.push_back( pfev.template get< Op::RefScalar<dim, Op::BulkDomain, spacedim>, dim >(fe) );
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
template<unsigned int dim, unsigned int spacedim = 3>
class OpVectorShape : public Op::Bulk::Base<dim, spacedim> {
public:
    /// Constructor
	OpVectorShape(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : Op::Bulk::Base<dim, spacedim>(pfev, {spacedim}, OpSizeType::pointOp, fe->n_dofs()), dispatch_op_(dispatch_op)
    {
        this->input_ops_.push_back( pfev.template get< Op::RefVector<dim, Op::BulkDomain, spacedim>, dim >(fe) );
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
        for (uint c=0; c<spacedim*n_dofs; ++c) {
            ref_shape_ovec(c) = ArenaOVec(ref_shape_vec(c));
        }

        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> result_ovec = elem_ovec * ref_shape_ovec;
        for (uint c=0; c<spacedim*n_dofs; ++c)
            result_vec(c) = result_ovec(c).get_vec();
    }

private:
    PatchOp<spacedim> &dispatch_op_;
};

// class OpVectorCovariantShape
// class OpVectorPiolaShape

/// Dispatch class of vector values
template<unsigned int dim, unsigned int spacedim = 3>
class DispatchVectorShape : public Op::Bulk::Base<dim, spacedim> {
public:
    /// Constructor
    DispatchVectorShape(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Bulk::Base<dim, spacedim>(pfev, {spacedim}, OpSizeType::pointOp, fe->n_dofs()), in_op_(nullptr)
    {
        switch (fe->fe_type()) {
            case FEVector:
            {
                in_op_ = new OpVectorShape<dim, spacedim>(pfev, fe, *this);
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpVectorCovariantShape<dim, spacedim>(pfev, fe, *this);
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpVectorPiolaShape<dim, spacedim>(pfev, fe, *this);
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
template<unsigned int dim, unsigned int spacedim = 3>
class OpGradScalarShape : public Op::Bulk::Base<dim, spacedim> {
public:
    /// Constructor
    OpGradScalarShape(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Bulk::Base<dim, spacedim>(pfev, {spacedim, 1}, OpSizeType::pointOp, fe->n_dofs())
    {
        ASSERT_EQ(fe->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape must be FEScalar!\n");
        this->input_ops_.push_back( pfev.template get< Op::InvJac<dim, Op::BulkDomain, spacedim>, dim >() );
        this->input_ops_.push_back( pfev.template get< Op::RefGradScalar<dim, Op::BulkDomain, spacedim>, dim >(fe) );
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
template<unsigned int dim, unsigned int spacedim = 3>
class OpGradVectorShape : public Op::Bulk::Base<dim, spacedim> {
public:
    /// Constructor
	OpGradVectorShape(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : Op::Bulk::Base<dim, spacedim>(pfev, {spacedim, spacedim}, OpSizeType::pointOp, fe->n_dofs()), dispatch_op_(dispatch_op)
    {
        this->input_ops_.push_back( pfev.template get< Op::InvJac<dim, Op::BulkDomain, spacedim>, dim >() );
        this->input_ops_.push_back( pfev.template get< Op::RefGradVector<dim, Op::BulkDomain, spacedim>, dim >(fe) );
	}

    void eval() override {
        auto inv_jac_vec = this->input_ops(0)->result_matrix();    // dim x spacedim
        auto ref_grads_vec = this->input_ops(1)->result_matrix();  // dim x spacedim
        auto result_vec = dispatch_op_.result_matrix();            // spacedim x spacedim

        uint n_dofs = this->n_dofs();

        Eigen::Matrix<ArenaOVec<double>, dim, 3> inv_jac_ovec;
        for (uint i=0; i<dim*spacedim; ++i) {
            inv_jac_ovec(i) = ArenaOVec(inv_jac_vec(i));
        }

        Eigen::Matrix<ArenaOVec<double>, dim, spacedim> ref_grads_ovec;
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i=0; i<spacedim*dim; ++i) {
                ref_grads_ovec(i) = ArenaOVec(ref_grads_vec(i_dof*3*dim + i));
            }

            Eigen::Matrix<ArenaOVec<double>, spacedim, spacedim> result_ovec = inv_jac_ovec.transpose() * ref_grads_ovec;
            for (uint i=0; i<spacedim; ++i) {
                for (uint j=0; j<spacedim; ++j) {
                   result_vec(j,i+spacedim*i_dof) = result_ovec(i,j).get_vec();
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
template<unsigned int dim, unsigned int spacedim = 3>
class DispatchGradVectorShape : public Op::Bulk::Base<dim, spacedim> {
public:
    /// Constructor
	DispatchGradVectorShape(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Bulk::Base<dim, spacedim>(pfev, {spacedim, spacedim}, OpSizeType::pointOp, fe->n_dofs()), in_op_(nullptr)
    {
        switch (fe->fe_type()) {
            case FEVector:
            {
                in_op_ = new OpGradVectorShape<dim, spacedim>(pfev, fe, *this);
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpGradVectorCovariantShape<dim, spacedim>(pfev, fe, *this);
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpGradVectorPiolaShape<dim, spacedim>(pfev, fe, *this);
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
template<unsigned int dim, unsigned int spacedim = 3>
class OpVectorSymGrad : public Op::VectorSymGradBase<dim, spacedim> {
public:
    /// Constructor
	OpVectorSymGrad(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::VectorSymGradBase<dim, spacedim>(pfev, {spacedim, spacedim}, OpSizeType::pointOp, fe->n_dofs())
    {
        this->bulk_side_ = 0;
        this->input_ops_.push_back( pfev.template get< DispatchGradVectorShape<dim, spacedim>, dim >(fe) );
    }
};

/// Evaluates vector divergence
template<unsigned int dim, unsigned int spacedim = 3>
class OpVectorDivergence : public Op::VectorDivergenceBase<dim, spacedim> {
public:
    /// Constructor
	OpVectorDivergence(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::VectorDivergenceBase<dim, spacedim>(pfev, {1}, OpSizeType::pointOp, fe->n_dofs())
    {
        this->bulk_side_ = 0;
        this->input_ops_.push_back( pfev.template get< DispatchGradVectorShape<dim, spacedim>, dim >(fe) );
    }
};

} // end of namespace Op::Bulk::Pt

} // end of namespace Op::Bulk

namespace Side {

/// Common ancestor of all side operations.
template<unsigned int dim, unsigned int spacedim = 3>
class Base : public PatchOp<spacedim> {
public:
    /// Constructor
	Base(PatchFEValues<spacedim> &pfev, std::initializer_list<uint> shape, OpSizeType size_type, uint n_dofs = 1)
    : PatchOp<spacedim>(dim, pfev, shape, size_type, n_dofs)
    {
        this->bulk_side_ = 1;
    }
};

namespace Pt {

/// Evaluates normal vector on quadrature points
template<unsigned int dim, unsigned int spacedim = 3>
class OpNormalVec : public Op::Side::Base<dim, spacedim> {
public:
    /// Constructor
    OpNormalVec(PatchFEValues<spacedim> &pfev)
    : Op::Side::Base<dim, spacedim>(pfev, {spacedim}, OpSizeType::elemOp)
    {
        this->input_ops_.push_back( pfev.template get< Op::InvJac<dim, Op::SideDomain, spacedim>, dim >() );
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

        for (uint i=0; i<spacedim; ++i) {
            normal_value(i) = normal_value(i) / norm_vec;
        }
    }
};

/// Evaluates scalar values
template<unsigned int dim, unsigned int spacedim = 3>
class OpScalarShape : public Op::Side::Base<dim, spacedim> {
public:
    /// Constructor
    OpScalarShape(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Side::Base<dim, spacedim>(pfev, {1}, OpSizeType::pointOp, fe->n_dofs())
    {
        ASSERT_EQ(fe->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape must be FEScalar!\n");
        this->input_ops_.push_back( pfev.template get< Op::RefScalar<dim, Op::SideDomain, spacedim>, dim >(fe) );
    }

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        this->allocate_result(ppv.n_points_, ppv.patch_arena());

        auto ref_vec = this->input_ops(0)->result_matrix();
        auto result_vec = this->result_matrix();

        uint n_dofs = this->n_dofs();
        uint n_sides = ppv.n_elems_;         // number of sides on patch
        uint n_patch_points = ppv.n_points_; // number of points on patch

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i_pt=0; i_pt<n_patch_points; ++i_pt) {
                result_vec(i_dof)(i_pt) = ref_vec(ppv.int_table_(4)(i_pt), i_dof)(i_pt / n_sides);
            }
        }
    }
};

/// Evaluates vector values (FEType == FEVector)
template<unsigned int dim, unsigned int spacedim = 3>
class OpVectorShape : public Op::Side::Base<dim, spacedim> {
public:
    /// Constructor
	OpVectorShape(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : Op::Side::Base<dim, spacedim>(pfev, {spacedim}, OpSizeType::pointOp, fe->n_dofs()), dispatch_op_(dispatch_op)
    {
        this->input_ops_.push_back( pfev.template get< Op::RefVector<dim, Op::SideDomain, spacedim>, dim >(fe) );
	}

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        dispatch_op_.allocate_result(ppv.n_points_, ppv.patch_arena());

        auto ref_shape_vec = this->input_ops(0)->result_matrix();  // dim+1 x spacedim
        auto result_vec = dispatch_op_.result_matrix();            // spacdim x 1

        uint n_dofs = this->n_dofs();
        uint n_sides = ppv.n_elems_;
        uint n_patch_points = ppv.n_points_;

        for (uint c=0; c<spacedim*n_dofs; c++)
        	result_vec(c) = ArenaVec<double>(n_patch_points, ppv.patch_arena());

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i_pt=0; i_pt<n_patch_points; ++i_pt)
                for (uint c=0; c<spacedim; c++)
                    result_vec(c,i_dof)(i_pt) = ref_shape_vec(ppv.int_table_(4)(i_pt),3*i_dof+c)(i_pt / n_sides);
        }
    }

private:
    PatchOp<spacedim> &dispatch_op_;
};

// class OpVectorCovariantShape
// class OpVectorPiolaShape

/// Dispatch class of vector values
template<unsigned int dim, unsigned int spacedim = 3>
class DispatchVectorShape : public Op::Side::Base<dim, spacedim> {
public:
    /// Constructor
    DispatchVectorShape(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Side::Base<dim, spacedim>(pfev, {spacedim}, OpSizeType::pointOp, fe->n_dofs()), in_op_(nullptr)
    {
        switch (fe->fe_type()) {
            case FEVector:
            {
                in_op_ = new OpVectorShape<dim, spacedim>(pfev, fe, *this);
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpVectorCovariantShape<dim, spacedim>(pfev, fe, *this);
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpVectorPiolaShape<dim, spacedim>(pfev, fe, *this);
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
template<unsigned int dim, unsigned int spacedim = 3>
class OpGradScalarShape : public Op::Side::Base<dim, spacedim> {
public:
    /// Constructor
    OpGradScalarShape(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Side::Base<dim, spacedim>(pfev, {spacedim, 1}, OpSizeType::pointOp, fe->n_dofs())
    {
        ASSERT_EQ(fe->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape must be FEScalar!\n");
        this->input_ops_.push_back( pfev.template get< Op::InvJac<dim, Op::SideDomain, spacedim>, dim >() );
        this->input_ops_.push_back( pfev.template get< Op::RefGradScalar<dim, Op::SideDomain, spacedim>, dim >(fe) );
    }

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        this->allocate_result(ppv.n_points_, ppv.patch_arena());

        auto ref_shape_grads = this->input_ops(1)->result_matrix();
        auto grad_scalar_shape_value = this->result_matrix();

        uint n_dofs = this->n_dofs();
        uint n_points = ref_shape_grads(0).data_size();
        uint n_sides = ppv.n_elems_;
        uint n_patch_points = ppv.n_points_;

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
template<unsigned int dim, unsigned int spacedim = 3>
class OpGradVectorShape : public Op::Side::Base<dim, spacedim> {
public:
    /// Constructor
	OpGradVectorShape(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : Op::Side::Base<dim, spacedim>(pfev, {spacedim, spacedim}, OpSizeType::pointOp, fe->n_dofs()), dispatch_op_(dispatch_op)
    {
        this->input_ops_.push_back( pfev.template get< Op::InvJac<dim, Op::SideDomain, spacedim>, dim >() );
        this->input_ops_.push_back( pfev.template get< Op::RefGradVector<dim, Op::SideDomain, spacedim>, dim >(fe) );
	}

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        auto inv_jac_value = this->input_ops(0)->result_matrix();    // dim x spacedim
        auto ref_vector_grad = this->input_ops(1)->result_matrix();  // n_sides*dim x spacedim

        uint n_dofs = this->n_dofs();
        uint n_points = ref_vector_grad(0).data_size();
        uint n_patch_sides = ppv.n_elems_;
        uint n_patch_points = ppv.n_points_;

        // Expands inverse jacobian to inv_jac_expd_value
        Eigen::Matrix<ArenaVec<double>, dim, 3> inv_jac_expd_value;
        for (uint i=0; i<dim*3; ++i) {
        	inv_jac_expd_value(i) = ArenaVec<double>( n_patch_points, ppv.patch_arena() );
        	for (uint j=0; j<n_patch_points; ++j)
        	    inv_jac_expd_value(i)(j) = inv_jac_value(i)(j%n_patch_sides);
        }

        // Fill ref shape gradients by q_point. DOF and side_idx
        Eigen::Matrix<ArenaVec<double>, dim, 3> ref_shape_grads_expd;
        for (uint i=0; i<spacedim*dim; ++i) {
            ref_shape_grads_expd(i) = ArenaVec<double>( n_patch_points, ppv.patch_arena() );
        }
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {

            for (uint i_pt=0; i_pt<n_points; ++i_pt) {
                uint i_begin = i_pt * n_patch_sides;
                for (uint i_sd=0; i_sd<n_patch_sides; ++i_sd) {
                    for (uint i_dim=0; i_dim<dim; ++i_dim) {
                        for (uint i_c=0; i_c<spacedim; ++i_c) {
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
template<unsigned int dim, unsigned int spacedim = 3>
class DispatchGradVectorShape : public Op::Side::Base<dim, spacedim> {
public:
    /// Constructor
	DispatchGradVectorShape(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::Side::Base<dim, spacedim>(pfev, {spacedim, spacedim}, OpSizeType::pointOp, fe->n_dofs()), in_op_(nullptr)
    {
        switch (fe->fe_type()) {
            case FEVector:
            {
                in_op_ = new OpGradVectorShape<dim, spacedim>(pfev, fe, *this);
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpGradVectorCovariantShape<dim, spacedim>(pfev, fe, *this);
                break;
            }
            case FEVectorPiola:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorPiola is not implemented yet!\n"); // temporary assert
                //in_op_ = new OpGradVectorPiolaShape<dim, spacedim>(pfev, fe, *this);
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
template<unsigned int dim, unsigned int spacedim = 3>
class OpVectorSymGrad : public Op::VectorSymGradBase<dim, spacedim> {
public:
    /// Constructor
	OpVectorSymGrad(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::VectorSymGradBase<dim, spacedim>(pfev, {spacedim, spacedim}, OpSizeType::pointOp, fe->n_dofs())
    {
        this->bulk_side_ = 1;
        this->input_ops_.push_back( pfev.template get< DispatchGradVectorShape<dim, spacedim>, dim >(fe) );
    }
};

/// Evaluates vector divergence
template<unsigned int dim, unsigned int spacedim = 3>
class OpVectorDivergence : public Op::VectorDivergenceBase<dim, spacedim> {
public:
    /// Constructor
	OpVectorDivergence(PatchFEValues<spacedim> &pfev, std::shared_ptr<FiniteElement<dim>> fe)
    : Op::VectorDivergenceBase<dim, spacedim>(pfev, {1}, OpSizeType::pointOp, fe->n_dofs())
    {
        this->bulk_side_ = 1;
        this->input_ops_.push_back( pfev.template get< DispatchGradVectorShape<dim, spacedim>, dim >(fe) );
    }
};

} // end of namespace Op::Side::Pt

} // end of namespace Side

} // end of namespace Op


#endif /* OP_FUNCTION_HH_ */
