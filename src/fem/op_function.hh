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
    static ElemDomain domain() {
        return ElemDomain::domain_bulk;
    }

    static inline constexpr uint n_nodes(uint dim) {
        return dim+1;
    }

    /// Return number of mesh entities (in this case elements) on patch
    static inline uint n_mesh_entities(PatchPointValues<3> &ppv) {
        return ppv.elem_dim_list_->size();
    }

    /// Return i_n-th node of i_elm-th element stored in PatchPointValues::elem_dim_list_
    static inline NodeAccessor<3> node(PatchPointValues<3> &ppv, unsigned int i_elm, unsigned int i_n) {
        return (*ppv.elem_dim_list_)[i_elm].node(i_n);
    }
};

/// Class used as template type for type resolution Bulk / Side
class SideDomain {
public:
    static ElemDomain domain() {
        return ElemDomain::domain_side;
    }

    static inline constexpr uint n_nodes(uint dim) {
        return dim;
    }

    /// Return number of mesh entities (in this case sides) on patch
    static inline uint n_mesh_entities(PatchPointValues<3> &ppv) {
        return ppv.side_list_.size();
    }

    /// Return i_n-th node of i_elm-th side stored in PatchPointValues::side_list_
    static inline NodeAccessor<3> node(PatchPointValues<3> &ppv, unsigned int i_elm, unsigned int i_n) {
        return ppv.side_list_[i_elm].node(i_n);
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
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class Coords : public PatchOp<spacedim> {
public:
    /// Constructor
    Coords(PatchFEValues<spacedim> &pfev, const Quadrature *quad)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim, Domain::n_nodes(dim)}) {
        this->domain_ = Domain::domain();
    }

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        uint n_elems = Domain::n_mesh_entities(ppv); // number of elements or sides on patch
        this->allocate_result( n_elems, this->patch_arena() );
        auto result = this->result_matrix();

        for (uint i_elm=0; i_elm<n_elems; ++i_elm)
            for (uint i_col=0; i_col<Domain::n_nodes(dim); ++i_col)
                for (uint i_row=0; i_row<spacedim; ++i_row) {
                    result(i_row, i_col)(i_elm) = ( *Domain::node(ppv, i_elm, i_col) )(i_row);
                }
    }

};

/// Evaluates Jacobians on Bulk (Element) / Side
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class Jac : public PatchOp<spacedim> {
public:
    /// Constructor
    Jac(PatchFEValues<spacedim> &pfev, const Quadrature *quad)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim, Domain::n_nodes(dim)-1})
    {
        this->domain_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< Op::Coords<dim, Domain, spacedim>, dim >(quad) );
    }

    void eval() override {
        auto jac_value = this->result_matrix();
        auto coords_value = this->input_ops(0)->result_matrix();
        for (unsigned int i=0; i<spacedim; i++)
            for (unsigned int j=0; j<Domain::n_nodes(dim)-1; j++)
                jac_value(i,j) = coords_value(i,j+1) - coords_value(i,0);
    }
};

/// Evaluates Jacobian determinants on Bulk (Element) / Side
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class JacDet : public PatchOp<spacedim> {
public:
    /// Constructor
	JacDet(PatchFEValues<spacedim> &pfev, const Quadrature *quad)
	: PatchOp<spacedim>(dim, pfev, quad, {1})
	{
        this->domain_ = Domain::domain();
	    this->input_ops_.push_back( pfev.template get< Op::Jac<dim, Domain, spacedim>, dim >(quad) );
	}

    void eval() override {
        auto jac_det_value = this->result_matrix();
        auto jac_value = this->input_ops(0)->result_matrix();
        jac_det_value(0) = eigen_arena_tools::determinant<spacedim, Domain::n_nodes(dim)-1>(jac_value).abs();
    }
};

/// Template specialization of previous: dim=1, domain=Side
template<>
class JacDet<1, Op::SideDomain, 3> : public PatchOp<3> {
public:
    /// Constructor
    JacDet(PatchFEValues<3> &pfev, const Quadrature *quad)
    : PatchOp<3>(1, pfev, quad, {1})
    {
        this->domain_ = Op::SideDomain::domain();
    }

    void eval() override {
        PatchPointValues<3> &ppv = this->ppv();
        uint n_sides = ppv.n_mesh_items();
        this->allocate_result( n_sides, this->patch_arena() );
        auto jac_det_value = this->result_matrix();
        for (uint i=0;i<n_sides; ++i) {
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
	InvJac(PatchFEValues<spacedim> &pfev, const Quadrature *quad)
    : PatchOp<spacedim>(dim, pfev, quad, {dim, spacedim})
    {
        this->domain_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< Op::Jac<dim, Domain, spacedim>, dim >(quad) );
    }

    void eval() override {
        auto inv_jac_value = this->result_matrix();
        auto jac_value = this->input_ops(0)->result_matrix();
        inv_jac_value = eigen_arena_tools::inverse<spacedim, dim>(jac_value);
    }
};


/**
 * Holds common functionality of patch operations.
 */
template<unsigned int spacedim = 3>
class FuncHelper {
public:
    /**
     * Copy reduced data from 'source' to 'target' ArenaVec. Mapping of reduced data is given by 'ppv' data.
     */
    static void fill_reduce_element_data_vec(PatchPointValues<spacedim> &ppv, ArenaVec<double> &source, ArenaVec<double> &target) {
        for (uint i_el=0; i_el<ppv.n_mesh_items(); ++i_el) {
            target( i_el ) = source( ppv.int_table_(patch_elem_on_domain)(i_el) );
        }
    }

    /**
     * Copy ref_side_on_quads data from 'source' to 'target' ArenaVec. Mapping of copied data is given by 'ppv' data.
     */
    static void copy_ref_side_on_quads_data(PatchOp<spacedim> &op, unsigned int i_comp,
            Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> &source,
	        Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> &target) {
        PatchPointValues<spacedim> &ppv = op.ppv();
        uint n_dofs = op.n_dofs();
        uint n_sides = ppv.n_mesh_items();               // number of sides on patch
        uint n_patch_points = n_sides * op.quad_size();  // number of points on patch
        uint n_comp = op.shape()[0];                     // number of components of result

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i_pt=0; i_pt<n_patch_points; ++i_pt)
                target(i_comp,i_dof)(i_pt) = source(ppv.int_table_(ref_side_on_quads)(i_pt),n_comp*i_dof+i_comp)(i_pt / n_sides);
        }
    }

    /**
     * Copy ref_side_on_sides vector data from 'source' to 'target' ArenaVec. Mapping of copied data is given by 'ppv' data.
     */
    static void copy_ref_side_on_sides_vector_data(PatchOp<spacedim> &op, unsigned int dim,
            Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> &source,
	        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> &target)
    {
        PatchPointValues<spacedim> &ppv = op.ppv();
        uint n_sides = ppv.n_mesh_items();
        for (uint i_dof=0; i_dof<op.n_dofs(); ++i_dof) {
            for (uint i_pt=0; i_pt<op.quad_size(); ++i_pt) {
                uint i_begin = i_pt * n_sides;
                for (uint i_sd=0; i_sd<n_sides; ++i_sd) {
                    for (uint i_c=0; i_c<dim; ++i_c) {
                        target(i_c, i_dof)(i_begin + i_sd) = source(ppv.int_table_(ref_side_on_sides)(i_sd), i_dof*dim+i_c)(i_pt);
                    }
                }
            }
        }
    }

    /**
     * Copy ref_side_on_sides tensor data from 'source' to 'target' ArenaVec. Mapping of copied data is given by 'ppv' data.
     */
    static void copy_ref_side_on_sides_tensor_data(PatchOp<spacedim> &op, unsigned int i_dof, unsigned int shape_m, unsigned int shape_n,
            Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> &source,
	        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> &target)
    {
        PatchPointValues<spacedim> &ppv = op.ppv();
        uint n_sides = ppv.n_mesh_items();
        uint n_points = source(0).data_size(); // number of ref points

        for (uint i_pt=0; i_pt<n_points; ++i_pt) {
            uint i_begin = i_pt * n_sides;
            for (uint i_sd=0; i_sd<n_sides; ++i_sd) {
                for (uint i_dim=0; i_dim<shape_m; ++i_dim) {
                    for (uint i_c=0; i_c<shape_n; ++i_c) {
                        target(i_dim, i_c)(i_begin + i_sd) = source(ppv.int_table_(ref_side_on_sides)(i_sd)*shape_m+i_dim, shape_n*i_dof+i_c)(i_pt);
                    }
                }
            }
        }
    }

    /**
     * Copy patch_elem_on_domain data from 'source' to 'target' ArenaVec. Mapping of copied data is given by 'ppv' data.
     */
    static void copy_patch_elem_on_domain_data(PatchOp<spacedim> &op, unsigned int dim,
            Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> &source,
	        Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> &target) {
        PatchPointValues<spacedim> &ppv = op.ppv();
        uint n_elems = ppv.n_mesh_items();          // number of elements on patch
        uint n_elem_points = op.quad_size();        // number of points on element

        for (uint i_c=0; i_c<spacedim*dim; ++i_c) {
            target(i_c) = ArenaVec<double>( n_elems*n_elem_points, op.patch_arena() );
            for (uint i_el=0; i_el<n_elems; ++i_el)
                for (uint i_pt=0; i_pt<n_elem_points; ++i_pt) {
                    target(i_c)( i_pt*n_elems + i_el ) = source(i_c)( ppv.int_table_(patch_elem_on_domain)(i_el) );
                }
        }
    }
private:
    /// Forbidden constructor
    FuncHelper() {}
};


/**
 * Reference data of PtCoords operation
 *
 * See note for PtCoords operation bellow.
 * TODO need specializations for SideDomain<dim> and SideDomain<1>
 */
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class RefBaryCoords : public PatchOp<spacedim> {
public:
    /// Constructor
	RefBaryCoords(PatchFEValues<spacedim> &pfev, const Quadrature *quad)
    : PatchOp<spacedim>(dim, pfev, quad, {dim+1})
    {
        this->domain_ = Domain::domain();
        double n_points = quad->size();
        this->allocate_result(n_points, pfev.asm_arena());
        for (uint i_p=0; i_p<n_points; ++i_p) {
            auto ref_bar_coords = RefElement<dim>::local_to_bary(quad->point<dim>(i_p));
            for (uint i_c=0; i_c<dim+1; ++i_c)
                this->result_(i_c)(i_p) = ref_bar_coords[i_c];
        }
    }

    void eval() override {}
};

/**
 * Evaluates coordinates of quadrature points
 *
 * Important note !!!
 * Usage of this operation is currently in L2 error calculation in DarcyFlow output.
 * It should not be used by any other operation. It should be unified with Field Coords in future.
 */
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class PtCoords : public PatchOp<spacedim> {
public:
    /// Constructor
    PtCoords(PatchFEValues<spacedim> &pfev, const Quadrature *quad)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim})
    {
        this->domain_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< Op::Coords<dim, Domain, spacedim>, dim >( pfev.element_quad(dim) ) );
        this->input_ops_.push_back( pfev.template get< Op::RefBaryCoords<dim, Domain, spacedim>, dim >(quad) );
    }

    void eval() override {
        auto coords_vec_elem = this->input_ops(0)->result_matrix();  // bulk: spacedim x (dim+1), side: spacedim x dim
        auto ref_bary_vec = this->input_ops(1)->result_matrix();     // bulk: dim+1, side: dim

        // Copy coords vector of elements registered on patch, convert to ovec
        PatchPointValues<spacedim> &ppv = this->ppv();
        uint n_elems = ppv.n_mesh_items();
        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> coords_vec(spacedim, Domain::n_nodes(dim));
        Eigen::Matrix<ArenaOVec<double>, spacedim, Domain::n_nodes(dim)> coords_ovec;
        for (uint i=0; i<spacedim*Domain::n_nodes(dim); ++i) {
            coords_vec(i) = ArenaVec<double>( n_elems, this->patch_arena() );
            FuncHelper<spacedim>::fill_reduce_element_data_vec( ppv, coords_vec_elem(i), coords_vec(i) );
            coords_ovec(i) = ArenaOVec<double>( coords_vec(i) );
        }

        // Convert ref_bary_vec to ovec
        Eigen::Vector<ArenaOVec<double>, Eigen::Dynamic> ref_bary_ovec(Domain::n_nodes(dim));
        for (uint i=0; i<Domain::n_nodes(dim); ++i) {
            ref_bary_ovec(i) = ArenaOVec<double>( ref_bary_vec(i) );
        }

        auto result_vec = this->result_matrix();
        Eigen::Vector<ArenaOVec<double>, Eigen::Dynamic> result_ovec = coords_ovec * ref_bary_ovec;
        for (uint i=0; i<spacedim; ++i) {
            result_vec(i) = result_ovec(i).get_vec();
        }
    }
};

/**
 * Fixed operation points weights
 * ElDomain (target) is equivalent with Domain (source)
 */
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class Weights : public PatchOp<spacedim> {
public:
    /// Constructor
    Weights(PatchFEValues<spacedim> &pfev, const Quadrature *quad)
    : PatchOp<spacedim>(dim, pfev, quad, {1})
    {
        this->domain_ = Domain::domain();
        // create result vector of weights operation in assembly arena
        const std::vector<double> &point_weights_vec = quad->get_weights();
        this->allocate_result(point_weights_vec.size(), pfev.asm_arena());
        for (uint i=0; i<point_weights_vec.size(); ++i)
            this->result_(0)(i) = point_weights_vec[i];
    }

    void eval() override {}
};


/**
 * Evaluates JxW on quadrature points
 */
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class JxW : public PatchOp<spacedim> {
public:
    /// Constructor
    JxW(PatchFEValues<spacedim> &pfev, const Quadrature *quad)
    : PatchOp<spacedim>(dim, pfev, quad, {1})
    {
        this->domain_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< Op::Weights<dim, Domain, spacedim>, dim >(quad) );
        this->input_ops_.push_back( pfev.template get< Op::JacDet<dim, Domain, spacedim>, dim >(pfev.element_quad(dim)) );
    }

    void eval() override {
        auto weights_value = this->input_ops(0)->result_matrix();
        auto jac_det_value_long = this->input_ops(1)->result_matrix();

        // Copy InvJac vector of elements registered on patch
        PatchPointValues<spacedim> &ppv = this->ppv();
        uint n_elems = ppv.n_mesh_items();
        ArenaVec<double> jac_det_value( n_elems, this->patch_arena() );
        FuncHelper<spacedim>::fill_reduce_element_data_vec(ppv, jac_det_value_long( 0 ), jac_det_value);

        ArenaOVec<double> weights_ovec( weights_value(0,0) );
        ArenaOVec<double> jac_det_ovec( jac_det_value );
        ArenaOVec<double> jxw_ovec = jac_det_ovec * weights_ovec;
        this->result_(0) = jxw_ovec.get_vec();
    }
};

template<unsigned int dim, unsigned int spacedim>
class JxW<dim, Op::SideDomain, spacedim> : public PatchOp<spacedim> {
public:
    /// Constructor
    JxW(PatchFEValues<spacedim> &pfev, const Quadrature *quad)
    : PatchOp<spacedim>(dim, pfev, quad, {1})
    {
        this->domain_ = Op::SideDomain::domain();
        this->input_ops_.push_back( pfev.template get< Op::Weights<dim, Op::SideDomain, spacedim>, dim >(quad) );
        this->input_ops_.push_back( pfev.template get< Op::JacDet<dim, Op::SideDomain, spacedim>, dim >(pfev.element_quad(dim)) );
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

/// Evaluates normal vector on side quadrature points
template<unsigned int dim, unsigned int spacedim = 3>
class NormalVec : public PatchOp<spacedim> {
public:
    /// Constructor
    NormalVec(PatchFEValues<spacedim> &pfev, const Quadrature *quad)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim})
    {
        this->domain_ = Op::SideDomain::domain();
        this->input_ops_.push_back( pfev.template get< Op::InvJac<dim, Op::BulkDomain, spacedim>, dim >(quad) );
    }

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        auto normal_value = this->result_matrix();
        auto inv_jac_value_elem = this->input_ops(0)->result_matrix(); // returns vector of inverse jacobians of all elements registered on patch

        // Copy InvJac vector of sides registered on patch
        uint n_sides = ppv.n_mesh_items();
        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> inv_jac_value(dim, spacedim);
        for (uint i=0; i<dim*spacedim; ++i) {
            inv_jac_value(i) = ArenaVec<double>( n_sides, this->patch_arena() );
        }
        for (uint i_c=0; i_c<dim*spacedim; ++i_c) {
            FuncHelper<spacedim>::fill_reduce_element_data_vec( ppv, inv_jac_value_elem(i_c), inv_jac_value(i_c) );
        }

        normal_value = inv_jac_value.transpose() * RefElement<dim>::normal_vector_array( ppv.int_table_(ref_side_on_sides) );

        ArenaVec<double> norm_vec( n_sides, this->patch_arena() );
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

/// Fixed operation of scalar shape reference values
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class RefScalar : public PatchOp<spacedim> {
public:
    /// Constructor
    RefScalar(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {1}, fe->n_dofs())
    {
        this->domain_ = Domain::domain();
        uint n_points = quad->size();

        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_scalar_value = this->result_matrix();
        for (unsigned int i_p = 0; i_p < n_points; i_p++)
            for (unsigned int i_dof = 0; i_dof < this->n_dofs_; i_dof++)
                ref_scalar_value(i_dof)(i_p) = fe->shape_value(i_dof, quad->point<dim>(i_p));
    }

    void eval() override {}
};

/// Template specialization of previous: Domain=SideDomain
template<unsigned int dim, unsigned int spacedim>
class RefScalar<dim, Op::SideDomain, spacedim> : public PatchOp<spacedim> {
public:
    /// Constructor
    RefScalar(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {dim+1}, fe->n_dofs())
    {
        this->domain_ = Op::SideDomain::domain();
        uint n_points = quad->size();

        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_scalar_value = this->result_matrix();
        for (unsigned int s=0; s<dim+1; ++s) {
            Quadrature side_q = quad->make_from_side<dim>(s);
            for (unsigned int i_p = 0; i_p < n_points; i_p++)
                for (unsigned int i_dof = 0; i_dof < this->n_dofs_; i_dof++)
                    ref_scalar_value(s, i_dof)(i_p) = fe->shape_value(i_dof, side_q.point<dim>(i_p));
        }
    }

    void eval() override {}
};

/// Fixed operation of vector shape reference values
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class RefVector : public PatchOp<spacedim> {
public:
    /// Constructor
    RefVector(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {fe->n_components()}, fe->n_dofs())
    {
        this->domain_ = Domain::domain();
        uint n_points = quad->size();

        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_vector_value = this->result_matrix();

        for (uint i_dof=0; i_dof<this->n_dofs_; ++i_dof)
            for (uint i_p=0; i_p<n_points; ++i_p)
                for (uint c=0; c<fe->n_components(); ++c)
                    ref_vector_value(c, i_dof)(i_p) = fe->shape_value(i_dof, quad->point<dim>(i_p), c);
    }

    void eval() override {}
};

/// Template specialization of previous: Domain=SideDomain
template<unsigned int dim, unsigned int spacedim>
class RefVector<dim, Op::SideDomain, spacedim> : public PatchOp<spacedim> {
public:
    /// Constructor
    RefVector(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {dim+1, fe->n_components()}, fe->n_dofs())
    {
        this->domain_ = Op::SideDomain::domain();
        uint n_points = quad->size();

        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_vector_value = this->result_matrix();

        for (unsigned int s=0; s<dim+1; ++s) {
            Quadrature side_q = quad->make_from_side<dim>(s);
            for (unsigned int i_p = 0; i_p < n_points; i_p++)
                for (unsigned int i_dof = 0; i_dof < this->n_dofs_; i_dof++)
                    for (uint c=0; c<fe->n_components(); ++c)
                        ref_vector_value(s,fe->n_components()*i_dof+c)(i_p) = fe->shape_value(i_dof, side_q.point<dim>(i_p), c);
        }
    }

    void eval() override {}
};

/// Fixed operation of vector shape reference values
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class RefTensor : public PatchOp<spacedim> {
public:
    /// Constructor
	RefTensor(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim, spacedim}, fe->n_dofs())
    {
        this->domain_ = Domain::domain();
        uint n_points = quad->size();
        uint n_comp = spacedim;

        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_tensor_value = this->result_matrix();

        for (uint i_col=0; i_col<n_comp; ++i_col) {
            for (uint i_row=0; i_row<n_comp; ++i_row)
                for (uint i_dof=0; i_dof<this->n_dofs_; ++i_dof)
                    for (uint i_p=0; i_p<n_points; ++i_p)
                        ref_tensor_value(i_row,n_comp*i_dof+i_col)(i_p) = fe->shape_value(i_dof, quad->point<dim>(i_p), i_col*n_comp+i_row);
        }
    }

    void eval() override {}
};

/// Fixed operation of gradient scalar reference values
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class RefGradScalar : public PatchOp<spacedim> {
public:
    /// Constructor
    RefGradScalar(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {dim, 1}, fe->n_dofs())
    {
        this->domain_ = Domain::domain();
        uint n_points = quad->size();

        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_scalar_value = this->result_matrix();
        for (uint i_row=0; i_row<ref_scalar_value.rows(); ++i_row)
            for (uint i_dof=0; i_dof<this->n_dofs_; ++i_dof)
                for (uint i_p=0; i_p<n_points; ++i_p)
                    ref_scalar_value(i_row, i_dof)(i_p) = fe->shape_grad(i_dof, quad->point<dim>(i_p))[i_row];
    }

    void eval() override {}
};

/// Template specialization of previous: Domain=SideDomain
template<unsigned int dim, unsigned int spacedim>
class RefGradScalar<dim, Op::SideDomain, spacedim> : public PatchOp<spacedim> {
public:
    /// Constructor
    RefGradScalar(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {dim+1, dim}, fe->n_dofs())
    {
        this->domain_ = Op::SideDomain::domain();
        uint n_points = quad->size();

        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_scalar_value = this->result_matrix();
        for (unsigned int s=0; s<dim+1; ++s) {
            Quadrature side_q = quad->make_from_side<dim>(s);
            for (uint i_dof=0; i_dof<this->n_dofs_; ++i_dof)
                for (uint i_p=0; i_p<n_points; ++i_p)
                    for (uint c=0; c<dim; ++c)
                        ref_scalar_value(s,dim*i_dof+c)(i_p) = fe->shape_grad(i_dof, side_q.point<dim>(i_p))[c];
        }
    }

    void eval() override {}
};

/// Fixed operation of gradient scalar reference values
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class RefGradVector : public PatchOp<spacedim> {
public:
    /// Constructor
	RefGradVector(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {dim, fe->n_components()}, fe->n_dofs())
    {
        this->domain_ = Domain::domain();
        uint n_points = quad->size();
        uint n_comp = fe->n_components();

        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_vector_value = this->result_matrix();
        for (uint i_c=0; i_c<n_comp; ++i_c) {
            for (uint i_dim=0; i_dim<dim; ++i_dim)
                for (uint i_dof=0; i_dof<this->n_dofs_; ++i_dof)
                    for (uint i_p=0; i_p<n_points; ++i_p)
                        ref_vector_value(i_dim,n_comp*i_dof+i_c)(i_p) = fe->shape_grad(i_dof, quad->point<dim>(i_p), i_c)[i_dim];
        }
    }

    void eval() override {}
};

/// Template specialization of previous: Domain=SideDomain
template<unsigned int dim, unsigned int spacedim>
class RefGradVector<dim, Op::SideDomain, spacedim> : public PatchOp<spacedim> {
public:
    /// Constructor
	RefGradVector(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {(dim+1)*dim, fe->n_components()}, fe->n_dofs())
    {
        this->domain_ = Op::SideDomain::domain();
        uint n_points = quad->size();
        uint n_sides = dim+1;
        uint n_comp = fe->n_components();

        this->allocate_result(n_points, pfev.asm_arena());
        auto ref_vector_value = this->result_matrix();
        for (uint i_sd=0; i_sd<n_sides; ++i_sd) {
            Quadrature side_q = quad->make_from_side<dim>(i_sd);
            for (uint i_c=0; i_c<n_comp; ++i_c)
                for (uint i_dim=0; i_dim<dim; ++i_dim)
                    for (uint i_dof=0; i_dof<this->n_dofs_; ++i_dof)
                        for (uint i_p=0; i_p<n_points; ++i_p)
                            ref_vector_value(i_sd*dim+i_dim, n_comp*i_dof+i_c)(i_p) = fe->shape_grad(i_dof, side_q.point<dim>(i_p), i_c)[i_dim];
        }
    }

    void eval() override {}
};

/// Evaluates scalar shape values
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class ScalarShape : public PatchOp<spacedim> {
public:
    /// Constructor
    ScalarShape(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {1}, fe->n_dofs())
    {
        ASSERT_EQ(fe->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape must be FEScalar!\n");
        this->domain_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< Op::RefScalar<dim, Domain, spacedim>, dim >(quad, fe) );
    }

    void eval() override {
        auto ref_vec = this->input_ops(0)->result_matrix();
        auto result_vec = this->result_matrix();

        uint n_dofs = this->n_dofs();
        uint n_elem = this->ppv().n_mesh_items();

        ArenaVec<double> elem_vec(n_elem, this->patch_arena());
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

/// Template specialization of previous: Domain=SideDomain
template<unsigned int dim, unsigned int spacedim>
class ScalarShape<dim, Op::SideDomain, spacedim> : public PatchOp<spacedim> {
public:
    /// Constructor
    ScalarShape(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {1}, fe->n_dofs())
    {
        ASSERT_EQ(fe->fe_type(), FEType::FEScalar).error("Type of FiniteElement of scalar_shape must be FEScalar!\n");
        this->domain_ = Op::SideDomain::domain();
        this->input_ops_.push_back( pfev.template get< Op::RefScalar<dim, Op::SideDomain, spacedim>, dim >(quad, fe) );
    }

    void eval() override {
        uint n_patch_points = this->ppv().n_mesh_items() * this->quad_size(); // number of points on patch
        this->allocate_result(n_patch_points, this->patch_arena());

        auto ref_vec = this->input_ops(0)->result_matrix();
        auto result_vec = this->result_matrix();

        FuncHelper<spacedim>::copy_ref_side_on_quads_data(*this, 0, ref_vec, result_vec);
    }
};

/// Evaluates vector values (FEType == FEVector)
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class VectorShape : public PatchOp<spacedim> {
public:
    /// Constructor
	VectorShape(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim}, fe->n_dofs()), dispatch_op_(dispatch_op)
    {
        this->domain_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< Op::RefVector<dim, Domain, spacedim>, dim >(quad, fe) );
	}

    void eval() override {
        auto ref_shape_vec = this->input_ops(0)->result_matrix();
        auto result_vec = dispatch_op_.result_matrix();

        uint n_dofs = this->n_dofs();
        uint n_elem = this->ppv().n_mesh_items();

        ArenaVec<double> elem_vec(n_elem, this->patch_arena());
        for (uint i=0; i<n_elem; ++i) {
            elem_vec(i) = 1.0;
        }
        ArenaOVec<double> elem_ovec(elem_vec);

        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> ref_shape_ovec(spacedim, n_dofs);
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

/// Template specialization of previous: Domain=SideDomain)
template<unsigned int dim, unsigned int spacedim>
class VectorShape<dim, Op::SideDomain, spacedim> : public PatchOp<spacedim> {
public:
    /// Constructor
	VectorShape(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim}, fe->n_dofs()), dispatch_op_(dispatch_op)
    {
        this->domain_ = Op::SideDomain::domain();
        this->input_ops_.push_back( pfev.template get< Op::RefVector<dim, Op::SideDomain, spacedim>, dim >(quad, fe) );
	}

    void eval() override {
        uint n_patch_points = this->ppv().n_mesh_items() * this->quad_size();
        dispatch_op_.allocate_result(n_patch_points, this->patch_arena());

        auto ref_shape_vec = this->input_ops(0)->result_matrix();  // dim+1 x spacedim
        auto result_vec = dispatch_op_.result_matrix();            // spacdim x 1

        for (uint c=0; c<spacedim; c++)
            FuncHelper<spacedim>::copy_ref_side_on_quads_data(*this, c, ref_shape_vec, result_vec);
    }

private:
    PatchOp<spacedim> &dispatch_op_;
};

/// Evaluates vector values (FEType == FEVectorPiola)
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class VectorPiolaShape : public PatchOp<spacedim> {
public:
    /// Constructor
	VectorPiolaShape(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim}, fe->n_dofs()), dispatch_op_(dispatch_op)
    {
        this->domain_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< Op::RefVector<dim, Domain, spacedim>, dim >(quad, fe) );            // input_ops_[0] ... RefVector dim x n_dofs
        this->input_ops_.push_back( pfev.template get< Op::Jac<dim, Domain, spacedim>, dim >(pfev.element_quad(dim)) );    // input_ops_[1] ... Jac       spacedim x dim
        this->input_ops_.push_back( pfev.template get< Op::JacDet<dim, Domain, spacedim>, dim >(pfev.element_quad(dim)) ); // input_ops_[2] ... JacDet    1
	}

    void eval() override {
        auto ref_shape_vec = this->input_ops(0)->result_matrix();
        auto jac_vec_elem = this->input_ops(1)->result_matrix();
        auto jac_det_vec_elem = this->input_ops(2)->result_matrix();
        auto result_vec = dispatch_op_.result_matrix();

        uint n_dofs = this->n_dofs();

        // Copy Jac and JacDet vector of elements registered on patch
        PatchPointValues<spacedim> &ppv = this->ppv();
        uint n_elems = ppv.n_mesh_items();
        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> jac_vec(spacedim, dim);
        for (uint i=0; i<spacedim*dim; ++i) {
            jac_vec(i) = ArenaVec<double>( n_elems, this->patch_arena() );
        }
        for (uint i_c=0; i_c<spacedim*dim; ++i_c) {
            FuncHelper<spacedim>::fill_reduce_element_data_vec( ppv, jac_vec_elem(i_c), jac_vec(i_c) );
        }
        ArenaVec<double> jac_det_vec( n_elems, this->patch_arena() );
        FuncHelper<spacedim>::fill_reduce_element_data_vec( ppv, jac_det_vec_elem(0), jac_det_vec );

        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> jac_div_det_vec = jac_vec / jac_det_vec;
        Eigen::Matrix<ArenaOVec<double>, spacedim, dim> jac_div_det_ovec;
        for (uint c=0; c<spacedim*dim; ++c) {
            jac_div_det_ovec(c) = ArenaOVec<double>(jac_div_det_vec(c));
        }

        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> ref_shape_ovec(dim, n_dofs);
        for (uint c=0; c<dim*n_dofs; ++c) {
            ref_shape_ovec(c) = ArenaOVec<double>(ref_shape_vec(c));
        }

        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> result_ovec = jac_div_det_ovec * ref_shape_ovec;
        for (uint c=0; c<spacedim*n_dofs; ++c) {
            result_vec(c) = result_ovec(c).get_vec();
        }
    }

private:
    PatchOp<spacedim> &dispatch_op_;
};

/// Template specialization of previous: Domain=SideDomain)
template<unsigned int dim, unsigned int spacedim>
class VectorPiolaShape<dim, Op::SideDomain, spacedim> : public PatchOp<spacedim> {
public:
    /// Constructor
	VectorPiolaShape(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim}, fe->n_dofs()), dispatch_op_(dispatch_op)
    {
        this->domain_ = Op::SideDomain::domain();
        this->input_ops_.push_back( pfev.template get< Op::RefVector<dim, Op::SideDomain, spacedim>, dim >(quad, fe) );            // input_ops_[0] ... dim x n_dofs
        this->input_ops_.push_back( pfev.template get< Op::Jac<dim, Op::BulkDomain, spacedim>, dim >(pfev.element_quad(dim)) );    // input_ops_[1] ... spacedim x dim
        this->input_ops_.push_back( pfev.template get< Op::JacDet<dim, Op::BulkDomain, spacedim>, dim >(pfev.element_quad(dim)) ); // input_ops_[2] ... 1
	}

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        uint n_dofs = this->n_dofs();
        uint n_sides = ppv.n_mesh_items();
        uint n_patch_points = n_sides * this->quad_size();
        uint n_points_per_side = this->quad_size();
        dispatch_op_.allocate_result(n_patch_points, this->patch_arena());

        auto ref_shape_vec = this->input_ops(0)->result_matrix();  // dim+1 x n_dofs
        auto jac_vec_elem = this->input_ops(1)->result_matrix();
        auto jac_det_vec_elem = this->input_ops(2)->result_matrix();
        auto result_vec = dispatch_op_.result_matrix();            // spacdim x 1

        // Copy Jac and JacDet vector of elements registered on patch
        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> jac_vec(spacedim, dim);
        for (uint i=0; i<spacedim*dim; ++i) {
            jac_vec(i) = ArenaVec<double>( n_sides, this->patch_arena() );
        }
        for (uint i_c=0; i_c<spacedim*dim; ++i_c) {
            FuncHelper<spacedim>::fill_reduce_element_data_vec( ppv, jac_vec_elem(i_c), jac_vec(i_c) );
        }
        ArenaVec<double> jac_det_vec( n_sides, this->patch_arena() );
        FuncHelper<spacedim>::fill_reduce_element_data_vec( ppv, jac_det_vec_elem(0), jac_det_vec );

        // Computes jacobian / determinant and converts result to ArenaOVec
        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> jac_div_det_vec = jac_vec / jac_det_vec;
        Eigen::Matrix<ArenaOVec<double>, spacedim, dim> jac_div_det_ovec;
        for (uint c=0; c<spacedim*dim; ++c) {
            jac_div_det_ovec(c) = ArenaOVec(jac_div_det_vec(c));
        }

        // Computes expand vector of previous result (jacobian / determinant)
        ArenaVec<double> side_points_vec(n_points_per_side, this->patch_arena());
        for (uint i=0; i<n_points_per_side; ++i) {
           side_points_vec(i) = 1.0;
        }
        ArenaOVec<double> side_points_ovec(side_points_vec);
        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> expand_jac_div_det_ovec = jac_div_det_ovec * side_points_ovec;
        Eigen::Matrix<ArenaVec<double>, spacedim, dim> expand_jac_div_det_vec;
        for (uint c=0; c<spacedim*this->dim_; ++c) {
            expand_jac_div_det_vec(c) = expand_jac_div_det_ovec(c).get_vec();
        }

        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> expand_ref_vec(dim, n_dofs);
        for (uint c=0; c<dim*n_dofs; c++)
            expand_ref_vec(c) = ArenaVec<double>(n_patch_points, this->patch_arena());

        FuncHelper<spacedim>::copy_ref_side_on_sides_vector_data(*this,dim, ref_shape_vec, expand_ref_vec);

        // computes operation result
        result_vec = expand_jac_div_det_vec * expand_ref_vec;
    }

private:
    PatchOp<spacedim> &dispatch_op_;
};

// class OpVectorCovariantShape
// class OpVectorPiolaShape

/// Dispatch class of vector values
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class DispatchVectorShape : public PatchOp<spacedim> {
public:
    /// Constructor
    DispatchVectorShape(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim}, fe->n_dofs()), in_op_(nullptr)
    {
        this->domain_ = Domain::domain();
        switch (fe->fe_type()) {
            case FEVector:
            {
                in_op_ = new VectorShape<dim, Domain, spacedim>(pfev, quad, fe, *this);
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                //in_op_ = new VectorCovariantShape<dim, Domain, spacedim>(pfev, quad, fe, *this);
                break;
            }
            case FEVectorPiola:
            {
                in_op_ = new VectorPiolaShape<dim, Domain, spacedim>(pfev, quad, fe, *this);
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

/// Evaluates vector values (FEType == FEVector)
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class TensorShape : public PatchOp<spacedim> {
public:
    /// Constructor
	TensorShape(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim, spacedim}, fe->n_dofs())
    {
        this->domain_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< Op::RefTensor<dim, Domain, spacedim>, dim >(quad, fe) );
	}

    void eval() override {
        auto ref_shape_vec = this->input_ops(0)->result_matrix();
        auto result_vec = this->result_matrix();

        uint n_dofs = this->n_dofs();
        uint n_elem = this->ppv().n_mesh_items();

        ArenaVec<double> elem_vec(n_elem, this->patch_arena());
        for (uint i=0; i<n_elem; ++i) {
            elem_vec(i) = 1.0;
        }
        ArenaOVec<double> elem_ovec(elem_vec);

        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> ref_shape_ovec(spacedim, spacedim);
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint c=0; c<spacedim*spacedim; ++c) {
                ref_shape_ovec(c) = ArenaOVec(ref_shape_vec(i_dof * spacedim*spacedim + c));
            }

            Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> result_ovec = elem_ovec * ref_shape_ovec;
            for (uint c=0; c<spacedim*spacedim; ++c)
                result_vec(c) = result_ovec(i_dof * spacedim*spacedim + c).get_vec();
        }
    }
};

/// Evaluates gradient scalar values
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class GradScalarShape : public PatchOp<spacedim> {
public:
    /// Constructor
    GradScalarShape(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim, 1}, fe->n_dofs())
    {
        ASSERT_EQ(fe->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape must be FEScalar!\n");
        this->domain_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< Op::InvJac<dim, Op::BulkDomain, spacedim>, dim >(pfev.element_quad(dim)) );
        this->input_ops_.push_back( pfev.template get< Op::RefGradScalar<dim, Domain, spacedim>, dim >(quad, fe) );
    }

    void eval() override {
        auto inv_jac_vec_elem = this->input_ops(0)->result_matrix();    // dim x spacedim=3
        auto ref_grads_vec = this->input_ops(1)->result_matrix();       // dim x n_dofs

        uint n_dofs = this->n_dofs();

        // Copy InvJac vector of elements registered on patch
        PatchPointValues<spacedim> &ppv = this->ppv();
        uint n_elems = ppv.n_mesh_items();
        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> inv_jac_vec(dim, spacedim);
        for (uint i=0; i<dim*spacedim; ++i) {
            inv_jac_vec(i) = ArenaVec<double>( n_elems, this->patch_arena() );
        }
        for (uint i_c=0; i_c<dim*spacedim; ++i_c) {
            FuncHelper<spacedim>::fill_reduce_element_data_vec( ppv, inv_jac_vec_elem(i_c), inv_jac_vec(i_c) );
        }

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

/// Template specialization of previous: Domain=SideDomain
template<unsigned int dim, unsigned int spacedim>
class GradScalarShape<dim, Op::SideDomain, spacedim> : public PatchOp<spacedim> {
public:
    /// Constructor
    GradScalarShape(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim, 1}, fe->n_dofs())
    {
        ASSERT_EQ(fe->fe_type(), FEType::FEScalar).error("Type of FiniteElement of grad_scalar_shape must be FEScalar!\n");
        this->domain_ = Op::SideDomain::domain();
        this->input_ops_.push_back( pfev.template get< Op::InvJac<dim, Op::BulkDomain, spacedim>, dim >(pfev.element_quad(dim)) );
        this->input_ops_.push_back( pfev.template get< Op::RefGradScalar<dim, Op::SideDomain, spacedim>, dim >(quad, fe) );
    }

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        auto ref_shape_grads = this->input_ops(1)->result_matrix();
        auto grad_scalar_shape_value = this->result_matrix();

        uint n_dofs = this->n_dofs();
        uint n_sides = ppv.n_mesh_items();
        uint n_patch_points = n_sides * this->quad_size();

        this->allocate_result(n_patch_points, this->patch_arena());

        // Expands inverse jacobian to inv_jac_expd_value
        auto inv_jac_value = this->input_ops(0)->result_matrix();
        Eigen::Matrix<ArenaVec<double>, dim, spacedim> inv_jac_expd_value;
        Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> inv_jac_expd_map(inv_jac_expd_value.data(), dim, spacedim);
        FuncHelper<spacedim>::copy_patch_elem_on_domain_data(*this, dim, inv_jac_value, inv_jac_expd_map);

        // Fill ref shape gradients by q_point. DOF and side_idx
        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> ref_shape_grads_expd(dim, n_dofs);
        for (uint i=0; i<dim*n_dofs; ++i) {
            ref_shape_grads_expd(i) = ArenaVec<double>( n_patch_points, this->patch_arena() );
        }
        FuncHelper<spacedim>::copy_ref_side_on_sides_vector_data(*this, dim, ref_shape_grads, ref_shape_grads_expd);

        // computes operation result
        grad_scalar_shape_value = inv_jac_expd_value.transpose() * ref_shape_grads_expd;
    }
};

/// Evaluates vector values (FEType == FEVector)
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class GradVectorShape : public PatchOp<spacedim> {
public:
    /// Constructor
	GradVectorShape(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim, spacedim}, fe->n_dofs()), dispatch_op_(dispatch_op)
    {
        this->domain_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< Op::InvJac<dim, Op::BulkDomain, spacedim>, dim >(pfev.element_quad(dim)) );
        this->input_ops_.push_back( pfev.template get< Op::RefGradVector<dim, Domain, spacedim>, dim >(quad, fe) );
	}

    void eval() override {
	    auto inv_jac_vec_elem = this->input_ops(0)->result_matrix();   // dim x spacedim
        auto ref_grads_vec = this->input_ops(1)->result_matrix();      // dim x spacedim
        auto result_vec = dispatch_op_.result_matrix();                // spacedim x spacedim

        uint n_dofs = this->n_dofs();

        // Copy InvJac vector of elements registered on patch
        PatchPointValues<spacedim> &ppv = this->ppv();
        uint n_elems = ppv.n_mesh_items();
        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> inv_jac_vec(dim, spacedim);
        for (uint i=0; i<dim*spacedim; ++i) {
            inv_jac_vec(i) = ArenaVec<double>( n_elems, this->patch_arena() );
        }
        for (uint i_c=0; i_c<dim*spacedim; ++i_c) {
            FuncHelper<spacedim>::fill_reduce_element_data_vec( ppv, inv_jac_vec_elem(i_c), inv_jac_vec(i_c) );
        }

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

/// Template specialization of previous: Domain=SideDomain)
template<unsigned int dim, unsigned int spacedim>
class GradVectorShape<dim, Op::SideDomain, spacedim> : public PatchOp<spacedim> {
public:
    /// Constructor
	GradVectorShape(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim, spacedim}, fe->n_dofs()), dispatch_op_(dispatch_op)
    {
        this->domain_ = Op::SideDomain::domain();
        this->input_ops_.push_back( pfev.template get< Op::InvJac<dim, Op::BulkDomain, spacedim>, dim >(pfev.element_quad(dim)) );
        this->input_ops_.push_back( pfev.template get< Op::RefGradVector<dim, Op::SideDomain, spacedim>, dim >(quad, fe) );
	}

    void eval() override {
        PatchPointValues<spacedim> &ppv = this->ppv();
        auto inv_jac_value = this->input_ops(0)->result_matrix();    // dim x spacedim
        auto ref_vector_grad = this->input_ops(1)->result_matrix();  // n_sides*dim x spacedim

        uint n_dofs = this->n_dofs();
        uint n_patch_sides = ppv.n_mesh_items();
        uint n_patch_points = n_patch_sides * this->quad_size();

        // Expands inverse jacobian to inv_jac_expd_value
        Eigen::Matrix<ArenaVec<double>, dim, spacedim> inv_jac_expd_value;
        Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> inv_jac_expd_map(inv_jac_expd_value.data(), dim, spacedim);
        FuncHelper<spacedim>::copy_patch_elem_on_domain_data(*this, dim, inv_jac_value, inv_jac_expd_map);

        // Fill ref shape gradients by q_point. DOF and side_idx
        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> ref_shape_grads_expd(dim, spacedim);
        for (uint i=0; i<spacedim*dim; ++i) {
            ref_shape_grads_expd(i) = ArenaVec<double>( n_patch_points, this->patch_arena() );
        }
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            // prepare copy of reference data by indices of sides on elements
            FuncHelper<spacedim>::copy_ref_side_on_sides_tensor_data(*this, i_dof, dim, spacedim, ref_vector_grad, ref_shape_grads_expd);

            // computes operation result
            Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > res_submat = dispatch_op_.result_sub_matrix(i_dof);
            res_submat = (inv_jac_expd_value.transpose() * ref_shape_grads_expd).transpose();
        }
    }

private:
    PatchOp<spacedim> &dispatch_op_;
};

/// Evaluates vector values (FEType == FEVectorPiola)
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class GradVectorPiolaShape : public PatchOp<spacedim> {
public:
    /// Constructor
	GradVectorPiolaShape(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim, spacedim}, fe->n_dofs()), dispatch_op_(dispatch_op)
    {
        this->domain_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< Op::InvJac<dim, Op::BulkDomain, spacedim>, dim >(pfev.element_quad(dim)) );
        this->input_ops_.push_back( pfev.template get< Op::RefGradVector<dim, Domain, spacedim>, dim >(quad, fe) );
        this->input_ops_.push_back( pfev.template get< Op::Jac<dim, Op::BulkDomain, spacedim>, dim >(pfev.element_quad(dim)) );
        this->input_ops_.push_back( pfev.template get< Op::JacDet<dim, Op::BulkDomain, spacedim>, dim >(pfev.element_quad(dim)) );
	}

    void eval() override {
        auto inv_jac_vec_elem = this->input_ops(0)->result_matrix();   // dim x spacedim
        auto ref_grads_vec = this->input_ops(1)->result_matrix();      // dim x dim
        auto jac_vec_elem = this->input_ops(2)->result_matrix();       // spacedim x dim
        auto jac_det_vec_elem = this->input_ops(3)->result_matrix();   // 1
        auto result_vec = dispatch_op_.result_matrix();                // spacedim x spacedim

        PatchPointValues<spacedim> &ppv = this->ppv();
        uint n_elems = ppv.n_mesh_items();
        uint n_dofs = this->n_dofs();

        // Copy InvJac and JacDet vector of elements registered on patch
        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> inv_jac_vec(dim, spacedim);
        for (uint i=0; i<dim*spacedim; ++i) {
            inv_jac_vec(i) = ArenaVec<double>( n_elems, this->patch_arena() );
        }
        for (uint i_c=0; i_c<spacedim*dim; ++i_c) {
            FuncHelper<spacedim>::fill_reduce_element_data_vec( ppv, inv_jac_vec_elem(i_c), inv_jac_vec(i_c) );
        }
        ArenaVec<double> jac_det_vec( n_elems, this->patch_arena() );
        FuncHelper<spacedim>::fill_reduce_element_data_vec( ppv, jac_det_vec_elem(0), jac_det_vec );

        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> inv_jac_div_det_vec = inv_jac_vec / jac_det_vec;
        Eigen::Matrix<ArenaOVec<double>, dim, spacedim> inv_jac_div_det_ovec;
        for (uint c=0; c<dim*spacedim; ++c) {
            inv_jac_div_det_ovec(c) = ArenaOVec<double>(inv_jac_div_det_vec(c));
        }

        Eigen::Matrix<ArenaVec<double>, spacedim, dim> jac_vec;
        Eigen::Map<Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic>> jac_vec_map(jac_vec.data(), spacedim, dim);
        FuncHelper<spacedim>::copy_patch_elem_on_domain_data(*this, dim, jac_vec_elem, jac_vec_map);

        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> ref_grads_ovec(dim, dim);
        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            for (uint i=0; i<dim*dim; ++i) {
                ref_grads_ovec(i) = ArenaOVec(ref_grads_vec(i_dof*dim*dim + i));
            }

            Eigen::Matrix<ArenaOVec<double>, spacedim, dim> inv_jac_multi_ref_ovec = inv_jac_div_det_ovec.transpose() * ref_grads_ovec;
            Eigen::Matrix<ArenaVec<double>, spacedim, dim> inv_jac_multi_ref_vec; // converts components ArenaOVec to ArenaVec
            for (uint i=0; i<spacedim*dim; ++i) {
                inv_jac_multi_ref_vec(i) = inv_jac_multi_ref_ovec(i);
            }

            Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> sub_result_vec = inv_jac_multi_ref_vec * jac_vec.transpose();
            for (uint i=0; i<spacedim; ++i) {
                for (uint j=0; j<spacedim; ++j) {
                    result_vec(j,i+spacedim*i_dof) = sub_result_vec(i,j);
                }
            }
        }
    }

private:
    PatchOp<spacedim> &dispatch_op_;
};

/// Template specialization of previous: Domain=SideDomain)
template<unsigned int dim, unsigned int spacedim>
class GradVectorPiolaShape<dim, Op::SideDomain, spacedim> : public PatchOp<spacedim> {
public:
    /// Constructor
	GradVectorPiolaShape(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe, PatchOp<spacedim> &dispatch_op)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim, spacedim}, fe->n_dofs()), dispatch_op_(dispatch_op)
    {
        this->domain_ = Op::SideDomain::domain();
        this->input_ops_.push_back( pfev.template get< Op::InvJac<dim, Op::BulkDomain, spacedim>, dim >(pfev.element_quad(dim)) );
        this->input_ops_.push_back( pfev.template get< Op::RefGradVector<dim, Op::SideDomain, spacedim>, dim >(quad, fe) );
        this->input_ops_.push_back( pfev.template get< Op::Jac<dim, Op::BulkDomain, spacedim>, dim >(pfev.element_quad(dim)) );
        this->input_ops_.push_back( pfev.template get< Op::JacDet<dim, Op::BulkDomain, spacedim>, dim >(pfev.element_quad(dim)) );
	}

    void eval() override {
        auto inv_jac_vec_elem = this->input_ops(0)->result_matrix();   // dim x spacedim
        auto ref_grads_vec = this->input_ops(1)->result_matrix();      // dim x dim
        auto jac_vec_elem = this->input_ops(2)->result_matrix();       // spacedim x dim
        auto jac_det_vec_elem = this->input_ops(3)->result_matrix();   // 1

        PatchPointValues<spacedim> &ppv = this->ppv();
        uint n_sides = ppv.n_mesh_items();
        uint n_dofs = this->n_dofs();
        uint n_patch_points = n_sides * this->quad_size();
        uint n_points_per_side = this->quad_size();

        // Copy InvJac, Jac and JacDet vector of elements registered on patch
        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> inv_jac_vec(dim, spacedim);
        for (uint i_c=0; i_c<dim*spacedim; ++i_c) {
            inv_jac_vec(i_c) = ArenaVec<double>( n_sides, this->patch_arena() );
            FuncHelper<spacedim>::fill_reduce_element_data_vec( ppv, inv_jac_vec_elem(i_c), inv_jac_vec(i_c) );
        }
        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> jac_vec(spacedim, dim);
        Eigen::Matrix<ArenaOVec<double>, spacedim, dim> jac_ovec;
        for (uint i_c=0; i_c<spacedim*dim; ++i_c) {
            jac_vec(i_c) = ArenaVec<double>( n_sides, this->patch_arena() );
            FuncHelper<spacedim>::fill_reduce_element_data_vec( ppv, jac_vec_elem(i_c), jac_vec(i_c) );
            jac_ovec(i_c) = ArenaOVec<double>(jac_vec(i_c));
        }
        ArenaVec<double> jac_det_vec( n_sides, this->patch_arena() );
        FuncHelper<spacedim>::fill_reduce_element_data_vec( ppv, jac_det_vec_elem(0), jac_det_vec );

        // Compute InvJac / Determinant, convert to ArenaOvec
        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> inv_jac_div_det_vec = inv_jac_vec / jac_det_vec;
        Eigen::Matrix<ArenaOVec<double>, dim, spacedim> inv_jac_div_det_ovec;
        for (uint c=0; c<dim*spacedim; ++c) {
            inv_jac_div_det_ovec(c) = ArenaOVec<double>(inv_jac_div_det_vec(c));
        }

        // Computes expand vector of previous result (inv_jac / determinant) and expand vector of Jacobian
        ArenaVec<double> side_points_vec(n_points_per_side, this->patch_arena());
        for (uint i=0; i<n_points_per_side; ++i) {
           side_points_vec(i) = 1.0;
        }
        ArenaOVec<double> side_points_ovec(side_points_vec);
        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> expand_inv_jac_div_det_ovec = inv_jac_div_det_ovec * side_points_ovec;
        Eigen::Matrix<ArenaVec<double>, dim, spacedim> expand_inv_jac_div_det_vec;
        for (uint c=0; c<spacedim*dim; ++c) {
            expand_inv_jac_div_det_vec(c) = expand_inv_jac_div_det_ovec(c).get_vec();
        }
        Eigen::Matrix<ArenaOVec<double>, Eigen::Dynamic, Eigen::Dynamic> expand_jac_ovec = jac_ovec * side_points_ovec;
        Eigen::Matrix<ArenaVec<double>, spacedim, dim> expand_jac_vec;
        for (uint c=0; c<spacedim*dim; ++c) {
            expand_jac_vec(c) = expand_jac_ovec(c).get_vec();
        }

        Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> expand_ref_vec(dim, dim);
        for (uint c=0; c<dim*dim; c++)
            expand_ref_vec(c) = ArenaVec<double>(n_patch_points, this->patch_arena());

        for (uint i_dof=0; i_dof<n_dofs; ++i_dof) {
            // prepare copy of reference data by indices of sides on elements
            FuncHelper<spacedim>::copy_ref_side_on_sides_tensor_data(*this, i_dof, dim, dim, ref_grads_vec, expand_ref_vec);

            // computes operation result
            Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > res_submat = dispatch_op_.result_sub_matrix(i_dof);
            res_submat = expand_inv_jac_div_det_vec.transpose() * expand_ref_vec * expand_jac_vec.transpose();
        }
    }

private:
    PatchOp<spacedim> &dispatch_op_;
};

// class OpGradVectorCovariantShape

/// Dispatch class of vector values
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class DispatchGradVectorShape : public PatchOp<spacedim> {
public:
    /// Constructor
	DispatchGradVectorShape(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim, spacedim}, fe->n_dofs()), in_op_(nullptr)
    {
        this->domain_ = Domain::domain();
        switch (fe->fe_type()) {
            case FEVector:
            {
                in_op_ = new GradVectorShape<dim, Domain, spacedim>(pfev, quad, fe, *this);
                break;
            }
            case FEVectorContravariant:
            {
                ASSERT_PERMANENT(false).error("Shape vector for FEVectorContravariant is not implemented yet!\n"); // temporary assert
                //in_op_ = new GradVectorCovariantShape<dim, Domain, spacedim>(pfev, quad, fe, *this);
                break;
            }
            case FEVectorPiola:
            {
                in_op_ = new GradVectorPiolaShape<dim, Domain, spacedim>(pfev, quad, fe, *this);
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
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class VectorSymGrad : public PatchOp<spacedim> {
public:
    /// Constructor
    VectorSymGrad(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim, spacedim}, fe->n_dofs())
    {
        this->domain_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< DispatchGradVectorShape<dim, Domain, spacedim>, dim >(quad, fe) );
    }

    void eval() override {
        for (uint i_dof=0; i_dof<this->n_dofs(); ++i_dof) {
            Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > grad_vector_dof = this->input_ops(0)->result_sub_matrix(i_dof);
            Eigen::Map< Eigen::Matrix<ArenaVec<double>, Eigen::Dynamic, Eigen::Dynamic> > sym_grad_dof = this->result_sub_matrix(i_dof);
            sym_grad_dof = 0.5 * (grad_vector_dof.transpose() + grad_vector_dof);
        }
    }
};

/// Evaluates vector vector divergence
template<unsigned int dim, class Domain, unsigned int spacedim = 3>
class VectorDivergence : public PatchOp<spacedim> {
public:
    /// Constructor
    VectorDivergence(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {1}, fe->n_dofs())
    {
        this->domain_ = Domain::domain();
        this->input_ops_.push_back( pfev.template get< DispatchGradVectorShape<dim, Domain, spacedim>, dim >(quad, fe) );
    }

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
	OpZero(PatchFEValues<spacedim> &pfev, const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe)
    : PatchOp<spacedim>(dim, pfev, quad, {spacedim, spacedim}, fe->n_dofs())
    {
        this->domain_ = Domain::domain();
        this->allocate_const_result( this->patch_fe_->patch_fe_data().zero_vec_ );
    }

    void eval() override {}
};

} // end of namespace Op


#endif /* OP_FUNCTION_HH_ */
