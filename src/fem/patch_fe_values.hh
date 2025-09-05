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
#include "fem/patch_op.hh"
#include "fem/op_accessors.hh"
#include "mesh/ref_element.hh"                // for RefElement
#include "mesh/accessors.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/arena_resource.hh"
#include "fem/arena_vec.hh"

template<unsigned int dim> class BulkValues;
template<unsigned int dim> class SideValues;
template<unsigned int dim> class JoinValues;



template<unsigned int spacedim = 3>
class PatchFEValues {
public:
    typedef typename PatchPointValues<spacedim>::PatchFeData PatchFeData;

    PatchFEValues()
    : patch_fe_data_(1024 * 1024, 256),
	  elems_dim_data_vec_(3),
      patch_point_vals_(2),
	  elements_map_(300, (uint)-1)
    {
        for (uint dim=1; dim<4; ++dim) {
            patch_point_vals_[0].push_back( PatchPointValues<spacedim>(&elems_dim_data_vec_[dim-1], bulk_domain) );
            patch_point_vals_[1].push_back( PatchPointValues<spacedim>(&elems_dim_data_vec_[dim-1], side_domain) );
        }
        used_domain_[bulk_domain] = false; used_domain_[side_domain] = false;
    }

    PatchFEValues(MixedPtr<FiniteElement> fe)
    : PatchFEValues<spacedim>()
    {
        fe_ = fe;

        // TODO move initialization zero_vec_ to patch_fe_data_ constructor when we will create separate ArenaVec of DOshape functions
        uint zero_vec_size = 300;
        patch_fe_data_.zero_vec_ = ArenaVec<double>(zero_vec_size, patch_fe_data_.asm_arena_);
        for (uint i=0; i<zero_vec_size; ++i) patch_fe_data_.zero_vec_(i) = 0.0;
    }


    /// Destructor
    ~PatchFEValues()
    {}

    /// Finalize initialization, creates child (patch) arena and passes it to PatchPointValue objects
    void init_finalize() {
        patch_fe_data_.patch_arena_ = patch_fe_data_.asm_arena_.get_child_arena();
    }

    /// Reset PatchpointValues structures
    void reset()
    {
        for (unsigned int i=0; i<spacedim; ++i) {
            if (used_domain_[bulk_domain]) patch_point_vals_[bulk_domain][i].reset();
            if (used_domain_[side_domain]) patch_point_vals_[side_domain][i].reset();
        }
        patch_fe_data_.patch_arena_->reset();
    }

    /// Reinit data.
    void reinit_patch()
    {
        for (auto * op : operations_) {
            op->eval();
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

    /**
     * @brief Returnd FiniteElement of \p component_idx for FESystem or \p fe for other types
     */
    template<unsigned int dim>
    std::shared_ptr<FiniteElement<dim>> fe_comp(std::shared_ptr< FiniteElement<dim> > fe, uint component_idx) {
        if (fe->fe_type() == FEMixedSystem) {
            FESystem<dim> *fe_sys = dynamic_cast<FESystem<dim>*>( fe.get() );
            return fe_sys->fe()[component_idx];
        } else {
            ASSERT_EQ(component_idx, 0).warning("Non-zero component_idx can only be used for FESystem.");
            return fe;
        }
    }

    /// Returns pointer to FiniteElement of given dimension.
    template<unsigned int dim>
    std::shared_ptr<FiniteElement<dim>> fe_dim() {
        return fe_[Dim<dim>{}];
    }

//    /// Return BulkValue object of dimension given by template parameter
//    template<unsigned int dim>
//    BulkValues<dim> bulk_values();
//
//    /// Return SideValue object of dimension given by template parameter
//    template<unsigned int dim>
//    SideValues<dim> side_values();
//
//    /// Return JoinValue object of dimension given by template parameter
//    template<unsigned int dim>
//    JoinValues<dim> join_values();

    /** Following methods are used during update of patch. **/

    /// Resize tables of patch_point_vals_
    void resize_tables() {
        for (uint i=0; i<spacedim; ++i) {
            if (used_domain_[bulk_domain]) patch_point_vals_[bulk_domain][i].resize_tables(*patch_fe_data_.patch_arena_);
            if (used_domain_[side_domain]) patch_point_vals_[side_domain][i].resize_tables(*patch_fe_data_.patch_arena_);
        }
        std::fill(elements_map_.begin(), elements_map_.end(), (uint)-1);
    }

    /// Register element to patch_point_vals_ table by dimension of element
    uint register_element(DHCellAccessor cell, uint element_patch_idx) {
        PatchPointValues<spacedim> &ppv = patch_point_vals_[bulk_domain][cell.dim()-1];
        if (elements_map_[element_patch_idx] != (uint)-1) {
    	    // Return index of element on patch if it is registered repeatedly
    	    return elements_map_[element_patch_idx];
    	}

        elements_map_[element_patch_idx] = ppv.elems_dim_data_->i_elem_;
        ppv.elems_dim_data_->elem_list_.push_back( cell.elm() );
        return ppv.elems_dim_data_->i_elem_++;
    }

    /// Register side to patch_point_vals_ table by dimension of side
    uint register_side(DHCellSide cell_side, uint element_patch_idx) {
        uint dim = cell_side.dim();
        uint elm_pos = register_element(cell_side.cell(), element_patch_idx);
        PatchPointValues<spacedim> &ppv = patch_point_vals_[side_domain][dim-1];

        ppv.int_table_(3)(ppv.i_side_) = cell_side.side_idx();
        ppv.int_table_(5)(ppv.i_side_) = elm_pos;
        ppv.side_list_.push_back( cell_side.side() );
        return ppv.i_side_++;
    }

    /// Register bulk point to patch_point_vals_ table by dimension of element
    uint register_bulk_point(DHCellAccessor cell, uint patch_elm_idx, uint elm_cache_map_idx, uint i_point_on_elem) {
        return patch_point_vals_[bulk_domain][cell.dim()-1].register_bulk_point(patch_elm_idx, elm_cache_map_idx, cell.elm_idx(), i_point_on_elem);
    }

    /// Register side point to patch_point_vals_ table by dimension of side
    uint register_side_point(DHCellSide cell_side, uint patch_side_idx, uint elm_cache_map_idx, uint i_point_on_side) {
        return patch_point_vals_[1][cell_side.dim()-1].register_side_point(patch_side_idx, elm_cache_map_idx, cell_side.elem_idx(),
                cell_side.side_idx(), i_point_on_side);
    }

    /// return reference to assembly arena
    inline AssemblyArena &asm_arena() {
    	return patch_fe_data_.asm_arena_;
    }

    /// same as previous but return constant reference
    inline const AssemblyArena &asm_arena() const {
    	return patch_fe_data_.asm_arena_;
    }

    /// return reference to patch arena
    inline PatchArena &patch_arena() const {
    	return *patch_fe_data_.patch_arena_;
    }

    /// Returns operation of given dim and OpType, creates it if doesn't exist
    template<class OpType, unsigned int dim>
    PatchOp<spacedim>* get(const Quadrature *quad) {
        std::string op_name = typeid(OpType).name();
        auto it = op_dependency_.find(op_name);
        if (it == op_dependency_.end()) {
            PatchOp<spacedim>* new_op = new OpType(*this, quad);
            op_dependency_.insert(std::make_pair(op_name, new_op));
            operations_.push_back(new_op);
            DebugOut().fmt("Create new operation '{}', dim: {}.\n", op_name, dim);
            return new_op;
        } else {
            return it->second;
        }
    }

    /// Returns operation of given dim and OpType, creates it if doesn't exist
    template<class OpType, unsigned int dim>
    PatchOp<spacedim>* get(const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe) {
        std::string op_name = typeid(OpType).name();
        auto it = op_dependency_.find(op_name);
        if (it == op_dependency_.end()) {
            PatchOp<spacedim>* new_op = new OpType(*this, quad, fe);
            op_dependency_.insert(std::make_pair(op_name, new_op));
            operations_.push_back(new_op);
            DebugOut().fmt("Create new operation '{}', dim: {}.\n", op_name, dim);
            return new_op;
        } else {
            return it->second;
        }
    }

    /// Print table of all used operations - development method
    void print_operations(ostream& stream) const {
        stream << endl << "Table of patch FE operations:" << endl;
        stream << std::setfill('-') << setw(160) << "" << endl;

        stream << std::setfill(' ') << " Operation" << std::setw(51) << "" << "Type" << std::setw(5) << "" << "Shape" << std::setw(2) << ""
                << "n DOFs" << std::setw(2) << "" << "Input operations" << std::endl;
        for (uint i=0; i<operations_.size(); ++i) {
            stream << " " << std::left << std::setw(60) << typeid(*operations_[i]).name() << "";
            stream << operations_[i]->dim_ << "D " << (operations_[i]->domain_ ? "side" : "bulk");
        	stream << "  " << std::setw(6) << operations_[i]->format_shape() << "" << " "
                << std::setw(7) << operations_[i]->n_dofs() << "" << " ";
            for (auto *i_o : operations_[i]->input_ops_) stream << typeid(*i_o).name() << "  ";
            stream << std::endl;
        }

        stream << std::setfill('=') << setw(160) << "" << endl;
    }

    /// Getter of patch_fe_data_
    PatchFeData &patch_fe_data() {
        return patch_fe_data_;
    }

    /// Mark domain (bulk or side) as used in assembly class
    inline void set_used_domain(fem_domain domain) {
        used_domain_[domain] = true;
    }

    /// Temporary method
    PatchPointValues<spacedim> &ppv(uint domain, uint dim) {
        ASSERT( domain<2 );
        ASSERT( (dim>0) && (dim<=3) )(dim);
    	return patch_point_vals_[domain][dim-1];
    }

    /// Temporary method
    void make_permanent_ppv_data() {
        for (uint i_dim=0; i_dim<3; ++i_dim)
            for (uint i_domain=0; i_domain<2; ++i_domain) {
                patch_point_vals_[i_domain][i_dim].n_points_.make_permanent();
                patch_point_vals_[i_domain][i_dim].n_mesh_items_.make_permanent();
            }
    }

private:
    PatchFeData patch_fe_data_;

    /// Sub objects of element data of dimensions 1,2,3
    std::vector< ElemsDimOata<spacedim> > elems_dim_data_vec_;

    /// Sub objects of bulk and side data of dimensions 1,2,3
    std::vector< std::vector<PatchPointValues<spacedim>> > patch_point_vals_;

    MixedPtr<FiniteElement> fe_;   ///< Mixed of shared pointers of FiniteElement object
    bool used_domain_[2];          ///< Pair of flags signs holds info if bulk and side quadratures are used

    std::vector< PatchOp<spacedim> *> operations_;
    std::unordered_map<std::string, PatchOp<spacedim> *> op_dependency_;
    //OperationSet< PatchOp<spacedim> > op_dependency_;

    std::vector<uint> elements_map_;    ///< Map of element patch indices to PatchOp::result_ and int_table_ tables

    friend class PatchOp<spacedim>;
};



#endif /* PATCH_FE_VALUES_HH_ */
