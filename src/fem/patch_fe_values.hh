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
 * @brief   Class PatchFEValues calculates finite element data on the actual
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

//template<unsigned int dim> class BulkValues;
//template<unsigned int dim> class SideValues;
//template<unsigned int dim> class JoinValues;



template<unsigned int spacedim = 3>
class PatchFEValues {
public:
    typedef typename PatchPointValues<spacedim>::PatchFeData PatchFeData;

    PatchFEValues()
    : patch_fe_data_(1024 * 1024, 256),
	  elem_dim_list_vec_(4),
      patch_point_vals_(2),
	  elements_map_(300, (uint)-1)
    {
        for (uint dim=0; dim<spacedim+1; ++dim) {
            patch_point_vals_[bulk_domain].push_back( PatchPointValues<spacedim>(&elem_dim_list_vec_[dim], bulk_domain) );
            patch_point_vals_[side_domain].push_back( PatchPointValues<spacedim>(&elem_dim_list_vec_[dim], side_domain) );
            element_quads_.push_back( QGauss(dim, 0) );
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
        for (unsigned int i=0; i<spacedim+1; ++i) {
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
     * @brief Returns the number of shape functions og higher dim element.
     */
    template<unsigned int dim>
    unsigned int n_dofs_high() const;

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

    /** Following methods are used during update of patch. **/

    /// Clear elements_map, set values to (-1)
    inline void prepare_new_patch(std::shared_ptr<EvalPoints> eval_points) {
        std::fill(elements_map_.begin(), elements_map_.end(), (uint)-1);
        if (used_domain_[bulk_domain]) patch_point_vals_[bulk_domain][0].resize_tables(eval_points->get_max_bulk_quad_size(0), *patch_fe_data_.patch_arena_);
        for (uint dim=1; dim<=3; ++dim) {
            if (used_domain_[bulk_domain]) patch_point_vals_[bulk_domain][dim].resize_tables(eval_points->get_max_bulk_quad_size(dim), *patch_fe_data_.patch_arena_);
            if (used_domain_[side_domain]) patch_point_vals_[side_domain][dim].resize_tables(eval_points->get_max_side_quad_size(dim), *patch_fe_data_.patch_arena_);
        }
    }

    /// Add elements, sides and quadrature points registered on patch
    template <unsigned int dim>
    inline void add_patch_points(const DimIntegrals<dim> &integrals, ElementCacheMap *element_cache_map) {
        // add bulk points
        for (auto integral_it : integrals.bulk_) {
            auto &patch_data = integral_it.second->patch_data();
            for (uint i_data=0; i_data<patch_data.permanent_size(); ++i_data) {
                uint element_patch_idx = element_cache_map->position_in_cache(patch_data[i_data].cell.elm_idx());
                uint elm_pos = this->register_element(patch_data[i_data].cell.elm(), element_patch_idx);
                uint i_point = 0;
                for (auto p : integral_it.second->points(element_patch_idx) ) {
                    this->register_bulk_point(patch_data[i_data].cell.elm(), elm_pos, p.value_cache_idx(), i_point++);
                }
            }
        }

    	// add boundary points
        for (auto integral_it : integrals.boundary_) {
            auto &patch_data = integral_it.second->patch_data();
            for (unsigned int i_data=0; i_data<patch_data.permanent_size(); ++i_data) {
            	ElementAccessor<spacedim> bdr_elem = patch_data[i_data].side.cond().element_accessor();
            	uint bdr_element_patch_idx = element_cache_map->position_in_cache(bdr_elem.idx(), true);
            	uint elem_pos = this->register_element(bdr_elem, bdr_element_patch_idx);
            	uint side_pos = this->register_side(patch_data[i_data].side, element_cache_map);
                uint i_point = 0;
                for (auto p : integral_it.second->points(patch_data[i_data].side) ) {
                    this->register_side_point(patch_data[i_data].side, side_pos, p.value_cache_idx(), i_point);
                    auto p_bdr = p.point_bdr( bdr_elem );
                    this->register_bulk_point(bdr_elem, elem_pos, p_bdr.value_cache_idx(), i_point++);
                }
            }
        }

    	// add edge points
        for (auto integral_it : integrals.edge_) {
            auto &patch_data = integral_it.second->patch_data();
            for (unsigned int i_data=0; i_data<patch_data.permanent_size(); ++i_data) {
            	auto range = patch_data[i_data].edge_side_range;
                for( DHCellSide edge_side : range )
                {
                	uint side_pos = this->register_side(edge_side, element_cache_map);
                    uint i_point = 0;
                    for (auto p : integral_it.second->points(edge_side) ) {
                        this->register_side_point(edge_side, side_pos, p.value_cache_idx(), i_point++);
                    }
                }
            }
        }

    	// add coupling points
        for (auto integral_it : integrals.coupling_) {
            auto &patch_data = integral_it.second->patch_data();
            uint side_pos, element_patch_idx, elm_pos=0;
            uint last_element_idx = -1;

            for (unsigned int i_data=0; i_data<patch_data.permanent_size(); ++i_data) {
                side_pos = this->register_side(patch_data[i_data].side, element_cache_map);
                if (patch_data[i_data].cell.elm_idx() != last_element_idx) {
                    element_patch_idx = element_cache_map->position_in_cache(patch_data[i_data].cell.elm_idx());
                    elm_pos = this->register_element(patch_data[i_data].cell.elm(), element_patch_idx);
                }

                uint i_bulk_point = 0, i_side_point = 0;
                for (auto p_high : integral_it.second->points(patch_data[i_data].side) )
                {
                    this->register_side_point(patch_data[i_data].side, side_pos, p_high.value_cache_idx(), i_side_point++);
                    if (patch_data[i_data].cell.elm_idx() != last_element_idx) {
                        auto p_low = p_high.lower_dim(patch_data[i_data].cell);
                        this->register_bulk_point(patch_data[i_data].cell.elm(), elm_pos, p_low.value_cache_idx(), i_bulk_point++);
                    }
                }
                last_element_idx = patch_data[i_data].cell.elm_idx();
            }
        }
    }

    /// Register element to patch_point_vals_ table by dimension of element
    uint register_element(ElementAccessor<spacedim> elem, uint element_patch_idx) {
        uint elem_pos = register_element_internal(elem, element_patch_idx);
        PatchPointValues<spacedim> &ppv = patch_point_vals_[bulk_domain][elem.dim()];
        auto map_it = ppv.n_elems_.insert( { (2*elem.idx() + uint(elem.region_idx().is_boundary()) ), ppv.i_mesh_item_} );
        bool is_elm_added = map_it.second;
        if (is_elm_added) {
            ppv.int_table_(patch_elem_on_domain)(ppv.i_mesh_item_) = elem_pos;
            ppv.i_mesh_item_++;
        }
        return map_it.first->second;
    }

    /// Register side to patch_point_vals_ table by dimension of side
    uint register_side(DHCellSide cell_side, ElementCacheMap *element_cache_map) {
        uint dim = cell_side.dim();
        uint element_patch_idx = element_cache_map->position_in_cache(cell_side.elem_idx());
        uint elm_pos = register_element_internal(cell_side.cell().elm(), element_patch_idx);
        PatchPointValues<spacedim> &ppv = patch_point_vals_[side_domain][dim];

        ppv.int_table_(patch_elem_on_domain)(ppv.i_mesh_item_) = elm_pos;
        ppv.int_table_(ref_side_on_sides)(ppv.i_mesh_item_) = cell_side.side_idx();
        ppv.side_list_.push_back( cell_side.side() );
        return ppv.i_mesh_item_++;
    }

    /// Register bulk point to patch_point_vals_ table by dimension of element
    uint register_bulk_point(ElementAccessor<spacedim> elem, uint patch_elm_idx, uint elm_cache_map_idx, uint i_point_on_elem) {
        return patch_point_vals_[bulk_domain][elem.dim()].register_bulk_point(patch_elm_idx, elm_cache_map_idx, elem.idx(), i_point_on_elem);
    }

    /// Register side point to patch_point_vals_ table by dimension of side
    uint register_side_point(DHCellSide cell_side, uint patch_side_idx, uint elm_cache_map_idx, uint i_point_on_side) {
        return patch_point_vals_[side_domain][cell_side.dim()].register_side_point(patch_side_idx, elm_cache_map_idx, cell_side.elem_idx(),
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
        auto tpl = OperationTplHash::op_tuple(op_name, quad->size());
        auto it = op_dependency_.find( tpl );
        if (it == op_dependency_.end()) {
            PatchOp<spacedim>* new_op = new OpType(*this, quad);
            op_dependency_[tpl] = new_op;
            operations_.push_back(new_op);
            DebugOut().fmt("Create new operation '{}', dim: {}, quad size: {}.\n", op_name, dim, quad->size());
            return new_op;
        } else {
            return it->second;
        }
    }

    /// Returns operation of given dim and OpType, creates it if doesn't exist
    template<class OpType, unsigned int dim>
    PatchOp<spacedim>* get_for_elem_quad() {
        return this->template get<OpType, dim>( this->element_quad(dim) );
    }

    /// Returns operation of given dim and OpType, creates it if doesn't exist
    template<class OpType, unsigned int dim>
    PatchOp<spacedim>* get(const Quadrature *quad, std::shared_ptr<FiniteElement<dim>> fe) {
        std::string op_name = typeid(OpType).name();
        auto tpl = OperationTplHash::op_tuple(op_name, quad->size());
        auto it = op_dependency_.find( tpl );
        if (it == op_dependency_.end()) {
            PatchOp<spacedim>* new_op = new OpType(*this, quad, fe);
            op_dependency_[tpl] = new_op;
            operations_.push_back(new_op);
            DebugOut().fmt("Create new operation '{}', dim: {}, quad size: {}.\n", op_name, dim, quad->size());
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
        ASSERT_LE(dim, spacedim);
    	return patch_point_vals_[domain][dim];
    }

    /// Marks data of last successfully added element to patch as permanent
    void make_permanent_ppv_data() {
        for (uint i_dim=0; i_dim<spacedim+1; ++i_dim)
            for (uint i_domain=0; i_domain<2; ++i_domain) {
                patch_point_vals_[i_domain][i_dim].n_mesh_items_.make_permanent();
            }
    }

    /// Return element quadrature (passed to element / side operations)
    const Quadrature* element_quad(unsigned int dim) const {
        ASSERT( dim <= 3 );
        return &element_quads_[dim];
    }

private:
    /// Register element to patch_point_vals_ table by dimension of element
    uint register_element_internal(ElementAccessor<spacedim> elem, uint element_patch_idx) {
        PatchPointValues<spacedim> &ppv = patch_point_vals_[bulk_domain][elem.dim()];
        if (elements_map_[element_patch_idx] != (uint)-1) {
    	    // Return index of element on patch if it is registered repeatedly
    	    return elements_map_[element_patch_idx];
    	}

        elements_map_[element_patch_idx] = ppv.elem_dim_list_->size();
        ppv.elem_dim_list_->push_back( elem );
        return ppv.elem_dim_list_->size() - 1;
    }

    PatchFeData patch_fe_data_;

    /// Sub objects of element data of dimensions 1,2,3
    std::vector< ElemDimList<spacedim> > elem_dim_list_vec_;

    /// Sub objects of bulk and side data of dimensions 1,2,3
    std::vector< std::vector<PatchPointValues<spacedim>> > patch_point_vals_;

    MixedPtr<FiniteElement> fe_;   ///< Mixed of shared pointers of FiniteElement object
    bool used_domain_[2];          ///< Pair of flags signs holds info if bulk and side quadratures are used

    std::vector< PatchOp<spacedim> *> operations_;
    OperationMap< PatchOp<spacedim> > op_dependency_;

    /**
     * Map of element patch indices to PatchOp::result_ and int_table_ tables
     *
     * TODO will be deleted after sorting elements in ElementCacheMap by dimension
     */
    std::vector<uint> elements_map_;

    /**
     * Array of element Quadratures of dim 0,1,2,3
     *
     * Items are used during construction of element/side operations. This solution solves duplicities
     * of these operations (with different quadrature sizes). Quadrature size has no effect on result
     * of these operations.
     */
    std::vector<Quadrature> element_quads_;

    friend class PatchOp<spacedim>;
};



#endif /* PATCH_FE_VALUES_HH_ */
