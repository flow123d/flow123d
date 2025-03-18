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
      patch_point_vals_(2)
    {
        for (uint dim=1; dim<4; ++dim) {
            patch_point_vals_[0].push_back( PatchPointValues(dim, 0, true, patch_fe_data_) );
            patch_point_vals_[1].push_back( PatchPointValues(dim, 0, false, patch_fe_data_) );
        }
        used_quads_[0] = false; used_quads_[1] = false;
    }

    PatchFEValues(unsigned int quad_order, MixedPtr<FiniteElement> fe)
    : patch_fe_data_(1024 * 1024, 256),
      patch_point_vals_(2),
      fe_(fe)
    {
        for (uint dim=1; dim<4; ++dim) {
            patch_point_vals_[0].push_back( PatchPointValues(dim, quad_order, true, patch_fe_data_) );
            patch_point_vals_[1].push_back( PatchPointValues(dim, quad_order, false, patch_fe_data_) );
        }
        used_quads_[0] = false; used_quads_[1] = false;

        // TODO move initialization zero_vec_ to patch_fe_data_ constructor when we will create separate ArenaVec of DOshape functions
        uint zero_vec_size = 300;
        patch_fe_data_.zero_vec_ = ArenaVec<double>(zero_vec_size, patch_fe_data_.asm_arena_);
        for (uint i=0; i<zero_vec_size; ++i) patch_fe_data_.zero_vec_(i) = 0.0;
    }


    /// Destructor
    ~PatchFEValues()
    {}

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
            patch_point_vals_[0][DIM-1].initialize(); // bulk
        } else {
            used_quads_[1] = true;
            patch_point_vals_[1][DIM-1].initialize(); // side
        }
    }

    /// Finalize initialization, creates child (patch) arena and passes it to PatchPointValue objects
    void init_finalize() {
        patch_fe_data_.patch_arena_ = patch_fe_data_.asm_arena_.get_child_arena();
    }

    /// Reset PatchpointValues structures
    void reset()
    {
        for (unsigned int i=0; i<spacedim; ++i) {
            if (used_quads_[0]) patch_point_vals_[0][i].reset();
            if (used_quads_[1]) patch_point_vals_[1][i].reset();
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

    /// Getter for bulk quadrature of given dimension
    Quadrature *get_bulk_quadrature(uint dim) const {
        ASSERT((dim>0) && (dim<=3))(dim).error("Dimension must be 1, 2 or 3.");
        return patch_point_vals_[0][dim-1].get_quadrature();
    }

    /// Getter for side quadrature of given dimension
    Quadrature *get_side_quadrature(uint dim) const {
        ASSERT((dim>0) && (dim<=3))(dim).error("Dimension must be 1, 2 or 3.");
        return patch_point_vals_[1][dim-1].get_quadrature();
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


    /// Return BulkValue object of dimension given by template parameter
    template<unsigned int dim>
    BulkValues<dim> bulk_values();

    /// Return SideValue object of dimension given by template parameter
    template<unsigned int dim>
    SideValues<dim> side_values();

    /// Return JoinValue object of dimension given by template parameter
    template<unsigned int dim>
    JoinValues<dim> join_values();

    /** Following methods are used during update of patch. **/

    /// Resize tables of patch_point_vals_
    void resize_tables(TableSizes table_sizes) {
        for (uint i=0; i<spacedim; ++i) {
            if (used_quads_[0]) patch_point_vals_[0][i].resize_tables(table_sizes.elem_sizes_[0][i], table_sizes.point_sizes_[0][i]);
            if (used_quads_[1]) patch_point_vals_[1][i].resize_tables(table_sizes.elem_sizes_[1][i], table_sizes.point_sizes_[1][i]);
        }
    }

    /// Register element to patch_point_vals_ table by dimension of element
    uint register_element(DHCellAccessor cell, uint element_patch_idx) {
        PatchPointValues<spacedim> &ppv = patch_point_vals_[0][cell.dim()-1];
    	if (ppv.elements_map_[element_patch_idx] != (uint)-1) {
    	    // Return index of element on patch if it is registered repeatedly
    	    return ppv.elements_map_[element_patch_idx];
    	}

        ppv.elements_map_[element_patch_idx] = ppv.i_elem_;
        ppv.elem_list_.push_back( cell.elm() );
        return ppv.i_elem_++;
    }

    /// Register side to patch_point_vals_ table by dimension of side
    uint register_side(DHCellSide cell_side) {
        uint dim = cell_side.dim();
        PatchPointValues<spacedim> &ppv = patch_point_vals_[1][dim-1];

        ppv.int_table_(3)(ppv.i_elem_) = cell_side.side_idx();
        ppv.elem_list_.push_back( cell_side.cell().elm() );
        ppv.side_list_.push_back( cell_side.side() );

        return ppv.i_elem_++;
    }

    /// Register bulk point to patch_point_vals_ table by dimension of element
    uint register_bulk_point(DHCellAccessor cell, uint elem_table_row, uint value_patch_idx, uint i_point_on_elem) {
        return patch_point_vals_[0][cell.dim()-1].register_bulk_point(elem_table_row, value_patch_idx, cell.elm_idx(), i_point_on_elem);
    }

    /// Register side point to patch_point_vals_ table by dimension of side
    uint register_side_point(DHCellSide cell_side, uint elem_table_row, uint value_patch_idx, uint i_point_on_side) {
        return patch_point_vals_[1][cell_side.dim()-1].register_side_point(elem_table_row, value_patch_idx, cell_side.elem_idx(),
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
    PatchOp<spacedim>* get() {
        std::string op_name = typeid(OpType).name();
        auto it = op_dependency_.find(op_name);
        if (it == op_dependency_.end()) {
            PatchOp<spacedim>* new_op = new OpType(*this);
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
    PatchOp<spacedim>* get(std::shared_ptr<FiniteElement<dim>> fe) {
        std::string op_name = typeid(OpType).name();
        auto it = op_dependency_.find(op_name);
        if (it == op_dependency_.end()) {
            PatchOp<spacedim>* new_op = new OpType(*this, fe);
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

private:
    PatchFeData patch_fe_data_;
    /// Sub objects of bulk and side data of dimensions 1,2,3
    std::vector< std::vector<PatchPointValues<spacedim>> > patch_point_vals_;

    MixedPtr<FiniteElement> fe_;   ///< Mixed of shared pointers of FiniteElement object
    bool used_quads_[2];           ///< Pair of flags signs holds info if bulk and side quadratures are used

    std::vector< PatchOp<spacedim> *> operations_;
    std::unordered_map<std::string, PatchOp<spacedim> *> op_dependency_;

    friend class PatchOp<spacedim>;
};



#endif /* PATCH_FE_VALUES_HH_ */
