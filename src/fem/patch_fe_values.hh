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
#include "fem/op_function.hh"
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
      patch_point_vals_bulk_{ {PatchPointValues(1, 0, true, patch_fe_data_),
                               PatchPointValues(2, 0, true, patch_fe_data_),
                               PatchPointValues(3, 0, true, patch_fe_data_)} },
      patch_point_vals_side_{ {PatchPointValues(1, 0, false, patch_fe_data_),
                               PatchPointValues(2, 0, false, patch_fe_data_),
                               PatchPointValues(3, 0, false, patch_fe_data_)} },
      operations_(3)
    {
        used_quads_[0] = false; used_quads_[1] = false;
    }

    PatchFEValues(unsigned int quad_order, MixedPtr<FiniteElement> fe)
    : patch_fe_data_(1024 * 1024, 256),
      patch_point_vals_bulk_{ {PatchPointValues(1, quad_order, true, patch_fe_data_),
    	                       PatchPointValues(2, quad_order, true, patch_fe_data_),
                               PatchPointValues(3, quad_order, true, patch_fe_data_)} },
      patch_point_vals_side_{ {PatchPointValues(1, quad_order, false, patch_fe_data_),
                               PatchPointValues(2, quad_order, false, patch_fe_data_),
                               PatchPointValues(3, quad_order, false, patch_fe_data_)} },
      fe_(fe),
	  operations_(3)
    {
        used_quads_[0] = false; used_quads_[1] = false;

        // TODO move initialization zero_vec_ to patch_fe_data_ constructor when we will create separate ArenaVec of DOshape functions
        uint zero_vec_size = 300;
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
//            if (used_quads_[0]) patch_point_vals_bulk_[i].reinit_patch();
//            if (used_quads_[1]) patch_point_vals_side_[i].reinit_patch();
            for (auto & it : operations_[i]) {
                it.second->eval();
            }
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
        return patch_point_vals_bulk_[dim-1].get_quadrature();
    }

    /// Getter for side quadrature of given dimension
    Quadrature *get_side_quadrature(uint dim) const {
        ASSERT((dim>0) && (dim<=3))(dim).error("Dimension must be 1, 2 or 3.");
        return patch_point_vals_side_[dim-1].get_quadrature();
    }

    /**
     * @brief Returns FiniteElement object of given dimension and component index.
     */
    template<unsigned int dim>
    std::shared_ptr<FiniteElement<dim>> fe_comp(uint component_idx) {
        ASSERT((dim>=0) && (dim<=3))(dim).error("Dimension must be 0, 1, 2 or 3.");
        std::shared_ptr<FiniteElement<dim>> fe = fe_[Dim<dim>{}];
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
        for (uint i=0; i<3; ++i) {
            for (auto & it : operations_[i]) {
                PatchOp<spacedim> *op = it.second;
                if (op->size_type() == elemOp) {
                    op->allocate_result( table_sizes.elem_sizes_[op->bulk_side_][i], *patch_fe_data_.patch_arena_ );
                } else if (op->size_type() == pointOp) {
                    op->allocate_result( table_sizes.point_sizes_[op->bulk_side_][i], *patch_fe_data_.patch_arena_ );
                }
            }

            if (used_quads_[0]) patch_point_vals_bulk_[i].resize_tables(table_sizes.elem_sizes_[0][i], table_sizes.point_sizes_[0][i]);
            if (used_quads_[1]) patch_point_vals_side_[i].resize_tables(table_sizes.elem_sizes_[1][i], table_sizes.point_sizes_[1][i]);
        }
    }

    /// Register element to patch_point_vals_ table by dimension of element
    uint register_element(DHCellAccessor cell, uint element_patch_idx) {
        ElementAccessor<spacedim> elm = cell.elm();
    	arma::mat coords(spacedim, cell.dim()+1);
        for (unsigned int i=0; i<cell.dim()+1; i++)
            coords.col(i) = *elm.node(i);

        PatchPointValues<spacedim> &ppv = patch_point_vals_bulk_[cell.dim()-1];
    	if (ppv.elements_map_[element_patch_idx] != (uint)-1) {
    	    // Return index of element on patch if it is registered repeatedly
    	    return ppv.elements_map_[element_patch_idx];
    	}

        auto coords_mat = ppv.op_el_coords_->result_matrix();
        std::size_t i_elem = ppv.i_elem_;
        for (uint i_col=0; i_col<coords.n_cols; ++i_col)
            for (uint i_row=0; i_row<coords.n_rows; ++i_row) {
                coords_mat(i_row, i_col)(i_elem) = coords(i_row, i_col);
            }

        ppv.elements_map_[element_patch_idx] = ppv.i_elem_;
        return ppv.i_elem_++;
    }

    /// Register side to patch_point_vals_ table by dimension of side
    uint register_side(DHCellSide cell_side) {
        ElementAccessor<spacedim> elm = cell_side.cell().elm();
        arma::mat elm_coords(spacedim, cell_side.dim()+1);
        for (unsigned int i=0; i<cell_side.dim()+1; i++)
            elm_coords.col(i) = *elm.node(i);

        arma::mat side_coords(spacedim, cell_side.dim());
        for (unsigned int n=0; n<cell_side.dim(); n++)
            for (unsigned int c=0; c<spacedim; c++)
                side_coords(c,n) = (*cell_side.side().node(n))[c];

        return patch_point_vals_side_[cell_side.dim()-1].register_side(elm_coords, side_coords, cell_side.side_idx());
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
    template<class OpType>
    PatchOp<spacedim>* get(uint dim) {
        std::string op_name = typeid(OpType).name();
        auto it = operations_[dim-1].find(op_name);
        if (it == operations_[dim-1].end()) {
            PatchOp<spacedim>* new_op = new OpType(dim, *this);
            operations_[dim-1].insert(std::make_pair(op_name, new_op));
            DebugOut().fmt("Create new operation '{}', dim: {}.\n", op_name, dim);
            return new_op;
        } else {
            return it->second;
        }
    }

    /// Returns operation of given dim and OpType, creates it if doesn't exist
    template<class OpType>
    PatchOp<spacedim>* get(uint dim, uint component_idx) {
        std::string op_name = typeid(OpType).name();
        auto it = operations_[dim-1].find(op_name);
        if (it == operations_[dim-1].end()) {
            PatchOp<spacedim>* new_op = new OpType(dim, *this, component_idx);
            operations_[dim-1].insert(std::make_pair(op_name, new_op));
            DebugOut().fmt("Create new operation '{}', dim: {}.\n", op_name, dim);
            return new_op;
        } else {
            return it->second;
        }
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
    std::array<PatchPointValues<spacedim>, 3> patch_point_vals_bulk_;  ///< Sub objects of bulk data of dimensions 1,2,3
    std::array<PatchPointValues<spacedim>, 3> patch_point_vals_side_;  ///< Sub objects of side data of dimensions 1,2,3

    MixedPtr<FiniteElement> fe_;   ///< Mixed of shared pointers of FiniteElement object
    bool used_quads_[2];           ///< Pair of flags signs holds info if bulk and side quadratures are used

    std::vector< std::map<std::string, PatchOp<spacedim> *> > operations_;

    template <class ValueType>
    friend class ElQ;
    template <class ValueType>
    friend class FeQ;
    friend class PatchOp<spacedim>;
    friend class Op::Bulk::El::OpCoords;
};



#endif /* PATCH_FE_VALUES_HH_ */
