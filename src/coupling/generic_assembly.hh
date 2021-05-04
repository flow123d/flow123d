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
 * @file    generic_assembly.hh
 * @brief
 */

#ifndef GENERIC_ASSEMBLY_HH_
#define GENERIC_ASSEMBLY_HH_

#include "quadrature/quadrature_lib.hh"
#include "fields/eval_subset.hh"
#include "fields/eval_points.hh"
#include "fields/field_value_cache.hh"
#include "coupling/balance.hh"
#include "tools/revertable_list.hh"



/// Allow set mask of active integrals.
enum ActiveIntegrals {
    none     =      0,
    bulk     = 0x0001,
    edge     = 0x0002,
    coupling = 0x0004,
    boundary = 0x0008
};


/// Set of all used integral necessary in assemblation
struct AssemblyIntegrals {
    std::array<std::shared_ptr<BulkIntegral>, 3> bulk_;          ///< Bulk integrals of elements of dimensions 1, 2, 3
    std::array<std::shared_ptr<EdgeIntegral>, 3> edge_;          ///< Edge integrals between elements of dimensions 1, 2, 3
    std::array<std::shared_ptr<CouplingIntegral>, 2> coupling_;  ///< Coupling integrals between elements of dimensions 1-2, 2-3
    std::array<std::shared_ptr<BoundaryIntegral>, 3> boundary_;  ///< Boundary integrals betwwen elements of dimensions 1, 2, 3 and boundaries
};


/**
 * @brief Generic class of assemblation.
 *
 * Class
 *  - holds assemblation structures (EvalPoints, Integral objects, Integral data table).
 *  - associates assemblation objects specified by dimension
 *  - provides general assemble method
 *  - provides methods that allow construction of element patches
 */
template < template<IntDim...> class DimAssembly>
class GenericAssembly
{
private:
	/**
	 * Helper structzre holds data of cell (bulk) integral
	 *
	 * Data is specified by cell and subset index in EvalPoint object
	 */
    struct BulkIntegralData {
    	/// Default constructor
        BulkIntegralData() {}

        /// Constructor with data mebers initialization
        BulkIntegralData(DHCellAccessor dhcell, unsigned int subset_idx)
        : cell(dhcell), subset_index(subset_idx) {}

        /// Copy constructor
        BulkIntegralData(const BulkIntegralData &other)
        : cell(other.cell), subset_index(other.subset_index) {}

        DHCellAccessor cell;          ///< Specified cell (element)
        unsigned int subset_index;    ///< Index (order) of subset in EvalPoints object
    };

	/**
	 * Helper structzre holds data of edge integral
	 *
	 * Data is specified by side and subset index in EvalPoint object
	 */
    struct EdgeIntegralData {
    	/// Default constructor
    	EdgeIntegralData()
    	: edge_side_range(make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide() ), make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide() )) {}

        /// Constructor with data mebers initialization
    	EdgeIntegralData(const EdgeIntegralData &other)
        : edge_side_range(other.edge_side_range), subset_index(other.subset_index) {}

        /// Copy constructor
    	EdgeIntegralData(RangeConvert<DHEdgeSide, DHCellSide> range, unsigned int subset_idx)
        : edge_side_range(range), subset_index(subset_idx) {}

    	RangeConvert<DHEdgeSide, DHCellSide> edge_side_range;   ///< Specified cell side (element)
        unsigned int subset_index;                              ///< Index (order) of subset in EvalPoints object
	};

	/**
	 * Helper structzre holds data of neighbour (coupling) integral
	 *
	 * Data is specified by cell, side and their subset indices in EvalPoint object
	 */
    struct CouplingIntegralData {
    	/// Default constructor
       	CouplingIntegralData() {}

        /// Constructor with data mebers initialization
       	CouplingIntegralData(DHCellAccessor dhcell, unsigned int bulk_idx, DHCellSide dhside, unsigned int side_idx)
        : cell(dhcell), bulk_subset_index(bulk_idx), side(dhside), side_subset_index(side_idx) {}

        /// Copy constructor
       	CouplingIntegralData(const CouplingIntegralData &other)
        : cell(other.cell), bulk_subset_index(other.bulk_subset_index), side(other.side), side_subset_index(other.side_subset_index) {}

        DHCellAccessor cell;
	    unsigned int bulk_subset_index;    ///< Index (order) of lower dim subset in EvalPoints object
        DHCellSide side;                   ///< Specified cell side (higher dim element)
	    unsigned int side_subset_index;    ///< Index (order) of higher dim subset in EvalPoints object
    };

	/**
	 * Helper structzre holds data of boundary integral
	 *
	 * Data is specified by side and subset indices of side and appropriate boundary element in EvalPoint object
	 */
    struct BoundaryIntegralData {
    	/// Default constructor
    	BoundaryIntegralData() {}

        /// Constructor with data mebers initialization
    	BoundaryIntegralData(unsigned int bdr_idx, DHCellSide dhside, unsigned int side_idx)
        : bdr_subset_index(bdr_idx), side(dhside), side_subset_index(side_idx) {}

        /// Copy constructor
    	BoundaryIntegralData(const BoundaryIntegralData &other)
        : bdr_subset_index(other.bdr_subset_index), side(other.side), side_subset_index(other.side_subset_index) {}

    	// We don't need hold ElementAccessor of boundary element, side.cond().element_accessor() provides it.
	    unsigned int bdr_subset_index;     ///< Index (order) of subset on boundary element in EvalPoints object
	    DHCellSide side;                   ///< Specified cell side (bulk element)
	    unsigned int side_subset_index;    ///< Index (order) of subset on side of bulk element in EvalPoints object
	};

public:

    /// Constructor
    GenericAssembly( typename DimAssembly<1>::EqData *eq_data, std::shared_ptr<Balance> balance, const DOFHandlerMultiDim * dh)
    : multidim_assembly_(eq_data),
	  bulk_integral_data_(20, 10),
	  edge_integral_data_(12, 6),
	  coupling_integral_data_(12, 6),
	  boundary_integral_data_(8, 4)
    {
        eval_points_ = std::make_shared<EvalPoints>();
        // first step - create integrals, then - initialize cache and initialize subobject of dimensions
        multidim_assembly_[1_d]->create_integrals(eval_points_, integrals_);
        multidim_assembly_[2_d]->create_integrals(eval_points_, integrals_);
        multidim_assembly_[3_d]->create_integrals(eval_points_, integrals_);
        element_cache_map_.init(eval_points_, dh);
        multidim_assembly_[1_d]->initialize(balance);
        multidim_assembly_[2_d]->initialize(balance);
        multidim_assembly_[3_d]->initialize(balance);
        active_integrals_ = multidim_assembly_[1_d]->n_active_integrals();
    }

    /// Getter to set of assembly objects
    inline MixedPtr<DimAssembly, 1> multidim_assembly() const {
        return multidim_assembly_;
    }

    /// Geter to EvalPoints object
    inline std::shared_ptr<EvalPoints> eval_points() const {
        return eval_points_;
    }

	/**
	 * @brief General assemble methods.
	 *
	 * Loops through local cells and calls assemble methods of assembly
	 * object of each cells over space dimension.
	 */
    void assemble(std::shared_ptr<DOFHandlerMultiDim> dh) {
        multidim_assembly_[1_d]->reallocate_cache(element_cache_map_);
        multidim_assembly_[1_d]->begin();

        bool add_into_patch = false; // control variable
        for(auto cell_it = dh->local_range().begin(); cell_it != dh->local_range().end(); )
        {

            if (!add_into_patch) {
        	    element_cache_map_.start_elements_update();
        	    add_into_patch = true;
            }

            START_TIMER("add_integrals_to_patch");
            this->add_integrals_of_computing_step(*cell_it);
            END_TIMER("add_integrals_to_patch");

            if (element_cache_map_.eval_point_data_.temporary_size() > CacheMapElementNumber::get()) {
                bulk_integral_data_.revert_temporary();
                edge_integral_data_.revert_temporary();
                coupling_integral_data_.revert_temporary();
                boundary_integral_data_.revert_temporary();
                element_cache_map_.eval_point_data_.revert_temporary();
                this->assemble_integrals();
                add_into_patch = false;
            } else {
                bulk_integral_data_.make_permanent();
                edge_integral_data_.make_permanent();
                coupling_integral_data_.make_permanent();
                boundary_integral_data_.make_permanent();
                element_cache_map_.eval_point_data_.make_permanent();
                if (element_cache_map_.eval_point_data_.temporary_size() == CacheMapElementNumber::get()) {
                    this->assemble_integrals();
                    add_into_patch = false;
                }
                ++cell_it;
            }
        }
        if (add_into_patch) {
            this->assemble_integrals();
        }

        multidim_assembly_[1_d]->end();
    }

    /// Return ElementCacheMap
    inline const ElementCacheMap &cache_map() const {
        return element_cache_map_;
    }

    /// Return BulkPoint range of appropriate dimension
    inline Range< BulkPoint > bulk_points(unsigned int element_patch_idx, unsigned int dim) const {
        ASSERT_DBG( dim > 0 ).error("Invalid cell dimension, must be 1, 2 or 3!\n");
        return integrals_.bulk_[dim-1]->points(element_patch_idx, &(element_cache_map_));
    }

    /// Return EdgePoint range of appropriate dimension
    inline Range< EdgePoint > edge_points(const DHCellSide &cell_side) const {
        ASSERT_DBG( cell_side.dim() > 0 ).error("Invalid cell dimension, must be 1, 2 or 3!\n");
	    return integrals_.edge_[cell_side.dim()-1]->points(cell_side, &(element_cache_map_));
    }

    /// Return CouplingPoint range of appropriate dimension
    inline Range< CouplingPoint > coupling_points(const DHCellSide &cell_side) const {
        ASSERT_DBG( cell_side.dim() > 1 ).error("Invalid cell dimension, must be 2 or 3!\n");
	    return integrals_.coupling_[cell_side.dim()-2]->points(cell_side, &(element_cache_map_));
    }

    /// Return BoundaryPoint range of appropriate dimension
    inline Range< BoundaryPoint > boundary_points(const DHCellSide &cell_side) const {
        ASSERT_DBG( cell_side.dim() > 0 ).error("Invalid cell dimension, must be 1, 2 or 3!\n");
	    return integrals_.boundary_[cell_side.dim()-1]->points(cell_side, &(element_cache_map_));
    }

private:
    /// Assembles the cell integrals for the given dimension.
    template<unsigned int dim>
    inline void assemble_cell_integrals() {
        // There is special solution of cell integrals solving that is faster than solution of other integrals.
    	// We can't use this solution for other integrals because it does not allow access to equivalent points
    	// throught edge, neigbour and boundary objects.
        for (unsigned int i=0; i<element_cache_map_.n_elements(); ++i) {
            unsigned int elm_start = element_cache_map_.element_chunk_begin_new(i);
            if (element_cache_map_.eval_point_data(elm_start).i_eval_point_ != 0) continue;
            multidim_assembly_[Dim<dim>{}]->cell_integral(i, element_cache_map_.eval_point_data(elm_start).dh_loc_idx_);
        }
    }

    /// Assembles the boundary side integrals for the given dimension.
    template<unsigned int dim>
    inline void assemble_boundary_side_integrals() {
        for (unsigned int i=0; i<boundary_integral_data_.permanent_size(); ++i) {
            if (boundary_integral_data_[i].side.dim() != dim) continue;
            multidim_assembly_[Dim<dim>{}]->boundary_side_integral(boundary_integral_data_[i].side);
        }
    }

    /// Assembles the edge integrals for the given dimension.
    template<unsigned int dim>
    inline void assemble_edge_integrals() {
        for (unsigned int i=0; i<edge_integral_data_.permanent_size(); ++i) {
        	auto range = edge_integral_data_[i].edge_side_range;
            if (range.begin()->dim() != dim) continue;
            multidim_assembly_[Dim<dim>{}]->edge_integral(edge_integral_data_[i].edge_side_range);
        }
    }

    /// Assembles the neighbours integrals for the given dimension.
    template<unsigned int dim>
    inline void assemble_neighbour_integrals() {
        for (unsigned int i=0; i<coupling_integral_data_.permanent_size(); ++i) {
            if (coupling_integral_data_[i].side.dim() != dim) continue;
            multidim_assembly_[Dim<dim>{}]->neigbour_integral(coupling_integral_data_[i].cell, coupling_integral_data_[i].side);
        }
    }

    /// Call assemblations when patch is filled
    void assemble_integrals() {
        START_TIMER("create_patch");
        element_cache_map_.create_patch();
        END_TIMER("create_patch");
        START_TIMER("cache_update");
        multidim_assembly_[1_d]->used_fields_.cache_update(element_cache_map_); // TODO replace with sub FieldSet
        END_TIMER("cache_update");
        element_cache_map_.finish_elements_update();

        {
            START_TIMER("assemble_volume_integrals");
            this->assemble_cell_integrals<1>();
            this->assemble_cell_integrals<2>();
            this->assemble_cell_integrals<3>();
            END_TIMER("assemble_volume_integrals");
        }

        {
            START_TIMER("assemble_fluxes_boundary");
            this->assemble_boundary_side_integrals<1>();
            this->assemble_boundary_side_integrals<2>();
            this->assemble_boundary_side_integrals<3>();
            END_TIMER("assemble_fluxes_boundary");
        }

        {
            START_TIMER("assemble_fluxes_elem_elem");
            this->assemble_edge_integrals<1>();
            this->assemble_edge_integrals<2>();
            this->assemble_edge_integrals<3>();
            END_TIMER("assemble_fluxes_elem_elem");
        }

        {
            START_TIMER("assemble_fluxes_elem_side");
            this->assemble_neighbour_integrals<2>();
            this->assemble_neighbour_integrals<3>();
            END_TIMER("assemble_fluxes_elem_side");
        }
        // clean integral data
        bulk_integral_data_.reset();
        edge_integral_data_.reset();
        coupling_integral_data_.reset();
        boundary_integral_data_.reset();
        element_cache_map_.clear_element_eval_points_map();
    }

    /**
     * Add data of integrals to appropriate structure and register elements to ElementCacheMap.
     *
     * Types of used integrals must be set in data member \p active_integrals_.
     */
    void add_integrals_of_computing_step(DHCellAccessor cell) {
        if (active_integrals_ & ActiveIntegrals::bulk)
    	    if (cell.is_own()) { // Not ghost
                this->add_volume_integral(cell);
    	    }

        for( DHCellSide cell_side : cell.side_range() ) {
            if (active_integrals_ & ActiveIntegrals::boundary)
                if (cell.is_own()) // Not ghost
                    if ( (cell_side.side().edge().n_sides() == 1) && (cell_side.side().is_boundary()) ) {
                        this->add_boundary_integral(cell_side);
                        continue;
                    }
            if (active_integrals_ & ActiveIntegrals::edge)
                if ( (cell_side.n_edge_sides() >= 2) && (cell_side.edge_sides().begin()->element().idx() == cell.elm_idx())) {
                    this->add_edge_integral(cell_side);
                }
        }

        if (active_integrals_ & ActiveIntegrals::coupling) {
            bool add_low = true;
        	for( DHCellSide neighb_side : cell.neighb_sides() ) { // cell -> elm lower dim, neighb_side -> elm higher dim
                if (cell.dim() != neighb_side.dim()-1) continue;
                this->add_coupling_integral(cell, neighb_side, add_low);
                add_low = false;
            }
        }
    }

    /// Add data of volume integral to appropriate data structure.
    inline void add_volume_integral(const DHCellAccessor &cell) {
        uint subset_idx = integrals_.bulk_[cell.dim()-1]->get_subset_idx();
        bulk_integral_data_.emplace_back(cell, subset_idx);

        unsigned int reg_idx = cell.elm().region_idx().idx();
        // Different access than in other integrals: We can't use range method CellIntegral::points
        // because it passes element_patch_idx as argument that is not known during patch construction.
        for (uint i=uint( eval_points_->subset_begin(cell.dim(), subset_idx) );
                  i<uint( eval_points_->subset_end(cell.dim(), subset_idx) ); ++i) {
            element_cache_map_.eval_point_data_.emplace_back(reg_idx, cell.elm_idx(), i, cell.local_idx());
        }
    }

    /// Add data of edge integral to appropriate data structure.
    inline void add_edge_integral(const DHCellSide &cell_side) {
        auto range = cell_side.edge_sides();
        edge_integral_data_.emplace_back(range, integrals_.edge_[range.begin()->dim()-1]->get_subset_idx());

        for( DHCellSide edge_side : range ) {
            unsigned int reg_idx = edge_side.element().region_idx().idx();
            for (auto p : integrals_.edge_[range.begin()->dim()-1]->points(edge_side, &element_cache_map_) ) {
                element_cache_map_.eval_point_data_.emplace_back(reg_idx, edge_side.elem_idx(), p.eval_point_idx(), edge_side.cell().local_idx());
            }
        }
    }

    /// Add data of coupling integral to appropriate data structure.
    inline void add_coupling_integral(const DHCellAccessor &cell, const DHCellSide &ngh_side, bool add_low) {
        coupling_integral_data_.emplace_back(cell, integrals_.coupling_[cell.dim()-1]->get_subset_low_idx(), ngh_side,
                integrals_.coupling_[cell.dim()-1]->get_subset_high_idx());

        unsigned int reg_idx_low = cell.elm().region_idx().idx();
        unsigned int reg_idx_high = ngh_side.element().region_idx().idx();
        for (auto p : integrals_.coupling_[cell.dim()-1]->points(ngh_side, &element_cache_map_) ) {
            element_cache_map_.eval_point_data_.emplace_back(reg_idx_high, ngh_side.elem_idx(), p.eval_point_idx(), ngh_side.cell().local_idx());

        	if (add_low) {
                auto p_low = p.lower_dim(cell); // equivalent point on low dim cell
                element_cache_map_.eval_point_data_.emplace_back(reg_idx_low, cell.elm_idx(), p_low.eval_point_idx(), cell.local_idx());
        	}
        }
    }

    /// Add data of boundary integral to appropriate data structure.
    inline void add_boundary_integral(const DHCellSide &bdr_side) {
        boundary_integral_data_.emplace_back(integrals_.boundary_[bdr_side.dim()-1]->get_subset_low_idx(), bdr_side,
                integrals_.boundary_[bdr_side.dim()-1]->get_subset_high_idx());

        unsigned int reg_idx = bdr_side.element().region_idx().idx();
        for (auto p : integrals_.boundary_[bdr_side.dim()-1]->points(bdr_side, &element_cache_map_) ) {
            element_cache_map_.eval_point_data_.emplace_back(reg_idx, bdr_side.elem_idx(), p.eval_point_idx(), bdr_side.cell().local_idx());

        	auto p_bdr = p.point_bdr(bdr_side.cond().element_accessor()); // equivalent point on boundary element
        	unsigned int bdr_reg = bdr_side.cond().element_accessor().region_idx().idx();
        	// invalid local_idx value, DHCellAccessor of boundary element doesn't exist
        	element_cache_map_.eval_point_data_.emplace_back(bdr_reg, bdr_side.cond().bc_ele_idx(), p_bdr.eval_point_idx(), -1);
        }
    }


    /// Assembly object
    MixedPtr<DimAssembly, 1> multidim_assembly_;

    /// Holds mask of active integrals.
    int active_integrals_;

    AssemblyIntegrals integrals_;                                 ///< Holds integral objects.
    std::shared_ptr<EvalPoints> eval_points_;                     ///< EvalPoints object shared by all integrals
    ElementCacheMap element_cache_map_;                           ///< ElementCacheMap according to EvalPoints

    // Following variables hold data of all integrals depending of actual computed element.
    // TODO sizes of arrays should be set dynamically, depend on number of elements in ElementCacheMap,
    RevertableList<BulkIntegralData>       bulk_integral_data_;      ///< Holds data for computing bulk integrals.
    RevertableList<EdgeIntegralData>       edge_integral_data_;      ///< Holds data for computing edge integrals.
    RevertableList<CouplingIntegralData>   coupling_integral_data_;  ///< Holds data for computing couplings integrals.
    RevertableList<BoundaryIntegralData>   boundary_integral_data_;  ///< Holds data for computing boundary integrals.
};


#endif /* GENERIC_ASSEMBLY_HH_ */
