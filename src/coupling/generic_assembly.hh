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
    struct BulkIntegralData {
        BulkIntegralData() {}

        DHCellAccessor cell;
        unsigned int subset_index;
    };

    struct EdgeIntegralData {
    	EdgeIntegralData()
    	: edge_side_range(make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide() ), make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide() )) {}

    	RangeConvert<DHEdgeSide, DHCellSide> edge_side_range;
        unsigned int subset_index;
	};

    struct CouplingIntegralData {
       	CouplingIntegralData() {}

        DHCellAccessor cell;
	    unsigned int bulk_subset_index;
        DHCellSide side;
	    unsigned int side_subset_index;
    };

    struct BoundaryIntegralData {
    	BoundaryIntegralData() {}

    	// We don't need hold ElementAccessor of boundary element, side.cond().element_accessor() provides it.
	    unsigned int bdr_subset_index; // index of subset on boundary element
	    DHCellSide side;
	    unsigned int side_subset_index;
	};

    /**
     * Temporary struct holds data od boundary element.
     *
     * It will be merged to BulkIntegralData after making implementation of DHCellAccessor on boundary elements.
     */
    struct BdrElementIntegralData {
    	BdrElementIntegralData() {}

        ElementAccessor<3> elm;
        unsigned int subset_index;
    };

public:

    /// Constructor
    GenericAssembly( typename DimAssembly<1>::EqDataDG *eq_data, std::shared_ptr<Balance> balance )
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
        element_cache_map_.init(eval_points_);
        multidim_assembly_[1_d]->initialize(balance);
        multidim_assembly_[2_d]->initialize(balance);
        multidim_assembly_[3_d]->initialize(balance);
        active_integrals_ = multidim_assembly_[1_d]->active_integrals_;
    }

    inline MixedPtr<DimAssembly, 1> multidim_assembly() const {
        return multidim_assembly_;
    }

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
                elm_idx_.clear();
                add_into_patch = false;
            } else {
                bulk_integral_data_.make_permanent();
                edge_integral_data_.make_permanent();
                coupling_integral_data_.make_permanent();
                boundary_integral_data_.make_permanent();
                element_cache_map_.eval_point_data_.make_permanent();
                if (element_cache_map_.eval_point_data_.temporary_size() == CacheMapElementNumber::get()) {
                    this->assemble_integrals();
                    elm_idx_.clear();
                    add_into_patch = false;
                }
                ++cell_it;
            }
        }
        if (add_into_patch) {
            this->assemble_integrals();
            elm_idx_.clear();
        }

        multidim_assembly_[1_d]->end();
    }

    /// Return ElementCacheMap
    inline const ElementCacheMap &cache_map() const {
        return element_cache_map_;
    }

    /// Return BulkPoint range of appropriate dimension
    Range< BulkPoint > bulk_points(const DHCellAccessor &cell) const {
        ASSERT_DBG( cell.dim() > 0 ).error("Invalid cell dimension, must be 1, 2 or 3!\n");
        return integrals_.bulk_[cell.dim()-1]->points(cell, &(element_cache_map_));
    }

    /// Return EdgePoint range of appropriate dimension
    Range< EdgePoint > edge_points(const DHCellSide &cell_side) const {
        ASSERT_DBG( cell_side.dim() > 0 ).error("Invalid cell dimension, must be 1, 2 or 3!\n");
	    return integrals_.edge_[cell_side.dim()-1]->points(cell_side, &(element_cache_map_));
    }

    /// Return CouplingPoint range of appropriate dimension
    Range< CouplingPoint > coupling_points(const DHCellSide &cell_side) const {
        ASSERT_DBG( cell_side.dim() > 1 ).error("Invalid cell dimension, must be 2 or 3!\n");
	    return integrals_.coupling_[cell_side.dim()-2]->points(cell_side, &(element_cache_map_));
    }

    /// Return BoundaryPoint range of appropriate dimension
    Range< BoundaryPoint > boundary_points(const DHCellSide &cell_side) const {
        ASSERT_DBG( cell_side.dim() > 0 ).error("Invalid cell dimension, must be 1, 2 or 3!\n");
	    return integrals_.boundary_[cell_side.dim()-1]->points(cell_side, &(element_cache_map_));
    }

private:
    /// Assembles the cell integrals for the given dimension.
    template<unsigned int dim>
    inline void assemble_cell_integrals() {
        for (unsigned int i=0; i<bulk_integral_data_.permanent_size(); ++i) {
            if (bulk_integral_data_[i].cell.dim() != dim) continue;
            multidim_assembly_[Dim<dim>{}]->cell_integral(bulk_integral_data_[i].cell);
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
        multidim_assembly_[1_d]->data_->cache_update(element_cache_map_);
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
    void add_volume_integral(const DHCellAccessor &cell) {
        BulkIntegralData data;
        data.cell = cell;
        data.subset_index = integrals_.bulk_[cell.dim()-1]->get_subset_idx();
        bulk_integral_data_.push_back(data);

        unsigned int reg_idx = cell.elm().region_idx().idx();
        for (auto p : integrals_.bulk_[cell.dim()-1]->points(cell, &element_cache_map_) ) {
            EvalPointData epd(reg_idx, cell.elm_idx(), p.eval_point_idx());
            element_cache_map_.eval_point_data_.push_back(epd);
        }
        elm_idx_.insert(cell.elm_idx());
    }

    /// Add data of edge integral to appropriate data structure.
    void add_edge_integral(const DHCellSide &cell_side) {
        EdgeIntegralData data;
        data.edge_side_range = cell_side.edge_sides();
        data.subset_index = integrals_.edge_[data.edge_side_range.begin()->dim()-1]->get_subset_idx();
        edge_integral_data_.push_back(data);

        for( DHCellSide edge_side : data.edge_side_range ) {
            unsigned int reg_idx = edge_side.element().region_idx().idx();
            for (auto p : integrals_.edge_[data.edge_side_range.begin()->dim()-1]->points(edge_side, &element_cache_map_) ) {
                EvalPointData epd(reg_idx, edge_side.elem_idx(), p.eval_point_idx());
                element_cache_map_.eval_point_data_.push_back(epd);
            }
            elm_idx_.insert(edge_side.elem_idx());
        }
    }

    /// Add data of coupling integral to appropriate data structure.
    void add_coupling_integral(const DHCellAccessor &cell, const DHCellSide &ngh_side, bool add_low) {
        CouplingIntegralData data;
        data.cell = cell;
        data.side = ngh_side;
        data.bulk_subset_index = integrals_.coupling_[cell.dim()-1]->get_subset_low_idx();
        data.side_subset_index = integrals_.coupling_[cell.dim()-1]->get_subset_high_idx();
        coupling_integral_data_.push_back(data);

        unsigned int reg_idx_low = cell.elm().region_idx().idx();
        unsigned int reg_idx_high = ngh_side.element().region_idx().idx();
        for (auto p : integrals_.coupling_[cell.dim()-1]->points(ngh_side, &element_cache_map_) ) {
            EvalPointData epd(reg_idx_high, ngh_side.elem_idx(), p.eval_point_idx());
            element_cache_map_.eval_point_data_.push_back(epd);

        	if (add_low) {
                auto p_low = p.lower_dim(cell); // equivalent point on low dim cell
               	EvalPointData epd_low(reg_idx_low, cell.elm_idx(), p_low.eval_point_idx());
                element_cache_map_.eval_point_data_.push_back(epd_low);
        	}
        }
        elm_idx_.insert(cell.elm_idx());
        elm_idx_.insert(ngh_side.elem_idx());
    }

    /// Add data of boundary integral to appropriate data structure.
    void add_boundary_integral(const DHCellSide &bdr_side) {
        BoundaryIntegralData data;
        data.side = bdr_side;
        data.bdr_subset_index = integrals_.boundary_[bdr_side.dim()-1]->get_subset_low_idx();
        data.side_subset_index = integrals_.boundary_[bdr_side.dim()-1]->get_subset_high_idx();
        boundary_integral_data_.push_back(data);

        unsigned int reg_idx = bdr_side.element().region_idx().idx();
        for (auto p : integrals_.boundary_[bdr_side.dim()-1]->points(bdr_side, &element_cache_map_) ) {
            EvalPointData epd(reg_idx, bdr_side.elem_idx(), p.eval_point_idx());
            element_cache_map_.eval_point_data_.push_back(epd);

        	auto p_bdr = p.point_bdr(bdr_side.cond().element_accessor()); // equivalent point on boundary element
        	unsigned int bdr_reg = bdr_side.cond().element_accessor().region_idx().idx();
        	EvalPointData epd_bdr(bdr_reg, bdr_side.cond().bc_ele_idx(), p_bdr.eval_point_idx());
        	element_cache_map_.eval_point_data_.push_back(epd_bdr);
        }
        elm_idx_.insert(bdr_side.elem_idx());
        auto bdr_elm_acc = bdr_side.cond().element_accessor();
        elm_idx_.insert(bdr_elm_acc.mesh_idx());
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

    /// Set of element idx used in patch, temporary data member, TODO it will be replaced
    std::set<unsigned int> elm_idx_;
};


#endif /* GENERIC_ASSEMBLY_HH_ */
