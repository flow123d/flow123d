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
#include "tools/tmp_size_list.hh"



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


/*
 * Specifies eval points by idx of region, element and eval point.
 *
 * TODO Add better description after finish implementation
 */
struct EvalPointData {
	EvalPointData() {}              ///< Default constructor
	/// Constructor sets all data members
	EvalPointData(unsigned int i_reg, unsigned int i_ele, unsigned int i_ep)
	: i_reg_(i_reg), i_element_(i_ele), i_eval_point_(i_ep) {}

	unsigned int i_reg_;            ///< region_idx of element
    unsigned int i_element_;        ///< mesh_idx of ElementAccessor appropriate to element
    unsigned int i_eval_point_;     ///< index of point in EvalPoint object
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
    /// Holds maximum number of neighbouring element adding in one step to patch (own element and its side, coupling, boundary neighbours).
    static const unsigned int maximal_neighbouring_elem = 9;

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
    GenericAssembly( typename DimAssembly<1>::EqDataDG *eq_data, std::shared_ptr<Balance> balance, int active_integrals )
    : multidim_assembly_(eq_data),
      active_integrals_(active_integrals),
	  bulk_integral_data_(12),
	  edge_integral_data_(12),
	  coupling_integral_data_(12),
	  boundary_integral_data_(8),
	  eval_point_data_(0)
    {
        eval_points_ = std::make_shared<EvalPoints>();
        // first step - create integrals, then - initialize cache and initialize subobject of dimensions
        multidim_assembly_[1_d]->create_integrals(eval_points_, integrals_, active_integrals_);
        multidim_assembly_[2_d]->create_integrals(eval_points_, integrals_, active_integrals_);
        multidim_assembly_[3_d]->create_integrals(eval_points_, integrals_, active_integrals_);
        element_cache_map_.init(eval_points_);
        multidim_assembly_[1_d]->initialize(balance);
        multidim_assembly_[2_d]->initialize(balance);
        multidim_assembly_[3_d]->initialize(balance);
        unsigned int ep_data_size = ElementCacheMap::n_cached_elements * eval_points_->max_size();
        eval_point_data_.resize(ep_data_size);
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
    void assemble(std::shared_ptr<DOFHandlerMultiDim> dh, const TimeStep &step) {
        multidim_assembly_[1_d]->reallocate_cache(element_cache_map_);
        multidim_assembly_[1_d]->begin();

        bool add_into_patch = false; // control variable
        for (auto cell : dh->local_range() )
        {
            if (!add_into_patch) {
        	    element_cache_map_.start_elements_update();
        	    add_into_patch = true;
            }

            this->add_integrals_of_computing_step(cell);
            element_cache_map_.cache_map_index(cell);
            eval_point_data_.finalize_tmp(); // check if there is space in FieldValueCache (if not call revert_tmp() and repeat iteration with same cell)

            unsigned int n_elements = 0; // TODO remove - temporary solution
            for (auto const& i : element_cache_map_.update_cache_data().region_cache_indices_map_) {
                n_elements += i.second.n_elements_;
            }
            if (n_elements > ElementCacheMap::n_cached_elements-GenericAssembly<DimAssembly>::maximal_neighbouring_elem) {
                this->assemble_integrals(step);
                add_into_patch = false;
            }
        }
        if (add_into_patch)
            this->assemble_integrals(step);

        multidim_assembly_[1_d]->end();
    }

    /// Return ElementCacheMap
    inline const ElementCacheMap &cache_map() const {
        return element_cache_map_;
    }

    /// Return BulkIntegral of appropriate dimension
    inline std::shared_ptr<BulkIntegral> bulk_integral(unsigned int dim) const {
        ASSERT_DBG( (dim>0) && (dim<=3) )(dim).error("Invalid dimension, must be 1, 2 or 3!\n");
	    return integrals_.bulk_[dim-1];
    }

    /// Return EdgeIntegral of appropriate dimension
    inline std::shared_ptr<EdgeIntegral> edge_integral(unsigned int dim) const {
        ASSERT_DBG( (dim>0) && (dim<=3) )(dim).error("Invalid dimension, must be 1, 2 or 3!\n");
	    return integrals_.edge_[dim-1];
    }

    /// Return CouplingIntegral between dimensions dim-1 and dim
    inline std::shared_ptr<CouplingIntegral> coupling_integral(unsigned int dim) const {
        ASSERT_DBG( (dim>1) && (dim<=3) )(dim).error("Invalid dimension, must be 2 or 3!\n");
	    return integrals_.coupling_[dim-2];
    }

    /// Return BoundaryIntegral of appropriate dimension
    inline std::shared_ptr<BoundaryIntegral> boundary_integral(unsigned int dim) const {
        ASSERT_DBG( (dim>0) && (dim<=3) )(dim).error("Invalid dimension, must be 1, 2 or 3!\n");
	    return integrals_.boundary_[dim-1];
    }

private:
    /// Call assemblations when patch is filled
    void assemble_integrals(const TimeStep &step) {
        unsigned int i;
        bulk_integral_data_.finalize_tmp();
        edge_integral_data_.finalize_tmp();
        coupling_integral_data_.finalize_tmp();
        boundary_integral_data_.finalize_tmp();
        element_cache_map_.prepare_elements_to_update();
        this->insert_eval_points_from_integral_data();
        element_cache_map_.create_elements_points_map();
        multidim_assembly_[1_d]->data_->cache_update(element_cache_map_);
        element_cache_map_.finish_elements_update();

        if (active_integrals_ & ActiveIntegrals::bulk) {
            START_TIMER("assemble_volume_integrals");
            for (i=0; i<bulk_integral_data_.final_size(); ++i) { // volume integral
                switch (bulk_integral_data_[i].cell.dim()) {
                case 1:
                    multidim_assembly_[1_d]->assemble_volume_integrals(bulk_integral_data_[i].cell);
                    break;
                case 2:
                    multidim_assembly_[2_d]->assemble_volume_integrals(bulk_integral_data_[i].cell);
                    break;
                case 3:
                    multidim_assembly_[3_d]->assemble_volume_integrals(bulk_integral_data_[i].cell);
                    break;
                }
            }
            END_TIMER("assemble_volume_integrals");
        }

        if (active_integrals_ & ActiveIntegrals::boundary) {
            START_TIMER("assemble_fluxes_boundary");
            for (i=0; i<boundary_integral_data_.final_size(); ++i) { // boundary integral
                switch (boundary_integral_data_[i].side.dim()) {
                case 1:
                    multidim_assembly_[1_d]->assemble_fluxes_boundary(boundary_integral_data_[i].side, step);
                    break;
                case 2:
                    multidim_assembly_[2_d]->assemble_fluxes_boundary(boundary_integral_data_[i].side, step);
                    break;
                case 3:
                    multidim_assembly_[3_d]->assemble_fluxes_boundary(boundary_integral_data_[i].side, step);
                    break;
                }
            }
            END_TIMER("assemble_fluxes_boundary");
        }

        if (active_integrals_ & ActiveIntegrals::edge) {
            START_TIMER("assemble_fluxes_elem_elem");
            for (i=0; i<edge_integral_data_.final_size(); ++i) { // edge integral
            	auto range = edge_integral_data_[i].edge_side_range;
                switch (range.begin()->dim()) {
                case 1:
                    multidim_assembly_[1_d]->assemble_fluxes_element_element(edge_integral_data_[i].edge_side_range);
                    break;
                case 2:
                    multidim_assembly_[2_d]->assemble_fluxes_element_element(edge_integral_data_[i].edge_side_range);
                    break;
                case 3:
                    multidim_assembly_[3_d]->assemble_fluxes_element_element(edge_integral_data_[i].edge_side_range);
                    break;
                }
            }
            END_TIMER("assemble_fluxes_elem_elem");
        }

        if (active_integrals_ & ActiveIntegrals::coupling) {
            START_TIMER("assemble_fluxes_elem_side");
            for (i=0; i<coupling_integral_data_.final_size(); ++i) { // coupling integral
                switch (coupling_integral_data_[i].side.dim()) {
                case 2:
                    multidim_assembly_[2_d]->assemble_fluxes_element_side(coupling_integral_data_[i].cell, coupling_integral_data_[i].side);
                    break;
                case 3:
                    multidim_assembly_[3_d]->assemble_fluxes_element_side(coupling_integral_data_[i].cell, coupling_integral_data_[i].side);
                    break;
                }
            }
            END_TIMER("assemble_fluxes_elem_side");
        }
        // clean integral data
        bulk_integral_data_.reset();
        edge_integral_data_.reset();
        coupling_integral_data_.reset();
        boundary_integral_data_.reset();
        eval_point_data_.reset();
    }

    /// Mark eval points in table of Element cache map.
    void insert_eval_points_from_integral_data() {
        for (unsigned int i=0; i<bulk_integral_data_.final_size(); ++i) {
            // add data to cache if there is free space, else return
        	unsigned int data_size = eval_points_->subset_size( bulk_integral_data_[i].cell.dim(), bulk_integral_data_[i].subset_index );
        	element_cache_map_.mark_used_eval_points(bulk_integral_data_[i].cell, bulk_integral_data_[i].subset_index, data_size);
        }

        for (unsigned int i=0; i<edge_integral_data_.final_size(); ++i) {
            // add data to cache if there is free space, else return
        	auto range = edge_integral_data_[i].edge_side_range;
        	for (DHCellSide edge_side : range) {
        	    unsigned int side_dim = edge_side.dim();
                unsigned int data_size = eval_points_->subset_size( side_dim, edge_integral_data_[i].subset_index ) / (side_dim+1);
                unsigned int start_point = data_size * edge_side.side_idx();
                element_cache_map_.mark_used_eval_points(edge_side.cell(), edge_integral_data_[i].subset_index, data_size, start_point);
        	}
        }

        for (unsigned int i=0; i<coupling_integral_data_.final_size(); ++i) {
            // add data to cache if there is free space, else return
            unsigned int bulk_data_size = eval_points_->subset_size( coupling_integral_data_[i].cell.dim(), coupling_integral_data_[i].bulk_subset_index );
            element_cache_map_.mark_used_eval_points(coupling_integral_data_[i].cell, coupling_integral_data_[i].bulk_subset_index, bulk_data_size);

            unsigned int side_dim = coupling_integral_data_[i].side.dim();
            unsigned int side_data_size = eval_points_->subset_size( side_dim, coupling_integral_data_[i].side_subset_index ) / (side_dim+1);
            unsigned int start_point = side_data_size * coupling_integral_data_[i].side.side_idx();
            element_cache_map_.mark_used_eval_points(coupling_integral_data_[i].side.cell(), coupling_integral_data_[i].side_subset_index, side_data_size, start_point);
        }

        for (unsigned int i=0; i<boundary_integral_data_.final_size(); ++i) {
            // add data to cache if there is free space, else return
            unsigned int bdr_data_size = eval_points_->subset_size( boundary_integral_data_[i].side.cond().element_accessor().dim(), boundary_integral_data_[i].bdr_subset_index );
            element_cache_map_.mark_used_eval_points(boundary_integral_data_[i].side.cond().element_accessor(), boundary_integral_data_[i].bdr_subset_index, bdr_data_size);

            unsigned int side_dim = boundary_integral_data_[i].side.dim();
            unsigned int data_size = eval_points_->subset_size( side_dim, boundary_integral_data_[i].side_subset_index ) / (side_dim+1);
            unsigned int start_point = data_size * boundary_integral_data_[i].side.side_idx();
            element_cache_map_.mark_used_eval_points(boundary_integral_data_[i].side.cell(), boundary_integral_data_[i].side_subset_index, data_size, start_point);
        }
    }

    /**
     * Add data of integrals to appropriate structure and register elements to ElementCacheMap.
     *
     * Types of used integrals must be set in data member \p active_integrals_.
     */
    void add_integrals_of_computing_step(DHCellAccessor cell) {
        element_cache_map_.add(cell);

        if (active_integrals_ & ActiveIntegrals::bulk)
    	    if (cell.is_own()) { // Not ghost
                this->add_volume_integral(cell);
    	    }

        for( DHCellSide cell_side : cell.side_range() ) {
            if (active_integrals_ & ActiveIntegrals::boundary)
                if (cell.is_own()) // Not ghost
                    if ( (cell_side.side().edge().n_sides() == 1) && (cell_side.side().is_boundary()) ) {
                        this->add_boundary_integral(cell_side);
                        element_cache_map_.add(cell_side);
                        auto bdr_elm_acc = cell_side.cond().element_accessor();
                        element_cache_map_.add(bdr_elm_acc);
                        continue;
                    }
            if (active_integrals_ & ActiveIntegrals::edge)
                if ( (cell_side.n_edge_sides() >= 2) && (cell_side.edge_sides().begin()->element().idx() == cell.elm_idx())) {
                    this->add_edge_integral(cell_side);
                	for( DHCellSide edge_side : cell_side.edge_sides() ) {
                		element_cache_map_.add(edge_side);
                    }
                }
        }

        if (active_integrals_ & ActiveIntegrals::coupling)
            for( DHCellSide neighb_side : cell.neighb_sides() ) { // cell -> elm lower dim, neighb_side -> elm higher dim
                if (cell.dim() != neighb_side.dim()-1) continue;
                this->add_coupling_integral(cell, neighb_side);
                element_cache_map_.add(cell);
                element_cache_map_.add(neighb_side);
            }
    }

    /// Add data of volume integral to appropriate data structure.
    void add_volume_integral(const DHCellAccessor &cell) {
        BulkIntegralData data;
        data.cell = cell;
        data.subset_index = integrals_.bulk_[cell.dim()-1]->get_subset_idx();
        bulk_integral_data_.add(data);

        unsigned int reg_idx = cell.elm().region_idx().idx();
        for (auto p : integrals_.bulk_[cell.dim()-1]->points(cell, &element_cache_map_) ) {
            EvalPointData epd(reg_idx, cell.elm_idx(), p.eval_point_idx());
        	eval_point_data_.add(epd);
        }
    }

    /// Add data of edge integral to appropriate data structure.
    void add_edge_integral(const DHCellSide &cell_side) {
        EdgeIntegralData data;
        data.edge_side_range = cell_side.edge_sides();
        data.subset_index = integrals_.edge_[data.edge_side_range.begin()->dim()-1]->get_subset_idx();
        edge_integral_data_.add(data);

        unsigned int reg_idx = cell_side.element().region_idx().idx();
        for (auto p : integrals_.edge_[data.edge_side_range.begin()->dim()-1]->points(cell_side, &element_cache_map_) ) {
            EvalPointData epd(reg_idx, cell_side.elem_idx(), p.eval_point_idx());
        	eval_point_data_.add(epd);

        	auto p_ghost = p.point_on(cell_side); // point on neghbouring side on one edge
        	unsigned int ghost_reg = p_ghost.dh_cell_side().element().region_idx().idx();
        	EvalPointData epd_ghost(ghost_reg, p_ghost.dh_cell_side().elem_idx(), p_ghost.eval_point_idx());
        	eval_point_data_.add(epd_ghost);
        }
    }

    /// Add data of coupling integral to appropriate data structure.
    void add_coupling_integral(const DHCellAccessor &cell, const DHCellSide &ngh_side) {
        CouplingIntegralData data;
        data.cell = cell;
        data.side = ngh_side;
        data.bulk_subset_index = integrals_.coupling_[cell.dim()-1]->get_subset_low_idx();
        data.side_subset_index = integrals_.coupling_[cell.dim()-1]->get_subset_high_idx();
        coupling_integral_data_.add(data);

        unsigned int reg_idx_low = cell.elm().region_idx().idx();
        unsigned int reg_idx_high = ngh_side.element().region_idx().idx();
        for (auto p : integrals_.coupling_[ngh_side.dim()-1]->points(ngh_side, &element_cache_map_) ) {
            EvalPointData epd(reg_idx_high, ngh_side.elem_idx(), p.eval_point_idx());
        	eval_point_data_.add(epd);

        	auto p_low = p.lower_dim(cell); // equivalent point on low dim cell
        	EvalPointData epd_low(reg_idx_low, cell.elm_idx(), p_low.eval_point_idx());
        	eval_point_data_.add(epd_low);
        }
    }

    /// Add data of boundary integral to appropriate data structure.
    void add_boundary_integral(const DHCellSide &bdr_side) {
        BoundaryIntegralData data;
        data.side = bdr_side;
        data.bdr_subset_index = integrals_.boundary_[bdr_side.dim()-1]->get_subset_low_idx();
        data.side_subset_index = integrals_.boundary_[bdr_side.dim()-1]->get_subset_high_idx();
        boundary_integral_data_.add(data);

        unsigned int reg_idx = bdr_side.element().region_idx().idx();
        for (auto p : integrals_.boundary_[bdr_side.dim()-1]->points(bdr_side, &element_cache_map_) ) {
            EvalPointData epd(reg_idx, bdr_side.elem_idx(), p.eval_point_idx());
        	eval_point_data_.add(epd);

        	auto p_bdr = p.point_bdr(bdr_side.cond().element_accessor()); // equivalent point on boundary element
        	unsigned int bdr_reg = p_bdr.elm_accessor().region_idx().idx();
        	EvalPointData epd_bdr(bdr_reg, p_bdr.elm_accessor().mesh_idx(), p_bdr.eval_point_idx());
        	eval_point_data_.add(epd_bdr);
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
    TmpSizeList<BulkIntegralData>       bulk_integral_data_;      ///< Holds data for computing bulk integrals.
    TmpSizeList<EdgeIntegralData>       edge_integral_data_;      ///< Holds data for computing edge integrals.
    TmpSizeList<CouplingIntegralData>   coupling_integral_data_;  ///< Holds data for computing couplings integrals.
    TmpSizeList<BoundaryIntegralData>   boundary_integral_data_;  ///< Holds data for computing boundary integrals.
    TmpSizeList<EvalPointData>          eval_point_data_;         ///< Holds data of evaluating points in patch.
};


#endif /* GENERIC_ASSEMBLY_HH_ */
