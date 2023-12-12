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
#include "fem/fe_values.hh"
#include "fem/patch_fe_values.hh"
#include "tools/revertable_list.hh"
#include "system/sys_profiler.hh"



/// Allow set mask of active integrals.
enum ActiveIntegrals {
    no_intg  =      0,
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
 * Common interface class for all Assembly classes.
 */
class GenericAssemblyBase
{
public:
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

        /// Copy constructor
    	EdgeIntegralData(const EdgeIntegralData &other)
        : edge_side_range(other.edge_side_range), subset_index(other.subset_index) {}

        /// Constructor with data mebers initialization
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

    GenericAssemblyBase() {}

    virtual ~GenericAssemblyBase(){}
    virtual void assemble(std::shared_ptr<DOFHandlerMultiDim> dh) = 0;

    /// Getter to EvalPoints object
    inline std::shared_ptr<EvalPoints> eval_points() const {
        return eval_points_;
    }

protected:
    AssemblyIntegrals integrals_;                                 ///< Holds integral objects.
    std::shared_ptr<EvalPoints> eval_points_;                     ///< EvalPoints object shared by all integrals
    ElementCacheMap element_cache_map_;                           ///< ElementCacheMap according to EvalPoints
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
class GenericAssembly : public GenericAssemblyBase
{
public:
    /// Constructor
    GenericAssembly( typename DimAssembly<1>::EqFields *eq_fields, typename DimAssembly<1>::EqData *eq_data)
    : multidim_assembly_(eq_fields, eq_data),
	  min_edge_sides_(2),
	  bulk_integral_data_(20, 10),
	  edge_integral_data_(12, 6),
	  coupling_integral_data_(12, 6),
	  boundary_integral_data_(8, 4)
    {
    	initialize();
    }

    /// Constructor
    GenericAssembly( typename DimAssembly<1>::EqFields *eq_fields, typename DimAssembly<1>::EqData *eq_data, DOFHandlerMultiDim* dh)
    : fe_values_(CacheMapElementNumber::get(), dh->ds()->fe()),
      multidim_assembly_(eq_fields, eq_data, &this->fe_values_),
      min_edge_sides_(2),
      bulk_integral_data_(20, 10),
      edge_integral_data_(12, 6),
      coupling_integral_data_(12, 6),
      boundary_integral_data_(8, 4)
    {
    	initialize();
    }

    /// Getter to set of assembly objects
    inline MixedPtr<DimAssembly, 1> multidim_assembly() const {
        return multidim_assembly_;
    }

    void set_min_edge_sides(unsigned int val) {
        min_edge_sides_ = val;
    }

	/**
	 * @brief General assemble methods.
	 *
	 * Loops through local cells and calls assemble methods of assembly
	 * object of each cells over space dimension.
	 *
	 * TODO:
	 * - make estimate of the cache fill for combination of (integral_type x element dimension)
	 * - add next cell to patch if current_patch_size + next_element_size <= fixed_cache_size
	 * - avoid reverting the integral data lists.
	 */
    void assemble(std::shared_ptr<DOFHandlerMultiDim> dh) override {
        START_TIMER( DimAssembly<1>::name() );
        this->reallocate_cache();
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

            if (element_cache_map_.get_simd_rounded_size() > CacheMapElementNumber::get()) {
                bulk_integral_data_.revert_temporary();
                edge_integral_data_.revert_temporary();
                coupling_integral_data_.revert_temporary();
                boundary_integral_data_.revert_temporary();
                element_cache_map_.eval_point_data_.revert_temporary();
                this->assemble_integrals(dh);
                add_into_patch = false;
            } else {
                bulk_integral_data_.make_permanent();
                edge_integral_data_.make_permanent();
                coupling_integral_data_.make_permanent();
                boundary_integral_data_.make_permanent();
                element_cache_map_.eval_point_data_.make_permanent();
                if (element_cache_map_.get_simd_rounded_size() == CacheMapElementNumber::get()) {
                    this->assemble_integrals(dh);
                    add_into_patch = false;
                }
                ++cell_it;
            }
        }
        if (add_into_patch) {
            this->assemble_integrals(dh);
        }

        multidim_assembly_[1_d]->end();
        END_TIMER( DimAssembly<1>::name() );
    }

    /// Return ElementCacheMap
    inline const ElementCacheMap &cache_map() const {
        return element_cache_map_;
    }

private:
    /// Common part of GenericAssemblz constructors.
    void initialize() {
        eval_points_ = std::make_shared<EvalPoints>();
        // first step - create integrals, then - initialize cache and initialize subobject of dimensions
        multidim_assembly_[1_d]->create_integrals(eval_points_, integrals_);
        multidim_assembly_[2_d]->create_integrals(eval_points_, integrals_);
        multidim_assembly_[3_d]->create_integrals(eval_points_, integrals_);
        element_cache_map_.init(eval_points_);
        multidim_assembly_[1_d]->initialize(&element_cache_map_);
        multidim_assembly_[2_d]->initialize(&element_cache_map_);
        multidim_assembly_[3_d]->initialize(&element_cache_map_);
        active_integrals_ = multidim_assembly_[1_d]->n_active_integrals();
    }

    /// Call assemblations when patch is filled
    void assemble_integrals(std::shared_ptr<DOFHandlerMultiDim> dh) {
        START_TIMER("create_patch");
        element_cache_map_.create_patch();
        END_TIMER("create_patch");
        START_TIMER("patch_reinit");
        patch_reinit(dh);
        END_TIMER("patch_reinit");
        START_TIMER("cache_update");
        multidim_assembly_[1_d]->eq_fields_->cache_update(element_cache_map_); // TODO replace with sub FieldSet
        END_TIMER("cache_update");
        element_cache_map_.finish_elements_update();

        {
            START_TIMER("assemble_volume_integrals");
            multidim_assembly_[1_d]->assemble_cell_integrals(bulk_integral_data_);
            multidim_assembly_[2_d]->assemble_cell_integrals(bulk_integral_data_);
            multidim_assembly_[3_d]->assemble_cell_integrals(bulk_integral_data_);
            END_TIMER("assemble_volume_integrals");
        }

        {
            START_TIMER("assemble_fluxes_boundary");
            multidim_assembly_[1_d]->assemble_boundary_side_integrals(boundary_integral_data_);
            multidim_assembly_[2_d]->assemble_boundary_side_integrals(boundary_integral_data_);
            multidim_assembly_[3_d]->assemble_boundary_side_integrals(boundary_integral_data_);
            END_TIMER("assemble_fluxes_boundary");
        }

        {
            START_TIMER("assemble_fluxes_elem_elem");
            multidim_assembly_[1_d]->assemble_edge_integrals(edge_integral_data_);
            multidim_assembly_[2_d]->assemble_edge_integrals(edge_integral_data_);
            multidim_assembly_[3_d]->assemble_edge_integrals(edge_integral_data_);
            END_TIMER("assemble_fluxes_elem_elem");
        }

        {
            START_TIMER("assemble_fluxes_elem_side");
            multidim_assembly_[2_d]->assemble_neighbour_integrals(coupling_integral_data_);
            multidim_assembly_[3_d]->assemble_neighbour_integrals(coupling_integral_data_);
            END_TIMER("assemble_fluxes_elem_side");
        }
        // clean integral data
        bulk_integral_data_.reset();
        edge_integral_data_.reset();
        coupling_integral_data_.reset();
        boundary_integral_data_.reset();
        element_cache_map_.clear_element_eval_points_map();
    }

    void patch_reinit(std::shared_ptr<DOFHandlerMultiDim> dh) {
        const std::vector<unsigned int> &elm_idx_vec = element_cache_map_.elm_idx_vec();
        std::array<PatchElementsList, 4> patch_elements;

        for (unsigned int i=0; i<elm_idx_vec.size(); ++i) {
            // Skip invalid element indices.
            if ( elm_idx_vec[i] == std::numeric_limits<unsigned int>::max() ) continue;

            ElementAccessor<3> elm(dh->mesh(), elm_idx_vec[i]);
            patch_elements[elm.dim()].push_back(std::make_pair(elm, i));
        }
        multidim_assembly_[1_d]->patch_reinit(patch_elements);
        multidim_assembly_[2_d]->patch_reinit(patch_elements);
        multidim_assembly_[3_d]->patch_reinit(patch_elements);
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
                if ( (cell_side.n_edge_sides() >= min_edge_sides_) && (cell_side.edge_sides().begin()->element().idx() == cell.elm_idx())) {
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
            element_cache_map_.add_eval_point(reg_idx, cell.elm_idx(), i, cell.local_idx());
        }
    }

    /// Add data of edge integral to appropriate data structure.
    inline void add_edge_integral(const DHCellSide &cell_side) {
        auto range = cell_side.edge_sides();
        edge_integral_data_.emplace_back(range, integrals_.edge_[range.begin()->dim()-1]->get_subset_idx());

        for( DHCellSide edge_side : range ) {
            unsigned int reg_idx = edge_side.element().region_idx().idx();
            for (auto p : integrals_.edge_[range.begin()->dim()-1]->points(edge_side, &element_cache_map_) ) {
                element_cache_map_.add_eval_point(reg_idx, edge_side.elem_idx(), p.eval_point_idx(), edge_side.cell().local_idx());
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
            element_cache_map_.add_eval_point(reg_idx_high, ngh_side.elem_idx(), p.eval_point_idx(), ngh_side.cell().local_idx());

        	if (add_low) {
                auto p_low = p.lower_dim(cell); // equivalent point on low dim cell
                element_cache_map_.add_eval_point(reg_idx_low, cell.elm_idx(), p_low.eval_point_idx(), cell.local_idx());
        	}
        }
    }

    /// Add data of boundary integral to appropriate data structure.
    inline void add_boundary_integral(const DHCellSide &bdr_side) {
        boundary_integral_data_.emplace_back(integrals_.boundary_[bdr_side.dim()-1]->get_subset_low_idx(), bdr_side,
                integrals_.boundary_[bdr_side.dim()-1]->get_subset_high_idx());

        unsigned int reg_idx = bdr_side.element().region_idx().idx();
        for (auto p : integrals_.boundary_[bdr_side.dim()-1]->points(bdr_side, &element_cache_map_) ) {
            element_cache_map_.add_eval_point(reg_idx, bdr_side.elem_idx(), p.eval_point_idx(), bdr_side.cell().local_idx());

        	BulkPoint p_bdr = p.point_bdr(bdr_side.cond().element_accessor()); // equivalent point on boundary element
        	unsigned int bdr_reg = bdr_side.cond().element_accessor().region_idx().idx();
        	// invalid local_idx value, DHCellAccessor of boundary element doesn't exist
        	element_cache_map_.add_eval_point(bdr_reg, bdr_side.cond().bc_ele_idx(), p_bdr.eval_point_idx(), -1);
        }
    }

    /// Calls cache_reallocate method on
    inline void reallocate_cache() {
        multidim_assembly_[1_d]->eq_fields_->cache_reallocate(this->element_cache_map_, multidim_assembly_[1_d]->used_fields_);
        // DebugOut() << "Order of evaluated fields (" << DimAssembly<1>::name() << "):" << multidim_assembly_[1_d]->eq_fields_->print_dependency();
    }


    PatchFEValues<3> fe_values_;                                     ///< Common FEValues object over all dimensions
    MixedPtr<DimAssembly, 1> multidim_assembly_;                     ///< Assembly object

    /// Holds mask of active integrals.
    int active_integrals_;

    /**
     * Minimal number of sides on edge.
     *
     * Edge integral is created and calculated if number of sides is greater or equal than this value. Default value
     * is 2 and can be changed
     */
    unsigned int min_edge_sides_;

    // Following variables hold data of all integrals depending of actual computed element.
    // TODO sizes of arrays should be set dynamically, depend on number of elements in ElementCacheMap,
    RevertableList<BulkIntegralData>       bulk_integral_data_;      ///< Holds data for computing bulk integrals.
    RevertableList<EdgeIntegralData>       edge_integral_data_;      ///< Holds data for computing edge integrals.
    RevertableList<CouplingIntegralData>   coupling_integral_data_;  ///< Holds data for computing couplings integrals.
    RevertableList<BoundaryIntegralData>   boundary_integral_data_;  ///< Holds data for computing boundary integrals.
};


#endif /* GENERIC_ASSEMBLY_HH_ */
