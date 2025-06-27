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
#include "fem/eval_subset.hh"
#include "fem/eval_points.hh"
#include "fem/element_cache_map.hh"
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


/**
 * Common interface class for all Assembly classes.
 */
class GenericAssemblyBase
{
public:
	/**
	 * Helper structure holds data of cell (bulk) integral
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
	 * Helper structure holds data of edge integral
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
	 * Helper structure holds data of neighbour (coupling) integral
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
	 * Helper structure holds data of boundary integral
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

    GenericAssemblyBase()
    {}

    virtual ~GenericAssemblyBase(){}
    virtual void assemble(std::shared_ptr<DOFHandlerMultiDim> dh) = 0;

    /// Getter to EvalPoints object
    inline std::shared_ptr<EvalPoints> eval_points() const {
        return eval_points_;
    }

protected:
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
    /**
     * Constructor.
     *
     * Used in equations working with 'old' FeValues objects in evaluation.
     * @param eq_fields   Descendant of FieldSet declared in equation
     * @param eq_data     Object defined in equation containing shared data of eqation and assembly class.
     */
    GenericAssembly( typename DimAssembly<1>::EqFields *eq_fields, typename DimAssembly<1>::EqData *eq_data)
    : GenericAssemblyBase(),
      use_patch_fe_values_(false),
      multidim_assembly_(eq_fields, eq_data)
    {
    	initialize();
    }

    /**
     * Constructor.
     *
     * Used in equations working with 'new' PatchFeValues objects in evaluation.
     * @param eq_fields   Descendant of FieldSet declared in equation
     * @param eq_data     Object defined in equation containing shared data of eqation and assembly class.
     * @param dh          DOF handler object
     */
    GenericAssembly( typename DimAssembly<1>::EqFields *eq_fields, typename DimAssembly<1>::EqData *eq_data, DOFHandlerMultiDim* dh)
    : GenericAssemblyBase(),
      fe_values_(eq_data->quad_order(), dh->ds()->fe()),
      use_patch_fe_values_(true),
      multidim_assembly_(eq_fields, eq_data, &this->fe_values_)
    {
    	initialize();
    }

    /// Getter to set of assembly objects
    inline MixedPtr<DimAssembly, 1> multidim_assembly() const {
        return multidim_assembly_;
    }

    /// Allows rewrite number of minimal edge sides.
    void set_min_edge_sides(unsigned int val) {
    	multidim_assembly_[1_d]->set_min_edge_sides(val);
    	multidim_assembly_[2_d]->set_min_edge_sides(val);
    	multidim_assembly_[3_d]->set_min_edge_sides(val);
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
            unsigned int cell_dim = cell_it->dim();
            if (!add_into_patch) {
        	    element_cache_map_.start_elements_update();
        	    add_into_patch = true;
            }

            START_TIMER("add_integrals_to_patch");
            bool is_patch_full = false;
            switch ( cell_dim ) {
            case 1:
                is_patch_full = multidim_assembly_[1_d]->add_integrals_of_computing_step(*cell_it, table_sizes_tmp_);
                break;
            case 2:
                is_patch_full = multidim_assembly_[2_d]->add_integrals_of_computing_step(*cell_it, table_sizes_tmp_);
                break;
            case 3:
                is_patch_full = multidim_assembly_[3_d]->add_integrals_of_computing_step(*cell_it, table_sizes_tmp_);
                break;
            default:
                ASSERT(false).error("Should not happen!");
            }
            END_TIMER("add_integrals_to_patch");

            if (is_patch_full) {
                element_cache_map_.eval_point_data_.revert_temporary();
                this->assemble_integrals();
                add_into_patch = false;
            } else {
                element_cache_map_.make_paermanent_eval_points();
                if (use_patch_fe_values_) {
                    table_sizes_.copy(table_sizes_tmp_);
                }
                if (element_cache_map_.get_simd_rounded_size() == CacheMapElementNumber::get()) {
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
        END_TIMER( DimAssembly<1>::name() );
    }

    /// Getter to ElementCacheMap
    inline const ElementCacheMap &cache_map() const {
        return element_cache_map_;
    }

private:
    /// Common part of GenericAssemblz constructors.
    void initialize() {
        eval_points_ = std::make_shared<EvalPoints>();
        // first step - create integrals, then - initialize cache and initialize subobject of dimensions
        eval_points_->create_integrals( {
            multidim_assembly_[1_d]->integrals(),
            multidim_assembly_[2_d]->integrals(),
		    multidim_assembly_[3_d]->integrals()
        } );
        multidim_assembly_[1_d]->set_eval_points(eval_points_);
        multidim_assembly_[2_d]->set_eval_points(eval_points_);
        multidim_assembly_[3_d]->set_eval_points(eval_points_);
        element_cache_map_.init(eval_points_);
        multidim_assembly_[1_d]->initialize(&element_cache_map_);
        multidim_assembly_[2_d]->initialize(&element_cache_map_);
        multidim_assembly_[3_d]->initialize(&element_cache_map_);
        if (use_patch_fe_values_) {
            fe_values_.init_finalize();
        }
        active_integrals_ = multidim_assembly_[1_d]->n_active_integrals();
    }

    /// Call assemblations when patch is filled
    void assemble_integrals() {
        START_TIMER("create_patch");
        element_cache_map_.create_patch();
        END_TIMER("create_patch");
        if (use_patch_fe_values_) {
            START_TIMER("patch_reinit");
            patch_reinit();
            END_TIMER("patch_reinit");
        }
        START_TIMER("cache_update");
        multidim_assembly_[1_d]->eq_fields_->cache_update(element_cache_map_); // TODO replace with sub FieldSet
        END_TIMER("cache_update");
        element_cache_map_.finish_elements_update();

        {
            START_TIMER("assemble_volume_integrals");
            multidim_assembly_[1_d]->assemble_cell_integrals();
            multidim_assembly_[2_d]->assemble_cell_integrals();
            multidim_assembly_[3_d]->assemble_cell_integrals();
            END_TIMER("assemble_volume_integrals");
        }

        {
            START_TIMER("assemble_fluxes_boundary");
            multidim_assembly_[1_d]->assemble_boundary_side_integrals();
            multidim_assembly_[2_d]->assemble_boundary_side_integrals();
            multidim_assembly_[3_d]->assemble_boundary_side_integrals();
            END_TIMER("assemble_fluxes_boundary");
        }

        {
            START_TIMER("assemble_fluxes_elem_elem");
            multidim_assembly_[1_d]->assemble_edge_integrals();
            multidim_assembly_[2_d]->assemble_edge_integrals();
            multidim_assembly_[3_d]->assemble_edge_integrals();
            END_TIMER("assemble_fluxes_elem_elem");
        }

        {
            START_TIMER("assemble_fluxes_elem_side");
            multidim_assembly_[1_d]->assemble_neighbour_integrals();
            multidim_assembly_[2_d]->assemble_neighbour_integrals();
            END_TIMER("assemble_fluxes_elem_side");
        }
        // clean integral data
        multidim_assembly_[1_d]->clean_integral_data();
        multidim_assembly_[2_d]->clean_integral_data();
        multidim_assembly_[3_d]->clean_integral_data();
        element_cache_map_.clear_element_eval_points_map();
        if (use_patch_fe_values_) {
            table_sizes_.reset();
            table_sizes_tmp_.reset();
            fe_values_.reset();
        }
    }

    /// Reinit PatchFeValues object during construction of patch
    void patch_reinit() {
    	fe_values_.resize_tables(table_sizes_);

        multidim_assembly_[1_d]->add_patch_bulk_points();
        multidim_assembly_[2_d]->add_patch_bulk_points();
        multidim_assembly_[3_d]->add_patch_bulk_points();

        multidim_assembly_[1_d]->add_patch_bdr_side_points();
        multidim_assembly_[2_d]->add_patch_bdr_side_points();
        multidim_assembly_[3_d]->add_patch_bdr_side_points();

        multidim_assembly_[1_d]->add_patch_edge_points();
        multidim_assembly_[2_d]->add_patch_edge_points();
        multidim_assembly_[3_d]->add_patch_edge_points();

        multidim_assembly_[1_d]->add_patch_coupling_integrals();
        multidim_assembly_[2_d]->add_patch_coupling_integrals();

        this->fe_values_.reinit_patch();
    }

    /// Calls cache_reallocate method on
    inline void reallocate_cache() {
        multidim_assembly_[1_d]->eq_fields_->cache_reallocate(this->element_cache_map_, multidim_assembly_[1_d]->used_fields_);
        // DebugOut() << "Order of evaluated fields (" << DimAssembly<1>::name() << "):" << multidim_assembly_[1_d]->eq_fields_->print_dependency();
    }


    PatchFEValues<3> fe_values_;                                     ///< Common FEValues object over all dimensions
    bool use_patch_fe_values_;                                       ///< Flag holds if common @p fe_values_ object is used in @p multidim_assembly_
    MixedPtr<DimAssembly, 1> multidim_assembly_;                     ///< Assembly object

    /// Holds mask of active integrals.
    int active_integrals_;

    /// Struct for pre-computing number of elements, sides, bulk points and side points on each dimension.
    PatchFEValues<3>::TableSizes table_sizes_;
    /// Same as previous but hold temporary values during adding elements, sides and points.
    PatchFEValues<3>::TableSizes table_sizes_tmp_;
};


#endif /* GENERIC_ASSEMBLY_HH_ */
