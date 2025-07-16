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


/// Holds common data shared between GenericAssemblz and Assembly<dim> classes.
struct AssemblyInternals {
public:
    AssemblyInternals() {}

    AssemblyInternals(MixedPtr<FiniteElement> fe) : fe_values_(fe)  {}

    std::shared_ptr<EvalPoints> eval_points_;                     ///< EvalPoints object shared by all integrals
    ElementCacheMap element_cache_map_;                           ///< ElementCacheMap according to EvalPoints
    PatchFEValues<3> fe_values_;                                  ///< Common FEValues object over all dimensions

    /// Struct for pre-computing number of elements, sides, bulk points and side points on each dimension.
    PatchFEValues<3>::TableSizes table_sizes_;
    /// Same as previous but hold temporary values during adding elements, sides and points.
    PatchFEValues<3>::TableSizes table_sizes_tmp_;
};


/**
 * Common interface class for all Assembly classes.
 */
class GenericAssemblyBase
{
public:
    GenericAssemblyBase()
    {}

    GenericAssemblyBase(MixedPtr<FiniteElement> fe)
    : asm_internals_(fe)
    {}

    virtual ~GenericAssemblyBase(){}
    virtual void assemble(std::shared_ptr<DOFHandlerMultiDim> dh) = 0;

    /// Getter to EvalPoints object
    inline std::shared_ptr<EvalPoints> eval_points() const {
        return asm_internals_.eval_points_;
    }

protected:
//    std::shared_ptr<EvalPoints> eval_points_;                     ///< EvalPoints object shared by all integrals
//    ElementCacheMap element_cache_map_;                           ///< ElementCacheMap according to EvalPoints
    AssemblyInternals asm_internals_;                             ///< Holds shared internals data
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
      multidim_assembly_(eq_fields, eq_data, &this->asm_internals_)
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
    : GenericAssemblyBase(dh->ds()->fe()),
      use_patch_fe_values_(true),
      multidim_assembly_(eq_fields, eq_data, &this->asm_internals_)
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
                asm_internals_.element_cache_map_.start_elements_update();
        	    add_into_patch = true;
            }

            START_TIMER("add_integrals_to_patch");
            bool is_patch_full = false;
            switch ( cell_dim ) {
            case 1:
                is_patch_full = multidim_assembly_[1_d]->add_integrals_of_computing_step(*cell_it);
                break;
            case 2:
                is_patch_full = multidim_assembly_[2_d]->add_integrals_of_computing_step(*cell_it);
                break;
            case 3:
                is_patch_full = multidim_assembly_[3_d]->add_integrals_of_computing_step(*cell_it);
                break;
            default:
                ASSERT(false).error("Should not happen!");
            }
            END_TIMER("add_integrals_to_patch");

            if (is_patch_full) {
                asm_internals_.element_cache_map_.eval_point_data_.revert_temporary();
                this->assemble_integrals();
                add_into_patch = false;
            } else {
                asm_internals_.element_cache_map_.make_paermanent_eval_points();
                if (use_patch_fe_values_) {
                    asm_internals_.table_sizes_.copy(asm_internals_.table_sizes_tmp_);
                }
                if (asm_internals_.element_cache_map_.get_simd_rounded_size() == CacheMapElementNumber::get()) {
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
        return asm_internals_.element_cache_map_;
    }

private:
    /// Common part of GenericAssemblz constructors.
    void initialize() {
        asm_internals_.eval_points_ = std::make_shared<EvalPoints>();
        // first step - create integrals, then - initialize cache and initialize subobject of dimensions
        asm_internals_.eval_points_->create_integrals( {
            multidim_assembly_[1_d]->integrals(),
            multidim_assembly_[2_d]->integrals(),
		    multidim_assembly_[3_d]->integrals()
        } );
        asm_internals_.element_cache_map_.init(asm_internals_.eval_points_);
        multidim_assembly_[1_d]->initialize();
        multidim_assembly_[2_d]->initialize();
        multidim_assembly_[3_d]->initialize();
        if (use_patch_fe_values_) {
            asm_internals_.fe_values_.init_finalize();
        }
    }

    /// Call assemblations when patch is filled
    void assemble_integrals() {
        START_TIMER("create_patch");
        asm_internals_.element_cache_map_.create_patch();
        END_TIMER("create_patch");
        if (use_patch_fe_values_) {
            START_TIMER("patch_reinit");
            patch_reinit();
            END_TIMER("patch_reinit");
        }
        START_TIMER("cache_update");
        multidim_assembly_[1_d]->eq_fields_->cache_update(asm_internals_.element_cache_map_); // TODO replace with sub FieldSet
        END_TIMER("cache_update");
        asm_internals_.element_cache_map_.finish_elements_update();

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
        asm_internals_.element_cache_map_.clear_element_eval_points_map();
        if (use_patch_fe_values_) {
            asm_internals_.table_sizes_.reset();
            asm_internals_.table_sizes_tmp_.reset();
            asm_internals_.fe_values_.reset();
        }
    }

    /// Reinit PatchFeValues object during construction of patch
    void patch_reinit() {
        asm_internals_.fe_values_.resize_tables(asm_internals_.table_sizes_);

        asm_internals_.fe_values_.add_patch_points(multidim_assembly_[1_d]->integrals(), multidim_assembly_[1_d]->integral_data(), &asm_internals_.element_cache_map_);
        asm_internals_.fe_values_.add_patch_points(multidim_assembly_[2_d]->integrals(), multidim_assembly_[2_d]->integral_data(), &asm_internals_.element_cache_map_);
        asm_internals_.fe_values_.add_patch_points(multidim_assembly_[3_d]->integrals(), multidim_assembly_[3_d]->integral_data(), &asm_internals_.element_cache_map_);

        asm_internals_.fe_values_.reinit_patch();
    }

    /// Calls cache_reallocate method on
    inline void reallocate_cache() {
        multidim_assembly_[1_d]->eq_fields_->cache_reallocate(asm_internals_.element_cache_map_, multidim_assembly_[1_d]->used_fields_);
        // DebugOut() << "Order of evaluated fields (" << DimAssembly<1>::name() << "):" << multidim_assembly_[1_d]->eq_fields_->print_dependency();
    }


//    PatchFEValues<3> fe_values_;                                     ///< Common FEValues object over all dimensions
    bool use_patch_fe_values_;                                       ///< Flag holds if common @p fe_values_ object is used in @p multidim_assembly_
    MixedPtr<DimAssembly, 1> multidim_assembly_;                     ///< Assembly object

//    /// Struct for pre-computing number of elements, sides, bulk points and side points on each dimension.
//    PatchFEValues<3>::TableSizes table_sizes_;
//    /// Same as previous but hold temporary values during adding elements, sides and points.
//    PatchFEValues<3>::TableSizes table_sizes_tmp_;
};


#endif /* GENERIC_ASSEMBLY_HH_ */
