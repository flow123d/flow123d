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
 * @file    assembly_base.hh
 * @brief
 */

#ifndef ASSEMBLY_BASE_HH_
#define ASSEMBLY_BASE_HH_


#include "coupling/generic_assembly.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/eval_points.hh"
#include "fem/element_cache_map.hh"
#include "fem/update_flags.hh"



/**
 * Base class define empty methods, these methods can be overwite in descendants.
 */
template <unsigned int dim>
class AssemblyBasePatch
{
public:
    /**
     * Constructor
     *
     * @param quad_order    Order of Quadrature objects.
     * @param patch_internals Holds shared data with GenericAssembly
     */
    AssemblyBasePatch(unsigned int quad_order, PatchInternals *patch_internals)
    : AssemblyBasePatch<dim>() {
    	patch_internals_ = patch_internals;
        quad_ = new QGauss(dim, 2*quad_order);
        quad_low_ = new QGauss(dim-1, 2*quad_order);
    }

	/// Destructor
    virtual ~AssemblyBasePatch() {
        delete quad_;
        delete quad_low_;
    }

    /**
     * Assembles the volume integrals on cell.
     *
     * Method can be overridden and implemented in descendant
     */
    virtual inline void cell_integral(FMT_UNUSED DHCellAccessor cell, FMT_UNUSED unsigned int element_patch_idx) {}

    /**
     * Assembles the fluxes on the boundary.
     *
     * Method can be overridden and implemented in descendant
     */
    virtual inline void boundary_side_integral(FMT_UNUSED DHCellSide cell_side) {}

    /**
     * Assembles the fluxes between sides on the edge.
     *
     * Method can be overwrite and implement in descendant
     */
    virtual inline void edge_integral(FMT_UNUSED RangeConvert<DHEdgeSide, DHCellSide> edge_side_range) {}

    /**
     * Assembles the fluxes between elements of different dimensions.
     *
     * Method can be overridden and implemented in descendant
     */
    virtual inline void dimjoin_intergral(FMT_UNUSED DHCellAccessor cell_lower_dim, FMT_UNUSED DHCellSide neighb_side) {}

    /**
     * Method prepares object before assemblation (e.g. balance, ...).
     *
     * Method can be overridden and implemented in descendant
     */
    virtual void begin() {}

    /**
     * Method finishes object after assemblation (e.g. balance, ...).
     *
     * Method can be overridden and implemented in descendant
     */
    virtual void end() {}

    /**
     * Create and return BulkIntegral accessor of given quadrature.
     *
     * Method is called from descendants during construction / initialization of assembly object.
     */
    std::shared_ptr<BulkIntegralAcc<dim>> create_bulk_integral(Quadrature *quad) {
        ASSERT_PERMANENT_EQ(quad->dim(), dim);
        std::tuple<uint, uint> tpl = IntegralTplHash::integral_tuple(dim, quad->size());
        auto result = integrals_.bulk_.insert({
                tpl,
                std::make_shared<BulkIntegralAcc<dim>>(*patch_internals_, quad)
            });
        return result.first->second;
    }

    /**
     * Create and return EdgeIntegral accessor of given quadrature.
     *
     * Method is called from descendants during construction / initialization of assembly object.
     */
    std::shared_ptr<EdgeIntegralAcc<dim>> create_edge_integral(Quadrature *quad) {
        ASSERT_PERMANENT_EQ(quad->dim()+1, dim);
        std::tuple<uint, uint> tpl = IntegralTplHash::integral_tuple(dim, quad->size());
        auto result = integrals_.edge_.insert({
                tpl,
                std::make_shared<EdgeIntegralAcc<dim>>(*patch_internals_, quad)
            });
        return result.first->second;
    }


    /**
     * Create and return CouplingIntegral accessor of given quadrature.
     *
     * Method is called from descendants during construction / initialization of assembly object.
     */
    std::shared_ptr<CouplingIntegralAcc<dim>> create_coupling_integral(Quadrature *quad) {
        if (dim == 3) return nullptr;

        ASSERT_PERMANENT_EQ(quad->dim(), dim);
        std::tuple<uint, uint> tpl = IntegralTplHash::integral_tuple(dim, quad->size());
        auto result = integrals_.coupling_.insert({
                tpl,
                std::make_shared<CouplingIntegralAcc<dim>>(*patch_internals_, quad)
            });
        return result.first->second;
    }


    /**
     * Create and return BoundaryIntegral accessor of given quadrature.
     *
     * Method is called from descendants during construction / initialization of assembly object.
     */
    std::shared_ptr<BoundaryIntegralAcc<dim>> create_boundary_integral(Quadrature *quad) {
        ASSERT_PERMANENT_EQ(quad->dim()+1, dim);
        std::tuple<uint, uint> tpl = IntegralTplHash::integral_tuple(dim, quad->size());
        auto result = integrals_.boundary_.insert({
                tpl,
                std::make_shared<BoundaryIntegralAcc<dim>>(*patch_internals_, quad)
            });
        return result.first->second;
    }


    /**
     * Add data of integrals to appropriate structure and register elements to ElementCacheMap.
     *
     * Return true if patch is full to its maximal capacity.
     * Method is called from GenericAssembly::assembly method.
     */
    bool add_integrals_of_computing_step(DHCellAccessor cell) {
        ASSERT_EQ(cell.dim(), dim);

        if (cell.is_own()) { // Not ghost
            if (integrals_.bulk_.size() > 0)
                ++( patch_internals_->fe_values_.ppv(bulk_domain, cell.dim()) ).n_mesh_items_;
            this->add_volume_integrals(cell);
        }

        for( DHCellSide cell_side : cell.side_range() ) {
            if (cell.is_own()) // Not ghost
                if ( (cell_side.side().edge().n_sides() == 1) && (cell_side.side().is_boundary()) ) {
                    this->add_boundary_integrals(cell_side);
                }
            if ( (cell_side.n_edge_sides() >= min_edge_sides_) && (cell_side.edge_sides().begin()->element().idx() == cell.elm_idx())) {
                this->add_edge_integrals(cell_side);
            }
        }

        add_coupling_integrals(cell);

        if (patch_internals_->element_cache_map_.get_simd_rounded_size() > CacheMapElementNumber::get()) {
            integrals_.revert_temporary();
            return true;
        } else {
            integrals_.make_permanent();
            return false;
        }
    }

    /**
     * Assembles the cell integrals for the given dimension.
     *
     * Method is called from GenericAssembly::assembly method.
     */
    virtual inline void assemble_cell_integrals() {
    	for (unsigned int i=0; i<integrals_.n_patch_cells(); ++i) {
            this->cell_integral(integrals_.bulk_.begin()->second->patch_data()[i].cell,
                                patch_internals_->element_cache_map_.position_in_cache(integrals_.bulk_.begin()->second->patch_data()[i].cell.elm_idx())
			                   );
    	}
    	// Possibly optimization but not so fast as we would assume (needs change interface of cell_integral)
        /*for (unsigned int i=0; i<element_cache_map_->n_elements(); ++i) {
            unsigned int elm_start = element_cache_map_->element_chunk_begin(i);
            if (element_cache_map_->eval_point_data(elm_start).i_eval_point_ != 0) continue;
            this->cell_integral(i, element_cache_map_->eval_point_data(elm_start).dh_loc_idx_);
        }*/
    }

    /**
     * Assembles the boundary side integrals for the given dimension.
     *
     * Method is called from GenericAssembly::assembly method.
     */
    inline void assemble_boundary_side_integrals() {
        for (unsigned int i=0; i<integrals_.n_patch_boundaries(); ++i) {
            this->boundary_side_integral( integrals_.boundary_.begin()->second->patch_data()[i].side );
        }
    }

    /**
     * Assembles the edge integrals for the given dimension.
     *
     * Method is called from GenericAssembly::assembly method.
     */
    inline void assemble_edge_integrals() {
        for (unsigned int i=0; i<integrals_.n_patch_edges(); ++i) {
            this->edge_integral(integrals_.edge_.begin()->second->patch_data()[i].edge_side_range);
        }
    }

    /**
     * Assembles the neighbours integrals for the given dimension.
     *
     * Method is called from GenericAssembly::assembly method.
     */
    inline void assemble_neighbour_integrals() {
        for (unsigned int i=0; i<integrals_.n_patch_neighbours(); ++i) {
            this->dimjoin_intergral(integrals_.coupling_.begin()->second->patch_data()[i].cell, integrals_.coupling_.begin()->second->patch_data()[i].side);
        }
    }

    /// Setter of min_edge_sides_
    void set_min_edge_sides(unsigned int val) {
        min_edge_sides_ = val;
    }

    /**
     * Clean all integral data structures
     *
     * Method is called from GenericAssembly::assembly method.
     */
    void clean_integral_data() {
        integrals_.reset();
    }

    /// Getter of integrals_
    const DimIntegrals<dim> &integrals() const {
    	return integrals_;
    }

protected:
	/**
	 * Default constructor.
	 *
	 * Be aware if you use this constructor. Quadrature objects must be initialized manually in descendant.
	 */
	AssemblyBasePatch()
	: quad_(nullptr), quad_low_(nullptr), patch_internals_(nullptr), min_edge_sides_(2) {}

    /**
     * Add data of volume integrals to appropriate data structure.
     *
     * Method is used internally in AssemblyBase
     */
    inline void add_volume_integrals(const DHCellAccessor &cell) {
        for (auto integral_it : integrals_.bulk_) {
            uint subset_idx = integral_it.second->get_subset_idx();
            integral_it.second->patch_data().emplace_back(cell);

            unsigned int reg_idx = cell.elm().region_idx().idx();
            // Different access than in other integrals: We can't use range method CellIntegral::points
            // because it passes element_patch_idx as argument that is not known during patch construction.
            for (uint i=uint( patch_internals_->eval_points_->subset_begin(dim, subset_idx) );
                      i<uint( patch_internals_->eval_points_->subset_end(dim, subset_idx) ); ++i) {
                patch_internals_->element_cache_map_.add_eval_point(reg_idx, cell.elm_idx(), i, cell.local_idx());
            }
        }
    }

    /**
     * Add data of edge integrals to appropriate data structure.
     *
     * Method is used internally in AssemblyBase
     */
    inline void add_edge_integrals(const DHCellSide &cell_side) {
	    auto range = cell_side.edge_sides();

        auto &ppv = patch_internals_->fe_values_.ppv(side_domain, cell_side.dim());
        for (auto integral_it : integrals_.edge_) {
            integral_it.second->patch_data().emplace_back(range);

            for( DHCellSide edge_side : range ) {
                add_side_points(integral_it.second, edge_side, ppv);
            }
        }
    }

    /**
     * Add data of boundary integrals to appropriate data structure.
     *
     * Method is used internally in AssemblyBase
     */
    inline void add_boundary_integrals(const DHCellSide &bdr_side) {
        auto &ppv_side = patch_internals_->fe_values_.ppv(side_domain, bdr_side.dim());
        auto &ppv_bdr = patch_internals_->fe_values_.ppv(bulk_domain, bdr_side.dim()-1);

        for (auto integral_it : integrals_.boundary_) {
            auto integral = integral_it.second;
            integral->patch_data().emplace_back(bdr_side);

            unsigned int reg_idx = bdr_side.element().region_idx().idx();
            ++ppv_side.n_mesh_items_;
            ++ppv_bdr.n_mesh_items_;
            for (auto p : integral->points(bdr_side) ) {
                patch_internals_->element_cache_map_.add_eval_point(reg_idx, bdr_side.elem_idx(), p.eval_point_idx(), bdr_side.cell().local_idx());

            	BulkPoint p_bdr = p.point_bdr(bdr_side.cond().element_accessor()); // equivalent point on boundary element
            	unsigned int bdr_reg = bdr_side.cond().element_accessor().region_idx().idx();
            	// invalid local_idx value, DHCellAccessor of boundary element doesn't exist
            	patch_internals_->element_cache_map_.add_eval_point(bdr_reg, bdr_side.cond().bc_ele_idx(), p_bdr.eval_point_idx(), -1);
            }
        }
    }

    /**
     * Add data of coupling integrals to appropriate data structure.
     *
     * Method is used internally in AssemblyBase
     */
    inline void add_coupling_integrals(const DHCellAccessor &cell) {
        if (dim==3) return;
        auto &ppv_low = patch_internals_->fe_values_.ppv(bulk_domain, cell.dim());
        auto &ppv_high = patch_internals_->fe_values_.ppv(side_domain, cell.dim()+1);
        uint i_int=0;
    	for (auto coupling_integral_it : integrals_.coupling_) {
    	    auto coupling_integral = coupling_integral_it.second;
            // Adds data of bulk points only if bulk point were not added during processing of bulk integral
            bool add_bulk_points = !( (integrals_.bulk_.size() > 0) & cell.is_own() );
            if (add_bulk_points) {
                // add points of low dim element only one time and only if they have not been added in BulkIntegral
                for( DHCellSide ngh_side : cell.neighb_sides() ) {
                    unsigned int reg_idx_low = cell.elm().region_idx().idx();
                    ++ppv_low.n_mesh_items_;
                    for (auto p : coupling_integral->points(ngh_side) ) {
                        auto p_low = p.lower_dim(cell); // equivalent point on low dim cell
                        patch_internals_->element_cache_map_.add_eval_point(reg_idx_low, cell.elm_idx(), p_low.eval_point_idx(), cell.local_idx());
                    }
                    break;
                }
            }
        	// Adds data of side points of all neighbour objects
        	for( DHCellSide ngh_side : cell.neighb_sides() ) { // cell -> elm lower dim, ngh_side -> elm higher dim
                coupling_integral->patch_data().emplace_back(cell, ngh_side);
                add_side_points(coupling_integral, ngh_side, ppv_high);
            }
            ++i_int;
    	}
    }

    /**
     * Common part of add_edge integrals and add_coupling_integrals methods
     *
     * Method is used internally in AssemblyBase
     */
    template <template <unsigned int> class IntegralAcc>
    inline void add_side_points(std::shared_ptr< IntegralAcc<dim> > &integral, DHCellSide cell_side, PatchPointValues<3> &ppv) {
        ++ppv.n_mesh_items_;
        unsigned int reg_idx = cell_side.element().region_idx().idx();
        for (auto p : integral->points(cell_side) ) {
            patch_internals_->element_cache_map_.add_eval_point(reg_idx, cell_side.elem_idx(), p.eval_point_idx(), cell_side.cell().local_idx());
        }
    }

    Quadrature *quad_;                                     ///< Quadrature used in assembling methods.
    Quadrature *quad_low_;                                 ///< Quadrature used in assembling methods (dim-1).
    DimIntegrals<dim> integrals_;                          ///< Set of used integrals.
    PatchInternals *patch_internals_;                     ///< Holds shared internals data with GeneriAssembly

    /**
     * Minimal number of sides on edge.
     *
     * Edge integral is created and calculated if number of sides is greater or equal than this value. Default value
     * is 2 and can be changed
     */
    unsigned int min_edge_sides_;
};

#endif /* ASSEMBLY_BASE_HH_ */
