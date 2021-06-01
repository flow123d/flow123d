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
#include "fields/eval_points.hh"
#include "fields/field_value_cache.hh"
#include "fem/update_flags.hh"



/**
 * Base class define empty methods, these methods can be overwite in descendants.
 */
template <unsigned int dim>
class AssemblyBase
{
public:
	/// Constructor
	AssemblyBase(unsigned int quad_order) {
        quad_ = new QGauss(dim, 2*quad_order);
        quad_low_ = new QGauss(dim-1, 2*quad_order);
	}

	// Destructor
    virtual ~AssemblyBase() {
        delete quad_;
        delete quad_low_;
    }

    /// Assembles the volume integrals on cell.
    inline void cell_integral(FMT_UNUSED DHCellAccessor cell, FMT_UNUSED unsigned int element_patch_idx) {}

    /// Assembles the fluxes on the boundary.
    inline void boundary_side_integral(FMT_UNUSED DHCellSide cell_side) {}

    /// Assembles the fluxes between sides on the edge.
    inline void edge_integral(FMT_UNUSED RangeConvert<DHEdgeSide, DHCellSide> edge_side_range) {}

    /// Assembles the fluxes between elements of different dimensions.
    inline void neigbour_integral(FMT_UNUSED DHCellAccessor cell_lower_dim, FMT_UNUSED DHCellSide neighb_side) {}

    /// Method prepares object before assemblation (e.g. balance, ...).
    virtual void begin() {}

    /// Method finishes object after assemblation (e.g. balance, ...).
    virtual void end() {}

    /// Getter of active_integrals.
    inline int n_active_integrals() const {
        return active_integrals_;
    }

    /// Create integrals according to dim of assembly object
    void create_integrals(std::shared_ptr<EvalPoints> eval_points, AssemblyIntegrals &integrals) {
    	if (active_integrals_ & ActiveIntegrals::bulk) {
    	    integrals_.bulk_ = eval_points->add_bulk<dim>(*quad_);
    		integrals.bulk_[dim-1] = integrals_.bulk_;
    	}
    	if (active_integrals_ & ActiveIntegrals::edge) {
    	    integrals_.edge_ = eval_points->add_edge<dim>(*quad_low_);
    	    integrals.edge_[dim-1] = integrals_.edge_;
    	}
       	if ((dim>1) && (active_integrals_ & ActiveIntegrals::coupling)) {
    	    integrals_.coupling_ = eval_points->add_coupling<dim>(*quad_low_);
       	    integrals.coupling_[dim-2] = integrals_.coupling_;
       	}
       	if (active_integrals_ & ActiveIntegrals::boundary) {
    	    integrals_.boundary_ = eval_points->add_boundary<dim>(*quad_low_);
       	    integrals.boundary_[dim-1] = integrals_.boundary_;
       	}
    }

    /// Return BulkPoint range of appropriate dimension
    inline Range< BulkPoint > bulk_points(unsigned int element_patch_idx) const {
        return integrals_.bulk_->points(element_patch_idx, element_cache_map_);
    }

    /// Return EdgePoint range of appropriate dimension
    inline Range< EdgePoint > edge_points(const DHCellSide &cell_side) const {
        ASSERT_DBG( cell_side.dim() > 0 ).error("Invalid cell dimension, must be 1, 2 or 3!\n");
	    return integrals_.edge_->points(cell_side, element_cache_map_);
    }

    /// Return CouplingPoint range of appropriate dimension
    inline Range< CouplingPoint > coupling_points(const DHCellSide &cell_side) const {
        ASSERT_DBG( cell_side.dim() > 1 ).error("Invalid cell dimension, must be 2 or 3!\n");
	    return integrals_.coupling_->points(cell_side, element_cache_map_);
    }

    /// Return BoundaryPoint range of appropriate dimension
    inline Range< BoundaryPoint > boundary_points(const DHCellSide &cell_side) const {
        ASSERT_DBG( cell_side.dim() > 0 ).error("Invalid cell dimension, must be 1, 2 or 3!\n");
	    return integrals_.boundary_->points(cell_side, element_cache_map_);
    }

protected:
    /// Set of integral of given dimension necessary in assemblation
    struct DimIntegrals {
        std::shared_ptr<BulkIntegral> bulk_;               ///< Bulk integrals of elements
        std::shared_ptr<EdgeIntegral> edge_;               ///< Edge integrals between elements of same dimensions
        std::shared_ptr<CouplingIntegral> coupling_;       ///< Coupling integrals between elements of dimensions dim and dim-1
        std::shared_ptr<BoundaryIntegral> boundary_;       ///< Boundary integrals betwwen side and boundary element of dim-1
    };

    /// Print update flags to string format.
    std::string print_update_flags(UpdateFlags u) const {
        std::stringstream s;
        s << u;
        return s.str();
    }

    Quadrature *quad_;                                     ///< Quadrature used in assembling methods.
    Quadrature *quad_low_;                                 ///< Quadrature used in assembling methods (dim-1).
    int active_integrals_;                                 ///< Holds mask of active integrals.
    DimIntegrals integrals_;                               ///< Set of used integrals.
    ElementCacheMap *element_cache_map_;                   ///< ElementCacheMap shared with GenericAssembly object.
};


#endif /* ASSEMBLY_BASE_HH_ */
