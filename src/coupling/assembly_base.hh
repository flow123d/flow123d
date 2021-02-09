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
    virtual void cell_integral(FMT_UNUSED DHCellAccessor cell) {}

    /// Assembles the fluxes on the boundary.
    virtual void boundary_side_integral(FMT_UNUSED DHCellSide cell_side, FMT_UNUSED const TimeStep &step) {}

    /// Assembles the fluxes between sides on the edge.
    virtual void edge_integral(FMT_UNUSED RangeConvert<DHEdgeSide, DHCellSide> edge_side_range) {}

    /// Assembles the fluxes between elements of different dimensions.
    virtual void neigbour_integral(FMT_UNUSED DHCellAccessor cell_lower_dim, FMT_UNUSED DHCellSide neighb_side) {}

    /// Method prepares object before assemblation (e.g. balance, ...).
    virtual void begin() {}

    /// Method finishes object after assemblation (e.g. balance, ...).
    virtual void end() {}

    /// Reallocate data caches of fields in equation.
    virtual void reallocate_cache(const ElementCacheMap &cache_map) =0;

    /// Create integrals according to dim of assembly object
    void create_integrals(std::shared_ptr<EvalPoints> eval_points, AssemblyIntegrals &integrals, int active_integrals) {
    	if (active_integrals & ActiveIntegrals::bulk)
    	    integrals.bulk_[dim-1] = eval_points->add_bulk<dim>(*quad_);
    	if (active_integrals & ActiveIntegrals::edge)
    	    integrals.edge_[dim-1] = eval_points->add_edge<dim>(*quad_low_);
       	if ((dim>1) && (active_integrals & ActiveIntegrals::coupling))
       	    integrals.coupling_[dim-2] = eval_points->add_coupling<dim>(*quad_low_);
       	if (active_integrals & ActiveIntegrals::boundary)
       	    integrals.boundary_[dim-1] = eval_points->add_boundary<dim>(*quad_low_);
    }

protected:
    Quadrature *quad_;                                     ///< Quadrature used in assembling methods.
    Quadrature *quad_low_;                                 ///< Quadrature used in assembling methods (dim-1).
};


#endif /* ASSEMBLY_BASE_HH_ */
