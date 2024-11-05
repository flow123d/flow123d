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
    typedef typename GenericAssemblyBase::BulkIntegralData BulkIntegralData;
    typedef typename GenericAssemblyBase::EdgeIntegralData EdgeIntegralData;
    typedef typename GenericAssemblyBase::CouplingIntegralData CouplingIntegralData;
    typedef typename GenericAssemblyBase::BoundaryIntegralData BoundaryIntegralData;

	/// Constructor
	AssemblyBase(unsigned int quad_order) {
        quad_ = new QGauss(dim, 2*quad_order);
        quad_low_ = new QGauss(dim-1, 2*quad_order);
	}

	/// Destructor
    virtual ~AssemblyBase() {
        delete quad_;
        delete quad_low_;
    }

    /// Assembles the volume integrals on cell.
    virtual inline void cell_integral(FMT_UNUSED DHCellAccessor cell, FMT_UNUSED unsigned int element_patch_idx) {}

    /// Assembles the fluxes on the boundary.
    virtual inline void boundary_side_integral(FMT_UNUSED DHCellSide cell_side) {}

    /// Assembles the fluxes between sides on the edge.
    virtual inline void edge_integral(FMT_UNUSED RangeConvert<DHEdgeSide, DHCellSide> edge_side_range) {}

    /// Assembles the fluxes between elements of different dimensions.
    virtual inline void dimjoin_intergral(FMT_UNUSED DHCellAccessor cell_lower_dim, FMT_UNUSED DHCellSide neighb_side) {}

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
    	    ASSERT_PERMANENT_PTR(quad_).error("Data member 'quad_' must be initialized if you use bulk integral!\n");
    	    integrals_.bulk_ = eval_points->add_bulk<dim>(*quad_);
    		integrals.bulk_[dim-1] = integrals_.bulk_;
    	}
    	if (active_integrals_ & ActiveIntegrals::edge) {
    	    ASSERT_PERMANENT_PTR(quad_low_).error("Data member 'quad_low_' must be initialized if you use edge integral!\n");
    	    integrals_.edge_ = eval_points->add_edge<dim>(*quad_low_);
    	    integrals.edge_[dim-1] = integrals_.edge_;
    	}
       	if ((dim>1) && (active_integrals_ & ActiveIntegrals::coupling)) {
    	    ASSERT_PERMANENT_PTR(quad_).error("Data member 'quad_' must be initialized if you use coupling integral!\n");
    	    ASSERT_PERMANENT_PTR(quad_low_).error("Data member 'quad_low_' must be initialized if you use coupling integral!\n");
    	    integrals_.coupling_ = eval_points->add_coupling<dim>(*quad_low_);
       	    integrals.coupling_[dim-2] = integrals_.coupling_;
       	}
       	if (active_integrals_ & ActiveIntegrals::boundary) {
    	    ASSERT_PERMANENT_PTR(quad_).error("Data member 'quad_' must be initialized if you use boundary integral!\n");
    	    ASSERT_PERMANENT_PTR(quad_low_).error("Data member 'quad_low_' must be initialized if you use boundary integral!\n");
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
        ASSERT( cell_side.dim() > 0 ).error("Invalid cell dimension, must be 1, 2 or 3!\n");
	    return integrals_.edge_->points(cell_side, element_cache_map_);
    }

    /// Return CouplingPoint range of appropriate dimension
    inline Range< CouplingPoint > coupling_points(const DHCellSide &cell_side) const {
        ASSERT( cell_side.dim() > 1 ).error("Invalid cell dimension, must be 2 or 3!\n");
	    return integrals_.coupling_->points(cell_side, element_cache_map_);
    }

    /// Return BoundaryPoint range of appropriate dimension
    inline Range< BoundaryPoint > boundary_points(const DHCellSide &cell_side) const {
        ASSERT( cell_side.dim() > 0 ).error("Invalid cell dimension, must be 1, 2 or 3!\n");
	    return integrals_.boundary_->points(cell_side, element_cache_map_);
    }

    /// Assembles the cell integrals for the given dimension.
    virtual inline void assemble_cell_integrals(const RevertableList<BulkIntegralData> &bulk_integral_data) {
    	for (unsigned int i=0; i<bulk_integral_data.permanent_size(); ++i) {
            if (bulk_integral_data[i].cell.dim() != dim) continue;
            this->cell_integral(bulk_integral_data[i].cell, element_cache_map_->position_in_cache(bulk_integral_data[i].cell.elm_idx()));
    	}
    	// Possibly optimization but not so fast as we would assume (needs change interface of cell_integral)
        /*for (unsigned int i=0; i<element_cache_map_->n_elements(); ++i) {
            unsigned int elm_start = element_cache_map_->element_chunk_begin(i);
            if (element_cache_map_->eval_point_data(elm_start).i_eval_point_ != 0) continue;
            this->cell_integral(i, element_cache_map_->eval_point_data(elm_start).dh_loc_idx_);
        }*/
    }

    /// Assembles the boundary side integrals for the given dimension.
    inline void assemble_boundary_side_integrals(const RevertableList<BoundaryIntegralData> &boundary_integral_data) {
        for (unsigned int i=0; i<boundary_integral_data.permanent_size(); ++i) {
            if (boundary_integral_data[i].side.dim() != dim) continue;
            this->boundary_side_integral(boundary_integral_data[i].side);
        }
    }

    /// Assembles the edge integrals for the given dimension.
    inline void assemble_edge_integrals(const RevertableList<EdgeIntegralData> &edge_integral_data) {
        for (unsigned int i=0; i<edge_integral_data.permanent_size(); ++i) {
        	auto range = edge_integral_data[i].edge_side_range;
            if (range.begin()->dim() != dim) continue;
            this->edge_integral(edge_integral_data[i].edge_side_range);
        }
    }

    /// Assembles the neighbours integrals for the given dimension.
    inline void assemble_neighbour_integrals(const RevertableList<CouplingIntegralData> &coupling_integral_data) {
        for (unsigned int i=0; i<coupling_integral_data.permanent_size(); ++i) {
            if (coupling_integral_data[i].side.dim() != dim) continue;
            this->dimjoin_intergral(coupling_integral_data[i].cell, coupling_integral_data[i].side);
        }
    }

    /// Register cell points of volume integral
    virtual inline void add_patch_bulk_points(FMT_UNUSED const RevertableList<BulkIntegralData> &bulk_integral_data) {}

    /// Register side points of boundary side integral
    virtual inline void add_patch_bdr_side_points(FMT_UNUSED const RevertableList<BoundaryIntegralData> &boundary_integral_data) {}

    /// Register side points of edge integral
    virtual inline void add_patch_edge_points(FMT_UNUSED const RevertableList<EdgeIntegralData> &edge_integral_data) {
    }

    /// Register bulk and side points of coupling integral
    virtual inline void add_patch_coupling_integrals(FMT_UNUSED const RevertableList<CouplingIntegralData> &coupling_integral_data) {}

protected:
    /// Set of integral of given dimension necessary in assemblation
    struct DimIntegrals {
        std::shared_ptr<BulkIntegral> bulk_;               ///< Bulk integrals of elements
        std::shared_ptr<EdgeIntegral> edge_;               ///< Edge integrals between elements of same dimensions
        std::shared_ptr<CouplingIntegral> coupling_;       ///< Coupling integrals between elements of dimensions dim and dim-1
        std::shared_ptr<BoundaryIntegral> boundary_;       ///< Boundary integrals betwwen side and boundary element of dim-1
    };

	/**
	 * Default constructor.
	 *
	 * Be aware if you use this constructor. Quadrature objects must be initialized manually in descendant.
	 */
	AssemblyBase()
	: quad_(nullptr), quad_low_(nullptr) {}

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


template <unsigned int dim>
class AssemblyBasePatch : public AssemblyBase<dim>
{
public:
    typedef typename GenericAssemblyBase::BulkIntegralData BulkIntegralData;
    typedef typename GenericAssemblyBase::EdgeIntegralData EdgeIntegralData;
    typedef typename GenericAssemblyBase::CouplingIntegralData CouplingIntegralData;
    typedef typename GenericAssemblyBase::BoundaryIntegralData BoundaryIntegralData;

	AssemblyBasePatch(PatchFEValues<3> *fe_values)
	: AssemblyBase<dim>(), fe_values_(fe_values) {
	    this->quad_ = fe_values_->get_quadrature(dim, true); // bulk quadrature
	    this->quad_low_  = fe_values_->get_quadrature(dim, false); // side quadrature
	}

    /// Register cell points of volume integral
    inline void add_patch_bulk_points(const RevertableList<BulkIntegralData> &bulk_integral_data) override {
        for (unsigned int i=0; i<bulk_integral_data.permanent_size(); ++i) {
            if (bulk_integral_data[i].cell.dim() != dim) continue;
            uint element_patch_idx = this->element_cache_map_->position_in_cache(bulk_integral_data[i].cell.elm_idx());
            uint elm_pos = fe_values_->register_element(bulk_integral_data[i].cell, element_patch_idx);
            uint i_point = 0;
            for (auto p : this->bulk_points(element_patch_idx) ) {
                unsigned int value_cache_idx = p.elm_cache_map()->element_eval_point(p.elem_patch_idx(), p.eval_point_idx());
                fe_values_->register_bulk_point(bulk_integral_data[i].cell, elm_pos, value_cache_idx, i_point++);
            }
        }
    }

    /// Register side points of boundary side integral
    inline void add_patch_bdr_side_points(const RevertableList<BoundaryIntegralData> &boundary_integral_data) override {
        for (unsigned int i=0; i<boundary_integral_data.permanent_size(); ++i) {
            if (boundary_integral_data[i].side.dim() != dim) continue;
        	uint side_pos = fe_values_->register_side(boundary_integral_data[i].side);
            uint i_point = 0;
            for (auto p : this->boundary_points(boundary_integral_data[i].side) ) {
        	    unsigned int value_cache_idx = p.elm_cache_map()->element_eval_point(p.elem_patch_idx(), p.eval_point_idx());
                fe_values_->register_side_point(boundary_integral_data[i].side, side_pos, value_cache_idx, i_point++);
            }
        }
    }

    /// Register side points of edge integral
    inline void add_patch_edge_points(const RevertableList<EdgeIntegralData> &edge_integral_data) override {
        for (unsigned int i=0; i<edge_integral_data.permanent_size(); ++i) {
        	auto range = edge_integral_data[i].edge_side_range;
            if (range.begin()->dim() != dim) continue;
            for( DHCellSide edge_side : range )
            {
            	uint side_pos = fe_values_->register_side(edge_side);
                uint i_point = 0;
                for (auto p : this->edge_points(edge_side) ) {
                    unsigned int value_cache_idx = p.elm_cache_map()->element_eval_point(p.elem_patch_idx(), p.eval_point_idx());
                    fe_values_->register_side_point(edge_side, side_pos, value_cache_idx, i_point++);
                }
            }
        }
    }

    /// Register bulk and side points of coupling integral
    inline void add_patch_coupling_integrals(const RevertableList<CouplingIntegralData> &coupling_integral_data) override {
        uint side_pos, element_patch_idx, elm_pos=0;
        uint last_element_idx = -1;

        for (unsigned int i=0; i<coupling_integral_data.permanent_size(); ++i) {
            if (coupling_integral_data[i].side.dim() != dim) continue;
            side_pos = fe_values_->register_side(coupling_integral_data[i].side);
            if (coupling_integral_data[i].cell.elm_idx() != last_element_idx) {
                element_patch_idx = this->element_cache_map_->position_in_cache(coupling_integral_data[i].cell.elm_idx());
                elm_pos = fe_values_->register_element(coupling_integral_data[i].cell, element_patch_idx);
            }

            uint i_bulk_point = 0, i_side_point = 0;
            for (auto p_high : this->coupling_points(coupling_integral_data[i].side) )
            {
                unsigned int value_cache_idx = p_high.elm_cache_map()->element_eval_point(p_high.elem_patch_idx(), p_high.eval_point_idx());
                fe_values_->register_side_point(coupling_integral_data[i].side, side_pos, value_cache_idx, i_side_point++);
                if (coupling_integral_data[i].cell.elm_idx() != last_element_idx) {
                    auto p_low = p_high.lower_dim(coupling_integral_data[i].cell);
            	    value_cache_idx = p_low.elm_cache_map()->element_eval_point(p_low.elem_patch_idx(), p_low.eval_point_idx());
                    fe_values_->register_bulk_point(coupling_integral_data[i].cell, elm_pos, value_cache_idx, i_bulk_point++);
                }
            }
            last_element_idx = coupling_integral_data[i].cell.elm_idx();
        }
    }

    /// Return BulkValues object
    inline unsigned int n_dofs() {
        return fe_values_->template n_dofs<dim>();
    }

    /// Return BulkValues object
    inline BulkValues<dim> bulk_values() {
        return fe_values_->template bulk_values<dim>();
    }

    /// Return SideValues object
    inline SideValues<dim> side_values() {
        return fe_values_->template side_values<dim>();
    }

    /// Return JoinValues object
    inline JoinValues<dim> join_values() {
        return fe_values_->template join_values<dim>();
    }

protected:
    PatchFEValues<3> *fe_values_;                          ///< Common FEValues object over all dimensions
};


#endif /* ASSEMBLY_BASE_HH_ */
