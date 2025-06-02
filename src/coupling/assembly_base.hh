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
//#include "fem/op_factory.hh"



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
    AssemblyBase(unsigned int quad_order)
    : min_edge_sides_(2),
	  bulk_integral_data_(20, 10),
      edge_integral_data_(12, 6),
      coupling_integral_data_(12, 6),
      boundary_integral_data_(8, 4) {
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

    /// Set shared_ptr to EvalPoints and create integral accessors
    void create_integrals(std::shared_ptr<EvalPoints> eval_points) {
        eval_points_ = eval_points;
        make_integrals();
    }

    /**
     * Add data of integrals to appropriate structure and register elements to ElementCacheMap.
     *
     * Types of used integrals must be set in data member \p active_integrals_.
     * Return true if patch is full to its maximal capacity.
     */
    bool add_integrals_of_computing_step(DHCellAccessor cell, PatchFEValues<3>::TableSizes &table_sizes_tmp) {
        if (integrals_.bulk_.size() > 0)
            if (cell.is_own()) { // Not ghost
                this->add_volume_integral(cell, table_sizes_tmp);
    	    }

        for( DHCellSide cell_side : cell.side_range() ) {
            if (integrals_.boundary_.size() > 0)
                if (cell.is_own()) // Not ghost
                    if ( (cell_side.side().edge().n_sides() == 1) && (cell_side.side().is_boundary()) ) {
                        this->add_boundary_integral(cell_side, table_sizes_tmp);
                        continue;
                    }
            if (integrals_.edge_.size() > 0)
                if ( (cell_side.n_edge_sides() >= min_edge_sides_) && (cell_side.edge_sides().begin()->element().idx() == cell.elm_idx())) {
                    this->add_edge_integral(cell_side, table_sizes_tmp);
                }
        }

        if (integrals_.coupling_.size() > 0) {
        	for (auto coupling_integral : integrals_.coupling_) {
                // Adds data of bulk points only if bulk point were not added during processing of bulk integral
                bool add_bulk_points = !( (integrals_.bulk_.size() > 0) & cell.is_own() );
                if (add_bulk_points) {
                    // add points of low dim element only one time and only if they have not been added in BulkIntegral
                    for( DHCellSide ngh_side : cell.neighb_sides() ) {
                        unsigned int reg_idx_low = cell.elm().region_idx().idx();
                        table_sizes_tmp.elem_sizes_[0][cell.dim()-1]++;
                        for (auto p : coupling_integral->points(ngh_side, element_cache_map_) ) {
                            auto p_low = p.lower_dim(cell); // equivalent point on low dim cell
                            element_cache_map_->add_eval_point(reg_idx_low, cell.elm_idx(), p_low.eval_point_idx(), cell.local_idx());
                            table_sizes_tmp.point_sizes_[0][cell.dim()-1]++;
                        }
                        break;
                    }
                }
            	// Adds data of side points of all neighbour objects
            	for( DHCellSide ngh_side : cell.neighb_sides() ) { // cell -> elm lower dim, ngh_side -> elm higher dim
                    coupling_integral_data_.emplace_back(cell, integrals_.coupling_[cell.dim()-1]->get_subset_low_idx(), ngh_side,
                            integrals_.coupling_[cell.dim()-1]->get_subset_high_idx());
                    table_sizes_tmp.elem_sizes_[1][cell.dim()]++;

                    unsigned int reg_idx_high = ngh_side.element().region_idx().idx();
                    for (auto p : coupling_integral->points(ngh_side, element_cache_map_) ) {
                        element_cache_map_->add_eval_point(reg_idx_high, ngh_side.elem_idx(), p.eval_point_idx(), ngh_side.cell().local_idx());
                        table_sizes_tmp.point_sizes_[1][cell.dim()]++;
                    }
                }
        	}
        }

        if (element_cache_map_->get_simd_rounded_size() > CacheMapElementNumber::get()) {
            bulk_integral_data_.revert_temporary();
            edge_integral_data_.revert_temporary();
            coupling_integral_data_.revert_temporary();
            boundary_integral_data_.revert_temporary();
            return true;
        } else {
            bulk_integral_data_.make_permanent();
            edge_integral_data_.make_permanent();
            coupling_integral_data_.make_permanent();
            boundary_integral_data_.make_permanent();
            return false;
        }

    }

    /// Add data of volume integral to appropriate data structure.
    inline void add_volume_integral(const DHCellAccessor &cell, PatchFEValues<3>::TableSizes &table_sizes_tmp) {
        ASSERT_EQ(cell.dim(), dim);

        for (auto integral_it : integrals_.bulk_) {
            uint subset_idx = integral_it->get_subset_idx();
            bulk_integral_data_.emplace_back(cell, subset_idx);

            unsigned int reg_idx = cell.elm().region_idx().idx();
            table_sizes_tmp.elem_sizes_[0][dim-1]++;
            // Different access than in other integrals: We can't use range method CellIntegral::points
            // because it passes element_patch_idx as argument that is not known during patch construction.
            for (uint i=uint( eval_points_->subset_begin(dim, subset_idx) );
                      i<uint( eval_points_->subset_end(dim, subset_idx) ); ++i) {
                element_cache_map_->add_eval_point(reg_idx, cell.elm_idx(), i, cell.local_idx());
                table_sizes_tmp.point_sizes_[0][dim-1]++;
            }
        }
    }

    /// Add data of edge integral to appropriate data structure.
    inline void add_edge_integral(const DHCellSide &cell_side, PatchFEValues<3>::TableSizes &table_sizes_tmp) {
	    auto range = cell_side.edge_sides();
        ASSERT_EQ(range.begin()->dim(), dim);

        for (auto integral_it : integrals_.edge_) {
            edge_integral_data_.emplace_back(range, integral_it->get_subset_idx());

            for( DHCellSide edge_side : range ) {
                unsigned int reg_idx = edge_side.element().region_idx().idx();
                table_sizes_tmp.elem_sizes_[1][dim-1]++;
                for (auto p : integral_it->points(edge_side, element_cache_map_) ) {
                    element_cache_map_->add_eval_point(reg_idx, edge_side.elem_idx(), p.eval_point_idx(), edge_side.cell().local_idx());
                    table_sizes_tmp.point_sizes_[1][dim-1]++;
                }
            }
        }
    }

    /// Add data of boundary integral to appropriate data structure.
    inline void add_boundary_integral(const DHCellSide &bdr_side, PatchFEValues<3>::TableSizes &table_sizes_tmp) {
        ASSERT_EQ(bdr_side.dim(), dim);

        for (auto integral_it : integrals_.boundary_) {
            boundary_integral_data_.emplace_back(integral_it->get_subset_low_idx(), bdr_side,
                    integral_it->get_subset_high_idx());

            unsigned int reg_idx = bdr_side.element().region_idx().idx();
            table_sizes_tmp.elem_sizes_[1][dim-1]++;
            for (auto p : integral_it->points(bdr_side, element_cache_map_) ) {
                element_cache_map_->add_eval_point(reg_idx, bdr_side.elem_idx(), p.eval_point_idx(), bdr_side.cell().local_idx());
                table_sizes_tmp.point_sizes_[1][dim-1]++;

            	BulkPoint p_bdr = p.point_bdr(bdr_side.cond().element_accessor()); // equivalent point on boundary element
            	unsigned int bdr_reg = bdr_side.cond().element_accessor().region_idx().idx();
            	// invalid local_idx value, DHCellAccessor of boundary element doesn't exist
            	element_cache_map_->add_eval_point(bdr_reg, bdr_side.cond().bc_ele_idx(), p_bdr.eval_point_idx(), -1);
            }
        }
    }

    /// Return BulkPoint range of appropriate dimension
    /// Obsolete method - must be removed
    inline Range< BulkPoint > bulk_points(unsigned int element_patch_idx) const {
        return integrals_.bulk_[0]->points(element_patch_idx, element_cache_map_);
    }

    /// Return EdgePoint range of appropriate dimension
    /// Obsolete method - must be removed
    inline Range< EdgePoint > edge_points(const DHCellSide &cell_side) const {
        ASSERT( cell_side.dim() > 0 ).error("Invalid cell dimension, must be 1, 2 or 3!\n");
	    return integrals_.edge_[0]->points(cell_side, element_cache_map_);
    }

    /// Return CouplingPoint range of appropriate dimension
    /// Obsolete method - must be removed
    inline Range< CouplingPoint > coupling_points(const DHCellSide &cell_side) const {
        ASSERT( cell_side.dim() > 1 ).error("Invalid cell dimension, must be 2 or 3!\n");
	    return integrals_.coupling_[0]->points(cell_side, element_cache_map_);
    }

    /// Return BoundaryPoint range of appropriate dimension
    /// Obsolete method - must be removed
    inline Range< BoundaryPoint > boundary_points(const DHCellSide &cell_side) const {
        ASSERT( cell_side.dim() > 0 ).error("Invalid cell dimension, must be 1, 2 or 3!\n");
	    return integrals_.boundary_[0]->points(cell_side, element_cache_map_);
    }

    /// Assembles the cell integrals for the given dimension.
    virtual inline void assemble_cell_integrals() {
    	for (unsigned int i=0; i<bulk_integral_data_.permanent_size(); ++i) {
            this->cell_integral(bulk_integral_data_[i].cell, element_cache_map_->position_in_cache(bulk_integral_data_[i].cell.elm_idx()));
    	}
    	// Possibly optimization but not so fast as we would assume (needs change interface of cell_integral)
        /*for (unsigned int i=0; i<element_cache_map_->n_elements(); ++i) {
            unsigned int elm_start = element_cache_map_->element_chunk_begin(i);
            if (element_cache_map_->eval_point_data(elm_start).i_eval_point_ != 0) continue;
            this->cell_integral(i, element_cache_map_->eval_point_data(elm_start).dh_loc_idx_);
        }*/
    }

    /// Assembles the boundary side integrals for the given dimension.
    inline void assemble_boundary_side_integrals() {
        for (unsigned int i=0; i<boundary_integral_data_.permanent_size(); ++i) {
            this->boundary_side_integral(boundary_integral_data_[i].side);
        }
    }

    /// Assembles the edge integrals for the given dimension.
    inline void assemble_edge_integrals() {
        for (unsigned int i=0; i<edge_integral_data_.permanent_size(); ++i) {
            this->edge_integral(edge_integral_data_[i].edge_side_range);
        }
    }

    /// Assembles the neighbours integrals for the given dimension.
    inline void assemble_neighbour_integrals() {
        for (unsigned int i=0; i<coupling_integral_data_.permanent_size(); ++i) {
            this->dimjoin_intergral(coupling_integral_data_[i].cell, coupling_integral_data_[i].side);
        }
    }

    /// Register cell points of volume integral
    virtual inline void add_patch_bulk_points() {}

    /// Register side points of boundary side integral
    virtual inline void add_patch_bdr_side_points() {}

    /// Register side points of edge integral
    virtual inline void add_patch_edge_points() {
    }

    /// Register bulk and side points of coupling integral
    virtual inline void add_patch_coupling_integrals() {}

    void set_min_edge_sides(unsigned int val) {
        min_edge_sides_ = val;
    }

    /// Clean all integral data structures
    void clean_integral_data() {
        bulk_integral_data_.reset();
        edge_integral_data_.reset();
        coupling_integral_data_.reset();
        boundary_integral_data_.reset();
    }

protected:
    /// Set of integral of given dimension necessary in assemblation
    struct DimIntegrals {
        std::vector< std::shared_ptr<BulkIntegral> > bulk_;          ///< Bulk integrals of elements
        std::vector< std::shared_ptr<EdgeIntegral> > edge_;          ///< Edge integrals between elements of same dimensions
        std::vector< std::shared_ptr<CouplingIntegral> > coupling_;  ///< Coupling integrals between elements of dimensions dim and dim-1
        std::vector< std::shared_ptr<BoundaryIntegral> > boundary_;  ///< Boundary integrals betwwen side and boundary element of dim-1
    };

    /**
     * Default constructor.
     *
     * Be aware if you use this constructor. Quadrature objects must be initialized manually in descendant.
     */
    AssemblyBase()
    : quad_(nullptr), quad_low_(nullptr),
      min_edge_sides_(2),
      bulk_integral_data_(20, 10),
      edge_integral_data_(12, 6),
      coupling_integral_data_(12, 6),
      boundary_integral_data_(8, 4) {}

    // Create integral accessors in descendants if accessors are needed
    virtual void make_integrals() {}

    /// Create and return BulkIntegral of given quadrature
    std::shared_ptr<BulkIntegral> create_bulk_integral(Quadrature *quad) {
        integrals_.bulk_.emplace_back( eval_points_->template add_bulk<dim>(*quad) );
        return integrals_.bulk_.back();
    }

    /// Create and return EdgeIntegral of given quadrature
    std::shared_ptr<EdgeIntegral> create_edge_integral(Quadrature *quad) {
        integrals_.edge_.emplace_back( eval_points_->template add_edge<dim>(*quad) );
        return integrals_.edge_.back();
    }

    /// Create and return CouplingIntegral of given quadrature
    std::shared_ptr<CouplingIntegral> create_coupling_integral(Quadrature *quad) {
        integrals_.coupling_.emplace_back( eval_points_->template add_coupling<dim>(*quad) );
        return integrals_.coupling_.back();
    }

    /// Create and return BoundaryIntegral of given quadrature
    std::shared_ptr<BoundaryIntegral> create_boundary_integral(Quadrature *quad) {
        integrals_.boundary_.emplace_back( eval_points_->template add_boundary<dim>(*quad) );
        return integrals_.boundary_.back();
    }

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
    std::shared_ptr<EvalPoints> eval_points_;              ///< EvalPoints shared with GenericAssembly object.

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
	    this->quad_ = fe_values_->get_bulk_quadrature(dim);
	    this->quad_low_  = fe_values_->get_side_quadrature(dim);
	}

    /// Register cell points of volume integral
    inline void add_patch_bulk_points() override {
        for (auto integral_it : this->integrals_.bulk_) {
            for (unsigned int i=0; i<this->bulk_integral_data_.permanent_size(); ++i) {
                if ( this->bulk_integral_data_[i].subset_index != (unsigned int)(integral_it->get_subset_idx()) ) continue;
                uint element_patch_idx = this->element_cache_map_->position_in_cache(this->bulk_integral_data_[i].cell.elm_idx());
                uint elm_pos = fe_values_->register_element(this->bulk_integral_data_[i].cell, element_patch_idx);
                uint i_point = 0;
                for (auto p : integral_it->points(element_patch_idx, this->element_cache_map_) ) {
                    fe_values_->register_bulk_point(this->bulk_integral_data_[i].cell, elm_pos, p.value_cache_idx(), i_point++);
                }
            }
        }
    }

    /// Register side points of boundary side integral
    inline void add_patch_bdr_side_points() override {
        for (auto integral_it : this->integrals_.boundary_) {
            for (unsigned int i=0; i<this->boundary_integral_data_.permanent_size(); ++i) {
                if ( this->boundary_integral_data_[i].bdr_subset_index != (unsigned int)(integral_it->get_subset_low_idx()) ) continue;
            	uint side_pos = fe_values_->register_side(this->boundary_integral_data_[i].side);
                uint i_point = 0;
                for (auto p : integral_it->points(this->boundary_integral_data_[i].side, this->element_cache_map_) ) {
                    fe_values_->register_side_point(this->boundary_integral_data_[i].side, side_pos, p.value_cache_idx(), i_point++);
                }
            }
        }
    }

    /// Register side points of edge integral
    inline void add_patch_edge_points() override {
        for (auto integral_it : this->integrals_.edge_) {
            for (unsigned int i=0; i<this->edge_integral_data_.permanent_size(); ++i) {
                if ( this->edge_integral_data_[i].subset_index != (unsigned int)(integral_it->get_subset_idx()) ) continue;
            	auto range = this->edge_integral_data_[i].edge_side_range;
                for( DHCellSide edge_side : range )
                {
                	uint side_pos = fe_values_->register_side(edge_side);
                    uint i_point = 0;
                    for (auto p : integral_it->points(edge_side, this->element_cache_map_) ) {
                        fe_values_->register_side_point(edge_side, side_pos, p.value_cache_idx(), i_point++);
                    }
                }
            }
        }
    }

    /// Register bulk and side points of coupling integral
    inline void add_patch_coupling_integrals() override {
        for (auto integral_it : this->integrals_.coupling_) {
            uint side_pos, element_patch_idx, elm_pos=0;
            uint last_element_idx = -1;

            for (unsigned int i=0; i<this->coupling_integral_data_.permanent_size(); ++i) {
                if ( this->coupling_integral_data_[i].bulk_subset_index != (unsigned int)(integral_it->get_subset_low_idx()) ) continue;
                side_pos = fe_values_->register_side(this->coupling_integral_data_[i].side);
                if (this->coupling_integral_data_[i].cell.elm_idx() != last_element_idx) {
                    element_patch_idx = this->element_cache_map_->position_in_cache(this->coupling_integral_data_[i].cell.elm_idx());
                    elm_pos = fe_values_->register_element(this->coupling_integral_data_[i].cell, element_patch_idx);
                }

                uint i_bulk_point = 0, i_side_point = 0;
                for (auto p_high : integral_it->points(this->coupling_integral_data_[i].side, this->element_cache_map_) )
                {
                    fe_values_->register_side_point(this->coupling_integral_data_[i].side, side_pos, p_high.value_cache_idx(), i_side_point++);
                    if (this->coupling_integral_data_[i].cell.elm_idx() != last_element_idx) {
                        auto p_low = p_high.lower_dim(this->coupling_integral_data_[i].cell);
                        fe_values_->register_bulk_point(this->coupling_integral_data_[i].cell, elm_pos, p_low.value_cache_idx(), i_bulk_point++);
                    }
                }
                last_element_idx = this->coupling_integral_data_[i].cell.elm_idx();
            }
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
