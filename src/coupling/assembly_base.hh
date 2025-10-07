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
class AssemblyBase
{
public:
//    typedef typename GenericAssemblyBase::BulkIntegralData BulkIntegralData;
//    typedef typename GenericAssemblyBase::EdgeIntegralData EdgeIntegralData;
//    typedef typename GenericAssemblyBase::CouplingIntegralData CouplingIntegralData;
//    typedef typename GenericAssemblyBase::BoundaryIntegralData BoundaryIntegralData;

    /**
     * Constructor
     *
     * @param quad_order    Order of Quadrature objects.
     * @param asm_internals Holds shared data with GenericAssembly
     */
    AssemblyBase(unsigned int quad_order, AssemblyInternals *asm_internals)
    : AssemblyBase<dim>() {
    	asm_internals_ = asm_internals;
        quad_ = new QGauss(dim, 2*quad_order);
        quad_low_ = new QGauss(dim-1, 2*quad_order);
    }

	/// Destructor
    virtual ~AssemblyBase() {
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

//    /// Create integrals according to dim of assembly object
//    void create_integrals(std::shared_ptr<EvalPoints> eval_points, AssemblyIntegrals &integrals) {
//    	if (active_integrals_ & ActiveIntegrals::bulk) {
//    	    ASSERT_PERMANENT_PTR(quad_).error("Data member 'quad_' must be initialized if you use bulk integral!\n");
//    	    integrals_.bulk_ = eval_points->add_bulk<dim>(*quad_);
//    		integrals.bulk_[dim-1] = integrals_.bulk_;
//    	}
//    	if (active_integrals_ & ActiveIntegrals::edge) {
//    	    ASSERT_PERMANENT_PTR(quad_low_).error("Data member 'quad_low_' must be initialized if you use edge integral!\n");
//    	    integrals_.edge_ = eval_points->add_edge<dim>(*quad_low_);
//    	    integrals.edge_[dim-1] = integrals_.edge_;
//    	}
//       	if ((dim>1) && (active_integrals_ & ActiveIntegrals::coupling)) {
//    	    ASSERT_PERMANENT_PTR(quad_).error("Data member 'quad_' must be initialized if you use coupling integral!\n");
//    	    ASSERT_PERMANENT_PTR(quad_low_).error("Data member 'quad_low_' must be initialized if you use coupling integral!\n");
//    	    integrals_.coupling_ = eval_points->add_coupling<dim>(*quad_low_);
//       	    integrals.coupling_[dim-2] = integrals_.coupling_;
//       	}
//       	if (active_integrals_ & ActiveIntegrals::boundary) {
//    	    ASSERT_PERMANENT_PTR(quad_).error("Data member 'quad_' must be initialized if you use boundary integral!\n");
//    	    ASSERT_PERMANENT_PTR(quad_low_).error("Data member 'quad_low_' must be initialized if you use boundary integral!\n");
//    	    integrals_.boundary_ = eval_points->add_boundary<dim>(*quad_low_);
//       	    integrals.boundary_[dim-1] = integrals_.boundary_;
//       	}
//    }

    /// Temporary method set pointers of integrals to GenericAssembly - IN DEVELOPMENT, OBSOLETE METHOD
//    void post_integrals_set(AssemblyIntegrals &integrals) {
//        if (integrals_.bulk_ != nullptr) {
//            integrals.bulk_[dim-1] = integrals_.bulk_;
//        }
//        if (integrals_.edge_ != nullptr) {
//            integrals.edge_[dim-1] = integrals_.edge_;
//     	  }
//        if ((dim>1) && (integrals_.coupling_ != nullptr)) {
//            integrals.coupling_[dim-2] = integrals_.coupling_;
//        }
//        if (integrals_.boundary_ != nullptr) {
//            integrals.boundary_[dim-1] = integrals_.boundary_;
//        }
//    }

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
                std::make_shared<BulkIntegralAcc<dim>>(asm_internals_->eval_points_, quad, &asm_internals_->fe_values_, &asm_internals_->element_cache_map_)
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
                std::make_shared<EdgeIntegralAcc<dim>>(asm_internals_->eval_points_, quad, &asm_internals_->fe_values_, &asm_internals_->element_cache_map_)
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
                std::make_shared<CouplingIntegralAcc<dim>>(asm_internals_->eval_points_, quad, &asm_internals_->fe_values_, &asm_internals_->element_cache_map_)
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
                std::make_shared<BoundaryIntegralAcc<dim>>(asm_internals_->eval_points_, quad, &asm_internals_->fe_values_, &asm_internals_->element_cache_map_)
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
                ++( asm_internals_->fe_values_.ppv(bulk_domain, cell.dim()) ).n_mesh_items_;
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

        add_coupling_integrala(cell);

        if (asm_internals_->element_cache_map_.get_simd_rounded_size() > CacheMapElementNumber::get()) {
            integral_data_.bulk_.revert_temporary();
            integral_data_.edge_.revert_temporary();
            integral_data_.coupling_.revert_temporary();
            integral_data_.boundary_.revert_temporary();
            return true;
        } else {
            integral_data_.bulk_.make_permanent();
            integral_data_.edge_.make_permanent();
            integral_data_.coupling_.make_permanent();
            integral_data_.boundary_.make_permanent();
            return false;
        }
    }

//    /// Return BulkPoint range of appropriate dimension
//    /// Obsolete method - will be replaced by 'points(integral, mesh_item)'
//    inline Range< BulkPoint > bulk_points(unsigned int element_patch_idx) const {
//        return integrals_.bulk_->points(element_patch_idx, &this->asm_internals_->element_cache_map_);
//    }
//
//    /// Return EdgePoint range of appropriate dimension
//    /// Obsolete method - will be replaced by 'points(integral, mesh_item)'
//    inline Range< EdgePoint > edge_points(const DHCellSide &cell_side) const {
//        ASSERT( cell_side.dim() > 0 ).error("Invalid cell dimension, must be 1, 2 or 3!\n");
//	    return integrals_.edge_->points(cell_side, &this->asm_internals_->element_cache_map_);
//    }
//
//    /// Return CouplingPoint range of appropriate dimension
//    /// Obsolete method - will be replaced by 'points(integral, mesh_item)'
//    inline Range< CouplingPoint > coupling_points(const DHCellSide &cell_side) const {
//        ASSERT( cell_side.dim() > 1 ).error("Invalid cell dimension, must be 2 or 3!\n");
//	    return integrals_.coupling_->points(cell_side, &this->asm_internals_->element_cache_map_);
//    }
//
//    /// Return BoundaryPoint range of appropriate dimension
//    /// Obsolete method - will be replaced by 'points(integral, mesh_item)'
//    inline Range< BoundaryPoint > boundary_points(const DHCellSide &cell_side) const {
//        ASSERT( cell_side.dim() > 0 ).error("Invalid cell dimension, must be 1, 2 or 3!\n");
//	    return integrals_.boundary_->points(cell_side, &this->asm_internals_->element_cache_map_);
//    }

    /// Assembles the cell integrals for the given dimension.
    virtual inline void assemble_cell_integrals() {
    	for (unsigned int i=0; i<integral_data_.bulk_.permanent_size(); ++i) {
            this->cell_integral(integral_data_.bulk_[i].cell, asm_internals_->element_cache_map_.position_in_cache(integral_data_.bulk_[i].cell.elm_idx()));
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
        for (unsigned int i=0; i<integral_data_.boundary_.permanent_size(); ++i) {
            this->boundary_side_integral(integral_data_.boundary_[i].side);
        }
    }

    /// Assembles the edge integrals for the given dimension.
    inline void assemble_edge_integrals() {
        for (unsigned int i=0; i<integral_data_.edge_.permanent_size(); ++i) {
            this->edge_integral(integral_data_.edge_[i].edge_side_range);
        }
    }

    /// Assembles the neighbours integrals for the given dimension.
    inline void assemble_neighbour_integrals() {
        for (unsigned int i=0; i<integral_data_.coupling_.permanent_size(); ++i) {
            this->dimjoin_intergral(integral_data_.coupling_[i].cell, integral_data_.coupling_[i].side);
        }
    }

    /// Setter of min_edge_sides_
    void set_min_edge_sides(unsigned int val) {
        min_edge_sides_ = val;
    }

//    /// Register cell points of volume integral
//    virtual inline void add_patch_bulk_points(FMT_UNUSED const RevertableList<BulkIntegralData> &bulk_integral_data) {}
//
//    /// Register side points of boundary side integral
//    virtual inline void add_patch_bdr_side_points(FMT_UNUSED const RevertableList<BoundaryIntegralData> &boundary_integral_data) {}
//
//    /// Register side points of edge integral
//    virtual inline void add_patch_edge_points(FMT_UNUSED const RevertableList<EdgeIntegralData> &edge_integral_data) {
//    }
//
//    /// Register bulk and side points of coupling integral
//    virtual inline void add_patch_coupling_integrals(FMT_UNUSED const RevertableList<CouplingIntegralData> &coupling_integral_data) {}

    /**
     * Clean all integral data structures
     *
     * Method is called from GenericAssembly::assembly method.
     */
    void clean_integral_data() {
        integral_data_.bulk_.reset();
        integral_data_.edge_.reset();
        integral_data_.coupling_.reset();
        integral_data_.boundary_.reset();
    }

    /// Getter of integrals_
    const DimIntegrals<dim> &integrals() const {
    	return integrals_;
    }

    /// Getter of integral_data_
    const IntegralData &integral_data() const {
        return integral_data_;
    }

protected:
//    /// Set of integral of given dimension necessary in assemblation
//    struct DimIntegrals {
//    	DimIntegrals() : bulk_(nullptr), edge_(nullptr), coupling_(nullptr), boundary_(nullptr) {}
//
//        std::shared_ptr<BulkIntegral> bulk_;               ///< Bulk integrals of elements
//        std::shared_ptr<EdgeIntegral> edge_;               ///< Edge integrals between elements of same dimensions
//        std::shared_ptr<CouplingIntegral> coupling_;       ///< Coupling integrals between elements of dimensions dim and dim-1
//        std::shared_ptr<BoundaryIntegral> boundary_;       ///< Boundary integrals betwwen side and boundary element of dim-1
//    };

	/**
	 * Default constructor.
	 *
	 * Be aware if you use this constructor. Quadrature objects must be initialized manually in descendant.
	 */
	AssemblyBase()
	: quad_(nullptr), quad_low_(nullptr), asm_internals_(nullptr), min_edge_sides_(2) {}

    /**
     * Add data of volume integrals to appropriate data structure.
     *
     * Method is used internally in AssemblyBase
     */
    inline void add_volume_integrals(const DHCellAccessor &cell) {
        auto &ppv = asm_internals_->fe_values_.ppv(bulk_domain, cell.dim());
        for (auto integral_it : integrals_.bulk_) {
            uint subset_idx = integral_it.second->get_subset_idx();
            integral_data_.bulk_.emplace_back(cell, subset_idx);

            unsigned int reg_idx = cell.elm().region_idx().idx();
            // Different access than in other integrals: We can't use range method CellIntegral::points
            // because it passes element_patch_idx as argument that is not known during patch construction.
            for (uint i=uint( asm_internals_->eval_points_->subset_begin(dim, subset_idx) );
                      i<uint( asm_internals_->eval_points_->subset_end(dim, subset_idx) ); ++i) {
                asm_internals_->element_cache_map_.add_eval_point(reg_idx, cell.elm_idx(), i, cell.local_idx());
                ++ppv.n_points_;
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

        auto &ppv = asm_internals_->fe_values_.ppv(side_domain, cell_side.dim());
        for (auto integral_it : integrals_.edge_) {
            integral_data_.edge_.emplace_back(range, integral_it.second->get_subset_idx());

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
        auto &ppv = asm_internals_->fe_values_.ppv(side_domain, bdr_side.dim());

        for (auto integral_it : integrals_.boundary_) {
            auto integral = integral_it.second;
            integral_data_.boundary_.emplace_back(integral->get_subset_low_idx(), bdr_side,
                    integral->get_subset_high_idx());

            unsigned int reg_idx = bdr_side.element().region_idx().idx();
            ++ppv.n_mesh_items_;
            for (auto p : integral->points(bdr_side) ) {
                asm_internals_->element_cache_map_.add_eval_point(reg_idx, bdr_side.elem_idx(), p.eval_point_idx(), bdr_side.cell().local_idx());
                ++ppv.n_points_;

            	BulkPoint p_bdr = p.point_bdr(bdr_side.cond().element_accessor()); // equivalent point on boundary element
            	unsigned int bdr_reg = bdr_side.cond().element_accessor().region_idx().idx();
            	// invalid local_idx value, DHCellAccessor of boundary element doesn't exist
            	asm_internals_->element_cache_map_.add_eval_point(bdr_reg, bdr_side.cond().bc_ele_idx(), p_bdr.eval_point_idx(), -1);
            }
        }
    }

    /**
     * Add data of coupling integrals to appropriate data structure.
     *
     * Method is used internally in AssemblyBase
     */
    inline void add_coupling_integrala(const DHCellAccessor &cell) {
        if (dim==3) return;
        auto &ppv_low = asm_internals_->fe_values_.ppv(bulk_domain, cell.dim());
        auto &ppv_high = asm_internals_->fe_values_.ppv(side_domain, cell.dim()+1);
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
                        asm_internals_->element_cache_map_.add_eval_point(reg_idx_low, cell.elm_idx(), p_low.eval_point_idx(), cell.local_idx());
                        ++ppv_low.n_points_;
                    }
                    break;
                }
            }
        	// Adds data of side points of all neighbour objects
        	for( DHCellSide ngh_side : cell.neighb_sides() ) { // cell -> elm lower dim, ngh_side -> elm higher dim
                integral_data_.coupling_.emplace_back(cell, coupling_integral->get_subset_low_idx(), ngh_side,
                        coupling_integral->get_subset_high_idx());
                add_side_points(coupling_integral, ngh_side, ppv_high);
            }
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
            asm_internals_->element_cache_map_.add_eval_point(reg_idx, cell_side.elem_idx(), p.eval_point_idx(), cell_side.cell().local_idx());
            ++ppv.n_points_;
        }
    }

    /// Print update flags to string format.
    std::string print_update_flags(UpdateFlags u) const {
        std::stringstream s;
        s << u;
        return s.str();
    }

    Quadrature *quad_;                                     ///< Quadrature used in assembling methods.
    Quadrature *quad_low_;                                 ///< Quadrature used in assembling methods (dim-1).
    DimIntegrals<dim> integrals_;                          ///< Set of used integrals.
    AssemblyInternals *asm_internals_;                     ///< Holds shared internals data with GeneriAssembly

    /**
     * Minimal number of sides on edge.
     *
     * Edge integral is created and calculated if number of sides is greater or equal than this value. Default value
     * is 2 and can be changed
     */
    unsigned int min_edge_sides_;

    IntegralData integral_data_;                           ///< Holds patch data for computing different types of integrals.
};


template <unsigned int dim>
class AssemblyBasePatch : public AssemblyBase<dim>
{
public:
//    typedef typename GenericAssemblyBase::BulkIntegralData BulkIntegralData;
//    typedef typename GenericAssemblyBase::EdgeIntegralData EdgeIntegralData;
//    typedef typename GenericAssemblyBase::CouplingIntegralData CouplingIntegralData;
//    typedef typename GenericAssemblyBase::BoundaryIntegralData BoundaryIntegralData;

    /**
     * Constructor.
     *
     * @param quad_order    Specification of Quadrature size.
     * @param asm_internals Holds shared data with GenericAssembly
     */
	AssemblyBasePatch(unsigned int quad_order, AssemblyInternals *asm_internals)
	: AssemblyBase<dim>(quad_order, asm_internals) {}

//    /// Register cell points of volume integral
//    inline void add_patch_bulk_points(const RevertableList<BulkIntegralData> &bulk_integral_data) override {
//        for (unsigned int i=0; i<bulk_integral_data.permanent_size(); ++i) {
//            if (bulk_integral_data[i].cell.dim() != dim) continue;
//            uint element_patch_idx = this->asm_internals_->element_cache_map_.position_in_cache(bulk_integral_data[i].cell.elm_idx());
//            uint elm_pos = this->asm_internals_->fe_values_.register_element(bulk_integral_data[i].cell, element_patch_idx);
//            uint i_point = 0;
//            for (auto p : this->bulk_points(element_patch_idx) ) {
//                this->asm_internals_->fe_values_.register_bulk_point(bulk_integral_data[i].cell, elm_pos, p.value_cache_idx(), i_point++);
//            }
//        }
//    }
//
//    /// Register side points of boundary side integral
//    inline void add_patch_bdr_side_points(const RevertableList<BoundaryIntegralData> &boundary_integral_data) override {
//        for (unsigned int i=0; i<boundary_integral_data.permanent_size(); ++i) {
//            if (boundary_integral_data[i].side.dim() != dim) continue;
//            uint element_patch_idx = this->asm_internals_->element_cache_map_.position_in_cache(boundary_integral_data[i].side.elem_idx());
//            uint side_pos = this->asm_internals_->fe_values_.register_side(boundary_integral_data[i].side, element_patch_idx);
//            uint i_point = 0;
//            for (auto p : this->boundary_points(boundary_integral_data[i].side) ) {
//                this->asm_internals_->fe_values_.register_side_point(boundary_integral_data[i].side, side_pos, p.value_cache_idx(), i_point++);
//            }
//        }
//    }
//
//    /// Register side points of edge integral
//    inline void add_patch_edge_points(const RevertableList<EdgeIntegralData> &edge_integral_data) override {
//        for (unsigned int i=0; i<edge_integral_data.permanent_size(); ++i) {
//        	auto range = edge_integral_data[i].edge_side_range;
//            if (range.begin()->dim() != dim) continue;
//            for( DHCellSide edge_side : range )
//            {
//                uint element_patch_idx = this->asm_internals_->element_cache_map_.position_in_cache(edge_side.elem_idx());
//                uint side_pos = this->asm_internals_->fe_values_.register_side(edge_side, element_patch_idx);
//                uint i_point = 0;
//                for (auto p : this->edge_points(edge_side) ) {
//                    this->asm_internals_->fe_values_.register_side_point(edge_side, side_pos, p.value_cache_idx(), i_point++);
//                }
//            }
//        }
//    }
//
//    /// Register bulk and side points of coupling integral
//    inline void add_patch_coupling_integrals(const RevertableList<CouplingIntegralData> &coupling_integral_data) override {
//        uint element_patch_idx, elm_pos=0;
//        uint last_element_idx = -1;
//
//        for (unsigned int i=0; i<coupling_integral_data.permanent_size(); ++i) {
//            if (coupling_integral_data[i].side.dim() != dim) continue;
//            element_patch_idx = this->asm_internals_->element_cache_map_.position_in_cache(coupling_integral_data[i].side.elem_idx());
//            uint side_pos = this->asm_internals_->fe_values_.register_side(coupling_integral_data[i].side, element_patch_idx);
//            if (coupling_integral_data[i].cell.elm_idx() != last_element_idx) {
//                element_patch_idx = this->asm_internals_->element_cache_map_.position_in_cache(coupling_integral_data[i].cell.elm_idx());
//                elm_pos = this->asm_internals_->fe_values_.register_element(coupling_integral_data[i].cell, element_patch_idx);
//            }
//
//            uint i_bulk_point = 0, i_side_point = 0;
//            for (auto p_high : this->coupling_points(coupling_integral_data[i].side) )
//            {
//                this->asm_internals_->fe_values_.register_side_point(coupling_integral_data[i].side, side_pos, p_high.value_cache_idx(), i_side_point++);
//                if (coupling_integral_data[i].cell.elm_idx() != last_element_idx) {
//                    auto p_low = p_high.lower_dim(coupling_integral_data[i].cell);
//                    this->asm_internals_->fe_values_.register_bulk_point(coupling_integral_data[i].cell, elm_pos, p_low.value_cache_idx(), i_bulk_point++);
//                }
//            }
//            last_element_idx = coupling_integral_data[i].cell.elm_idx();
//        }
//    }

    /// Return BulkValues object
    inline unsigned int n_dofs() {
        return this->asm_internals_->fe_values_.template n_dofs<dim>();
    }

    /// Return number of DOFs of higher dim element
    inline unsigned int n_dofs_high() {
        return this->asm_internals_->fe_values_.template n_dofs_high<dim>();
    }

//    /// Return BulkValues object
//    inline BulkValues<dim> bulk_values() {
//        return this->asm_internals_->fe_values_.template bulk_values<dim>();
//    }
//
//    /// Return SideValues object
//    inline SideValues<dim> side_values() {
//        return this->asm_internals_->fe_values_.template side_values<dim>();
//    }
//
//    /// Return JoinValues object
//    inline JoinValues<dim> join_values() {
//        return this->asm_internals_->fe_values_.template join_values<dim>();
//    }

};


#endif /* ASSEMBLY_BASE_HH_ */
