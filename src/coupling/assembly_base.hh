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
    /**
     * Constructor
     *
     * @param quad_order   Order of Quadrature objects.
     */
    AssemblyBase(unsigned int quad_order)
    : AssemblyBase<dim>() {
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

    /**
     * Set shared_ptr to EvalPoints and create integral accessors
     *
     * Method is called from GenericAssembly during initialization of assembly object
     */
    void set_eval_points(std::shared_ptr<EvalPoints> eval_points) {
        eval_points_ = eval_points;
    }

    /**
     * Add data of integrals to appropriate structure and register elements to ElementCacheMap.
     *
     * Return true if patch is full to its maximal capacity.
     * Method is called from GenericAssembly::assembly method.
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
                    integral_data_.coupling_.emplace_back(cell, coupling_integral->get_subset_low_idx(), ngh_side,
                            coupling_integral->get_subset_high_idx());
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

    /**
     * Create point range of different type of integral accessors.
     *
	 * @param integral  Integral whose range is created.
	 * @param mesh_item Parameter of range method specialized by type of integral:
	 *                  unsigned int - index of element on patch in case of BulkIntegral
	 *                  DHCellSide - accessor to cell side in case of other integral types
	 * @return          Point range of appropriate integral.
     */
    template <class QIntegral>
    Range< typename QIntegral::PointType > points(std::shared_ptr<QIntegral> integral, typename QIntegral::MeshItem mesh_item) const {
    	return integral->points(mesh_item, element_cache_map_);
    }

    /**
     * Assembles the cell integrals for the given dimension.
     *
     * Method is called from GenericAssembly::assembly method.
     */
    virtual inline void assemble_cell_integrals() {
    	for (unsigned int i=0; i<integral_data_.bulk_.permanent_size(); ++i) {
            this->cell_integral(integral_data_.bulk_[i].cell, element_cache_map_->position_in_cache(integral_data_.bulk_[i].cell.elm_idx()));
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
        for (unsigned int i=0; i<integral_data_.boundary_.permanent_size(); ++i) {
            this->boundary_side_integral(integral_data_.boundary_[i].side);
        }
    }

    /**
     * Assembles the edge integrals for the given dimension.
     *
     * Method is called from GenericAssembly::assembly method.
     */
    inline void assemble_edge_integrals() {
        for (unsigned int i=0; i<integral_data_.edge_.permanent_size(); ++i) {
            this->edge_integral(integral_data_.edge_[i].edge_side_range);
        }
    }

    /**
     * Assembles the neighbours integrals for the given dimension.
     *
     * Method is called from GenericAssembly::assembly method.
     */
    inline void assemble_neighbour_integrals() {
        for (unsigned int i=0; i<integral_data_.coupling_.permanent_size(); ++i) {
            this->dimjoin_intergral(integral_data_.coupling_[i].cell, integral_data_.coupling_[i].side);
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
        integral_data_.bulk_.reset();
        integral_data_.edge_.reset();
        integral_data_.coupling_.reset();
        integral_data_.boundary_.reset();
    }

    /// Getter of integrals_
    const DimIntegrals &integrals() const {
    	return integrals_;
    }

    /// Getter of integral_data_
    const IntegralData &integral_data() const {
    	return integral_data_;
    }

protected:
    /**
     * Default constructor.
     *
     * Be aware if you use this constructor. Quadrature objects must be initialized manually in descendant.
     */
    AssemblyBase()
    : quad_(nullptr), quad_low_(nullptr),
      min_edge_sides_(2) {
        this->set_integral_data_lists();
    }

    /**
     * Allow to set base parameters (reserved_size and enlarged_by) of RevertibleLists
     * in \p integral_data_ in descendant classes.
     *
     * Example of definition in descendant class:
     @code
     void set_integral_data_lists() override {
         this->integral_data_.bulk_.reinit_default_list( 20 * (dim+1), dim+1 );
         this->integral_data_.edge_.reinit_default_list( 20 * dim+, dim );
         // in case of need add call of 'reinit_default_list' on other lists of 'integral_data_' struct
     }
     @endcode
     */
    virtual void set_integral_data_lists() {}

    /**
     * Add data of volume integral to appropriate data structure.
     *
     * Method is used internally in AssemblyBase
     */
    inline void add_volume_integral(const DHCellAccessor &cell, PatchFEValues<3>::TableSizes &table_sizes_tmp) {
        ASSERT_EQ(cell.dim(), dim);

        table_sizes_tmp.elem_sizes_[0][dim-1]++;
        for (auto integral_it : integrals_.bulk_) {
            uint subset_idx = integral_it->get_subset_idx();
            integral_data_.bulk_.emplace_back(cell, subset_idx);

            unsigned int reg_idx = cell.elm().region_idx().idx();
            // Different access than in other integrals: We can't use range method CellIntegral::points
            // because it passes element_patch_idx as argument that is not known during patch construction.
            for (uint i=uint( eval_points_->subset_begin(dim, subset_idx) );
                      i<uint( eval_points_->subset_end(dim, subset_idx) ); ++i) {
                element_cache_map_->add_eval_point(reg_idx, cell.elm_idx(), i, cell.local_idx());
                table_sizes_tmp.point_sizes_[0][dim-1]++;
            }
        }
    }

    /**
     * Add data of edge integral to appropriate data structure.
     *
     * Method is used internally in AssemblyBase
     */
    inline void add_edge_integral(const DHCellSide &cell_side, PatchFEValues<3>::TableSizes &table_sizes_tmp) {
	    auto range = cell_side.edge_sides();
        ASSERT_EQ(range.begin()->dim(), dim);

        for (auto integral_it : integrals_.edge_) {
            integral_data_.edge_.emplace_back(range, integral_it->get_subset_idx());

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

    /**
     * Add data of boundary integral to appropriate data structure.
     *
     * Method is used internally in AssemblyBase
     */
    inline void add_boundary_integral(const DHCellSide &bdr_side, PatchFEValues<3>::TableSizes &table_sizes_tmp) {
        ASSERT_EQ(bdr_side.dim(), dim);

        for (auto integral_it : integrals_.boundary_) {
            integral_data_.boundary_.emplace_back(integral_it->get_subset_low_idx(), bdr_side,
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

    /**
     * Create and return BulkIntegral of given quadrature.
     *
     * Method is called from descendants during construction / initialization of assembly object.
     */
    std::shared_ptr<BulkIntegral> create_bulk_integral(Quadrature *quad) {
        auto result = integrals_.bulk_.insert( std::make_shared<BulkIntegral>(quad, quad->dim()) );
	    return *result.first;
    }

    /**
     * Create and return EdgeIntegral of given quadrature.
     *
     * Method is called from descendants during construction / initialization of assembly object.
     */
    std::shared_ptr<EdgeIntegral> create_edge_integral(Quadrature *quad) {
        auto result = integrals_.edge_.insert( std::make_shared<EdgeIntegral>(quad, quad->dim()+1) );
	    return *result.first;
    }

    /**
     * Create and return CouplingIntegral of given quadrature.
     *
     * Method is called from descendants during construction / initialization of assembly object.
     */
    std::shared_ptr<CouplingIntegral> create_coupling_integral(Quadrature *quad) {
        if (dim==3) return nullptr;

        auto result = integrals_.coupling_.insert( std::make_shared<CouplingIntegral>(quad, quad->dim()) );
	    return *result.first;
    }

    /**
     * Create and return BoundaryIntegral of given quadrature.
     *
     * Method is called from descendants during construction / initialization of assembly object.
     */
    std::shared_ptr<BoundaryIntegral> create_boundary_integral(Quadrature *quad) {
        auto result = integrals_.boundary_.insert( std::make_shared<BoundaryIntegral>(quad, quad->dim()+1) );
	    return *result.first;
    }

    /// Print update flags to string format.
    std::string print_update_flags(UpdateFlags u) const {
        std::stringstream s;
        s << u;
        return s.str();
    }

    Quadrature *quad_;                                     ///< Quadrature used in assembling methods.
    Quadrature *quad_low_;                                 ///< Quadrature used in assembling methods (dim-1).
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

    /// Holds data for computing different types of integrals.
    IntegralData integral_data_;

};


template <unsigned int dim>
class AssemblyBasePatch : public AssemblyBase<dim>
{
public:
    /**
     * Constructor.
     *
     * @param fe_values   PatchFeValues object shared between assembly objects of different dimensions.
     */
	AssemblyBasePatch(PatchFEValues<3> *fe_values)
	: AssemblyBase<dim>(), fe_values_(fe_values) {
	    this->quad_ = fe_values_->get_bulk_quadrature(dim);
	    this->quad_low_  = fe_values_->get_side_quadrature(dim);
	}

    /// Return number of DOFs
    inline unsigned int n_dofs() {
        return fe_values_->template n_dofs<dim>();
    }

    /// Return number of DOFs of higher dim element
    inline unsigned int n_dofs_high() {
        return fe_values_->template n_dofs_high<dim>();
    }

    /// Return BulkValues object
    inline BulkValues<dim> bulk_values(const Quadrature *quad) {
        ASSERT_EQ(quad->dim(), dim);
	    return fe_values_->template bulk_values<dim>(quad);
    }

    /// Same as previous but pass default Quadrature argument
    inline BulkValues<dim> bulk_values() {
        return fe_values_->template bulk_values<dim>(this->quad_);
    }

    /// Return SideValues object
    inline SideValues<dim> side_values(const Quadrature *quad) {
        ASSERT_EQ(quad->dim(), dim-1);
        return fe_values_->template side_values<dim>(quad);
    }

    /// Same as previous but pass default Quadrature argument
    inline SideValues<dim> side_values() {
        return fe_values_->template side_values<dim>(this->quad_low_);
    }

    /// Return SideValues object
    inline SideValues<dim+1> side_values_high_dim(const Quadrature *quad) {
        ASSERT_EQ(quad->dim(), dim);
        return fe_values_->template side_values<dim+1>(quad);
    }

    /// Same as previous but pass default Quadrature argument
    inline SideValues<dim+1> side_values_high_dim() {
        return fe_values_->template side_values<dim+1>(this->quad_);
    }

    /// Return JoinValues object
    inline JoinValues<dim> join_values(const Quadrature *quad, const Quadrature *quad_low) {
        return fe_values_->template join_values<dim>(quad, quad_low);
    }

    /// Same as previous but pass default Quadrature arguments
    inline JoinValues<dim> join_values() {
        return fe_values_->template join_values<dim>(this->quad_, this->quad_low_);
    }

protected:
    PatchFEValues<3> *fe_values_;                          ///< Common FEValues object over all dimensions
};


#endif /* ASSEMBLY_BASE_HH_ */
