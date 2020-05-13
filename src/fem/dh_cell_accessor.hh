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
 * @file    dh_cell_accessor.hh
 * @brief
 * @author  David Flanderka
 */

#ifndef DH_CELL_ACCESSOR_HH_
#define DH_CELL_ACCESSOR_HH_

#include <armadillo>
#include "mesh/accessors.hh"
#include "mesh/neighbours.h"
#include "fem/finite_element.hh"
#include "fem/dofhandler.hh"
#include "system/index_types.hh"

class DHCellSide;
class DHNeighbSide;
class DHEdgeSide;
template <int spacedim> class ElementAccessor;

/**
 * @brief Cell accessor allow iterate over DOF handler cells.
 *
 * Iterating is possible over different ranges of local and ghost elements.
 *
 * Iterator is defined by:
 *  - DOF handler
 *  - local index of DOF cell (iterated value)
 */
class DHCellAccessor {
public:
    /// Default invalid accessor.
	DHCellAccessor()
    : dof_handler_(NULL)
    {}

    /**
     * DOF cell accessor.
     */
	DHCellAccessor(const DOFHandlerMultiDim *dof_handler, unsigned int loc_idx)
    : dof_handler_(dof_handler), loc_ele_idx_(loc_idx)
    {}

    /// Return local index to element (index of DOF handler).
    inline unsigned int local_idx() const {
        ASSERT_LT_DBG(loc_ele_idx_, dof_handler_->global_to_local_el_idx_.size()).error("Local element index is out of range!\n");
        return loc_ele_idx_;
    }

    /// Return serial idx to element of loc_ele_idx_.
    inline unsigned int elm_idx() const {
        unsigned int ds_lsize = dof_handler_->el_ds_->lsize();
        if (local_idx()<ds_lsize) return dof_handler_->mesh()->get_el_4_loc()[loc_ele_idx_]; //own elements
        else return dof_handler_->ghost_4_loc[loc_ele_idx_-ds_lsize]; //ghost elements
    }

    /// Return ElementAccessor to element of loc_ele_idx_.
    inline const ElementAccessor<3> elm() const {
    	return dof_handler_->mesh()->element_accessor( elm_idx() );
    }

    /**
     * @brief Fill vector of the global indices of dofs associated to the cell.
     *
     * @param indices Vector of dof indices on the cell.
     */
    unsigned int get_dof_indices(std::vector<LongIdx> &indices) const
    { return dof_handler_->get_dof_indices( *this, indices ); }

    /**
     * @brief Returns the local indices of dofs associated to the cell on the local process.
     *
     * @param indices Array of dof indices on the cell.
     */
    LocDofVec get_loc_dof_indices() const
    { return dof_handler_->get_loc_dof_indices(loc_ele_idx_); }

    /// Return number of dofs on given cell.
    unsigned int n_dofs() const;

    /**
     * @brief Return dof on a given cell.
     * @param idof Number of dof on the cell.
     */
    const Dof &cell_dof(unsigned int idof) const;

    /// Return dimension of element appropriate to cell.
    inline unsigned int dim() const {
    	return elm().dim();
    }

    /// Return DOF handler
    inline const DOFHandlerMultiDim *dh() const{
        return dof_handler_;
    }

    /**
     * @brief Returns finite element object for given space dimension.
     */
    template<unsigned int dim>
    FEPtr<dim> fe() const {
        ElementAccessor<3> elm_acc = this->elm();
        return dof_handler_->ds_->fe(elm_acc).get<dim>();
    }

    /// Check validity of accessor (see default constructor)
    inline bool is_valid() const {
        return dof_handler_ != NULL;
    }

    /// Getter of elm_cache_index_.
    inline unsigned int element_cache_index() const {
        return elm_cache_index_;
    }

    /// Setter of elm_cache_index_.
    inline void set_element_cache_index(unsigned int idx) const {
        elm_cache_index_ = idx;
    }

    /// Returns range of cell sides
    Range<DHCellSide> side_range() const;

    /// Returns range of neighbour cell of lower dimension corresponding to cell of higher dimension
    RangeConvert<DHNeighbSide, DHCellSide> neighb_sides() const;

    /// Return true if accessor represents own element (false for ghost element)
    inline bool is_own() const {
    	return (loc_ele_idx_ < dof_handler_->el_ds_->lsize());
    }

    /// Create new accessor with same local idx and given DOF handler. Actual and given DOF handler must be create on same Mesh.
    DHCellAccessor cell_with_other_dh(const DOFHandlerMultiDim * dh) const{
    	ASSERT( (dh->mesh()->n_nodes() == dof_handler_->mesh()->n_nodes()) && (dh->mesh()->n_elements() == dof_handler_->mesh()->n_elements()) )
    			.error("Incompatible DOF handlers!");
    	return DHCellAccessor(dh, loc_ele_idx_);
    }

    /// Iterates to next local element.
    inline void inc() {
        loc_ele_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const DHCellAccessor& other) const {
    	return (loc_ele_idx_ == other.loc_ele_idx_);
    }

    /// Comparison of accessors.
    bool operator!=(const DHCellAccessor& other) const {
    	return (loc_ele_idx_ != other.loc_ele_idx_);
    }


private:
    /// Pointer to the DOF handler owning the element.
    const DOFHandlerMultiDim * dof_handler_;
    /// Index into DOFHandler::el_4_loc array.
    unsigned int loc_ele_idx_;

    /// Optional member used in field evaluation, holds index of cell in field data cache.
    mutable unsigned int elm_cache_index_;

    friend class DHCellSide;
    friend class DHEdgeSide;
    friend class DHNeighbSide;
};


/**
 * @brief Side accessor allows to iterate over sides of DOF handler cell.
 *
 * Iterator is defined by:
 *  - DOF handler cell
 *  - index of Side in DOF cell (iterated value)
 */
class DHCellSide {
public:

    /**
     * @brief Default invalid accessor.
     *
     * Create invalid \p dh_cell_accessor_.
     */
	DHCellSide() {}

    /**
     * DOF cell side accessor.
     */
	DHCellSide(const DHCellAccessor &dh_cell_accessor, unsigned int side_idx)
    : dh_cell_accessor_(dh_cell_accessor), side_idx_(side_idx) {}

	/// Check validity of accessor (see default constructor)
    inline virtual bool is_valid() const {
        return dh_cell_accessor_.is_valid();
    }

    /// Return Side of given cell and side_idx.
    inline Side side() const {
    	ASSERT( this->is_valid() );
   		return Side(dh_cell_accessor_.dof_handler_->mesh(), dh_cell_accessor_.elm_idx(), side_idx_ );
    }

    /// Return DHCellAccessor appropriate to the side.
    inline const DHCellAccessor &cell() const {
    	return dh_cell_accessor_;
    }

    /// Return DHCellAccessor appropriate to the side.
    inline DHCellAccessor &cell() {
    	return dh_cell_accessor_;
    }

    /// Return dimension of element appropriate to the side.
    inline unsigned int dim() const {
    	return cell().dim();
    }

    /// Side centre.
    inline arma::vec3 centre() const {
    	return side().centre();
    }

    inline ElementAccessor<3> element() const {
    	return side().element();
    }

    inline unsigned int elem_idx() const {
    	return side().elem_idx();
    }

    inline Boundary cond() const {
        return side().cond();
    }

	inline unsigned int side_idx() const {
	   return side_idx_;
	}

	inline double measure() const {
	   return side().measure();
	}

	inline double diameter() const {
	   return side().diameter();
	}


    /// Returns range of all sides looped over common Edge.
    RangeConvert<DHEdgeSide, DHCellSide> edge_sides() const;

    /**
     * Returns total number of sides appropriate to Edge that owns actual cell side.
     *
     * return empty range if no element connected to Edge is local
     */
    unsigned int n_edge_sides() const;

    /// Iterates to next local element.
    inline virtual void inc() {
        side_idx_++;
    }

    /// Comparison of accessors.
	inline bool operator ==(const DHCellSide &other) {
		return this->elem_idx() == other.elem_idx() && side_idx_ == other.side_idx_;
	}

	inline bool operator !=(const DHCellSide &other) const {
		return this->elem_idx() != other.elem_idx() || side_idx_ != other.side_idx_;
	}

private:
    /// Appropriate DHCellAccessor.
    DHCellAccessor dh_cell_accessor_;
    /// Index of side.
    unsigned int side_idx_;

    friend class DHEdgeSide;
};


/**
 * @brief Class allows to iterate over sides of edge.
 *
 * Iterator is defined by:
 *  - DOF handler
 *  - global index of Edge (constant value)
 *  - index of Side in Edge (iterated value)
 *
 * Note: Class is used only internally. Appropriate range method (DHCellSide::edge_sides) uses convertible iterators
 * and returns corresponding DHCellSide.
 */
class DHEdgeSide {
public:
    /// Default invalid accessor.
	DHEdgeSide() : dof_handler_(nullptr), edge_idx_(0) {}

    /**
     * Valid accessor allows iterate over sides.
     */
	DHEdgeSide(const DHCellSide &cell_side, unsigned int side_idx)
    : dof_handler_(cell_side.dh_cell_accessor_.dof_handler_),
	  edge_idx_(cell_side.dh_cell_accessor_.elm()->edge_idx(cell_side.side_idx_)),
	  side_idx_(side_idx)
    {}

	/// Check validity of accessor (see default constructor)
    inline bool is_valid() const {
        return (dof_handler_!=nullptr);
    }

    /// Iterates to next edge side.
    inline void inc() {
        side_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const DHEdgeSide& other) {
    	return (edge_idx_ == other.edge_idx_) && (side_idx_ == other.side_idx_);
    }

    /// This class is implicitly convertible to DHCellSide.
    operator DHCellSide() const {
    	SideIter side = dof_handler_->mesh()->edge(edge_idx_).side(side_idx_);
    	DHCellAccessor cell = dof_handler_->cell_accessor_from_element( side->elem_idx() );
        return DHCellSide(cell, side->side_idx());
    }

private:
    /// Pointer to the DOF handler owning the element.
    const DOFHandlerMultiDim * dof_handler_;
    /// Global index of Edge.
    unsigned int edge_idx_;
    /// Index of side owned by Edge.
    unsigned int side_idx_;
};


/**
 * @brief Class allows to iterate over sides of neighbour.
 *
 * Class returns only local cells (owns + ghosts), non-local cells are skipped.
 *
 * Iterator is defined by:
 *  - DOF handler cell accessor
 *  - index of Neighbour on cell (iterated value)
 *  - maximal index of Neighbour (allow to skip non-local cells)
 *
 * Note: Class is used only internally. Appropriate range method (DHCellAccessor::neighb_sides) uses convertible
 * iterators and returns corresponding DHCellSide.
 */
class DHNeighbSide {
public:
    /// Default invalid accessor.
	DHNeighbSide() {}

    /**
     * @brief Valid accessor allows iterate over neighbor sides.
     *
     * @param dh_cell    Element of lower dim.
     * @param neighb_idx Index of neighbour.
     * @param max_idx    Maximal index of neighbour, method inc() doesn't set neighb_idx_ on higher value.
     */
	DHNeighbSide(const DHCellAccessor &dh_cell, unsigned int neighb_idx, unsigned int max_idx)
    : dh_cell_(dh_cell), neighb_idx_(neighb_idx), max_idx_(max_idx)
    {
	    // Skip non-local cells
	    while ( (neighb_idx_<max_idx_) && not_local_cell() ) {
        	neighb_idx_++;
        }
    }

	/// Check validity of accessor (see default constructor of DHCellAccessor)
    inline bool is_valid() const {
        return dh_cell_.is_valid();
    }

    /// Iterates to next neighbour side.
    inline void inc() {
    	// Skip non-local cells
        do {
            neighb_idx_++;
        	if (neighb_idx_>=max_idx_) break; //stop condition at the end item of range
        } while ( not_local_cell() );
    }

    /// Comparison of accessors.
    bool operator==(const DHNeighbSide& other) {
    	return (neighb_idx_ == other.neighb_idx_);
    }

    /// This class is implicitly convertible to DHCellSide.
    operator DHCellSide() const {
        SideIter side = dh_cell_.elm()->neigh_vb[neighb_idx_]->side();
        DHCellAccessor cell = dh_cell_.dof_handler_->cell_accessor_from_element( side->elem_idx() );
        return DHCellSide(cell, side->side_idx());
    }

private:
    /// Check if cell side of neighbour is not local (allow skip invalid accessors).
    inline bool not_local_cell() {
        return ( dh_cell_.dof_handler_->global_to_local_el_idx_.end() ==
            dh_cell_.dof_handler_->global_to_local_el_idx_.find((LongIdx)dh_cell_.elm()->neigh_vb[neighb_idx_]->side()->elem_idx()) );
    }

    /// Appropriate cell accessor.
    DHCellAccessor dh_cell_;
    /// Index into neigh_vb array
    unsigned int neighb_idx_;
    /// Maximal index into neigh_vb array
    unsigned int max_idx_;
};


/*************************************************************************************
 * Implementation of inlined methods.
 */


inline unsigned int DHCellAccessor::n_dofs() const
{
    switch (this->dim()) {
        case 1:
            return fe<1>()->n_dofs();
            break;
        case 2:
            return fe<2>()->n_dofs();
            break;
        case 3:
            return fe<3>()->n_dofs();
            break;
    }
    return 0; // only fix compiler warning
}


inline const Dof &DHCellAccessor::cell_dof(unsigned int idof) const
{
    switch (this->dim())
    {
        case 1:
            return fe<1>()->dof(idof);
            break;
        case 2:
            return fe<2>()->dof(idof);
            break;
        case 3:
            return fe<3>()->dof(idof);
            break;
    }

    ASSERT(0)(this->dim()).error("Unsupported FE dimension.");
    // cannot be reached:
    return fe<1>()->dof(idof);;
}


inline Range<DHCellSide> DHCellAccessor::side_range() const {
	auto bgn_it = make_iter<DHCellSide>( DHCellSide(*this, 0) );
	auto end_it = make_iter<DHCellSide>( DHCellSide(*this, dim()+1) );
	return Range<DHCellSide>(bgn_it, end_it);
}


inline RangeConvert<DHNeighbSide, DHCellSide> DHCellAccessor::neighb_sides() const {
	unsigned int upper_bound = this->elm()->n_neighs_vb();
	auto bgn_it = make_iter<DHNeighbSide, DHCellSide>( DHNeighbSide(*this, 0, upper_bound) );
	auto end_it = make_iter<DHNeighbSide, DHCellSide>( DHNeighbSide(*this, upper_bound, upper_bound) );
	return RangeConvert<DHNeighbSide, DHCellSide>(bgn_it, end_it);
}


inline RangeConvert<DHEdgeSide, DHCellSide> DHCellSide::edge_sides() const {
	return RangeConvert<DHEdgeSide, DHCellSide>(make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide( *this, 0) ),
	                                            make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide( *this, n_edge_sides()) ));
}


inline unsigned int DHCellSide::n_edge_sides() const {
    unsigned int edge_idx = dh_cell_accessor_.elm()->edge_idx(side_idx_);
    Edge edg = dh_cell_accessor_.dof_handler_->mesh()->edge(edge_idx);
    for (uint sid=0; sid<edg.n_sides(); sid++)
        if ( dh_cell_accessor_.dof_handler_->el_is_local(edg.side(sid)->element().idx()) ) return edg.n_sides();
    return 0;
}



#endif /* DH_CELL_ACCESSOR_HH_ */
