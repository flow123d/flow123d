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

#include "mesh/accessors.hh"
#include "mesh/sides.h"
#include "mesh/neighbours.h"
#include "fem/dofhandler.hh"

class DHCellSide;
class DHNeighbSide;
class DHEdgeSide;

/**
 * Cell accessor allow iterate over DOF handler cells.
 *
 * Iterating is possible over different ranges of local and ghost elements.
 */
class DHCellAccessor {
public:
    /**
     * Default invalid accessor.
     */
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
    	ASSERT_LT_DBG(loc_ele_idx_, dof_handler_->el_ds_->lsize()+dof_handler_->ghost_4_loc.size()).error("Local element index is out of range!\n");
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
    unsigned int get_dof_indices(std::vector<int> &indices) const;

    /**
     * @brief Returns the indices of dofs associated to the cell on the local process.
     *
     * @param indices Array of dof indices on the cell.
     */
    unsigned int get_loc_dof_indices(std::vector<LongIdx> &indices) const;

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

    /**
     * @brief Returns finite element object for given space dimension.
     */
    template<unsigned int dim>
    FiniteElement<dim> *fe() const {
        ElementAccessor<3> elm_acc = this->elm();
        return dof_handler_->ds_->fe<dim>(elm_acc);
    }

    /// Check validity of accessor (see default constructor)
    inline bool is_valid() const {
        return dof_handler_ != NULL;
    }

    /// Returns range of cell sides
    Range<DHCellSide> side_range() const;

    /// Returns range of neighbour cell of lower dimension corresponding to cell of higher dimension
    Range<DHNeighbSide> neighb_sides() const;

    /// Return true if accessor represents own element (false for ghost element)
    inline bool is_own() const {
    	return (loc_ele_idx_ < dof_handler_->el_ds_->lsize());
    }

    /// Iterates to next local element.
    inline void inc() {
        loc_ele_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const DHCellAccessor& other) {
    	return (loc_ele_idx_ == other.loc_ele_idx_);
    }

private:
    /// Pointer to the DOF handler owning the element.
    const DOFHandlerMultiDim * dof_handler_;
    /// Index into DOFHandler::el_4_loc array.
    unsigned int loc_ele_idx_;

    friend class DHCellSide;
    friend class DHEdgeSide;
};


/**
 * Side accessor allows to iterate over sides of DOF handler cell.
 */
class DHCellSide {
public:

    /**
     * Default invalid accessor.
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
    inline const Side * side() const {
    	ASSERT( this->is_valid() );
   		return new Side( const_cast<const Mesh*>(dh_cell_accessor_.dof_handler_->mesh()), dh_cell_accessor_.elm_idx(), side_idx_ );
    }

    /// Return DHCellAccessor appropriate to the side.
    inline const DHCellAccessor cell() const {
    	return dh_cell_accessor_;
    }

    /// Return dimension of element appropriate to the side.
    inline unsigned int dim() const {
    	return cell().dim();
    }

    /// Returns range of all sides looped over common Edge.
    RangeConvert<DHEdgeSide, DHCellSide> edge_sides() const;

    /// Iterates to next local element.
    inline virtual void inc() {
        side_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const DHCellSide& other) {
    	return (side_idx_ == other.side_idx_);
    }

private:
    /// Appropriate DHCellAccessor.
    DHCellAccessor dh_cell_accessor_;
    /// Index of side.
    unsigned int side_idx_;

    friend class DHEdgeSide;
};


/**
 * Class allows to iterate over sides of edge.
 *
 * Iterator is defined by:
 *  - DOF handler
 *  - global index of Edge (constant value)
 *  - index of Side in Edge (iterated value)
 */
class DHEdgeSide {
public:
    /**
     * Default invalid accessor.
     */
	DHEdgeSide() : dof_handler_(nullptr), edge_idx_(0) {}

    /**
     * Valid accessor allows iterate over sides.
     */
	DHEdgeSide(const DHCellSide &cell_side, unsigned int side_idx)
    : dof_handler_(cell_side.dh_cell_accessor_.dof_handler_), edge_idx_(cell_side.dh_cell_accessor_.elm()->edge_idx(cell_side.side_idx_)), side_idx_(side_idx)
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
    	SideIter side = dof_handler_->mesh()->edges[edge_idx_].side(side_idx_);
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
 * Class allows to iterate over sides of neighbour.
 */
class DHNeighbSide {
public:
    /**
     * Default invalid accessor.
     */
	DHNeighbSide() {}

    /**
     * Valid accessor allows iterate over neighbor sides.
     *
     * @param dh_cell    Element of lower dim.
     * @param neighb_idx Index of neighbour.
     */
	DHNeighbSide(const DHCellAccessor &dh_cell, unsigned int neighb_idx)
    : dh_cell_(dh_cell), neighb_idx_(neighb_idx)
    {}

	/// Check validity of accessor (see default constructor)
    inline bool is_valid() const {
        return dh_cell_.is_valid();
    }

    /// Return DHCellSide according to this object.
    inline DHCellSide cell_side() const {
    	unsigned int side_idx = neighbour()->side()->side_idx();
    	return DHCellSide(dh_cell_, side_idx);
    }

    /// Return Neighbour object according to this object.
    inline Neighbour * neighbour() const {
    	ASSERT( this->is_valid() );
    	return dh_cell_.elm()->neigh_vb[neighb_idx_];
    }

    /// Iterates to next edge side.
    inline void inc() {
    	neighb_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const DHNeighbSide& other) {
    	return (neighb_idx_ == other.neighb_idx_);
    }

private:
    /// Appropriate cell accessor.
    DHCellAccessor dh_cell_;
    /// Index into neigh_vb array
    unsigned int neighb_idx_;
};


/*************************************************************************************
 * Implementation of inlined methods.
 */

inline unsigned int DHCellAccessor::get_dof_indices(std::vector<int> &indices) const
{
  ASSERT_LT( loc_ele_idx_+1, dof_handler_->cell_starts.size() )(loc_ele_idx_)(dof_handler_->cell_starts.size());
  unsigned int ndofs = 0;
  ndofs = dof_handler_->cell_starts[loc_ele_idx_+1]-dof_handler_->cell_starts[loc_ele_idx_];
  for (unsigned int k=0; k<ndofs; k++)
    indices[k] = dof_handler_->dof_indices[dof_handler_->cell_starts[loc_ele_idx_]+k];

  return ndofs;
}


inline unsigned int DHCellAccessor::get_loc_dof_indices(std::vector<LongIdx> &indices) const
{
  unsigned int ndofs = 0;
  ndofs = dof_handler_->cell_starts[loc_ele_idx_+1]-dof_handler_->cell_starts[loc_ele_idx_];
  for (unsigned int k=0; k<ndofs; k++)
    indices[k] = dof_handler_->cell_starts[loc_ele_idx_]+k;

  return ndofs;
}


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
}


inline Range<DHCellSide> DHCellAccessor::side_range() const {
	auto bgn_it = make_iter<DHCellSide>( DHCellSide(*this, 0) );
	auto end_it = make_iter<DHCellSide>( DHCellSide(*this, dim()+1) );
	return Range<DHCellSide>(bgn_it, end_it);
}


inline Range<DHNeighbSide> DHCellAccessor::neighb_sides() const {
	auto bgn_it = make_iter<DHNeighbSide>( DHNeighbSide(*this, 0) );
	auto end_it = make_iter<DHNeighbSide>( DHNeighbSide(*this, this->elm()->n_neighs_vb()) );
	return Range<DHNeighbSide>(bgn_it, end_it);
}


inline RangeConvert<DHEdgeSide, DHCellSide> DHCellSide::edge_sides() const {
	unsigned int edge_idx = dh_cell_accessor_.elm()->edge_idx(side_idx_);
	Edge *edg = &dh_cell_accessor_.dof_handler_->mesh()->edges[edge_idx];
	unsigned int upper_bound = 0; // return empty range if no element connected to Edge is local
	for (int sid=0; sid<edg->n_sides; sid++)
	    if ( dh_cell_accessor_.dof_handler_->el_is_local(edg->side(sid)->element().idx()) )
	    {
	    	upper_bound = edg->n_sides;
	    	break;
	    }

	return RangeConvert<DHEdgeSide, DHCellSide>(make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide( *this, 0) ),
	                                            make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide( *this, upper_bound) ));
}


#endif /* DH_CELL_ACCESSOR_HH_ */
