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
    : dof_handler_(dof_handler), loc_ele_idx_(loc_idx), neighb_idx_(Mesh::undef_idx)
    {}

    /// Return local index to element (index of DOF handler).
    inline unsigned int local_idx() const {
    	ASSERT_LT_DBG(loc_ele_idx_, dof_handler_->el_ds_->lsize()).error("Method 'local_idx()' can't be used for ghost cells!\n");
        return loc_ele_idx_;
    }

    /// Return serial idx of appropriate element of DH cell accessor.
    inline unsigned int elm_idx() const {
    	//different for loop over local cells and neighbour elements
    	if (neighb_idx_==Mesh::undef_idx) return this->local_elm_idx();
    	else return dof_handler_->mesh()->element_accessor( this->local_elm_idx() )->neigh_vb[neighb_idx_]->side()->element().idx();
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

    /// Returns range of neighbour cells of higher dimension
    Range<DHCellAccessor> neighb_cells() const;

    /// Iterates to next local element.
    inline void inc() {
    	if (neighb_idx_==Mesh::undef_idx) loc_ele_idx_++;
    	else neighb_idx_++;
    }

    /// Comparison of accessors.
    bool operator==(const DHCellAccessor& other) {
    	return (loc_ele_idx_ == other.loc_ele_idx_) && (loc_ele_idx_ == other.loc_ele_idx_);
    }

private:
    /**
     * DOF cell accessor allows iterate over neighbour cells of higher dimension.
     */
	DHCellAccessor(const DOFHandlerMultiDim *dof_handler, unsigned int loc_idx, unsigned int neighb_idx)
    : dof_handler_(dof_handler), loc_ele_idx_(loc_idx), neighb_idx_(neighb_idx)
    {}

    /// Return serial idx to element of loc_ele_idx_.
    inline unsigned int local_elm_idx() const {
    	unsigned int ds_lsize = dof_handler_->el_ds_->lsize();
        if (loc_ele_idx_<ds_lsize) return dof_handler_->el_index(loc_ele_idx_); //own elements
        else return dof_handler_->ghost_4_loc[loc_ele_idx_-ds_lsize]; //ghost elements
    }

    /// Pointer to the DOF handler owning the element.
    const DOFHandlerMultiDim * dof_handler_;
    /// Index into DOFHandler::el_4_loc array.
    unsigned int loc_ele_idx_;
    /// Index into neigh_vb array (data member is used only for iterate over neighbour cells of higher dimension)
    unsigned int neighb_idx_;

    friend class DHCellSide;
};


/**
 * Side accessor allow iterate over sides of DOF handler cell.
 *
 * TODO: complete description of iterate over different ranges
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
    inline virtual const Side * side() const {
    	ASSERT( this->is_valid() );
   		return new Side( const_cast<const Mesh*>(dh_cell_accessor_.dof_handler_->mesh()), dh_cell_accessor_.elm_idx(), side_idx_ );
    }

    /// Returns range of all sides looped over common Edge.
    Range<DHCellSide> edge_sides() const;

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
};


/**
 *
 */
class DHEdgeSide : public DHCellSide {
public:
    /**
     * Default invalid accessor.
     */
	DHEdgeSide()
	: dof_handler_(nullptr) {}

    /**
     * Valid accessor allows iterate over sides.
     */
	DHEdgeSide(const DOFHandlerMultiDim * dof_handler, unsigned int edge_idx, unsigned int side_idx)
    : dof_handler_(dof_handler), edge_idx_(edge_idx), side_idx_(side_idx)
    {}

	/// Check validity of accessor (see default constructor)
    inline bool is_valid() const {
        return (dof_handler_ != nullptr);
    }

    /// Return Side of given cell and side_idx.
    inline const Side * side() const {
    	ASSERT( this->is_valid() );
   		return &(*dof_handler_->mesh()->edges[edge_idx_].side(side_idx_));
    }

    /// Iterates to next edge side.
    inline void inc() {
        side_idx_++;
    }

private:
    /// Pointer to the DOF handler owning the element.
    const DOFHandlerMultiDim * dof_handler_;
    /// Index of cell / edge (different for any AccessorTypes) through that is iterated.
    unsigned int edge_idx_;
    /// Index of side owned by Edge.
    unsigned int side_idx_;
};


/*************************************************************************************
 * Implementation of inlined methods.
 */

inline unsigned int DHCellAccessor::get_dof_indices(std::vector<int> &indices) const
{
  unsigned int elem_idx = this->elm_idx();
  ASSERT_LT( dof_handler_->row_4_el[elem_idx]+1, dof_handler_->cell_starts.size() )(dof_handler_->row_4_el[elem_idx])(dof_handler_->cell_starts.size());
  unsigned int ndofs = 0;
  ndofs = dof_handler_->cell_starts[dof_handler_->row_4_el[elem_idx]+1]-dof_handler_->cell_starts[dof_handler_->row_4_el[elem_idx]];
  for (unsigned int k=0; k<ndofs; k++)
    indices[k] = dof_handler_->dof_indices[dof_handler_->cell_starts[dof_handler_->row_4_el[elem_idx]]+k];

  return ndofs;
}


inline unsigned int DHCellAccessor::get_loc_dof_indices(std::vector<LongIdx> &indices) const
{
  unsigned int elem_idx = this->elm_idx();
  unsigned int ndofs = 0;
  ndofs = dof_handler_->cell_starts[dof_handler_->row_4_el[elem_idx]+1]-dof_handler_->cell_starts[dof_handler_->row_4_el[elem_idx]];
  for (unsigned int k=0; k<ndofs; k++)
    indices[k] = dof_handler_->cell_starts[dof_handler_->row_4_el[elem_idx]]+k;

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


inline Range<DHCellAccessor> DHCellAccessor::neighb_cells() const {
	auto bgn_it = make_iter<DHCellAccessor>( DHCellAccessor(dof_handler_, loc_ele_idx_, 0) );
	auto end_it = make_iter<DHCellAccessor>( DHCellAccessor(dof_handler_, loc_ele_idx_, this->elm()->n_neighs_vb()) );
	return Range<DHCellAccessor>(bgn_it, end_it);
}


inline Range<DHCellSide> DHCellSide::edge_sides() const {
	unsigned int edge_idx = dh_cell_accessor_.elm()->edge_idx(side_idx_);
	auto dh = dh_cell_accessor_.dof_handler_;
	return Range<DHCellSide>(make_iter<DHCellSide>( DHEdgeSide( dh, edge_idx, 0) ),
	                         make_iter<DHCellSide>( DHEdgeSide( dh, edge_idx, dh->mesh()->edges[edge_idx].n_sides) ));
}


#endif /* DH_CELL_ACCESSOR_HH_ */
