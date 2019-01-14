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
#include "fem/dofhandler.hh"

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
};


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


#endif /* DH_CELL_ACCESSOR_HH_ */
