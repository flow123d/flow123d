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
 * @file    dof_element_accessors.hh
 * @brief
 * @author  David Flanderka
 */

#ifndef DOF_ELEMENT_ACCESSORS_HH_
#define DOF_ELEMENT_ACCESSORS_HH_

#include "mesh/accessors.hh"

class DofElementAccessor {
public:
    /**
     * Default invalid accessor.
     */
    DofElementAccessor()
    : dof_handler_(NULL)
    {}

    /**
     * DOF element accessor.
     */
    DofElementAccessor(const DOFHandlerMultiDim *dof_handler, unsigned int loc_idx)
    : dof_handler_(dof_handler), loc_ele_idx_(loc_idx)
    {}

    /// Return local index to element (index of DOF handler).
    inline unsigned int local_idx() const {
        return loc_ele_idx_;
    }

    /// Return serial idx to element of loc_ele_idx_.
    inline unsigned int element_idx() const {
        return dof_handler_->el_index(loc_ele_idx_);
    }

    /// Return ElementAccessor to element of loc_ele_idx_.
    inline const ElementAccessor<3> element_accessor() const {
    	return dof_handler_->mesh()->element_accessor(loc_ele_idx_);
    }

    /// Iterates to next local element.
    inline void inc() {
        loc_ele_idx_++;
    }

    bool operator==(const DofElementAccessor& other) {
    	return (loc_ele_idx_ == other.loc_ele_idx_);
    }

    /**
     * -> dereference operator
     *
     * Return ElementAccessor to element of loc_ele_idx_. Allow to simplify code:
 @code
     DofElementAccessor dof_ac(dh, loc_index);
     unsigned int dim;
     dim = dof_ac.element_accessor().dim();  // full format of access to element
     dim = dof_ac->dim();                    // short format with dereference operator
 @endcode
     */
    inline const ElementAccessor<3> operator ->() const {
    	return dof_handler_->mesh()->element_accessor(loc_ele_idx_);
    }

private:
    /// Pointer to the DOF handler owning the element.
    const DOFHandlerMultiDim * dof_handler_;
    /// Index into DOFHandler::el_4_loc array.
    unsigned int loc_ele_idx_;
};


#endif /* DOF_ELEMENT_ACCESSORS_HH_ */
