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

template <int spacedim>
class DofElementAccessor : public ElementAccessor<spacedim> {
public:
    /**
     * Default invalid accessor.
     */
    DofElementAccessor()
    : ElementAccessor<spacedim>()
    {}

    /**
     * Element accessor.
     */
    DofElementAccessor(const DOFHandlerMultiDim *dof_handler, unsigned int loc_idx)
    : dof_handler_(dof_handler), loc_ele_idx_(loc_idx)
    {
        this->mesh_ = dof_handler_->mesh();
    	set_data();
    }

    /// Iterates to next local element.
    inline void inc() override {
        loc_ele_idx_++;
		set_data();
    }

private:
    /// Set data of parent class ElementAccessor
    inline void set_data() {
    	this->element_idx_ = dof_handler_->el_index(loc_ele_idx_);
    	this->boundary_ = (this->element_idx_ >= this->mesh_->n_elements());
    	if ( dof_handler_->el_index(loc_ele_idx_) < this->mesh_->n_elements() ) {
    		this->r_idx_ = this->element()->region_idx();
    		this->dim_ = this->element()->dim();
    	}
    }

    /// Pointer to the DOF handler owning the element.
    const DOFHandlerMultiDim * dof_handler_;
    /// Index into DOFHandler::el_4_loc array.
    unsigned int loc_ele_idx_;
};


#endif /* DOF_ELEMENT_ACCESSORS_HH_ */
