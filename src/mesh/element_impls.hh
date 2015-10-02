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
 * @file    element_impls.hh
 * @brief   
 */

#ifndef ELEMENT_IMPLS_HH_
#define ELEMENT_IMPLS_HH_

#include "elements.h"
#include "sides.h"
#include "side_impl.hh"

inline unsigned int Element::dim() const {
    return dim_;
}


inline unsigned int Element::index() const {
    return mesh_->element.index( this );
}


inline unsigned int Element::n_nodes() const {
    return dim()+1;
}



inline unsigned int Element::n_sides() const {
    return dim()+1;
}



inline SideIter Element::side(const unsigned int loc_index) {
    return SideIter( Side(this, loc_index) );
}

inline const SideIter Element::side(const unsigned int loc_index) const {
    return SideIter( Side(this, loc_index) );
}
#endif /* ELEMENT_IMPLS_HH_ */
