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
 * @file    side_impl.h
 * @brief   
 */

#ifndef SIDE_IMPL_HH_
#define SIDE_IMPL_HH_

#include "mesh/mesh.h"

///////////////////////////
// SIDE inline implementation
///////////////////////////

    inline Side::Side(const Mesh * mesh, unsigned int elem_idx, unsigned int set_lnum)
    : mesh_(mesh), elem_idx_(elem_idx), side_idx_(set_lnum)
    {
        mesh_->check_element_size(elem_idx);
    }


    inline unsigned int Side::n_nodes() const {
        return dim()+1;
    }

    inline const Mesh * Side::mesh() const {
        return this->mesh_;
    }

    inline unsigned int Side::side_idx() const {
        return side_idx_;
    }


    inline unsigned int Side::elem_idx() const {
        return elem_idx_;
    }


    inline bool Side::valid() const {
        return mesh_!= NULL;
    }


    inline void Side::inc() {
    	side_idx_++;
    }



///////////////////////////
// SIDE_ITER inline implementation
///////////////////////////

    inline SideIter::SideIter(const Side &side)
    : side_(side)
    {}

    inline bool SideIter::operator==(const SideIter &other) {
        return (side_.mesh() == other.side_.mesh() ) && ( side_.elem_idx() == other.side_.elem_idx() )
        		&& ( side_.side_idx() == other.side_.side_idx() );
    }

    inline bool SideIter::operator!=(const SideIter &other) {
        return !( *this == other);
    }

    inline const Side & SideIter::operator *() const
            { return side_; }

    inline const Side * SideIter::operator ->() const
            { return &side_; }

    inline SideIter & SideIter::operator ++ () {
        side_.inc();
        return (*this);
    }
#endif /* SIDE_IMPL_HH_ */
