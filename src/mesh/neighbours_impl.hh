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
 * @file    neighbours_impl.hh
 * @brief   
 */

#ifndef NEIGBOURS_IMPL_HH_
#define NEIGBOURS_IMPL_HH_

#include "mesh/sides.h"
#include "mesh/edges.h"


    // side of the edge in higher dim. mesh
    inline SideIter Neighbour::side() {
        ASSERT( edge()->n_sides == 1 , "VB neighbouring with %d sides.\n", edge()->n_sides);
        //DBGMSG("VB neighbouring with %d sides.\n", edge_->n_sides);
        return edge()->side(0);
    }

    inline unsigned int Neighbour::edge_idx() {
        return edge_idx_;
    }

    // edge of lower dimensional mesh in VB neigh.
    inline Edge *Neighbour::edge() {
        return &( element_->mesh_->edges[ edge_idx_] );
    }

    // element of higher dimension mesh in VB neigh.
    inline ElementIter Neighbour::element() {
        return element_;
    }

#endif /* NEIGBOURS_IMPL_HH_ */
