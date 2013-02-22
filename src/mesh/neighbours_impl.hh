/*
 * neigbours_impl.hh
 *
 *  Created on: Jun 28, 2012
 *      Author: jb
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
