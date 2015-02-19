/*
 * side_impl.hh
 *
 *  Created on: Jun 17, 2012
 *      Author: jb
 */

#ifndef SIDE_IMPL_HH_
#define SIDE_IMPL_HH_

#include "mesh/mesh.h"
#include "mesh/edges.h"


inline Side::Side(const Element * ele, unsigned int set_lnum)
: element_(ele), el_idx_(set_lnum)
{
    ASSERT(mesh()->element.full_iter( const_cast<Element *>(element_) ), "Wrong initialization of the Side.\n");
}


    inline unsigned int Side::n_nodes() const {
        return dim()+1;
    }

    inline unsigned int Side::dim() const {
        return element_->dim()-1;
    }

    // returns true for all sides either on boundary or connected to vb neigboring
    inline bool Side::is_external() const {
        return edge()->n_sides == 1;
    }

    inline const Node * Side::node(unsigned int i) const {
        int i_n = mesh()->side_nodes[dim()][el_idx_][i];

        return element_->node[ i_n ];
    }

    inline ElementFullIter Side::element() const {
        ASSERT( valid(), "Wrong use of uninitialized accessor.\n");
        return mesh()->element.full_iter( const_cast<Element *>(element_) );
    }

    inline Mesh * Side::mesh() const {
        return element_->mesh_;
    }

    inline unsigned int Side::edge_idx() const {
        return element_->edge_idx_[el_idx()];
    }

    inline Edge * Side::edge() const {
        if (edge_idx() == Mesh::undef_idx) return NULL;
        else return &( mesh()->edges[ edge_idx() ] );
    }

    inline Boundary * Side::cond() const {
            if (cond_idx() == Mesh::undef_idx) return NULL;
            else return &( mesh()->boundary_[ cond_idx() ] );
    }

    inline unsigned int Side::cond_idx() const {
         if (element_->boundary_idx_ == NULL) return Mesh::undef_idx;
         else return element_->boundary_idx_[el_idx()];
    }


    inline unsigned int Side::el_idx() const {
        return el_idx_;
    }


    inline bool Side::valid() const {
        return element_!= NULL;
    }


    inline void Side::inc() {
        el_idx_++;
    }


    inline void *Side::make_ptr() const {
        return (void *)((long int) element_ << (2 + el_idx_) );
    }
#endif /* SIDE_IMPL_HH_ */
