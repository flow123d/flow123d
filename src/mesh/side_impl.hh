/*
 * side_impl.hh
 *
 *  Created on: Jun 17, 2012
 *      Author: jb
 */

#ifndef SIDE_IMPL_HH_
#define SIDE_IMPL_HH_

#include "mesh/mesh.h"



    inline unsigned int Side::n_nodes() const {
        return dim()+1;
    }

    inline unsigned int Side::dim() const {
        return element_->dim-1;
    }

    // returns true for all sides either on boundary or connected to vb neigboring
    inline bool Side::is_external() const {
        return edge_->n_sides == 1;
    }

    inline const Node * Side::node(unsigned int i) const {
        // cout << "sn dim: " << dim() << "side: " << lnum << "node: "<< i << endl;
        int i_n = mesh_->side_nodes[dim()][el_idx_][i];
        // cout << "el node: "<< i_n << "nn: " << element->n_nodes << endl;

        return element_->node[ i_n ];
    }

    inline ElementFullIter Side::element()  {
        return mesh_->element.full_iter( element_ );
    }

    inline const Mesh * Side::mesh() const {
        return mesh_;
    }

    inline Edge * Side::edge() {
        return edge_;
    }

    inline unsigned int Side::el_idx() const {
        return el_idx_;
    }


#endif /* SIDE_IMPL_HH_ */
