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
        return element->dim-1;
    }

    // returns true for all sides either on boundary or connected to vb neigboring
    inline bool Side::is_external() const {
        return edge->n_sides == 1;
    }

    inline const Node * Side::node(unsigned int i) const {
        // cout << "sn dim: " << dim() << "side: " << lnum << "node: "<< i << endl;
        int i_n = mesh->side_nodes[dim()][lnum][i];
        // cout << "el node: "<< i_n << "nn: " << element->n_nodes << endl;

        return element->node[ i_n ];
    }
/*
    inline ElementFullIter Side::element()  {
        return mesh->element.full_iter( element );
    }

    inline const Mesh * Side::mesh() const {
        return mesh;
    }

    inline Edge * Side::edge() {
        return edge;
    }

    inline unsigned int Side::el_idx() const {
        return lnum;
    }
*/

#endif /* SIDE_IMPL_HH_ */
