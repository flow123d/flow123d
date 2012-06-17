/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @ingroup mesh
 * @brief  Mesh construction
 *
 */

#include <unistd.h>


#include "system/system.hh"
#include "mesh/mesh.h"

// think about following dependencies
#include "mesh/boundaries.h"


//TODO: sources, concentrations, initial condition  and similarly boundary conditions should be
// instances of a Element valued field
// concentrations is in fact reimplemented in transport REMOVE it HERE

// After removing non-geometrical things from mesh, this should be part of mash initializing.
#include "mesh/topology.cc"
#include "mesh/msh_reader.h"
#include "mesh/msh_gmshreader.h"

void count_element_types(Mesh*);
void read_node_list(Mesh*);

Mesh::Mesh() {
    xprintf(Msg, " - Mesh()     - version with node_vector\n");

    n_materials = NDEF;

    n_insides = NDEF;
    n_exsides = NDEF;
    n_neighs = NDEF;
    neighbour = NULL;
    l_neighbour = NULL;

    // Hashes
    n_lines = 0;
    n_triangles = 0;
    n_tetrahedras = 0;

}


//=============================================================================
// COUNT ELEMENT TYPES
//=============================================================================

void Mesh::count_element_types() {
    F_ENTRY;

    FOR_ELEMENTS(this, elm)
    switch (elm->dim) {
        case 1:
            n_lines++;
            break;
        case 2:
            n_triangles++;
            break;
        case 3:
            n_tetrahedras++;
            break;
    }
}


void Mesh::setup_topology() {
    F_ENTRY;
    Mesh *mesh=this;

    /// initialize mesh topology (should be handled inside mesh object)
    read_neighbour_list(mesh);


    count_element_types();

    // topology
    node_to_element(mesh);
    element_to_side_both(mesh);

    neigh_vv_to_element(mesh);
    //element_to_neigh_vv(mesh);
    neigh_vb_to_element_and_side(mesh);
    neigh_bv_to_side(mesh);
    element_to_neigh_vb(mesh);

    //side_shape_specific(mesh);
    side_to_node(mesh);
    neigh_bb_topology(mesh);

    make_edge_list(mesh);

    neigh_bb_to_edge_both(mesh);

    edge_to_side_both(mesh);
    neigh_vb_to_edge_both(mesh);
    side_types(mesh);
    count_side_types(mesh);
    xprintf(MsgVerb, "Topology O.K.\n")/*orig verb 4*/;

    read_boundary(mesh);

}

void Mesh::setup_materials( MaterialDatabase &base)
{
    xprintf( MsgVerb, "   Element to material... ")/*orig verb 5*/;
    FOR_ELEMENTS(this, ele ) {
        ele->material=base.find_id(ele->mid);
        INPUT_CHECK( ele->material != base.end(),
                "Reference to undefined material %d in element %d\n", ele->mid, ele.id() );
    }
    xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}


//-----------------------------------------------------------------------------
// vim: set cindent:
