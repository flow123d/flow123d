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
    n_sides_ = NDEF;

    // number of element of particular dimension
    n_lines = 0;
    n_triangles = 0;
    n_tetrahedras = 0;

    // indices of side nodes in element node array
    // Currently this is made ad libitum
    // with some ordering here we can get sides with correct orientation.
    // This speedup normal calculation.

    side_nodes.resize(3); // three side dimensions
    for(int i=0; i < 3; i++) {
        side_nodes[i].resize(i+2); // number of sides
        for(int j=0; j < i+2; j++)
            side_nodes[i][j].resize(i+1);
    }

    side_nodes[0][0][0] = 0;
    side_nodes[0][1][0] = 1;


    side_nodes[1][0][0] = 0;
    side_nodes[1][0][1] = 1;

    side_nodes[1][1][0] = 1;
    side_nodes[1][1][1] = 2;

    side_nodes[1][2][0] = 2;
    side_nodes[1][2][1] = 0;


    side_nodes[2][0][0] = 1;
    side_nodes[2][0][1] = 2;
    side_nodes[2][0][2] = 3;

    side_nodes[2][1][0] = 0;
    side_nodes[2][1][1] = 2;
    side_nodes[2][1][2] = 3;

    side_nodes[2][2][0] = 0;
    side_nodes[2][2][1] = 1;
    side_nodes[2][2][2] = 3;

    side_nodes[2][3][0] = 0;
    side_nodes[2][3][1] = 1;
    side_nodes[2][3][2] = 2;


}


unsigned int Mesh::n_sides()
{
    if (n_sides_ == NDEF) {
        n_sides_=0;
        FOR_ELEMENTS(this, ele) n_sides_ += ele->n_sides();
    }
    return n_sides_;
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

    count_element_types();

    // topology
    node_to_element(mesh);

    read_neighbours();
    edge_to_side(mesh);

    neigh_vb_to_element_and_side(mesh);
    element_to_neigh_vb(mesh);
    neigh_vb_to_edge_both(mesh);

    count_side_types(mesh);
    xprintf(MsgVerb, "Topology O.K.\n")/*orig verb 4*/;

    read_boundary(mesh);

}


void Mesh::read_neighbours() {
    FILE    *in;   // input file
    char     line[ LINE_SIZE ];   // line of data file
    unsigned int id;

    xprintf( Msg, "Reading neighbours...");
    const std::string& file_name = IONameHandler::get_instance()->get_input_file_name(OptGetStr( "Input", "Neighbouring", "\\" ));
    in = xfopen( file_name, "rt" );
    skip_to( in, "$Neighbours" );
    xfgets( line, LINE_SIZE - 2, in );

    unsigned int n_neighs = atoi( xstrtok( line) );
    INPUT_CHECK( n_neighs > 0 ,"Number of neighbours  < 1 in read_neighbour_list()\n");
    vb_neighbours_.resize( n_neighs );


    FOR_NEIGHBOURS(this,  ngh ) {
        xfgets( line, LINE_SIZE - 2, in );

        id              = atoi( xstrtok( line) );
        ngh->type            = atoi( xstrtok( NULL) );

        switch( ngh->type ) {
            case BB_E:
                xprintf(UsrErr, "Not supported - Neighboring of type (10) - of elements of same dimension without local side number!\n");
                break;
            case BB_EL:
                ngh->n_elements = atoi( xstrtok( NULL) );
                INPUT_CHECK(!( ngh->n_elements < 2 ),"Neighbour %d has bad number of elements: %d\n", id, ngh->n_elements );
                ngh->n_sides = ngh->n_elements;

                ngh->eid = (int*) xmalloc( ngh->n_elements * sizeof( int ) );
                ngh->sid = (int*) xmalloc( ngh->n_elements * sizeof( int ) );
                ngh->element_ = (ElementIter *) xmalloc( ngh->n_elements *
                        sizeof( ElementIter  ) );
                ngh->side_ = new SideIter [ ngh->n_elements ];


                for( int i = 0; i < ngh->n_sides; i++) {
                    ngh->element_[ i ] = NULL;
                    ngh->eid[ i ] = atoi( xstrtok( NULL) );
                    ngh->sid[ i ] = atoi( xstrtok( NULL) );
                }

                break;
            case VB_ES:
                ngh->n_elements = 2;
                ngh->n_sides = 2;
                ngh->eid = (int*) xmalloc( ngh->n_elements * sizeof( int ) );
                ngh->element_ = (ElementIter *) xmalloc( ngh->n_elements * sizeof( ElementIter  ) );
                ngh->sid = (int*) xmalloc( ngh->n_sides * sizeof( int ) );
                ngh->side_ = new SideIter [ ngh->n_sides ];

                ngh->eid[ 0 ] = atoi( xstrtok( NULL) );
                ngh->eid[ 1 ] = atoi( xstrtok( NULL) );
                ngh->sid[ 0 ] = NDEF;
                ngh->sid[ 1 ] = atoi( xstrtok( NULL) );
                ngh->sigma = atof( xstrtok( NULL) );
                //ngh->element[ 0 ] = NULL;
                //ngh->element[ 1 ] = NULL;

                // sid[ 0 ] (and side[ 0 ]) doesn't have defined value. I use sid[ 1 ] (and
                // side[ 1 ]) instead to correspond with elm1. Using sid[ 0 ] for describing
                // side of elm1 would be confusing and error-prone.
                break;
            case VV_2E:
                xprintf(UsrErr, "Not supported - Neighboring of type (30) - Noncompatible only elements!\n");
                break;
            default:
                xprintf(UsrErr,"Neighbour %d is of the unsupported type %d\n", id, ngh->type );
                break;
        }
    }

    //xprintf( Msg, " %d neighbours readed. ", n_vb_neighbours() );
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
