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
#include "input/input_type.hh"

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


Input::Type::Record BoundarySegment::get_input_type() {
    using namespace Input::Type;

    static Record rec("BoundarySegment","Record with specification of boundary segments,\n"
            "i.e. subsets of the domain boundary where we prescribe one boundary condition.\n"
            "Currently very GMSH oriented. (NOT IMPLEMENTED YET)");
    if (!rec.is_finished()) {
        rec.declare_key("physical_domains", Array(Integer(0)),
                "Numbers of physical domains (submeshes) that forms a segment and that "
                "will be removed from the computational mesh.");
        rec.declare_key("elements", Array(Integer(0)),
                        "Numbers of elements that forms a segment and that "
                        "will be removed from the computational mesh.");
        rec.declare_key("sides", Array( Array(Integer(0), 2,2) ),
                        "Pairs [ element, local_side] specifying sides that are part of the boundary segment."
                        "Sides are NOT removed from computation.");
        rec.finish();
    }
    return rec;
}

Input::Type::Record Mesh::get_input_type() {
    using namespace Input::Type;

    static Record rec("Mesh","Record with mesh related data." );
    if (!rec.is_finished()) {
        rec.declare_key("mesh_file", FileName::input(), Default::obligatory(),
                "Input file with mesh description.");
        rec.declare_key("boundary_segmants", Array( BoundarySegment::get_input_type() ),
                "Array with specification of boundary segments");
        rec.declare_key("neighbouring", FileName::input(), Default::obligatory(),
                "File with mesh connectivity data.");
        rec.finish();
    }
    return rec;
}

Mesh::Mesh(Input::Record in_record)
: in_record_(in_record) {

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
    switch (elm->dim()) {
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

/**
 *  Setup whole topology for read mesh.
 */
void Mesh::setup_topology() {
    F_ENTRY;
    Mesh *mesh=this;


    count_element_types();

    // topology
    //node_to_element();

    read_neighbours();
    edge_to_side();

    neigh_vb_to_element_and_side();
    element_to_neigh_vb();

    count_side_types();

    //read_boundary(mesh);

    xprintf(MsgVerb, "Topology O.K.\n")/*orig verb 4*/;


    // cleanup
    {
        vector<Neighbour_both> empty_vec;
        neighbours_.swap(empty_vec);
    }

}


/**
 *   Creates back references from nodes to elements.
 *
 *   TODO: This is not necessary after the topology setup so
 *   we should make that as an independent structure which can be easily deleted.
 */
void Mesh::node_to_element()
{
    F_ENTRY;
/*
    int li;
    NodeIter nod;
    ElementIter ele;

    xprintf( MsgVerb, "   Node to element... ");

    // Set counter of elements in node to zero
    FOR_NODES(this,  nod )
        nod->n_elements = 0;
    // Count elements
    FOR_ELEMENTS(this,  ele )
        FOR_ELEMENT_NODES( ele, li ) {
            nod = ele->node[ li ];
            (nod->n_elements)++;
        }
    // Allocate arrays
    FOR_NODES(this,  nod ) {
                if (nod->n_elements == 0)
                        continue;
            nod->element = (ElementIter *) xmalloc( nod->n_elements * sizeof( ElementIter ) );
        nod->aux = 0;
    }
    // Set poiners in arrays
    FOR_ELEMENTS(this,  ele )
        FOR_ELEMENT_NODES( ele, li ) {
            nod = ele->node[ li ];
            nod->element[ nod->aux ] = ele;
            (nod->aux)++;
        }
    xprintf( MsgVerb, "O.K.\n");*/
}

//=============================================================================
//
//=============================================================================
void Mesh::count_side_types()
{
    struct Side *sde;

    n_insides = 0;
    n_exsides = 0;
    FOR_SIDES(this,  sde )
        if (sde->is_external()) n_exsides++;
        else n_insides++;
}



void Mesh::read_neighbours() {
    FILE    *in;   // input file
    char     line[ LINE_SIZE ];   // line of data file
    unsigned int id;

    xprintf( Msg, "Reading neighbours...A\n");
    in = xfopen(  in_record_.val<FilePath>("neighbouring"), "rt" );
    skip_to( in, "$Neighbours" );
    xfgets( line, LINE_SIZE - 2, in );

    unsigned int n_neighs = atoi( xstrtok( line) );
    INPUT_CHECK( n_neighs > 0 ,"Number of neighbours  < 1 in read_neighbour_list()\n");
    neighbours_.resize( n_neighs );

    n_bb_neigh = 0;
    n_vb_neigh = 0;

    for(vector<Neighbour_both>::iterator ngh= neighbours_.begin();
            ngh != neighbours_.end(); ++ngh ) {
        xfgets( line, LINE_SIZE - 2, in );

        id              = atoi( xstrtok( line) );
        ngh->type            = atoi( xstrtok( NULL) );

        switch( ngh->type ) {
            case BB_E:
                xprintf(UsrErr, "Not supported - Neighboring of type (10) - of elements of same dimension without local side number!\n");
                break;
            case BB_EL:
                n_bb_neigh++;
                ngh->n_sides = atoi( xstrtok( NULL) );
                INPUT_CHECK(!( ngh->n_sides < 2 ),"Neighbour %d has bad number of elements: %d\n", id, ngh->n_sides );

                ngh->eid = new int [ngh->n_sides];
                ngh->sid = new int [ngh->n_sides];

                for( int i = 0; i < ngh->n_sides; i++) {
                    ngh->eid[ i ] = atoi( xstrtok( NULL) );
                    ngh->sid[ i ] = atoi( xstrtok( NULL) );
                }

                break;
            case VB_ES:
                n_vb_neigh++;
                ngh->n_sides = 2;
                ngh->eid = new int [ngh->n_sides];
                ngh->sid = new int [ngh->n_sides];

                ngh->eid[ 0 ] = atoi( xstrtok( NULL) );
                ngh->eid[ 1 ] = atoi( xstrtok( NULL) );
                ngh->sid[ 0 ] = NDEF;
                ngh->sid[ 1 ] = atoi( xstrtok( NULL) );

                ngh->sigma = atof( xstrtok( NULL) );
                break;
            case VV_2E:
                xprintf(UsrErr, "Not supported - Neighboring of type (30) - Noncompatible only elements!\n");
                break;
            default:
                xprintf(UsrErr,"Neighbour %d is of the unsupported type %d\n", id, ngh->type );
                break;
        }
    }

    xprintf( Msg, " %d VB neighbours %d BB neigs. readed. ", n_vb_neigh, n_bb_neigh );
}



void Mesh::edge_to_side()
{
    F_ENTRY;

    struct Edge *edg;
    Element *ele;

    xprintf( MsgVerb, "   Edge to side and back... \n");

    // count edges
    unsigned int n_edges = n_sides();
    for(vector<Neighbour_both>::iterator it= neighbours_.begin();
        it != neighbours_.end(); ++it )
        if ( it->type == BB_EL ) n_edges -= ( it->n_sides - 1 );

    // create edge vector
    edge.resize(n_edges);
    xprintf( Msg, "Created  %d edges.\n.", n_edges );

    // set edge, side connections
    unsigned int i_edge=0;
    for(vector<Neighbour_both>::iterator it= neighbours_.begin();
            it != neighbours_.end(); ++it ) {

        if ( it->type != BB_EL ) continue;

        edg = &( edge[i_edge++] );

        // init edge (can init all its data)
        edg->n_sides = it->n_sides;
        edg->side_ = new SideIter [edg->n_sides];

        for(int si=0; si < it->n_sides; si++) {
            ele = element.find_id( it->eid[si] );
            edg->side_[ si ] = ele->side( it->sid[ si ] );
            ele->edges_[ it->sid[ si ] ] = edg;
        }
    }

    // now the external ones ( pair all remaining edges with external sides)
    FOR_SIDES(this, sde) {
        if ( sde->edge() == NULL ) {
            edg = &( edge[i_edge++] );

            // make external edges and edges on neighborings.
            edg->n_sides = 1;
            edg->side_ = new SideIter [ edg->n_sides ];
            edg->side_[ 0 ] = sde;

            sde->element()->edges_[sde->el_idx()] = edg;
        }
    }
    ASSERT(i_edge == n_edges, "Actual number of edges %d do not match size %d of its array.\n", i_edge, n_edges);

    //FOR_SIDES(mesh, side) ASSERT(side->edge != NULL, "Empty side %d !\n", side->id);
}



/**
 * Make
 */
void Mesh::neigh_vb_to_element_and_side()
{

    vb_neighbours_.resize( n_vb_neigh );

    ElementIter ele_lower, ele_higher;
    Edge *edg;

    xprintf( MsgVerb, "   Creating %d VB neighbours... ", n_vb_neigh);

    vector<Neighbour>::iterator new_ngh = vb_neighbours_.begin();

    for(vector<Neighbour_both>::iterator ngh= neighbours_.begin(); ngh != neighbours_.end(); ++ngh ) {

        if ( ngh->type != VB_ES ) continue;


        ele_lower = element.find_id( ngh->eid[0]);
        ele_higher = element.find_id( ngh->eid[1] );
        edg = ele_higher->side( ngh->sid[ 1 ] )->edge();

        ASSERT(edg->n_sides == 1, "Edge with %d\n", edg->n_sides);
        ASSERT( ele_higher == edg->side(0)->element(),"Diff els.\n");
        new_ngh->reinit(  ele_lower, edg , ngh->sigma);


        //DBGMSG(" %d %d -> %d %d\n", ngh->eid[0], ngh->eid[1],
        //        new_ngh->element()->index(),
        //        new_ngh->side()->element()->index());

        ++new_ngh;
    }

    ASSERT( new_ngh == vb_neighbours_.end(), "Some VB neigbourings wasn't set.\n");


    xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}




//=============================================================================
//
//=============================================================================
void Mesh::element_to_neigh_vb()
{

    xprintf( MsgVerb, "   Element to neighbours of vb2 type... ")/*orig verb 5*/;

    FOR_ELEMENTS(this,ele) ele->n_neighs_vb =0;

    // count vb neighs per element
    FOR_NEIGHBOURS(this,  ngh )  ngh->element_->n_neighs_vb++;

    // Allocation of the array per element
    FOR_ELEMENTS(this,  ele )
        if( ele->n_neighs_vb > 0 ) {
            ele->neigh_vb = new struct Neighbour* [ele->n_neighs_vb];
            ele->n_neighs_vb=0;
        }

    // fill
    ElementIter ele;
    FOR_NEIGHBOURS(this,  ngh ) {
        ele = ngh->element();
        ele->neigh_vb[ ele->n_neighs_vb++ ] = &( *ngh );
    }

    xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
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
