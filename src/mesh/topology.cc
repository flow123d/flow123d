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
 *
 * @file
 * @ingroup mesh
 * @brief Functions for construction of all pointers in the Mesh.
 *
 */


// remove dependency on following:
#include "materials.hh"
#include "mesh/boundaries.h"


static void element_to_material(Mesh*, MaterialDatabase &);
static void node_to_element(Mesh*);
static void element_to_side_both(Mesh*);
static void neigh_vv_to_element(Mesh*);
static void element_to_neigh_vv(Mesh*);
static void neigh_vb_to_element_and_side(Mesh*);
//static void neigh_bv_to_side(Mesh *);
static void element_to_neigh_vb(Mesh*);
static void side_to_node(Mesh*);
//static void neigh_bb_topology(Mesh*);
//static void neigh_bb_to_edge_both(Mesh*);
static void edge_to_side_both(Mesh*);
static void neigh_vb_to_edge_both(Mesh*);
static void side_types(Mesh*);
static void count_side_types(Mesh*);
//static void boundary_to_side_both(Mesh*);
static void side_to_node_line(ElementFullIter );
static void side_to_node_triangle(ElementFullIter );
static void side_to_node_tetrahedron(ElementFullIter );
static void neigh_bb_to_element(struct Neighbour*,Mesh*);
static void neigh_bb_el_to_side(struct Neighbour*);
//static void neigh_bb_e_to_side(Mesh *mesh, struct Neighbour*);
static int elements_common_sides(ElementFullIter ,ElementFullIter ,int[]);
static int elements_common_sides_1D(ElementFullIter ,ElementFullIter ,int[]);
static int elements_common_sides_2D(ElementFullIter ,ElementFullIter ,int[]);
static int elements_common_sides_3D(ElementFullIter ,ElementFullIter ,int[]);



//=============================================================================
//
//=============================================================================

//=============================================================================
//
//=============================================================================
void element_to_side_both(Mesh* mesh)
{
    F_ENTRY;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

	// count sides
    {
        int n_sides = 0;
        FOR_ELEMENTS(mesh, ele) n_sides+=ele->n_sides;
        mesh->sides.resize(n_sides);
        xprintf( Msg, "   Creating %d sides ...\n", n_sides);
    }

	{
	    int i_side=0;

	    FOR_ELEMENTS(mesh, ele ) {
	        for(int i_lside=0; i_lside< ele->n_sides; i_lside++) {
	            mesh->sides[i_side].reinit(mesh, ele,  i_side , i_lside);
	            ele->side[i_lside]=&( mesh->sides[i_side] );
	            i_side++;
	        }
	    }
	}
}
//=============================================================================
//
//=============================================================================
void neigh_vv_to_element(Mesh* mesh)
{
	int li, aux;
	struct Neighbour *ngh;
	ElementIter el;

	xprintf( MsgVerb, "   Neighbour of vv2 type to element... ")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

    FOR_NEIGHBOURS(mesh, ngh ) {
		if( ngh->type != VV_2E )
			continue;
		if( ngh->eid[ 0 ] > ngh->eid[ 1 ] ) {
			aux           = ngh->eid[ 1 ];
			ngh->eid[ 1 ] = ngh->eid[ 0 ];
			ngh->eid[ 0 ] = aux;
		}
		for( li = 0; li < 2; li++ ) {
		    el = mesh->element.find_id( ngh->eid[li] );
		    INPUT_CHECK( NONULL(el), "Reference to undefined element %d in neighbour %d\n", ngh->eid[ li ], ngh->id );
			ngh->element[ li ] = el;
		}
	}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}
//=============================================================================
//
//=============================================================================
/*
void element_to_neigh_vv(Mesh* mesh)
{
	struct Neighbour *ngh;
    ElementFullIter ele = ELEMENT_FULL_ITER_NULL(mesh );

	xprintf( MsgVerb, "   Element to neighbours of vv2 type... ");
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

    // Counting the neighbours using the aux variable
	FOR_ELEMENTS(mesh, ele )
		ele->aux = 0;
	FOR_NEIGHBOURS(mesh, ngh ) {
		if( ngh->type != VV_2E )
			continue;
		ele = ELEMENT_FULL_ITER(mesh, ngh->element[ 0 ]);
		ele->aux++;
		ele = ELEMENT_FULL_ITER(mesh, ngh->element[ 1 ]);
		ele->aux++;
	}
	// Allocation of the array
	FOR_ELEMENTS(mesh,  ele ) {
		ele->n_neighs_vv = ele->aux;
		if( ele->n_neighs_vv > 0 )
			ele->neigh_vv = (struct Neighbour**) xmalloc(
				ele->n_neighs_vv * sizeof( struct Neighbour* ) );
	}
	// Fill the array
	FOR_ELEMENTS(mesh,  ele )
		ele->aux = 0;
	FOR_NEIGHBOURS(mesh,  ngh ) {
		if( ngh->type != VV_2E )
			continue;
		ele = ELEMENT_FULL_ITER(mesh, ngh->element[ 0 ]);
		ele->neigh_vv[ ele->aux ] = ngh;
                if (ele->aux > 0)
                        if( ele->neigh_vv[ ele->aux - 1 ]->element[ 1 ] >
                            ele->neigh_vv[ ele->aux ]->element[ 1 ] )
                                xprintf(UsrErr,"The neighbouring elements of the element %d, ",
                                        "have to be given in order from the lowest id to the highest one. ",
                                        "Check the input file of the neighbours.", ele.id());
		ele->aux++;
		ele = ELEMENT_FULL_ITER(mesh, ngh->element[ 1 ]);
		ele->neigh_vv[ ele->aux ] = ngh;
                if (ele->aux > 0)
                        if( ele->neigh_vv[ ele->aux - 1 ]->element[ 1 ] >
                            ele->neigh_vv[ ele->aux ]->element[ 1 ] )
                                xprintf(UsrErr,"The neighbouring elements of the element %d, ",
                                        "have to be given in order from the lowest id to the highest one. ",
                                        "Check the input file of the neighbours.", ele.id());
		ele->aux++;
	}
	xprintf( MsgVerb, "O.K.\n");
}*/
//=============================================================================
//
//=============================================================================
void neigh_vb_to_element_and_side(Mesh* mesh)
{
	int li;
	struct Neighbour *ngh;
	ElementIter ele;

	xprintf( MsgVerb, "   Neighbour of vb2 type to element... ")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

    FOR_NEIGHBOURS(mesh,  ngh ) {
		if( ngh->type != VB_ES )
			continue;
		for( li = 0; li < 2; li++ ) {
            ele = mesh->element.find_id( ngh->eid[li] );
            INPUT_CHECK( NONULL(ele), "Reference to undefined element %d in neighbour %d\n", ngh->eid[ li ], ngh->id );
            ngh->element[ li ] = ele;
		}
		ele = ngh->element[ 1 ];
		INPUT_CHECK(!( (ngh->sid[ 1 ] < 0) || (ngh->sid[ 1 ] >= ele->n_sides) ),
		        "Neighbour %d has nonsensual reference to side %d\n",
				ngh->id, ngh->sid[ 1 ] );
		ngh->side[ 1 ] = ele->side[ ngh->sid[ 1 ] ];
	}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}

//=============================================================================
//
//=============================================================================
void element_to_neigh_vb(Mesh* mesh)
{
	struct Neighbour *ngh;
	ElementIter ele;

	xprintf( MsgVerb, "   Element to neighbours of vb2 type... ")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

    FOR_ELEMENTS(mesh,ele) ele->n_neighs_vb =0;
    // count vb neighs
    FOR_NEIGHBOURS(mesh,  ngh )
        if( ngh->type == VB_ES ) {
            ele = ngh->element[ 0 ];
            ele->n_neighs_vb++;
        }
	// Allocation of the array
	FOR_ELEMENTS(mesh,  ele )
		if( ele->n_neighs_vb > 0 ) {
			ele->neigh_vb = new struct Neighbour* [ele->n_neighs_vb];
			ele->n_neighs_vb=0;
		}
	// fill
	FOR_NEIGHBOURS(mesh,  ngh )
		if( ngh->type == VB_ES ) {
		    ele = ngh->element[ 0 ];
		    ele->neigh_vb[ ele->n_neighs_vb++ ] = ngh;
		}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}


//=============================================================================
//
//=============================================================================
void node_to_element(Mesh* mesh)
{
    F_ENTRY;

	int li;
	NodeIter nod;
	ElementIter ele;

	xprintf( MsgVerb, "   Node to element... ")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

	// Set counter of elements in node to zero
	FOR_NODES(mesh,  nod )
		nod->n_elements = 0;
	// Count elements
	FOR_ELEMENTS(mesh,  ele )
		FOR_ELEMENT_NODES( ele, li ) {
			nod = ele->node[ li ];
			(nod->n_elements)++;
		}
	// Allocate arrays
	FOR_NODES(mesh,  nod ) {
                if (nod->n_elements == 0)
                        continue;
       		nod->element = (ElementIter *) xmalloc( nod->n_elements * sizeof( ElementIter ) );
		nod->aux = 0;
	}
	// Set poiners in arrays
	FOR_ELEMENTS(mesh,  ele )
		FOR_ELEMENT_NODES( ele, li ) {
			nod = ele->node[ li ];
			nod->element[ nod->aux ] = ele;
			(nod->aux)++;
		}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}

//=============================================================================
//
//=============================================================================
void count_side_types(Mesh* mesh)
{
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");
    struct Side *sde;

	mesh->n_insides = 0;
	mesh->n_exsides = 0;
	FOR_SIDES(mesh,  sde )
		if (sde->is_external()) mesh->n_exsides++;
		else mesh->n_insides++;
}
/*
//=============================================================================
//
//=============================================================================
void neigh_bb_topology(Mesh* mesh)
{
	struct Neighbour *ngh;

	xprintf( MsgVerb, "   Neighbour BB to element and side... ");
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

    FOR_NEIGHBOURS(mesh,  ngh ) {
		if( ngh->type != BB_EL )
			continue;
		neigh_bb_to_element( ngh, mesh );
		neigh_bb_el_to_side( ngh );
	}
	xprintf( MsgVerb, "O.K.\n");
}

//=============================================================================
//
//=============================================================================
void neigh_bb_el_to_side( struct Neighbour *ngh )
{
	int li;
	ElementIter ele;

    ASSERT(!( ngh == NULL ),"Neighbour is NULL\n");

    FOR_NEIGH_ELEMENTS( ngh, li ) {
		ele = ngh->element[ li ];
		INPUT_CHECK(!( (ngh->sid[ li ] < 0) || (ngh->sid[ li ] >= ele->n_sides) ),
		        "Neighbor %d has incorrecst reference to side %d\n", ngh->id, ngh->sid[ li ] );
		ngh->side[ li ] = ele->side[ ngh->sid[ li ] ];
	}
}
*/

//=============================================================================
//
//=============================================================================
void edge_to_side_both(Mesh* mesh)
{
	struct Side *sde;
	struct Edge *edg;
	struct Neighbour *ngh;
	Element *ele;
    int si;

	xprintf( MsgVerb, "   Edge to side and back... \n")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

    // count edges
    unsigned int n_edges = mesh->n_sides();
    FOR_NEIGHBOURS(mesh,  ngh )
        if ( ngh->type == BB_EL ) n_edges -= ( ngh->n_elements - 1 );

    // create edge vector
    mesh->edge.resize(n_edges);
    xprintf( MsgLog, "Created  %d edges.\n.", mesh->n_edges() );

    // set edge, side connections
    unsigned int i_edge=0;
	FOR_NEIGHBOURS(mesh,  ngh ) if(  ngh->type == BB_EL ) {
	    edg = &( mesh->edge[i_edge++] );

        // init edge (can init all its data)
		edg->n_sides = ngh->n_sides;
		edg->side = (struct Side**) xmalloc( edg->n_sides *
				sizeof( struct Side* ) );
		FOR_NEIGH_SIDES( ngh, si ) {
		    ele = mesh->element.find_id( ngh->eid[si] );
			sde = ele->side[ ngh->sid[ si ] ];
			edg->side[ si ] = sde;
			sde->edge_ = edg;
		}
	}

	// now the external ones ( pair all remaining edges with external sides)
	FOR_SIDES(mesh, sde)
	    if ( sde->edge() == NULL ) {
	        edg = &( mesh->edge[i_edge++] );

	        // make external edges and edges on neighborings.
	        edg->n_sides = 1;
	        edg->side = (struct Side**) xmalloc( edg->n_sides *
	                    sizeof( struct Side* ) );

	        sde->edge_ = edg;
	        edg->side[ 0 ] = &( *sde );

	    }

	ASSERT(i_edge == mesh->n_edges(), "Actual number of edges %d do not match size of its array.\n", i_edge);

    //FOR_SIDES(mesh, side) ASSERT(side->edge != NULL, "Empty side %d !\n", side->id);
}
//=============================================================================
//
//=============================================================================
void neigh_vb_to_edge_both(Mesh* mesh)
{
	struct Side *sde;
	struct Edge *edg;
	struct Neighbour *ngh;

	xprintf( MsgVerb, "   Neighbour VB to edge and back... ")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

	FOR_NEIGHBOURS(mesh,  ngh ) {
		if( ngh->type != VB_ES )
			continue;
		sde = ngh->side[ 1 ];
		edg = sde->edge();
		edg->neigh_vb = ngh;
		ngh->edge = edg;
	}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}

//-----------------------------------------------------------------------------
// vim: set cindent:
