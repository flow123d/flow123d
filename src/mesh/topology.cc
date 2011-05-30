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
static void neigh_bv_to_side(Mesh *);
static void element_to_neigh_vb(Mesh*);
static void side_to_node(Mesh*);
static void neigh_bb_topology(Mesh*);
static void neigh_bb_to_edge_both(Mesh*);
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
static void neigh_bb_e_to_side(Mesh *mesh, struct Neighbour*);
//static void source_to_element_both(Mesh*);
//static void concentration_to_element(Mesh*);
//static void transport_bcd_to_boundary(Mesh*);
//static void initial_to_element(Mesh*);
static int elements_common_sides(ElementFullIter ,ElementFullIter ,int[]);
static int elements_common_sides_1D(ElementFullIter ,ElementFullIter ,int[]);
static int elements_common_sides_2D(ElementFullIter ,ElementFullIter ,int[]);
static int elements_common_sides_3D(ElementFullIter ,ElementFullIter ,int[]);



//=============================================================================
//
//=============================================================================
void element_to_material(Mesh* mesh, MaterialDatabase &base)
{

	xprintf( MsgVerb, "   Element to material... ")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");
	FOR_ELEMENTS( ele ) {
	    ele->material=base.find_id(ele->mid);
		INPUT_CHECK( ele->material != base.end(),
		        "Reference to undefined material %d in element %d\n", ele->mid, ele.id() );
	}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}
//=============================================================================
//
//=============================================================================
void element_to_side_both(Mesh* mesh)
{
	int li;

	struct Side *sde;

	xprintf( MsgVerb, "   Element to side and back... ")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");
	sde = mesh->side;
	FOR_ELEMENTS( ele ) {
		FOR_ELEMENT_SIDES( ele, li ) {
			ASSERT(!( sde == NULL ),"Inconsistency between number of elements and number of sides\n");
			ele->side[ li ] = sde;
			sde->element = ele;
			sde->lnum    = li;
			sde = sde->next;
		}
	}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
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

    FOR_NEIGHBOURS( ngh ) {
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
void element_to_neigh_vv(Mesh* mesh)
{
	struct Neighbour *ngh;
    ElementFullIter ele = ELEMENT_FULL_ITER_NULL;

	xprintf( MsgVerb, "   Element to neighbours of vv2 type... ")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

    // Counting the neighbours using the aux variable
	FOR_ELEMENTS( ele )
		ele->aux = 0;
	FOR_NEIGHBOURS( ngh ) {
		if( ngh->type != VV_2E )
			continue;
		ele = ELEMENT_FULL_ITER(ngh->element[ 0 ]);
		ele->aux++;
		ele = ELEMENT_FULL_ITER(ngh->element[ 1 ]);
		ele->aux++;
	}
	// Allocation of the array
	FOR_ELEMENTS( ele ) {
		ele->n_neighs_vv = ele->aux;
		if( ele->n_neighs_vv > 0 )
			ele->neigh_vv = (struct Neighbour**) xmalloc(
				ele->n_neighs_vv * sizeof( struct Neighbour* ) );
	}
	// Fill the array
	FOR_ELEMENTS( ele )
		ele->aux = 0;
	FOR_NEIGHBOURS( ngh ) {
		if( ngh->type != VV_2E )
			continue;
		ele = ELEMENT_FULL_ITER(ngh->element[ 0 ]);
		ele->neigh_vv[ ele->aux ] = ngh;
                if (ele->aux > 0)
                        if( ele->neigh_vv[ ele->aux - 1 ]->element[ 1 ] >
                            ele->neigh_vv[ ele->aux ]->element[ 1 ] )
                                xprintf(UsrErr,"The neighbouring elements of the element %d, ",
                                        "have to be given in order from the lowest id to the highest one. ",
                                        "Check the input file of the neighbours.", ele.id());
		ele->aux++;
		ele = ELEMENT_FULL_ITER(ngh->element[ 1 ]);
		ele->neigh_vv[ ele->aux ] = ngh;
                if (ele->aux > 0)
                        if( ele->neigh_vv[ ele->aux - 1 ]->element[ 1 ] >
                            ele->neigh_vv[ ele->aux ]->element[ 1 ] )
                                xprintf(UsrErr,"The neighbouring elements of the element %d, ",
                                        "have to be given in order from the lowest id to the highest one. ",
                                        "Check the input file of the neighbours.", ele.id());
		ele->aux++;
	}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}
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

    FOR_NEIGHBOURS( ngh ) {
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
void neigh_bv_to_side(Mesh* mesh)
{
	struct Neighbour *ngh;

    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

	FOR_NEIGHBOURS( ngh ) {
		if( ngh->type != VB_ES )
			continue;
                ngh->side[1]->neigh_bv = ngh;
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

    // Counting the neighbours using the aux variable
	FOR_ELEMENTS( ele )
		ele->aux = 0;
	FOR_NEIGHBOURS( ngh ) {
		if( ngh->type != VB_ES )
			continue;
		ele = ngh->element[ 0 ];
		ele->aux++;
	}
	// Allocation of the array
	FOR_ELEMENTS( ele ) {
		ele->n_neighs_vb = ele->aux;
		if( ele->n_neighs_vb > 0 )
			ele->neigh_vb = (struct Neighbour**) xmalloc(
				ele->n_neighs_vb * sizeof( struct Neighbour* ) );
	}
	// Fill the array
	FOR_ELEMENTS( ele )
		ele->aux = 0;
	FOR_NEIGHBOURS( ngh ) {
		if( ngh->type != VB_ES )
			continue;
		ele = ngh->element[ 0 ];
		ele->neigh_vb[ ele->aux ] = ngh;
		ele->aux++;
	}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}
//=============================================================================
//
//=============================================================================
void side_to_node(Mesh* mesh)
{
	ElementIter ele;

	xprintf( MsgVerb, "   Side to node... ")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

    FOR_ELEMENTS( ele )
		switch( ele->type ) {
			case LINE:
				side_to_node_line( ele );
				break;
			case TRIANGLE:
				side_to_node_triangle( ele );
				break;
			case TETRAHEDRON:
				side_to_node_tetrahedron( ele );
				break;
		}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}
//=============================================================================
//
//=============================================================================
void side_to_node_line( ElementFullIter ele )
{
	struct Side *sde;

    ASSERT(!( ele == NULL ),"element is NULL\n");

	sde = ele->side[ 0 ];
	sde->node[ 0 ] = ele->node[ 0 ];
	sde = ele->side[ 1 ];
	sde->node[ 0 ] = ele->node[ 1 ];
}
//=============================================================================
//
//=============================================================================
void side_to_node_triangle( ElementFullIter ele )
{
	struct Side *sde;

    ASSERT(!( ele == NULL ),"element is NULL\n");

    sde = ele->side[ 0 ];
	sde->node[ 0 ] = ele->node[ 0 ];
	sde->node[ 1 ] = ele->node[ 1 ];
	sde = ele->side[ 1 ];
	sde->node[ 0 ] = ele->node[ 1 ];
	sde->node[ 1 ] = ele->node[ 2 ];
	sde = ele->side[ 2 ];
	sde->node[ 0 ] = ele->node[ 2 ];
	sde->node[ 1 ] = ele->node[ 0 ];
}
//=============================================================================
//
//=============================================================================
void side_to_node_tetrahedron( ElementFullIter ele )
{
	struct Side *sde;

    ASSERT(!( ele == NULL ),"element is NULL\n");

    sde = ele->side[ 0 ];
	sde->node[ 0 ] = ele->node[ 1 ];
	sde->node[ 1 ] = ele->node[ 2 ];
	sde->node[ 2 ] = ele->node[ 3 ];
	sde = ele->side[ 1 ];
	sde->node[ 0 ] = ele->node[ 0 ];
	sde->node[ 1 ] = ele->node[ 2 ];
	sde->node[ 2 ] = ele->node[ 3 ];
	sde = ele->side[ 2 ];
	sde->node[ 0 ] = ele->node[ 0 ];
	sde->node[ 1 ] = ele->node[ 1 ];
	sde->node[ 2 ] = ele->node[ 3 ];
	sde = ele->side[ 3 ];
	sde->node[ 0 ] = ele->node[ 0 ];
	sde->node[ 1 ] = ele->node[ 1 ];
	sde->node[ 2 ] = ele->node[ 2 ];
}


//=============================================================================
//
//=============================================================================
void node_to_element(Mesh* mesh)
{
	int li;
	NodeIter nod;
	ElementIter ele;

	xprintf( MsgVerb, "   Node to element... ")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

	// Set counter of elements in node to zero
	FOR_NODES( nod )
		nod->n_elements = 0;
	// Count elements
	FOR_ELEMENTS( ele )
		FOR_ELEMENT_NODES( ele, li ) {
			nod = ele->node[ li ];
			(nod->n_elements)++;
		}
	// Allocate arrays
	FOR_NODES( nod ) {
                if (nod->n_elements == 0)
                        continue;
       		nod->element = (ElementIter *) xmalloc( nod->n_elements * sizeof( ElementIter ) );
		nod->aux = 0;
	}
	// Set poiners in arrays
	FOR_ELEMENTS( ele )
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
void side_types(Mesh* mesh)
{
	struct Side *sde;

	xprintf( MsgVerb, "   Side types... ")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

	FOR_SIDES( sde )
		sde->type = ( sde->edge->n_sides == 1 ? EXTERNAL : INTERNAL );
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}
//=============================================================================
//
//=============================================================================
void count_side_types(Mesh* mesh)
{
	struct Side *sde;

	xprintf( MsgVerb, "   Counting side types... ")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

	mesh->n_insides = 0;
	mesh->n_exsides = 0;
	FOR_SIDES( sde )
		switch( sde->type ) {
			case INTERNAL:
				mesh->n_insides++;
				break;
			case EXTERNAL:
				mesh->n_exsides++;
				break;
			default:
				xprintf(PrgErr,"Unknown type of side - %d in side %d\n",
					   sde->type, sde->id );
				break;
		}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}
//=============================================================================
//
//=============================================================================
void neigh_bb_topology(Mesh* mesh)
{
	struct Neighbour *ngh;

	xprintf( MsgVerb, "   Neighbour BB to element and side... ")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

    FOR_NEIGHBOURS( ngh ) {
		if( ngh->type != BB_E && ngh->type != BB_EL )
			continue;
		neigh_bb_to_element( ngh, mesh );
		if( ngh->type == BB_EL )
			neigh_bb_el_to_side( ngh );
		else
			neigh_bb_e_to_side( mesh, ngh );
	}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}
//=============================================================================
//
//=============================================================================
void neigh_bb_to_element(struct Neighbour* ngh, Mesh* mesh)
{
	int li;
	ElementIter ele;

    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");
    ASSERT(!( ngh == NULL ),"Neighbor is NULL\n");

	FOR_NEIGH_ELEMENTS( ngh, li ) {
	    ele = mesh->element.find_id( ngh->eid[li] );
	    INPUT_CHECK( NONULL(ele), "Reference to undefined element %d neighbor %d\n", ngh->eid[ li ], ngh->id );
		ngh->element[ li ] = ele;
	}
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
//=============================================================================
//
//=============================================================================
void neigh_bb_e_to_side(Mesh *mesh, struct Neighbour *ngh )
{
	int li, s[ 2 ];
	ElementFullIter e0 = ELEMENT_FULL_ITER_NULL;
    ElementFullIter e1 = ELEMENT_FULL_ITER_NULL;

    ASSERT(!( ngh == NULL ),"Neighbour is NULL\n");

    FOR_NEIGH_ELEMENTS( ngh, li ) {
		if( li == 0 ) {
			e0 = ELEMENT_FULL_ITER(ngh->element[ 0 ]);
			continue;
		}
		e1 = ELEMENT_FULL_ITER(ngh->element[ li ]);
		if( elements_common_sides( e0, e1, s ) == false )
		       xprintf(UsrErr,"In neighbour %d elements %d and %d do not have common side\n",
				  ngh->id, e0.id(), e1.id() );
		if( li == 1 ) {
			e0->aux = s[ 0 ];
			ngh->side[ 0 ] = e0->side[ s [ 0 ] ];
		} else {
			INPUT_CHECK( s[ 0 ] == e0->aux ,"Cannot find free side for edge %d\n", ngh->id); // TODO message ??
		}
		ngh->side[ li ] = e1->side[ s [ 1 ] ];
	}
}
//=============================================================================
//
//=============================================================================
int elements_common_sides( ElementFullIter e0, ElementFullIter e1, int s[] )
{
	int rc;

    ASSERT(!( (e0 == NULL) || (e1 == NULL) || (s == NULL) ),"NULL argument in elements_common_sides()\n");

	s[ 0 ] = NDEF;
	s[ 1 ] = NDEF;
	if( e0->dim != e1->dim )
		return false;
	switch( e0->dim ) {
		case 1:
			rc = elements_common_sides_1D( e0, e1, s );
			break;
		case 2:
			rc = elements_common_sides_2D( e0, e1, s );
			break;
		case 3:
			rc = elements_common_sides_3D( e0, e1, s );
			break;
	}
	return rc;
}
//=============================================================================
//
//=============================================================================
int elements_common_sides_1D( ElementFullIter e0, ElementFullIter e1, int s[] )
{
    ASSERT(!( (e0 == NULL) || (e1 == NULL) || (s == NULL) ),"NULL argument!\n");

    if( e0->node[ 0 ] == e1->node[ 0 ] ) {
		s[ 0 ] = 0;
		s[ 1 ] = 0;
		return true;
	}
	if( e0->node[ 0 ] == e1->node[ 1 ] ) {
		s[ 0 ] = 0;
		s[ 1 ] = 1;
		return true;
	}
	if( e0->node[ 1 ] == e1->node[ 0 ] ) {
		s[ 0 ] = 1;
		s[ 1 ] = 0;
		return true;
	}
	if( e0->node[ 1 ] == e1->node[ 1 ] ) {
		s[ 0 ] = 1;
		s[ 1 ] = 1;
		return true;
	}
	return false;
}
//=============================================================================
//
//=============================================================================
int elements_common_sides_2D( ElementFullIter e0, ElementFullIter e1, int s[] )
{
	int i, j, cnt, ieq[ 2 ], jeq[ 2 ], a;

    ASSERT(!( (e0 == NULL) || (e1 == NULL) || (s == NULL) ),"NULL argument!\n");

    cnt = 0;
	for( i = 0; i < 3; i++ )
		for( j = 0; j < 3; j++ )
			if( e0->node[ i ] == e1->node[ j ] ) {
				if( cnt < 2 ) {
					ieq[ cnt ] = i;
					jeq[ cnt ] = j;
				}
				cnt++;
			}
	if( cnt != 2 )
		return false;
	// Order ieq[] and jeq[]
	if( ieq[ 0 ] > ieq[ 1 ] ) {
		a        = ieq[ 1 ];
		ieq[ 1 ] = ieq[ 0 ];
		ieq[ 0 ] = a;
	}
	if( jeq[ 0 ] > jeq[ 1 ] ) {
		a        = jeq[ 1 ];
		jeq[ 1 ] = jeq[ 0 ];
		jeq[ 0 ] = a;
	}
	// This is little tricky, so short explanation:
	// ieq[] is ordered, so there are three possible states:
	// ieq[ 0 ]      ieq[ 1 ]      side in triangle
	//   0              1		   0
	//   0              2		   2
	//   1              2		   1
	// The same holds for jeq[], of course.
	s[ 0 ] = ieq[ 0 ] == 1 ? 1 : ( ieq[ 1 ] == 1 ? 0 : 2 );
	s[ 1 ] = jeq[ 0 ] == 1 ? 1 : ( jeq[ 1 ] == 1 ? 0 : 2 );
	return true;
}
//=============================================================================
//
//=============================================================================
int elements_common_sides_3D( ElementFullIter e0, ElementFullIter e1, int s[] )
{
	int i, j, cnt, ieq[ 3 ], jeq[ 3 ], a[ 4 ];

    ASSERT(!( (e0 == NULL) || (e1 == NULL) || (s == NULL) ),"NULL argument!\n");

    cnt = 0;
	for( i = 0; i < 4; i++ )
		for( j = 0; j < 4; j++ )
			if( e0->node[ i ] == e1->node[ j ] ) {
				if( cnt < 3 ) {
					ieq[ cnt ] = i;
					jeq[ cnt ] = j;
				}
				cnt++;
			}
	if( cnt != 3 )
		return false;
	// Number of side is the missing number from set {0,1,2,3}
	// in arrays ieq[], jeq[]...
	// First element
	for( i = 0; i < 4; i++ )
		a[ i ] = 1;
	for( i = 0; i < 3; i++ )
		a[ ieq[ i ] ] = 0;
	for( i = 0; i < 4; i++ )
		if( a[ i ] == 1 )
			s[ 0 ] = i;
	// Second element
	for( i = 0; i < 4; i++ )
		a[ i ] = 1;
	for( i = 0; i < 3; i++ )
		a[ jeq[ i ] ] = 0;
	for( i = 0; i < 4; i++ )
		if( a[ i ] == 1 )
			s[ 1 ] = i;
	return true;
}
//=============================================================================
//
//=============================================================================
void neigh_bb_to_edge_both(Mesh* mesh)
{
	struct Neighbour *ngh;

	xprintf( MsgVerb, "   Neighbour BB to edge and back... ")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

    EdgeFullIter edg = mesh->edge.begin();
	FOR_NEIGHBOURS( ngh ) {
		if( ngh->type != BB_E && ngh->type != BB_EL )
			continue;
		ngh->edge = edg;
		edg->neigh_bb = ngh;
		++edg;
		ASSERT( edg != mesh->edge.end() ,"Inconsistency between number of neighbours and number of edges\n");
	}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}
//=============================================================================
//
//=============================================================================
void edge_to_side_both(Mesh* mesh)
{
	struct Side *sde;
	struct Edge *edg;
	struct Neighbour *ngh;
    int si;

	xprintf( MsgVerb, "   Edge to side and back... ")/*orig verb 5*/;
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

    // first, the internal sides
	FOR_NEIGHBOURS( ngh ) {
		if( ngh->type != BB_E && ngh->type != BB_EL )
			continue;
		edg = ngh->edge;
		edg->n_sides = ngh->n_sides;
		edg->side = (struct Side**) xmalloc( edg->n_sides *
				sizeof( struct Side* ) );
		FOR_NEIGH_SIDES( ngh, si ) {
			sde = ngh->side[ si ];
			edg->side[ si ] = sde;
			sde->edge = edg;
		}
	}
	// now the external ones
	sde = mesh->side;
	FOR_EDGES( edg ) {
		if( edg->n_sides != NDEF )
			continue;
		edg->n_sides = 1;
		edg->side = (struct Side**) xmalloc( edg->n_sides *
					sizeof( struct Side* ) );
		while( sde->edge != NULL ) {
			sde = sde->next;
			INPUT_CHECK(!( sde == NULL ),"No next side during external edge initialization!");
		}
		sde->edge = edg;
		edg->side[ 0 ] = sde;
	}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;

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

	FOR_NEIGHBOURS( ngh ) {
		if( ngh->type != VB_ES )
			continue;
		sde = ngh->side[ 1 ];
		edg = sde->edge;
		edg->neigh_vb = ngh;
		ngh->edge = edg;
	}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}
/*
//=============================================================================
//
//=============================================================================
void concentration_to_element(Mesh* mesh)
{
	struct Concentration *con;
	ElementFullIter ele = ELEMENT_FULL_ITER_NULL;

	xprintf( MsgVerb, "   Concentration to element... ");
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

    FOR_CONCENTRATIONS( con ) {
        ele = mesh->element.find_id( con->eid );
		INPUT_CHECK( ele.is_valid(), "Reference to undefined element in concentration %d\n", con->id );
		ele->start_conc = con;
	}
	xprintf( MsgVerb, "O.K.\n");
}

void transport_bcd_to_boundary(Mesh* mesh)
{
	struct Transport_bcd *tbc;
	BoundaryFullIter bcd = BOUNDARY_NULL;

	xprintf( MsgVerb, "   Transport boundary condition to boundary condition... ");
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");

    FOR_TRANSPORT_BCDS( tbc ) {
		bcd = mesh->boundary.find_id( tbc->bid );
        INPUT_CHECK( bcd.is_valid(), "Reference to undefined boundary condition in transport BCDs %d\n", tbc->id );
		bcd->transport_bcd = tbc;
	}
	xprintf( MsgVerb, "O.K.\n");
}*/
//=============================================================================
//
//=============================================================================
/*
void initial_to_element(Mesh* mesh)
{
	struct Initial *icd;
	ElementFullIter ele=ELEMENT_FULL_ITER_NULL;

	xprintf( MsgVerb, "   Initial conditions to element... ")
    ASSERT(!( mesh == NULL ),"Mesh is NULL\n");
    INPUT_CHECK(!(mesh->n_initials != mesh->n_elements()),"Different number of initial conditions and elements\n");

    FOR_INITIALS( icd ) {
        ele = mesh->element.find_id( icd->id );
        INPUT_CHECK( NONULL(ele), "Reference to undefined element in initial condition %d\n", icd->id );

        INPUT_CHECK(!(ele->n_sides != icd->n_sides),
                "Different number of sides in the mesh and in the initial condition. Element %d\n", ele.id());
        INPUT_CHECK(!(ele->n_neighs_vv != icd->n_neighs_vv),
                "Different number of VV neighbours in the mesh and in the initial condition. Element %d\n", ele.id());
        INPUT_CHECK(!(ele->initial != NULL),
                "Element %d has more than one initial conditions\n", ele.id());
        ele->initial = icd;
    }
	xprintf( MsgVerb, "O.K.\n")
}*/
//-----------------------------------------------------------------------------
// vim: set cindent:
