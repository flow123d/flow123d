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
 * @file   hashes.cc
 * @brief  Creation of pseudo-hash structures
 */

#include "system.hh"
#include "hashes.h"
#include "mesh.h"
#include "problem.h"
//#include "materials.hh"
#include "boundaries.h"
#include "neighbours.h"

#include "concentrations.h"
#include "transport_bcd.h"
#include "edges.h"
#include "sides.h"

#include "constantdb.h"
#include "mesh/ini_constants_mesh.hh"

static void create_boundary_hash(Mesh*);
static void create_neighbour_hash(Mesh*);
static void create_source_hash(Mesh*);
static void create_concentration_hash(Mesh*);
static void create_transport_bcd_hash(Mesh*);

//=============================================================================
//
//=============================================================================
void make_hashes( struct Problem *problem )
{
	xprintf( MsgVerb, "Generating hash tables...\n")/*orig verb 5*/;

        Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

	{
	Edge *edg;
	Side *side;

	// make edge hash
	mesh->max_edg_id=mesh->n_edges-1;
	mesh->edge_hash=(Edge **)malloc(mesh->n_edges*sizeof(Edge *));
	FOR_EDGES(edg) mesh->edge_hash[edg->id]=edg;
	// make side hash
	mesh->max_side_id=mesh->n_sides-1;
	mesh->side_hash=(Side **)malloc(mesh->n_sides*sizeof(Side *));
	FOR_SIDES(side) mesh->side_hash[side->id]=side;
	}

//	create_boundary_hash( mesh );
    create_neighbour_hash( mesh );
//	if( mesh->source != NULL )
//		create_source_hash( mesh );
    if( mesh->concentration != NULL) {
        create_concentration_hash( mesh );
        create_transport_bcd_hash( mesh );
    }
	xprintf( MsgVerb, "Hashes generated O.K.\n")/*orig verb 5*/;
}
#if 0
//=============================================================================
//
//=============================================================================
void create_boundary_hash(Mesh* mesh)
{
	int hi;
	struct Boundary *bcd;

	xprintf( MsgVerb, "   Boundary hash... ")/*orig verb 5*/;
	// Get maximal id
	mesh->max_bou_id = NDEF;
	FOR_BOUNDARIES( bcd )
		if( bcd->id > mesh->max_bou_id )
			mesh->max_bou_id = bcd->id;
	ASSERT(!( mesh->max_bou_id == NDEF ),"Maximal boundary index < 0\n");
	// Alloc and clear
	mesh->boundary_hash = (struct Boundary**) xmalloc(
			  ( mesh->max_bou_id + 1 ) * sizeof( struct Boundary* ) );
	for( hi = 0; hi <= mesh->max_bou_id; hi++ )
		mesh->boundary_hash[ hi ] = NULL;
	// Fill hash
	FOR_BOUNDARIES( bcd ) {
		INPUT_CHECK(!( mesh->boundary_hash[ bcd->id ] != NULL ),"Multiple definition of boundary %d\n", bcd->id );
		mesh->boundary_hash[ bcd->id ] = bcd;
	}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}
#endif
//=============================================================================
//
//=============================================================================
void create_neighbour_hash(Mesh* mesh)
{
	int hi;
	struct Neighbour *ngh;

	xprintf( MsgVerb, "   Neighbour hash... ")/*orig verb 5*/;
	// Get maximal id
	mesh->max_ngh_id = NDEF;
	FOR_NEIGHBOURS( ngh )
		if( ngh->id > mesh->max_ngh_id )
			mesh->max_ngh_id = ngh->id;
	ASSERT(!( mesh->max_ngh_id == NDEF ),"Maximal neighbour index < 0\n");
	// Alloc and clear
	mesh->neighbour_hash = (struct Neighbour**) xmalloc(
			  ( mesh->max_ngh_id + 1 ) * sizeof( struct Neighbour* ) );
	for( hi = 0; hi <= mesh->max_ngh_id; hi++ )
		mesh->neighbour_hash[ hi ] = NULL;
	// Fill hash
	FOR_NEIGHBOURS( ngh ) {
		INPUT_CHECK(!( mesh->neighbour_hash[ ngh->id ] != NULL ),"Multiple definition of neighbour %d\n", ngh->id );
		mesh->neighbour_hash[ ngh->id ] = ngh;
	}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}
//=============================================================================
//
//=============================================================================
/*
void create_source_hash(Mesh* mesh)
{
	int hi;
	struct Source *src;

	xprintf( MsgVerb, "   Neighbour hash... ")
	// Get maximal id
	mesh->max_src_id = NDEF;
	FOR_SOURCES( src )
		if( src->id > mesh->max_src_id )
			mesh->max_src_id = src->id;
	ASSERT(!( mesh->max_src_id == NDEF ),"Maximal source index < 0\n");
	// Alloc and clear
	mesh->source_hash = (struct Source**) xmalloc(
			  ( mesh->max_src_id + 1 ) * sizeof( struct Source* ) );
	for( hi = 0; hi <= mesh->max_src_id; hi++ )
		mesh->source_hash[ hi ] = NULL;
	// Fill hash
	FOR_SOURCES( src ) {
		INPUT_CHECK(!( mesh->source_hash[ src->id ] != NULL ),"Multiple definition of source %d\n", src->id );
		mesh->source_hash[ src->id ] = src;
	}
	xprintf( MsgVerb, "O.K.\n");
}*/
//=============================================================================
//
//=============================================================================
void create_concentration_hash(Mesh* mesh)
{
	int hi;
	struct Concentration *con;

	xprintf( MsgVerb, "   Concentration hash... ")/*orig verb 5*/;
	// Get maximal id
	mesh->max_con_id = NDEF;
	FOR_CONCENTRATIONS( con )
		if( con->id > mesh->max_con_id )
			mesh->max_con_id = con->id;
	ASSERT(!( mesh->max_con_id == NDEF ),"Maximal concentration index < 0\n");
	// Alloc and clear
	mesh->concentration_hash = (struct Concentration**) xmalloc(
			  ( mesh->max_con_id + 1 ) * sizeof( struct Concentration* ) );
	for( hi = 0; hi <= mesh->max_con_id; hi++ )
		mesh->concentration_hash[ hi ] = NULL;
	// Fill hash
	FOR_CONCENTRATIONS( con ) {
		INPUT_CHECK(!( mesh->concentration_hash[ con->id ] != NULL ),"Multiple definition of concentration %d\n", con->id );
		mesh->concentration_hash[ con->id ] = con;
	}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}
//=============================================================================
//
//=============================================================================
void create_transport_bcd_hash(Mesh* mesh)
{
	int hi;
	struct Transport_bcd *tbc;

	xprintf( MsgVerb, "   Concentration hash... ")/*orig verb 5*/;
	// Get maximal id
	mesh->max_tbc_id = NDEF;
	FOR_TRANSPORT_BCDS( tbc )
		if( tbc->id > mesh->max_tbc_id )
			mesh->max_tbc_id = tbc->id;
	ASSERT(!( mesh->max_tbc_id == NDEF ),"Maximal transport boundary index < 0\n");
	// Alloc and clear
	mesh->transport_bcd_hash = (struct Transport_bcd**) xmalloc(
			  ( mesh->max_tbc_id + 1 ) * sizeof( struct Transport_bcd* ) );
	for( hi = 0; hi <= mesh->max_tbc_id; hi++ )
		mesh->transport_bcd_hash[ hi ] = NULL;
	// Fill hash
	FOR_TRANSPORT_BCDS( tbc ) {
		INPUT_CHECK((mesh->transport_bcd_hash[ tbc->id ]==NULL),"Multiple definition of transport_bcd %d\n", tbc->id );
		mesh->transport_bcd_hash[ tbc->id ] = tbc;
	}
	xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}
//-----------------------------------------------------------------------------
// vim: set cindent:
