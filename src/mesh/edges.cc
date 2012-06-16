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
 * @brief ???
 * @ingroup mesh
 *
 */

#include "system/system.hh"
#include "mesh/mesh.h"

static struct Edge *new_edge(void);
static void add_to_edge_list(Mesh*,struct Edge*);
static void init_edge(struct Edge*);
static int count_edges(Mesh*);

static int number_of_common_nodes_ss(struct Side*,struct Side*);
static int number_of_common_nodes_se(struct Side*,ElementIter );


Edge::Edge()
: n_sides(NDEF),
  side(NULL),
  neigh_vb(NULL),
  neigh_bb(NULL),
  c_row(0),
  f_val(0),
  f_rhs(0)
{

}


//=============================================================================
// CREATE AND PREFILL LIST OF EDGES
//=============================================================================
void make_edge_list(Mesh* mesh)
{
    F_ENTRY;
	int edi;
	struct Edge *edg;

	ASSERT(!( mesh == NULL ),"NULL as argument of function make_edge_list()\n");
	xprintf( Msg, "Creating edges... ")/*orig verb 2*/;
	mesh->edge.resize(count_edges( mesh ));

	xprintf( MsgVerb, " O.K. %d edges created.", mesh->n_edges() )/*orig verb 4*/;
	xprintf( Msg, "O.K.\n")/*orig verb 2*/;
}
/*
//=============================================================================
//
//=============================================================================
int number_of_common_nodes_ss( struct Side *s0, struct Side *s1 )
{
	int i, j;
	int ncn;

	ASSERT(!( (s0 == NULL) || (s1 == NULL) ),"NULL argument Side\n");
	ncn = 0;
	FOR_SIDE_NODES( s0, i )
		FOR_SIDE_NODES( s1, j )
			if( s0->node[ i ] == s1->node[ j ] )
				ncn++;
	return ncn;
}
//=============================================================================
//
//=============================================================================
int number_of_common_nodes_se( struct Side *sde, ElementIter ele )
{
	int i, j;
	int ncn;

	INPUT_CHECK(!( (sde == NULL) || (ele == NULL) ),"NULL argument\n");
	ncn = 0;
	FOR_SIDE_NODES( sde, i )
		FOR_ELEMENT_NODES( ele, j )
			if( sde->node[ i ] == ele->node[ j ] )
				ncn++;
	return ncn;
}*/
//=============================================================================
//
//=============================================================================
int count_edges(Mesh* mesh)
{
	int rc;
	struct Neighbour *ngh;

	rc = mesh->n_sides();
	//cout << "n_edges: " << rc << endl;
	FOR_NEIGHBOURS(mesh,  ngh ) {
		if( ngh->type == BB_E || ngh->type == BB_EL )
			rc -= ( ngh->n_elements - 1 );
	}
	//cout << "n_edges: " << rc << endl;
	return rc;
}


//=============================================================================
// CALCULATE PROPERTIES OF ALL EDGES OF THE MESH
//=============================================================================
void edge_calculation_mh(Mesh* mesh)
{
	int edi;
	struct Edge *edg;

	xprintf( Msg, "Calculating properties of edges... ")/*orig verb 2*/;
	ASSERT(!( mesh == NULL ),"NULL as argument of function edge_calculation_mh()\n");
	edi = 0;
	FOR_EDGES(mesh,  edg ) {
		edg->c_row = mesh->n_sides() + mesh->n_elements() + edi;
		edg->f_rhs=0.0;
		if( edg->neigh_vb == NULL )
			edg->f_val = 0.0;
		else
			edg->f_val = -1.0 * edg->neigh_vb->sigma * edg->side[0]->metric();
		edi++;
	}
	xprintf( Msg, "O.K.\n")/*orig verb 2*/;
}
//-----------------------------------------------------------------------------
// vim: set cindent:
