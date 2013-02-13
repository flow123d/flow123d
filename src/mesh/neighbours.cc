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
 * @brief  Initialize neighbouring
 *
 */

#include "system/system.hh"
#include "neighbours.h"
#include "mesh/mesh.h"
/*
static void parse_neighbour_line(struct Neighbour*,char*);
static char supported_neighbour_type(int);
static void neighbour_type_specific(struct Neighbour* , char * line );
//static void neighbour_specs_bb_e(struct Neighbour*);
static void neighbour_specs_bb_el(struct Neighbour*, char * line );
static void neighbour_specs_vb_es(struct Neighbour*, char * line );
static void neighbour_specs_vv_2e(struct Neighbour*);
*/

//=============================================================================
// READ DATA OF ALL NEIGHBOURS
//=============================================================================
void read_neighbour_list(Mesh* mesh)
{

}

//=============================================================================
// INIT DATA OF NEIGHBOUR
//=============================================================================
/*
Neighbour_both::Neighbour_both()
{
	type        = NDEF;
	n_sides     = NDEF;
	sid  	     = NULL;
	eid         = NULL;

	sigma        = 0.0;
}
*/

Neighbour::Neighbour()
: edge_idx_(-1)
{}

void Neighbour::reinit(ElementIter ele, unsigned int edg_idx, double sigma_in)
{
    element_=ele;
    edge_idx_=edg_idx;
    sigma = sigma_in;
}



/**
 *  boundary - boundary neigbouring, given by list of neigbouring elements,
 *             sides are given implicitely by shared nodes
 */
/*
void neighbour_specs_bb_e( struct Neighbour *ngh )
{
	int ei;

	xstrtok( ngh->line);
	xstrtok( NULL);
	ngh->n_elements = atoi( xstrtok( NULL) );
	INPUT_CHECK(!( ngh->n_elements < 2 ),"Neighbour %d has bad number of elements: %d\n", ngh->id, ngh->n_elements );
	ngh->n_sides = ngh->n_elements;
	ngh->eid = (int*) xmalloc( ngh->n_elements * sizeof( int ) );
	ngh->element = (ElementIter *) xmalloc( ngh->n_elements * sizeof( ElementIter  ) );
	ngh->side = (struct Side**) xmalloc( ngh->n_elements * sizeof( struct Side* ) );
	FOR_NEIGH_ELEMENTS( ngh, ei ) {
		ngh->element[ ei ] = NULL;
		ngh->side[ ei ] = NULL;
	}
	FOR_NEIGH_ELEMENTS( ngh, ei )
		ngh->eid[ ei ] = atoi( xstrtok( NULL) );
	xfree( ngh->line );
	ngh->line = NULL;
}*/
/**
 * volume - volume neighbouring of different dimensions (non-compatible)
 */
/*
void neighbour_specs_vv_2e( struct Neighbour *ngh )
{
	xstrtok( ngh->line);
	xstrtok( NULL);
	ngh->n_elements = 2;
	ngh->element = (ElementIter *) xmalloc( ngh->n_elements * sizeof( ElementIter  ) );
	ngh->eid = (int*) xmalloc( ngh->n_elements * sizeof( int ) );
	ngh->eid[ 0 ] = atoi( xstrtok( NULL) );
	ngh->eid[ 1 ] = atoi( xstrtok( NULL) );
	ngh->geom_factor = atof( xstrtok( NULL) ); // fraction of lower dim. element measure
    ngh->sigma= atof( xstrtok( NULL) ); // transition coefficient
	ngh->element[ 0 ] = NULL;
	ngh->element[ 1 ] = NULL;
	xfree( ngh->line );
	ngh->line = NULL;
}*/
//-----------------------------------------------------------------------------
// vim: set cindent:
