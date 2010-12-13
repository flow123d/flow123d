/*!
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief  Initialize neighbouring
 *
 */

#include <strings.h>

#include "constantdb.h"
#include "system.hh"
#include "xio.h"
#include "neighbours.h"
#include "mesh.h"

static struct Neighbour *new_neighbour(void);
static void add_to_neighbour_list(Mesh*, struct Neighbour*);
static void init_neighbour_list(Mesh*);
static void init_neighbour(struct Neighbour*);
static void parse_neighbour_line(struct Neighbour*,char*);
static char supported_neighbour_type(int);
static void neighbour_type_specific(struct Neighbour*);
static void neighbour_specs_bb_e(struct Neighbour*);
static void neighbour_specs_bb_el(struct Neighbour*);
static void neighbour_specs_vb_es(struct Neighbour*);
static void neighbour_specs_vv_2e(struct Neighbour*);


//=============================================================================
// READ DATA OF ALL NEIGHBOURS
//=============================================================================
void read_neighbour_list(Mesh* mesh)
{
	FILE	*in;   // input file
	char     line[ LINE_SIZE ];   // line of data file
	struct Neighbour *ngh;

	ASSERT(!( mesh == NULL ),"NULL as argument of function read_neighbour_list()\n");
	xprintf( Msg, "Reading neighbours...")/*orig verb 2*/;
	in = xfopen( OptGetStr( "Input", "Neighbouring", "\\" ), "rt" );
	skip_to( in, "$Neighbours" );
	xfgets( line, LINE_SIZE - 2, in );
	mesh->n_neighs = atoi( xstrtok( line) );
	INPUT_CHECK(!( mesh->n_neighs < 1 ),"Number of neighbours  < 1 in read_neighbour_list()\n");
	init_neighbour_list( mesh );
	FOR_NEIGHBOURS( ngh ) {
		xfgets( line, LINE_SIZE - 2, in );
		parse_neighbour_line( ngh, line );
		neighbour_type_specific( ngh );
	}
	xprintf( MsgVerb, " %d neighbours readed. ", mesh->n_neighs )/*orig verb 4*/;
	xprintf( Msg, "O.K.\n")/*orig verb 2*/;
}
//=============================================================================
// INIT LIST OF NEIGHBOURS
//=============================================================================
void init_neighbour_list(Mesh* mesh)
{
	int ngi;
	struct Neighbour *ngh;

	ASSERT(NONULL( mesh ),"NULL as argument of function init_neighbours_list()\n");
    ASSERT(!( mesh->n_neighs < 1 ),"Number of neighbours < 1 in function init_neighbours_list()\n");

	for( ngi = 0; ngi < mesh->n_neighs; ngi++ ) {
		ngh = new_neighbour();
		ASSERT(!( ngh == NULL ),"Cannot create neighbour %d\n", ngi );
		add_to_neighbour_list( mesh, ngh );
	}
}
//=============================================================================
// CREATE NEW NEIGHBOUR
//=============================================================================
struct Neighbour *new_neighbour( void )
{
	struct Neighbour *ngh;

	ngh = (struct Neighbour*) xmalloc( sizeof( struct Neighbour ) );
	init_neighbour( ngh );
	return ngh;
}
//=============================================================================
// INIT DATA OF NEIGHBOUR
//=============================================================================
void init_neighbour( struct Neighbour *ngh )
{
	ASSERT(!( ngh == NULL ),"NULL as argument of function init_neighbour()\n");
	ngh->id          = NDEF;
	ngh->type        = NDEF;
	ngh->n_sides     = NDEF;
	ngh->n_elements  = NDEF;
	ngh->sid  	     = NULL;
	ngh->eid         = NULL;
    ngh->flux        = 0.0;
	ngh->sigma        = 0.0;
	ngh->geom_factor    = 0.0;
	ngh->edge        = NULL;
	ngh->side        = NULL;
	ngh->element     = NULL;
	ngh->line	 = NULL;
	ngh->prev	 = NULL;
	ngh->next	 = NULL;
	ngh->aux	 = NDEF;
	ngh->faux	 = 0.0;
}
//=============================================================================
//
//=============================================================================
void add_to_neighbour_list(Mesh* mesh, struct Neighbour* ngh)
{
	ASSERT(!( (mesh == NULL) || (ngh == NULL) ),"NULL as an argument of function add_to_neighbour_list()\n");
	// First element in the list
	if( mesh->neighbour == NULL && mesh->l_neighbour == NULL ) {
		mesh->neighbour = ngh;
		mesh->l_neighbour = ngh;
		ngh->prev = NULL;
		ngh->next = NULL;
		return;
	}
	// If something is wrong with the list
	ASSERT(!( (mesh->neighbour == NULL) || (mesh->l_neighbour == NULL) ),"Inconsistency in the neighbour list\n");
	// Add after last neighbour
	ngh->next = NULL;
	ngh->prev = mesh->l_neighbour;
	mesh->l_neighbour->next = ngh;
	mesh->l_neighbour = ngh;
}
//=============================================================================
// PARSE LINE
//=============================================================================
void parse_neighbour_line( struct Neighbour *ngh, char *line )
{
	int id, type;

	ASSERT(!( (ngh == NULL) || (line == NULL) ),"NULL as argument of function parse_neighbour_line()\n");
	ngh->line	= xstrcpy( line );
	id              = atoi( xstrtok( line) );
	//TODO: co kdyz id=0 ?
	INPUT_CHECK(!( id < 0 ),"Id number of neighbour must be > 0\n");
	type            = atoi( xstrtok( NULL) );
	if( supported_neighbour_type( type ) == false )
		xprintf(UsrErr,"Neighbour %d is of the unsupported type %d\n", id, type );
	ngh->id        = id;
	ngh->type      = type;
}
//=============================================================================
//
//=============================================================================
char supported_neighbour_type( int type )
{
	switch( type ) {
		case BB_E:
		case BB_EL:
		case VB_ES:
		case VV_2E:
			return true;
	}
	return false;
}
//=============================================================================
//
//=============================================================================
void neighbour_type_specific( struct Neighbour *ngh )
{
	switch( ngh->type ) {
		case BB_E:
			neighbour_specs_bb_e( ngh );
			break;
		case BB_EL:
			neighbour_specs_bb_el( ngh );
			break;
		case VB_ES:
			neighbour_specs_vb_es( ngh );
			break;
		case VV_2E:
			neighbour_specs_vv_2e( ngh );
			break;
	}
}
/**
 *  boundary - boundary neigbouring, given by list of neigbouring elements,
 *             sides are given implicitely by shared nodes
 */
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
}
/**
 *  boundary - boundary neigbouring,
 *  given by list of neigbouring sides given by
 *  their elements and local side numbers,
 */
void neighbour_specs_bb_el( struct Neighbour *ngh )
{
	int ei;

	xstrtok( ngh->line);
	xstrtok( NULL);
	ngh->n_elements = atoi( xstrtok( NULL) );
	INPUT_CHECK(!( ngh->n_elements < 2 ),"Neighbour %d has bad number of elements: %d\n", ngh->id, ngh->n_elements );
	ngh->n_sides = ngh->n_elements;
	ngh->eid = (int*) xmalloc( ngh->n_elements * sizeof( int ) );
	ngh->sid = (int*) xmalloc( ngh->n_elements * sizeof( int ) );
	ngh->element = (ElementIter *) xmalloc( ngh->n_elements *
			sizeof( ElementIter  ) );
	ngh->side = (struct Side**) xmalloc( ngh->n_elements *
			sizeof( struct Side* ) );
	FOR_NEIGH_ELEMENTS( ngh, ei ) {
		ngh->element[ ei ] = NULL;
		ngh->side[ ei ] = NULL;
	}
	FOR_NEIGH_ELEMENTS( ngh, ei ) {
		ngh->eid[ ei ] = atoi( xstrtok( NULL) );
		ngh->sid[ ei ] = atoi( xstrtok( NULL) );
	}
	xfree( ngh->line );
	ngh->line = NULL;
}
/**
 *  boundary - volume neighbouring of different dimensions (compatible)
 */
void neighbour_specs_vb_es( struct Neighbour *ngh )
{
	xstrtok( ngh->line);
	xstrtok( NULL);
	ngh->n_elements = 2;
	ngh->n_sides = 2;
	ngh->eid = (int*) xmalloc( ngh->n_elements * sizeof( int ) );
	ngh->element = (ElementIter *) xmalloc( ngh->n_elements * sizeof( ElementIter  ) );
	ngh->sid = (int*) xmalloc( ngh->n_sides * sizeof( int ) );
	ngh->side = (struct Side**) xmalloc( ngh->n_sides *	sizeof( struct Side* ) );
	ngh->eid[ 0 ] = atoi( xstrtok( NULL) );
	ngh->eid[ 1 ] = atoi( xstrtok( NULL) );
	ngh->sid[ 0 ] = NDEF;
	ngh->sid[ 1 ] = atoi( xstrtok( NULL) );
	ngh->sigma = atof( xstrtok( NULL) );
	ngh->element[ 0 ] = NULL;
	ngh->element[ 1 ] = NULL;
	ngh->side[ 0 ] = NULL;
	ngh->side[ 1 ] = NULL;
	xfree( ngh->line );
	ngh->line = NULL;
	// sid[ 0 ] (and side[ 0 ]) doesn't have defined value. I use sid[ 1 ] (and
	// side[ 1 ]) instead to correspond with elm1. Using sid[ 0 ] for describing
	// side of elm1 would be confusing and error-prone.
}
/**
 * volume - volume neighbouring of different dimensions (non-compatible)
 */
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
}
//-----------------------------------------------------------------------------
// vim: set cindent:
