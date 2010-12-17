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
 * @brief  Read and initialize sources
 *
 */

#include <strings.h>

#include "system.hh"
#include "xio.h"
#include "sources.h"
#include "mesh.h"


#if __NOTHING__

static struct Source *new_source(void);
static void add_to_source_list(Mesh*,struct Source*);
static void init_source_list(Mesh*);
static void init_source(struct Source*);
static void parse_source_line(struct Source*,char*);

//=============================================================================
// READ DATA OF ALL SOURCES
//=============================================================================
void read_source_list(Mesh* mesh)
{
	FILE	*in;   // input file
	char     line[ LINE_SIZE ];   // line of data file
	struct Source *src;

	ASSERT(!( mesh == NULL ),"NULL as argument of function read_source_list()\n");
	xprintf( Msg, "Reading sources... ")/*orig verb 2*/;
	if( mesh->sources_fname == NULL ) {
		mesh->n_sources = NDEF;
		mesh->source = NULL;
		mesh->l_source = NULL;
		xprintf( Msg, "File with sources not defined, skipping.\n ")/*orig verb 2*/;
		xprintf( Msg, "O.K.\n")/*orig verb 2*/;
		return;
	}
	in = xfopen( mesh->sources_fname, "rt" );
	skip_to( in, "$Sources" );
	xfgets( line, LINE_SIZE - 2, in );
	mesh->n_sources = atoi( xstrtok( line) );
	ASSERT(!( mesh->n_sources < 1 ),"Number of sources < 1 in function read_source_list()\n");
	init_source_list( mesh );
	FOR_SOURCES( src ) {
		xfgets( line, LINE_SIZE - 2, in );
		parse_source_line( src, line );
	}
	xfclose( in );
	xprintf( MsgVerb, " %d sources readed. ", mesh->n_sources )/*orig verb 4*/;
	xprintf( Msg, "O.K.\n")/*orig verb 2*/;
}
//=============================================================================
// INIT LIST OF SOURCES
//=============================================================================
void init_source_list(Mesh* mesh)
{
	int si;
	struct Source *src;

	ASSERT(!( mesh == NULL ),"NULL as argument of function init_source_list()\n");
	for( si = 0; si < mesh->n_sources; si++ ) {
		src = new_source();
		ASSERT(!( src == NULL ),"Cannot create source %d\n", si );
		add_to_source_list( mesh, src );
	}
}
//=============================================================================
// CREATE NEW SOURCE
//=============================================================================
struct Source *new_source( void )
{
	struct Source *src;

	src = (struct Source*) xmalloc( sizeof( struct Source ) );
	init_source( src );
	return src;
}
//=============================================================================
// INIT DATA OF SOURCE
//=============================================================================
void init_source( struct Source *src )
{
	ASSERT(!( src == NULL ),"NULL as argument of function init_source()\n");
	src->id          = NDEF;
	src->type        = NDEF;
	src->eid         = NDEF;
	src->density     = 0.0;
	src->element     = NULL;
	src->prev	 = NULL;
	src->next	 = NULL;
	src->aux         = NDEF;
	src->faux        = 0.0;
}
//=============================================================================
//
//=============================================================================
void add_to_source_list(Mesh* mesh, struct Source* src)
{
	ASSERT(!( (mesh == NULL) || (src == NULL) ),"NULL as an argument of function add_to_source_list()\n");
	// First source in the list
	if( mesh->source == NULL && mesh->l_source == NULL ) {
		mesh->source = src;
		mesh->l_source = src;
		src->prev = NULL;
		src->next = NULL;
		return;
	}
	// If something is wrong with the list
	ASSERT(!( (mesh->source == NULL) || (mesh->l_source == NULL) ),"Inconsistency in the source list\n");
	// Add after last source
	src->next = NULL;
	src->prev = mesh->l_source;
	mesh->l_source->next = src;
	mesh->l_source = src;
}
//=============================================================================
// PARSE LINE
//=============================================================================
void parse_source_line( struct Source *src, char *line )
{
	ASSERT(!( (src == NULL) || (line == NULL) ),"NULL as argument of function parse_source_line()\n");
	src->id = atoi( xstrtok( line) );
	//TODO: co kdyz id=0 ?
	INPUT_CHECK(!( src->id < 0 ),"Id number of source must be > 0\n");
	src->type     = atoi( xstrtok( NULL) );
	src->eid      = atoi( xstrtok( NULL) );
	src->density  = atof( xstrtok( NULL) );
}
//-----------------------------------------------------------------------------
// vim: set cindent:

#endif
