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
 *
 */

#include <strings.h>
#include <string.h>

#include "system.hh"
#include "xio.h"
#include "mesh.h"
#include "concentrations.h"
#include "constantdb.h"

static struct Concentration *new_concentration(int);
static void add_to_concentration_list(Mesh*, struct Concentration*);
static void init_concentration_list(Mesh*);
static void init_concentration(struct Concentration*,int);
static void parse_concentration_line(struct Concentration*,char*);

//=============================================================================
// READ DATA OF CONCENTRATIONS
//=============================================================================
void read_concentration_list(Mesh* mesh)
{
	FILE	*in;		  // input file
	char     line[ LINE_SIZE ]; // line of data file
	struct Concentration *con;

	ASSERT(!( mesh == NULL ),"NULL as argument of function read_concentration_list()\n");
	xprintf( Msg, "Reading concentrations...")/*orig verb 2*/;
	in = xfopen( ConstantDB::getInstance()->getChar("Concentration_fname"), "rt" );
	skip_to( in, "$Concentrations" );
	xfgets( line, LINE_SIZE - 2, in );
	mesh->n_concentrations = atoi( xstrtok( line) );
	INPUT_CHECK(!( mesh->n_concentrations < 1 ),"Number of concentrations < 1 in function read_concentration_list()\n");
    INPUT_CHECK(!( mesh->n_elements() != mesh->n_concentrations),"Different number of elements and concentrations\n");
	init_concentration_list( mesh );
	FOR_CONCENTRATIONS( con ) {
		xfgets( line, LINE_SIZE - 2, in );
		parse_concentration_line( con, line );
	}

        // FOR_ELEMENTS(ele) {
        // xprintf(MsgVerb, " %f " , ele->start_conc->conc[0]);

	xfclose( in );
	xprintf( MsgVerb, " %d concentrations readed. ", mesh->n_concentrations )/*orig verb 4*/;
	xprintf( Msg, "O.K.\n")/*orig verb 2*/;
}
//=============================================================================
// INIT LIST OF CONCENTRATIONS
//=============================================================================
void init_concentration_list(Mesh* mesh)
{
	int ci;
	struct Concentration *con;

	ASSERT(NONULL(mesh),"NULL as argument of function init_concentration_list()\n");
	for( ci = 0; ci < mesh->n_concentrations; ci++ ) {
		con = new_concentration( mesh->n_substances );
		ASSERT(NONULL(con),"Cannot create concentration %d\n",ci);
		add_to_concentration_list( mesh, con );
	}
}
//=============================================================================
// CREATE NEW CONCENTRATION
//=============================================================================
struct Concentration *new_concentration( int n_subst )
{
	struct Concentration *con;

	con = (struct Concentration*) xmalloc( sizeof( struct Concentration ) );
	init_concentration( con, n_subst );
	return con;
}
//=============================================================================
// INIT DATA OF CONCENTRATION
//=============================================================================
void init_concentration( struct Concentration *con, int n_subst )
{
	int sbi;

	ASSERT(!( con == NULL ),"Cannot create concentrations, NULL parameter.\n");
    ASSERT(!( n_subst < 1 ),"Number of substances must be positive.\n");
	con->id       = NDEF;
	con->eid      = NDEF;
	con->n_subst  = n_subst;
	con->conc     = (double*) xmalloc( n_subst * sizeof( double ) );
	for( sbi = 0; sbi < con->n_subst; sbi++ )
		con->conc[ sbi ] = 0.0;
	con->prev     = NULL;
	con->next     = NULL;
	con->aux      = NDEF;
	con->faux     = 0.0;
}
//=============================================================================
//
//=============================================================================
void add_to_concentration_list(Mesh* mesh, struct Concentration* con)
{
	ASSERT(!( (mesh == NULL) || (con == NULL) ),"NULL as an argument of function add_to_concentration_list()\n");
	// First concentration in the list
	if( mesh->concentration == NULL && mesh->l_concentration == NULL ) {
		mesh->concentration = con;
		mesh->l_concentration = con;
		con->prev = NULL;
		con->next = NULL;
		return;
	}
	// If something is wrong with the list
	ASSERT(!( (mesh->concentration == NULL) || (mesh->l_concentration == NULL) ),"Inconsistency in the concentration list\n");
	// Add after last boundary
	con->next = NULL;
	con->prev = mesh->l_concentration;
	mesh->l_concentration->next = con;
	mesh->l_concentration = con;
}
//=============================================================================
// PARSE LINE
//=============================================================================
void parse_concentration_line( struct Concentration *con, char *line )
{
	int sbi;

	ASSERT(!( (con == NULL) || (line == NULL) ),"NULL as argument of function parse_concentration_line()\n");
	con->id    = atoi( xstrtok( line) );
	// TODO: id musi byt >0 nebo >=0 ???
	INPUT_CHECK(!( con->id < 0 ),"Id number of concentration must be > 0\n");
	con->eid    = atoi( xstrtok( NULL) );
	for( sbi = 0; sbi < con->n_subst; sbi++ )
	        con->conc[ sbi ] = atof( xstrtok( NULL) );
}
//-----------------------------------------------------------------------------
// vim: set cindent:
