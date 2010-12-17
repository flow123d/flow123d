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
 * @brief  Initial conditions setup
 *
 */

#include <strings.h>

#include "constantdb.h"
#include "system.hh"
#include "xio.h"
#include "initials.h"
#include "mesh.h"

static struct Initial *new_initial(void);
static void add_to_initial_list(Mesh*, struct Initial*);
static void init_initial_list(Mesh*);
static void init_initial(struct Initial*);
static void parse_initial_line(struct Initial*,char*);

//=============================================================================
// READ DATA OF initial CONDITIONS
//=============================================================================
void read_initial_list(Mesh* mesh)
{
	FILE	*in;		  // input file
	char     line[ LINE_SIZE ]; // line of data file
	struct Initial *icd;
        double time;

	ASSERT(!( mesh == NULL ),"NULL as argument of function read_initial_list()\n");
	xprintf( Msg, "Reading initial conditions...")/*orig verb 2*/;
	in = xfopen( OptGetStr( "Input", "Initial", "\\" ), "rt" );
	skip_to( in, "$FlowField" );
	xfgets( line, LINE_SIZE - 2, in );
	time = atof( xstrtok( line) );
	xfgets( line, LINE_SIZE - 2, in );
	mesh->n_initials = atoi( xstrtok( line) );
	INPUT_CHECK(!( mesh->n_initials < 1 ),"Number of initials < 1 in function read_initial_list()\n");
	init_initial_list( mesh );
	FOR_INITIALS( icd ) {
		xfgets( line, LINE_SIZE - 2, in );
		parse_initial_line( icd, line );
	}
	xfclose( in );
	xprintf( MsgVerb, " %d conditions readed. ", mesh->n_initials )/*orig verb 4*/;
	xprintf( Msg, "O.K.\n")/*orig verb 2*/;
}
//=============================================================================
// INIT LIST OF initial CONDITIONS
//=============================================================================
void init_initial_list(Mesh* mesh)
{
	int bi;
	struct Initial *icd;

	ASSERT(NONULL( mesh ),"NULL as argument of function init_initial_list()\n");
    ASSERT(!( mesh->n_initials < 1 ),"Number of conditions < 1 in function init_initial_list()\n");

	for( bi = 0; bi < mesh->n_initials; bi++ ) {
		icd = new_initial();
		ASSERT(!( icd == NULL ),"Cannot create initial %d\n", bi );
		add_to_initial_list( mesh, icd );
	}
}
//=============================================================================
// CREATE NEW initial CONDITION
//=============================================================================
struct Initial *new_initial( void )
{
	struct Initial *icd;

	icd = (struct Initial*) xmalloc( sizeof( struct Initial ) );
	init_initial( icd );
	return icd;
}
//=============================================================================
// INIT DATA OF initial CONDITION
//=============================================================================
void init_initial( struct Initial *icd )
{
	ASSERT(!( icd == NULL ),"NULL as argument of function init_initial()\n");
	icd->id       = NDEF;
	icd->eid      = NDEF;
	icd->epress   = 0.0;
	icd->n_sides  = NDEF;
	icd->spress   = NULL;
	icd->sflux    = NULL;
	icd->n_neighs_vv = NDEF;
	icd->hid      = NULL;
	icd->prev     = NULL;
	icd->next     = NULL;
	icd->aux      = NDEF;
	icd->faux     = 0.0;
}
//=============================================================================
//
//=============================================================================
void add_to_initial_list(Mesh* mesh, struct Initial *icd )
{
	ASSERT(!( (mesh == NULL) || (icd == NULL) ),"NULL as an argument of function add_to_initial_list()\n");
	// First initial in the list
	if( mesh->initial == NULL && mesh->l_initial == NULL ) {
		mesh->initial = icd;
		mesh->l_initial = icd;
		icd->prev = NULL;
		icd->next = NULL;
		return;
	}
	// If something is wrong with the list
	ASSERT(!( (mesh->initial == NULL) || (mesh->l_initial == NULL) ),"Inconsistency in the initial list\n");
	// Add after last initial
	icd->next = NULL;
	icd->prev = mesh->l_initial;
	mesh->l_initial->next = icd;
	mesh->l_initial = icd;
}
//=============================================================================
// PARSE LINE
//=============================================================================
void parse_initial_line( struct Initial *icd, char *line )
{
	int i;

	ASSERT(!( (icd == NULL) || (line == NULL) ),"NULL as argument of function parse_initial_line()\n");
	icd->id    = atoi( xstrtok( line) );
	// TODO: muze byt ID i nulove?
	INPUT_CHECK(!( icd->id < 0 ),"Id number of initial must be > 0\n");
	icd->eid    = atoi( xstrtok( NULL) );
    // TODO: muze byt ID i nulove?
	INPUT_CHECK(!( icd->eid < 0 ),"Id number of element must be > 0\n");
	icd->epress  = atof( xstrtok( NULL) );
	icd->n_sides  = atoi( xstrtok( NULL) );
    icd->spress = (double *) xmalloc(icd->n_sides * sizeof(double));
    icd->sflux = (double *) xmalloc(icd->n_sides * sizeof(double));
    for (i = 0; i < icd->n_sides; i++)
        icd->spress[i] = atof(xstrtok(NULL));
    for (i = 0; i < icd->n_sides; i++)
        icd->sflux[i] = atof(xstrtok(NULL));
    icd->n_neighs_vv = atoi(xstrtok(NULL));
    if (icd->n_neighs_vv > 0) {
        icd->hid = (int *) xmalloc(icd->n_neighs_vv * sizeof(int));
        for (i = 0; i < icd->n_neighs_vv; i++)
            icd->hid[i] = atoi(xstrtok(NULL));
    }
}
//-----------------------------------------------------------------------------
// vim: set cindent:

