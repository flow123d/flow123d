/*!
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief  Read and init transport boundary conditions list
 *
 */

#include <strings.h>

#include "system.hh"
#include "xio.h"
#include "mesh.h"
#include "transport_bcd.h"
#include "constantdb.h"

static struct Transport_bcd *new_transport_bcd(int);
static void add_to_transport_bcd_list(Mesh*, struct Transport_bcd*);
static void init_transport_bcd_list(Mesh*);
static void init_transport_bcd(struct Transport_bcd*,int);
static void parse_transport_bcd_line(struct Transport_bcd*,char*);

//=============================================================================
// READ DATA OF TRANSPORT BOUNDARY CONDITIONS
//=============================================================================
void read_transport_bcd_list(Mesh* mesh)
{
	FILE	*in;		  // input file
	char     line[ LINE_SIZE ]; // line of data file
	struct Transport_bcd *tbcd;

	ASSERT(!( mesh == NULL ),"NULL as argument of function read_transport_bcd_list()\n");
	xprintf( Msg, "Reading transport boundary conditions...")/*orig verb 2*/;
	in = xfopen( ConstantDB::getInstance()->getChar("Transport_bcd_fname"), "rt" );
	skip_to( in, "$Transport_BCD" );
	xfgets( line, LINE_SIZE - 2, in );
	mesh->n_transport_bcd = atoi( xstrtok( line) );
	INPUT_CHECK(!( mesh->n_transport_bcd < 1 ),"Number of transport boundaries < 1 in function read_transport_bcd_list()\n");
    INPUT_CHECK( mesh->n_boundaries() == mesh->n_transport_bcd,
            "Different number of boundary conditions %d and transport boundary conditions %d\n",mesh->n_boundaries(), mesh->n_transport_bcd);
	init_transport_bcd_list( mesh );
	FOR_TRANSPORT_BCDS( tbcd ) {
		xfgets( line, LINE_SIZE - 2, in );
		parse_transport_bcd_line( tbcd, line );
	}
	xfclose( in );
	xprintf( MsgVerb, " %d transport conditions readed. ", mesh->n_transport_bcd )/*orig verb 4*/;
	xprintf( Msg, "O.K.\n")/*orig verb 2*/;
}
//=============================================================================
// INIT LIST OF TRANSPORT BOUNDARY CONDITIONS
//=============================================================================
void init_transport_bcd_list(Mesh* mesh)
{
	int tbi;
	struct Transport_bcd *tbcd;

	ASSERT(!( mesh == NULL ),"NULL as argument of function init_transport_bcd_list()\n");
	ASSERT(!( mesh->n_transport_bcd < 1 ),"Number of transport conditions < 1 in function init_transport_bcd_list()\n");
	for( tbi = 0; tbi < mesh->n_transport_bcd; tbi++ ) {
		tbcd = new_transport_bcd( mesh->n_substances );
		ASSERT(!( tbcd == NULL ),"Cannot create transport boundary %d\n", tbi );
		add_to_transport_bcd_list( mesh, tbcd );
	}
}
//=============================================================================
// CREATE NEW TRANSPORT BOUNDARY CONDITION
//=============================================================================
struct Transport_bcd *new_transport_bcd( int n_subst )
{
	struct Transport_bcd *tbcd;

	tbcd = (struct Transport_bcd*) xmalloc( sizeof( struct Transport_bcd ) );
	init_transport_bcd( tbcd, n_subst );
	return tbcd;
}
//=============================================================================
// INIT DATA OF TRANSPORT BOUNDARY CONDITION
//=============================================================================
void init_transport_bcd( struct Transport_bcd *tbcd, int n_subst )
{
	int sbi;

	ASSERT(!( tbcd == NULL ),"NULL as argument of function init_transport_bcd()\n");
	ASSERT(!( n_subst < 1 ),"Number of substances must be positive\n");
	tbcd->id       = NDEF;
	tbcd->bid      = NDEF;
	tbcd->n_subst  = n_subst;
	tbcd->conc     = (double*) xmalloc( n_subst * sizeof( double ) );
	for( sbi = 0; sbi < tbcd->n_subst; sbi++ )
		tbcd->conc[ sbi ] = 0.0;
	tbcd->prev     = NULL;
	tbcd->next     = NULL;
    tbcd->aux      = NDEF;
	tbcd->faux     = 0.0;
}
//=============================================================================
//
//=============================================================================
void add_to_transport_bcd_list(Mesh* mesh, struct Transport_bcd* tbcd)
{
	ASSERT(!( (mesh == NULL) || (tbcd == NULL) ),"NULL as an argument of function add_to_transport_bcd_list()\n");
	// First transport boundary in the list
	if( (mesh->transport_bcd == NULL) && (mesh->l_transport_bcd == NULL) ) {
		mesh->transport_bcd = tbcd;
		mesh->l_transport_bcd = tbcd;
		tbcd->prev = NULL;
		tbcd->next = NULL;
		return;
	}
	// If something is wrong with the list
	ASSERT(!( (mesh->transport_bcd == NULL) || (mesh->l_transport_bcd == NULL) ),"Inconsistency in the transport boundary list\n");
	// Add after last boundary
	tbcd->next = NULL;
	tbcd->prev = mesh->l_transport_bcd;
	mesh->l_transport_bcd->next = tbcd;
	mesh->l_transport_bcd = tbcd;
}
//=============================================================================
// PARSE LINE
//=============================================================================
void parse_transport_bcd_line( struct Transport_bcd *tbcd, char *line )
{
	int sbi;

	ASSERT(!( (tbcd == NULL) || (line == NULL) ),"NULL as argument of function parse_transport_bcd_line()\n");
	tbcd->id    = atoi( xstrtok( line) );
	//TODO: muze byt id=0 ?
	INPUT_CHECK(!( tbcd->id < 0 ),"Id number of transport boundary must be > 0\n");
	tbcd->bid  = atoi( xstrtok( NULL) );
	for( sbi = 0; sbi < tbcd->n_subst; sbi++ )
	        tbcd->conc[ sbi ] = atof( xstrtok( NULL) );
}
//-----------------------------------------------------------------------------
// vim: set cindent:
