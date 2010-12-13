#include "system.hh"
#include "mesh.h"

static struct Edge *new_edge(void);
static void add_to_edge_list(Mesh*,struct Edge*);
static void init_edge(struct Edge*);
static int count_edges(Mesh*);

static int number_of_common_nodes_ss(struct Side*,struct Side*);
static int number_of_common_nodes_se(struct Side*,ElementIter );

//=============================================================================
// CREATE AND PREFILL LIST OF EDGES
//=============================================================================
void make_edge_list(Mesh* mesh)
{
	int edi;
	struct Edge *edg;

	ASSERT(!( mesh == NULL ),"NULL as argument of function make_edge_list()\n");
	xprintf( Msg, "Creating edges... ")/*orig verb 2*/;
	mesh->n_edges = count_edges( mesh );
	for( edi = 0; edi < mesh->n_edges; edi++ ) {
		edg = new_edge();
		ASSERT(!( edg == NULL ),"Cannot create edge %d\n", edi );
		add_to_edge_list( mesh, edg );
		edg->id = edi;
	}
	xprintf( MsgVerb, " O.K. %d edges created.", mesh->n_edges )/*orig verb 4*/;
	xprintf( Msg, "O.K.\n")/*orig verb 2*/;
}

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
}
//=============================================================================
//
//=============================================================================
int count_edges(Mesh* mesh)
{
	int rc;
	struct Neighbour *ngh;

	rc = mesh->n_sides;
	FOR_NEIGHBOURS( ngh ) {
		if( ngh->type == BB_E || ngh->type == BB_EL )
			rc -= ( ngh->n_elements - 1 );
	}
	return rc;
}
//=============================================================================
// CREATE NEW EDGE
//=============================================================================
struct Edge *new_edge( void )
{
	struct Edge *edg;

	edg = (struct Edge*) xmalloc( sizeof( struct Edge ) );
	init_edge( edg );
	return edg;
}
//=============================================================================
// INIT DATA OF PARTICULAR EDGE
//=============================================================================
void init_edge( struct Edge *edg )
{
	ASSERT(!( edg == NULL ),"NULL as argument of function init_edge()\n");
	edg->id          = NDEF;
	edg->n_sides     = NDEF;
	edg->side        = NULL;
	edg->neigh_vb    = NULL;
	edg->neigh_bb    = NULL;
	edg->prev        = NULL;
	edg->next        = NULL;
	edg->c_row       = NDEF;
	edg->f_val	  = 0.0;
	edg->aux         = NDEF;
	edg->faux        = 0.0;
}
//=============================================================================
//
//=============================================================================
void add_to_edge_list(Mesh* mesh, struct Edge* edg)
{
	ASSERT(!( (mesh == NULL) || (edg == NULL) ),"NULL as an argument of function add_to_edge_list()\n");
	// First edge in the list
	if( mesh->edge == NULL && mesh->l_edge == NULL ) {
		mesh->edge = edg;
		mesh->l_edge = edg;
		edg->prev = NULL;
		edg->next = NULL;
		return;
	}
	// If something is wrong with the list
	ASSERT(!( (mesh->edge == NULL) || (mesh->l_edge == NULL) ),"Inconsistency in the edge list\n");
	// Add after last edge
	edg->next = NULL;
	edg->prev = mesh->l_edge;
	mesh->l_edge->next = edg;
	mesh->l_edge = edg;
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
	FOR_EDGES( edg ) {
		edg->c_row = mesh->n_sides + mesh->n_elements() + edi;
		edg->f_rhs=0.0;
		if( edg->neigh_vb == NULL )
			edg->f_val = 0.0;
		else
			edg->f_val = -1.0 * edg->neigh_vb->sigma * edg->side[0]->metrics;
		edi++;
	}
	xprintf( Msg, "O.K.\n")/*orig verb 2*/;
}
//-----------------------------------------------------------------------------
// vim: set cindent:
