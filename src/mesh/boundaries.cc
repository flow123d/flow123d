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
 * @brief Boundary conditions
 * @ingroup mesh
 *
 */

#include "system/system.hh"
#include "mesh/boundaries.h"
#include "mesh/mesh.h"

static struct Boundary *new_boundary(void);
static void add_to_boundary_list(struct Mesh*,struct Boundary*);
static void init_boundary_list(struct Mesh*);
static void parse_boundary_line(struct Boundary*,char*);

//=============================================================================
// READ DATA OF BOUNDARY CONDITIONS
//=============================================================================
void read_boundary( struct Mesh *mesh )
{
	FILE	*in;		  // input file
	char     line[ LINE_SIZE ]; // line of data file
//	int where;
	int bcd_id, n_tags;
	BoundaryFullIter bcd = BOUNDARY_NULL(mesh);
	ElementFullIter ele = ELEMENT_FULL_ITER_NULL(mesh);

	ASSERT(!( mesh == NULL ),"NULL as argument of function read_boundary_list()\n");
	xprintf( Msg, "Reading boundary conditions...")/*orig verb 2*/;
	const std::string& fname = IONameHandler::get_instance()->get_input_file_name(OptGetStr( "Input", "Boundary", "\\" ));
	in = xfopen( fname, "rt" );
	skip_to( in, "$BoundaryConditions" );
	xfgets( line, LINE_SIZE - 2, in );

	int n_boundaries = atoi( xstrtok( line) );
	INPUT_CHECK( n_boundaries >= 1 ,"Number of bounaries < 1 in function read_boundary_list()\n");
	mesh->boundary.reserve(n_boundaries);

	int group_number=0;

	for(int i_bcd=0; i_bcd < n_boundaries; i_bcd++) {
	    // Read one line
        xfgets( line, LINE_SIZE - 2, in );
        // Parse the line
        bcd_id    = atoi( xstrtok( line) );
        bcd = mesh->boundary.add_item(bcd_id);
        // DBGMSG("boundary id: %d \n",bcd_id);

        bcd->type  = atoi( xstrtok( NULL) );

        // physical data - should be moved to water_linsys
        switch( bcd->type ) {
            case DIRICHLET:
                bcd->scalar = atof( xstrtok( NULL) );
                break;
            case NEUMANN:
                bcd->flux   = atof( xstrtok( NULL) );
                break;
            case NEWTON:
                bcd->scalar = atof( xstrtok( NULL) );
                bcd->sigma  = atof( xstrtok( NULL) );
                            break;
            default :
                xprintf(UsrErr,"Unknown type of boundary condition - cond # %d, type %c\n", bcd_id, bcd->type );
        }

        unsigned int where  = atoi( xstrtok( NULL) );
        int eid, sid, n_exterior;
        SideIter sde;

        switch( where ) {
            case 2: // SIDE_EL - BC given by element and its local side number
                eid = atoi( xstrtok( NULL) );
                sid = atoi( xstrtok( NULL) );

                // find and set the side
                ele = mesh->element.find_id( eid );
                if( sid < 0 || sid >= ele->n_sides() )
                     xprintf(UsrErr,"Boundary %d has incorrect reference to side %d\n", bcd_id, sid );
                sde = ele->side( sid );
                ele->boundaries_[ sid ] = bcd;
                bcd->side=sde;

                break;
            case 3: // SIDE_E - BC given only by element, apply to all its sides
                eid = atoi( xstrtok( NULL) );

                // find and set all exterior sides, possibly add more boundaries
                ele = mesh->element.find_id( eid );
                n_exterior=0;
                FOR_ELEMENT_SIDES(ele, li) {
                    sde = ele->side( li );
                    if ( sde->is_external() ) {
                        if (n_exterior > 0) {
                            xprintf(UsrErr, "Implicitly setting BC %d on more then one exterior sides of the element %d.\n", bcd_id, eid);
                            //BoundaryFullIter new_bcd = mesh->boundary.add_item();
                            //*new_bcd = *bcd;
                            //bcd=new_bcd;
                        }
                        ele->boundaries_[ sid ] = bcd;
                        bcd->side=sde;
                        n_exterior++;
                    }
                }

                break;
            default:
                xprintf(UsrErr,"Unknown entity for boundary condition - cond # %d, ent. %c\n", bcd_id, where );
                break;
        }
        //TODO: if group is necessary set it for all bcd in case where == SIDE_E
        n_tags  = atoi( xstrtok( NULL) );
        if( n_tags > 0 ) {
            int group_id = atoi( xstrtok( NULL) );
            flow::VectorId<int>::FullIter group_iter( mesh->bcd_group_id.find_id(group_id) );

            if ( group_iter == mesh->bcd_group_id.end() ) {
                // not found -> create new group
                group_iter = mesh->bcd_group_id.add_item(group_id);
            }
            bcd->group = group_iter.index();   // in fact we do not use integres stored in the vector, but we use index
        }


	}

	xfclose( in );
	xprintf( MsgVerb, " %d conditions readed. ", mesh->n_boundaries() )/*orig verb 4*/;
	xprintf( Msg, "O.K.\n")/*orig verb 2*/;
}





//=============================================================================
// FILL EDGE PART OF MH MATRIX WHICH BELONGS TO BOUNDARIES
//=============================================================================
// TODO: be more sure that BC are applied well
/*
void boundary_calculation_mh( struct Mesh *mesh )
{
	struct Boundary *bcd;
	struct Edge 	*edg;

	xprintf( Msg, "Calculating properties of boundaries... ");
	ASSERT( NONULL(mesh) ,"NULL as argument of function boundary_calculation_mh()\n");
	FOR_BOUNDARIES(mesh,  bcd ) {
		edg=bcd->side->edge();
		// following code may not work if a BC is applied to an edge with neighbouring
		// in such a case we need to add to the f_val, nevertheless the original code
		// does not do it as well
		//ASSERT(edg->n_sides == 1,"Boundary condition %d on an internal edge %d!\n",bcd->id,edg->id);
		//ASSERT(edg->neigh_vb == NULL,"Boundary condition %d on an neighbouring edge %d!",bcd->id,edg->id);
		switch( bcd->type ) {
			// for the dirichlet condition it should be checked, that the corresponding
		    // row in the MH matrix is zero except the diagonal
			case DIRICHLET:
				edg->f_val = -1.0;
				edg->f_rhs = -1.0 * bcd->scalar;
				break;
			case NEUMANN:
				edg->f_val =  0.0;
				edg->f_rhs = bcd->flux * bcd->side->metric();
				break;
			case NEWTON:
				edg->f_val = -1.0 * bcd->side->metric()*bcd->sigma;
				edg->f_rhs = edg->f_val * bcd->scalar;
				break;
			default:
				xprintf(PrgErr,"Unknown type of boundary condition\n");
		}
	}
	xprintf( Msg, "O.K.\n");
}*/
//-----------------------------------------------------------------------------
// vim: set cindent:

