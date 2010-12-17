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
 * @brief Element preprocessing
 *
 */

#include "constantdb.h"
#include "mesh/ini_constants_mesh.hh"

#include "transport.h"

#include "system.hh"
#include "elements.h"
#include "problem.h"
#include "mesh.h"
#include "initials.h"
#include "preprocess.h"

static void make_ele_initials(struct Problem*);

//=============================================================================
//SET VALUES FROM INITIAL CONDITIONS ETC...
//=============================================================================
void preprocess( struct Problem *problem )
{
    F_ENTRY;

	xprintf( Msg, "Preprocessing...")/*orig verb 2*/;
    if (ConstantDB::getInstance()->getInt("Problem_type") == UNSTEADY_SATURATED) // presunout sem
    	make_ele_initials( problem );            // pp koncentraci

                /*init_concentration_list( mesh );
	FOR_CONCENTRATIONS( con ) {
		xfgets( line, LINE_SIZE - 2, in );
		parse_concentration_line( con, line );
	}
	xfclose( in );*/

	//settings for variable density


}
//=============================================================================
//SET ELEMENT'S VALUES FROM INITIAL CONDITIONS
//=============================================================================
void make_ele_initials(struct Problem *problem)
{
	xprintf( Msg, "   Setting element's values from initial conditions...")/*orig verb 2*/;

        Mesh* mesh = (Mesh*)ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);
        ElementIter ele;
        FOR_ELEMENTS( ele ){
                ele->pscalar = ele->initial->epress;
        }
	xprintf( Msg, "O.K.\n")/*orig verb 2*/;
        return;
}


