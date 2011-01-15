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
 * @brief Reactions
 * Read input data and perform reaction computation
 *
 */

#include "constantdb.h"
#include "mesh/ini_constants_mesh.hh"

#include "transport.h"

#include "system.hh"
#include "xio.h"
#include "mesh.h"
#include "reaction.h"
#include "materials.hh"
#include "elements.h"

static int reaction_type_specific_coefficients( int st );
static char supported_reaction_type( int st );
void parse_reaction_line( struct Transport *transport, int i, char *line);

//=============================================================================
//      TRANSPORT REACTION
/*! @brief MAIN REACTION COMPUTATION (zero order)
 *
 * input{Transport,elm_pos,sbi}
 *      perform reaction computation on element with elm_pos position
 *
 */
//=============================================================================
void transport_reaction(struct Transport *transport, int elm_pos, MaterialDatabase::Iter material, int sbi)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    struct Reaction *rct;
        int i;
        double ***conc=transport->conc,***pconc=transport->pconc;

        // TODO: remove epos_id, find_id
        for(i=0;i < transport->n_reaction;i++){
				rct = &transport->reaction[i];
                switch(rct->type){
                        case 0:
                         if(rct->sbi != sbi)
                                break;
                         if(rct->coef[0] != -1.0){
                        	 conc[MOBILE][rct->sbi][elm_pos] += rct->coef[0] * transport->time_step;
                        	 pconc[MOBILE][rct->sbi][elm_pos] = conc[MOBILE][rct->sbi][elm_pos];
                      //  	 printf("%f\t%f\n",rct->coef[0],transport->time_step);
                      //  	 getchar();
                         }
                         else{
                        	 conc[MOBILE][rct->sbi][elm_pos] += 1/(material->size);
                        	 pconc[MOBILE][rct->sbi][elm_pos] = conc[MOBILE][rct->sbi][elm_pos];
                         }
                         break;
                }
        }
}
//=============================================================================
// PARSE REACTION LINE
/*! @brief PARSE REACTION LINE
 *
 * input{Transport,reaction_number,line}
 *
 */
//=============================================================================
void parse_reaction_line( struct Transport *transport, int i, char *line)
{
	int ci,pc;     // si - substance index
	struct Reaction *rct;

	rct = &transport->reaction[i];

	rct->type      = atoi( xstrtok( line) );
    if(supported_reaction_type(rct->type) == true)
        pc =  reaction_type_specific_coefficients(rct->type);
    else
        xprintf(UsrErr,"Unsupported reaction type #%d.\n", rct->type);

    rct->sbi       = atoi( xstrtok( NULL) );
    INPUT_CHECK(!(rct->sbi > transport->n_substances - 1),"Unknown substance ID#%d.\n", rct->sbi);

    rct->coef = (double*)xmalloc(pc * sizeof(double));
    for(ci = 0 ; ci < pc ; ci++)
    	rct->coef[ci] = atof( xstrtok( NULL) );
}
//=============================================================================
// SUPPORTED REACTION TYPE
/*! @brief SUPPORTED REACTION TYPE
 *
 * input{reaction_type}
 *
 */
//=============================================================================
char supported_reaction_type( int st )
{
	switch (st) {
    case 0:
        return true;
    }
    return false;
}
//=============================================================================
// REACTION TYPE SPECIFIC COEFFICIENTS
/*! @brief REACTION TYPE SPECIFIC COEFFICIENTS
 *
 * input{reaction_type}
 *
 */
//=============================================================================
int reaction_type_specific_coefficients( int st )
{
    switch (st) {
    case 0:
        return 1;
    }
    return -1;
}
//=============================================================================
// READ REACTIONS LIST
/*! @brief READ REACTIONS LIST FROM *.MTR FILE
 *
 * input{Transport}
 *
 */
//=============================================================================
void read_reaction_list( struct Transport *transport )
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    FILE* in;   // input file
    char  line[ LINE_SIZE ],string[ LINE_SIZE ];   // line of data file
    bool  found;
    int   ki,r,i;

		ASSERT(!( mesh == NULL ),"NULL as argument of function read_reaction_list()\n");
		xprintf( Msg, "Reading reactions... ")/*orig verb 2*/;
        in = xfopen( ConstantDB::getInstance()->getChar("Material_fname"), "rt" );
        found=skip_to( in, "$Reactions" );
        ASSERT( found , "Can not find section: $Reactions." );
        xfgets( line, LINE_SIZE - 2, in );
        sscanf( line, "%s", string ); //test
        r = 0;

			while ((strcmp(string,"$EndReactions")) != 0){
				r++;
                xfgets( line, LINE_SIZE - 2, in );
                sscanf( line, "%s", string ); //test
            }

            transport->n_reaction = r;
            ASSERT(!( r == 0 ),"No reaction haven't defined\n");
            transport->reaction =(struct Reaction*)xmalloc(transport->n_reaction * sizeof(struct Reaction));
            xfclose( in );

            in = xfopen( ConstantDB::getInstance()->getChar("Material_fname"), "rt" );
            skip_to( in, "$Reactions" );
			for(i=0;i<transport->n_reaction;i++){
				xfgets( line, LINE_SIZE - 2, in );
				parse_reaction_line(transport, i ,line );
			}
            xfclose( in );
            xprintf( Msg, "O.K.\n")/*orig verb 2*/;
}
