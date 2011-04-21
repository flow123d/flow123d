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
 * @ingroup io
 *
 * @brief Break through curve.
 * @todo Rewrite to be applyable to general fields and produce specified output. documentation!!!
 *
 *
 */

#include "constantdb.h"
#include "mesh/ini_constants_mesh.hh"
#include "transport.h"

#include "read_ini.h"
#include "xio.h"
#include "output.h"
#include "btc.h"
#include "system/system.hh"
#include "mesh/mesh.h"
#include "problem.h"

static int *BTC_elm_list( int n_btc, char *line );
static int count_BTC_elms( char *line );
static void btc_init(struct BTC *btc);

//=============================================================================

//=============================================================================
void btc_init(struct BTC *btc){

	char *Btc;

	Btc					= IONameHandler::get_instance()->get_output_file_name(OptGetStr( "Output", "BTC_elms", "-9999" )).c_str();
	btc->n_BTC_elms		= count_BTC_elms( Btc );
	btc->BTC_elm		= BTC_elm_list( btc->n_BTC_elms, Btc );

}

/*
//=============================================================================

//=============================================================================
void btc_check(struct Transport *transport) {
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    struct BTC *btc;
    int i;

    transport->btc = (struct BTC*) xmalloc(sizeof(struct BTC));
    btc = transport->btc;

    btc_init(btc);

    // TODO: remove find_id
    if (btc->BTC_elm != NULL) {
        if (btc->BTC_elm[0] != -9999) {
            for (i = 0; i < btc->n_BTC_elms; i++) {
                INPUT_CHECK( NONULL(mesh->element.find_id(btc->BTC_elm[i])),
                    "Unknown element #%d in BTC.\n", btc->BTC_elm[i]);
            }
            output_transport_init_BTC(transport);
            output_transport_time_BTC(transport, 0.0);
        } else {
            xfree(transport->btc);
            transport->btc = NULL;
        }
    }
}*/
//=============================================================================

//=============================================================================
int count_BTC_elms( char *line )
{
	char *l;
	int rc;
	ASSERT(!( line == NULL ),"NULL as an argument of function count_BTC_elms()\n");
	l = xstrcpy( line );
	rc = 0;
	if( strtok( l, " \t,;" ) == NULL )
		return rc;
	do
		rc++;
	while( strtok( NULL, " \t,;" ) != NULL );
	return rc;
}
//=============================================================================
//
//=============================================================================
int *BTC_elm_list( int n_btc, char *line )
{
	int *cp;
	int i;

	ASSERT(!( line == NULL ),"NULL as an argument of function BTC_elm_list()\n");
	ASSERT(!( n_btc < 0 ),"Number of BTC cannot be negative in BTC_elm_list()\n");
	if( n_btc == 0 )
		return NULL;
	cp = (int*) xmalloc( n_btc * sizeof( int ) );
	for( i = 0; i < n_btc; i++ )
		cp[ i ] = atoi( strtok( i == 0 ? line : NULL , " \t,;" ) );
	return cp;
}
//==============================================================================
// INITIALIZE TRANSPORT OUTPUT FILE FOR BTC
//==============================================================================
void output_transport_init_BTC(struct Transport *transport)
{
        FILE **out;
        int i,iel,sbi;
        out = open_temp_files(transport, "%s.btc", "wt" );
        for(i=0; i < 4; i++)
        	if(out[i] == NULL)
        		continue;
        	else{
        		for(sbi=0;sbi<transport->n_substances;sbi++)
        		for(iel=0;iel<transport->btc->n_BTC_elms;iel++)
        			xfprintf(out[i],"\t%d",transport->btc->BTC_elm[iel]);
        		xfprintf(out[i],"\n");
        		xfclose(out[i]);
				}
        xfree(out);
        }

/*
//==============================================================================
// TRANSPORT OUTPUT IN TIME (breakthrough curve in single element in all zones)
//==============================================================================
void output_transport_time_BTC(struct Transport *transport, double time)
{
        FILE **out;
//        ElementIter ele;
//        TNode* nod;
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    if(transport->btc == NULL)
    	return;


	char dbl_fmt[ 16 ];
	int sbi,iel,i;
	int n_subst,el;
	n_subst = transport->n_substances;
 	sprintf( dbl_fmt, "%%.%dg\t", ConstantDB::getInstance()->getInt("Out_digit"));
        out = open_temp_files(transport, "%s.btc", "at" );
        for(i=0; i < 4; i++){
        if(out[i] == NULL) continue;
        xfprintf( out[i], dbl_fmt, time);

        for(el=0;el < mesh->n_elements();el++){
		for (iel=0;iel<transport->btc->n_BTC_elms;iel++)
		if (mesh->epos_id[el] == transport->btc->BTC_elm[iel])
		{
			for( sbi = 0; sbi < n_subst; sbi++ )
                         switch(i)
                         {
                         case MOBILE:
                         xfprintf(out[0],dbl_fmt,transport->conc[sbi][MOBILE][el]);
                         break;
                         case IMMOBILE:
                         xfprintf(out[1],dbl_fmt,transport->conc[sbi][IMMOBILE][el]);
                         break;
                         case MOBILE_SORB:
                         xfprintf(out[2],dbl_fmt,transport->conc[sbi][MOBILE_SORB][el]);
                         break;
                         case IMMOBILE_SORB:
                         xfprintf(out[3],dbl_fmt,transport->conc[sbi][IMMOBILE_SORB][el]);
                         break;
                         }
                }
        }
        xfprintf( out[i], "\n" );
        xfclose( out[i] );
        }
        xfree(out);
}
*/
