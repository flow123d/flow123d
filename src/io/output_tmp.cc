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
 * @file    output_tmp.cc
 * @brief   The functions for outputs to EVIL file formats.
 *
 */

#if 0

#include "io/output.h"
#include "xio.h"
#include "mesh.h"
//#include "problem.h"

/**
 * \brief This function does some initialization of evil file format.
 */
void output_flow_field_init(char *fname)
{
    FILE *out;

    if( OptGetBool("Output", "Write_output_file", "no") == false )
        return;

    xprintf( Msg, "%s: Writing output files... %s ", __func__, fname);

    out = xfopen( fname, "wt" );
    xfprintf( out, "$DataFormat\n" );
    xfprintf( out, "1.0 0 %d\n", sizeof( double ) );
    xfprintf( out, "$EndDataFormat\n" );
    xfclose( out );

    xprintf( Msg, "O.K.\n");
}

/**
 * \brief This function write data to evil file format.
 */
void output_flow_field_in_time(struct Problem *problem, double time)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    int i,cit;
    ElementIter ele;
    FILE *out;
    char dbl_fmt[ 16 ];

    ASSERT(!( problem == NULL ),"NULL as argument of function output_flow_field_in_time()\n");
    if( OptGetBool("Output", "Write_output_file", "no") == false )
        return;

    xprintf( Msg, "%s: Writing output file %s ... ", __func__, problem->out_fname_2);

    out = xfopen( problem->out_fname_2, "at" );
    sprintf( dbl_fmt, "%%.%dg ", ConstantDB::getInstance()->getInt("Out_digit"));
    xfprintf( out, "$FlowField\n" );
    xfprintf( out, "T = ");
    xfprintf( out, dbl_fmt, time);
    xfprintf( out, "\n%d\n", mesh->n_elements() );
    cit = 0;

    FOR_ELEMENTS( ele ) {
        xfprintf( out, "%d ", cit);
        xfprintf( out, "%d ", ele.id());
        xfprintf( out, dbl_fmt, ele->scalar);
        xfprintf( out, " %d ", ele->n_sides);
        for (i = 0; i < ele->n_sides; i++)
            xfprintf( out, dbl_fmt, ele->side[i]->scalar);
        xfprintf( out, "\t");
        for (i = 0; i < ele->n_sides; i++)
            xfprintf( out, dbl_fmt, ele->side[i]->flux);
        xfprintf( out, "\t");
        xfprintf( out, "%d ", ele->n_neighs_vv);
        for (i = 0; i < ele->n_neighs_vv; i++)
            xfprintf( out, "%d ", ele->neigh_vv[i]->id);
        xfprintf( out, "\n" );
        cit ++;
    }
    xfprintf( out, "$EndFlowField\n\n" );
    xfclose( out );

    xprintf( Msg, "O.K.\n");
}

#endif
