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
 * @brief  Setup of the Problem data
 *
 */

#include "mesh/ini_constants_mesh.hh"
#include "constantdb.h"

#include <stdlib.h>    // atof()
#include <strings.h>
#include <string.h>    // strcmpi()

#include "system.hh"
#include "xio.h"
#include "mesh.h"
#include "materials.hh"
#include "darcy_flow_mh.hh"
#include "solve.h"
#include "output.h"
#include "problem.h"
#include "read_ini.h"
#include "constantdb.h"

static void check_ini_values(struct Problem*);
static struct Pos_view_params *pos_view_par_read( char *view_char );
static int pos_format_ID(char * format);
/*static void init_log_file(struct Problem*);*/

//=============================================================================
// Initialize a Problem - by values from opt/ini file; check values
//=============================================================================
void problem_init(struct Problem *problem)
{
    char *view_char;
    F_ENTRY;

    // [Global]
    ConstantDB::getInstance()->setInt("Problem_type", OptGetInt("Global", "Problem_type", NULL));

    // problem -> dens         = OptGetBool( "Global", "Density_on", "no" );


    // [Output]
    ConstantDB::getInstance()->setInt("Out_digit", OptGetInt("Output", "Output_digits", "6"));

    if(ConstantDB::getInstance()->getInt("Goal") == COMPUTE_MH) {
        problem->out_fname_2 = OptGetFileName("Output", "Output_file_2", NULL);
    }

    if(problem->out_fname_2[0] == '\\')
        problem->out_fname_2 = OptGetFileName("Output", "Output_file", NULL);

    ConstantDB::getInstance()->setInt("Out_file_type", OptGetInt("Output", "Output_file_type", "1"));

    if(ConstantDB::getInstance()->getInt("Out_file_type") == 1)
        problem->out_fname_2 = OptGetFileName("Output", "Output_file", NULL);

    ConstantDB::getInstance()->setChar("Pos_format", OptGetStr("Output", "POS_format", "ASCII"));
    ConstantDB::getInstance()->setInt("Pos_format_id", pos_format_ID( (char*)ConstantDB::getInstance()->getChar("Pos_format") ));

    // TODO: Proper implementation of cross section
#if 0
    problem->ftrans_out       = get_b( "Output", "Write_ftrans_out", false );
    problem->cross_section    = get_b( "Output", "Cross_section", false );         //jh
    problem->cs_params        = get_s( "Output", "Cs_params", "0 0 0 0 0 0 0" );        //jh
//    problem->res_run          = get_b( "Output", "Cs_results_run", false );           //jh
//    problem->res_fin          = get_b( "Output", "Cs_results_final", false );
    problem->specify_elm_output =  get_b( "Output", "Specify_elm_type", false );   //jh temp
    problem->output_elm_type  = get_i( "Output", "Output_elm_type", 1 );        //jh temp
    problem->fsec_params       = get_s( "Output", "FCs_params", "0 0 0 0 0" );
//    problem->CF_params         = get_s( "Output", "ConfFlow_params", "0");
#endif

    view_char = OptGetStr("Output", "POS_view_params", "0 0 0 0 1 1 1 0 0");
    problem->pos_view_params = pos_view_par_read(view_char);

    // [Constants]
    ConstantDB::getInstance()->setDouble("G", OptGetDbl("Constants", "g", "1.0"));
    ConstantDB::getInstance()->setDouble("Rho", OptGetDbl("Constants", "rho", "1.0"));

    // [Transport]
    if(ConstantDB::getInstance()->getInt("Goal") == COMPUTE_MH) {
        ConstantDB::getInstance()->setChar("Concentration_fname", OptGetFileName("Transport", "Concentration", "\\"));
        ConstantDB::getInstance()->setChar("Transport_bcd_fname", OptGetFileName("Transport", "Transport_BCD", "\\"));
    }

    // Material Database
    problem->material_database = new MaterialDatabase(OptGetStr( "Input", "Material", "\\" ));



    // Initialize sub structures by NULL
    problem->water = NULL;
    problem->transport = NULL;

    check_ini_values(problem);
}
//=============================================================================
//
//=============================================================================
// TODO : this should be removed - each test should be done where the variable is used
void check_ini_values( struct Problem *problem )
{
	ASSERT(NONULL( problem ),"NULL as argument of function check_ini_values()\n");

        int type = ConstantDB::getInstance()->getInt("Problem_type");
	INPUT_CHECK(
                 ( type == STEADY_SATURATED ) ||
                 ( type == UNSTEADY_SATURATED ) ||
                 ( type == PROBLEM_DENSITY),
                 "Unsupported type of problem: %d\n", type );
	if( strcmpi( OptGetStr("Input", "Mesh", "\\"), "\\"  ) == 0 )
		xprintf(UsrErr,"Name of mesh file must be defined\n");
	if( strcmpi( OptGetStr( "Input", "Material", "\\" ), "\\"  ) == 0 )
		xprintf(UsrErr,"Name of material properties file must be defined\n");
	if( strcmpi( OptGetStr( "Input", "Boundary", "\\" ), "\\"  ) == 0 )
		xprintf(UsrErr,"Name of boundary condition file must be defined\n");
	if( strcmpi( OptGetStr( "Input", "Neighbouring", "\\" ), "\\"  ) == 0 )
		xprintf(UsrErr,"Name of file describing neighbouring must be defined\n");
	// if( OptGetStr( "Input", "Sources", "\\" ), "\\"  ) == 0 )
	//	problem->sources_fname = NULL;
	if( strcmpi( OptGetFileName("Output", "Output_file", "\\"), "\\"  ) == 0 )
		xprintf(UsrErr,"Name of output file must be defined\n");
	INPUT_CHECK(!( (ConstantDB::getInstance()->getInt("Out_digit") < 1)
                || (ConstantDB::getInstance()->getInt("Out_digit") > 16) ), "Number of digits of output must be between %d and %d\n", 1, 16);
	INPUT_CHECK( DBL_GT(ConstantDB::getInstance()->getDouble("G"), 0.0), "Gravitotional acceleration has to be greater than ZERO\n");
	INPUT_CHECK( DBL_GT(ConstantDB::getInstance()->getDouble("Rho"), 0.0), "Density of fluid has to be greater than ZERO\n");


    //TODO: proc je tohle zakomentovane? vypada to smysluplne...
    //INPUT_CHECK(!( problem->transport_on == true && problem->n_substances < 1 ),"Number of substances must be positive\n");

    if( (OptGetBool("Transport", "Transport_on", "no") != true) && (type == PROBLEM_DENSITY )) //|| problem->dens == true ) )
		xprintf(UsrErr,"Transport must be ON for variable-density calculation\n");


/*	if( (problem->dens_eps < ZERO) || (problem->dens_eps > 1e2) )
		xprintf(UsrErr,"Accuracy of solver be between %f and %f\n");*/
	// i will add new error msg later
}
//=============================================================================
//
//=============================================================================
/*void init_log_file( struct Problem *problem )
{
	ASSERT(!( problem == NULL ),"NULL as argument of function init_log_file()\n");
	ASSERT(!( problem->log_fname == NULL ),"NULL as name of log file in function init_log_file()\n");
	xremove( problem->log_fname );
	problem->log_init = true;
}*/
//=============================================================================
//
//=============================================================================
struct Pos_view_params *pos_view_par_read( char *view_char )
{
//  Bug: Formlorn message O.K. in output
//  DF add specification wher the O.K. is writen.
    xprintf( Msg, "problem.cc  (struct Pos_view_params *pos_view_par_read( char *view_char )): " );

	struct Pos_view_params* par = (struct Pos_view_params*) xmalloc( sizeof( struct Pos_view_params ) );

	if( sscanf( view_char, "%lf %lf %lf %lf %lf %lf %lf %lf",
                &par->x_ang,
                &par->y_ang,
                &par->z_ang,
                &par->x_sca,
                &par->y_sca,
                &par->z_sca,
                &par->x_tra,
                &par->y_tra ) == 8 )
        	xprintf( Msg, "O.K.\n")/*orig verb 2*/;


	return par;
}
//=============================================================================
//
//=============================================================================
int pos_format_ID(char* format){

        if(strcmp(format,"ASCII") == 0)
            return POS_ASCII;
        if(strcmp(format,"BIN") == 0)
            return POS_BIN;
        if(strcmp(format, "VTK_SERIAL_ASCII") == 0)
            return VTK_SERIAL_ASCII;
        if(strcmp(format, "VTK_PARALLEL_ASCII") == 0)
            return VTK_PARALLEL_ASCII;
        xprintf(UsrErr,"Unknown file format: %s.\n", format );
        return (-1);
}


//-----------------------------------------------------------------------------
// vim: set cindent:
