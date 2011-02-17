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

#ifndef MAKE_PROBLEM_H
#define MAKE_PROBLEM_H

// types of goals
#define COMPUTE_MH         1
#define CONVERT_TO_POS     2
//#define COMPUTE_MH_NEW     3

// problem types
#define STEADY_SATURATED   1
#define UNSTEADY_SATURATED 2
#define PROBLEM_DENSITY    3

struct Solver;
struct DarcyFlowMH;
class Mesh;
struct Reaction;
struct Transport;
class MaterialDatabase;

//=============================================================================
// STRUCTURE OF THE SOLVED TASK
//=============================================================================

struct Problem {
    // Global
    //double stop_time; // Number of time steps
    //double save_step; // Step for outputing results
    //double time_step; // Time step for computation
    //	int              dens;            // Density Yes/NO
    //	bool              dens_step;            //

    // Output
    char *out_fname_2; // Name of output file of type 2
    // Chemistry
    bool semchemie_on; //Enable to compute chemistry, YES/NO, NO defalut

    struct Pos_view_params *pos_view_params;

    /* tohle se musi okomentovat !!!
            // NEW
            int		ftrans_out;       // velocity field as coeficients for FTRANS code
            bool              cross_section;   // Include cross_section YES/NO     //jh
            char            *cs_params;       // Params for cross_section         //jh
            struct Section  *section;         // Geometric data of cross-section  //jh
            int              res_run;         // Whether make an output file with running results  //jh
            int              res_fin;         // Whether make an output file with final results    //jh
            int              specify_elm_output;   // Yes if want to choose only one type of elements  //jh temp
            int              output_elm_type; // Output elm type - 1 = 1D, 2 = 2D, 3 = 3D //jh temp

            struct FSection  *fsec;
            char            *fsec_params;

          //  struct ConcFlow  *concflow;
          //  char            *CF_params;
     */

    //Transport
    struct Transport *transport;

    // Reaction
    struct Reaction *react;
    struct Reaction **reaction_hash;

    struct DarcyFlowMH *water; // Global MHsystem of the system

    MaterialDatabase *material_database;

    struct Read_ini *ini;
};

//=============================================================================
// STRUCTURE OF THE POS view params             // jh
//=============================================================================

struct Pos_view_params {
    double x_ang, y_ang, z_ang; // rotation angles
    double x_sca, y_sca, z_sca; // scaling
    double x_tra, y_tra; //        translation
};

void problem_init(struct Problem* problem);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
