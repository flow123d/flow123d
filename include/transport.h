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

#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include <petscmat.h>

#include "par_distribution.hh"
#include "mesh.h"

struct Problem;
struct TMatrix;
struct BTC;
struct Reaction;
//=============================================================================
// TRANSPORT
//=============================================================================
#define MOBILE 0
#define IMMOBILE 1
#define MOBILE_SORB 2
#define IMMOBILE_SORB 3

#define MAX_PHASES 4

struct Transport
{


//		double           stop_time;       // Number of time steps
//		double           save_step;       // Step for outputing results
		double           time_step;       // Time step for computation
		double			 max_step;		// bounded by CFL

        struct TMatrix* tmatrix;

        // only local part
        double ***conc;
        double ***pconc;

        // global
        double ***out_conc;
        //double ****node_conc;	// zatim nepouzite


    	//Density
        bool density;			// Density Yes/NO
		int              dens_step;            //
		double update_dens_time;
        double ***prev_conc;
        double *scalar_it;

    	int		max_dens_it;	// maximum number of iterations in the variable-density iteration cycle
    	bool		dens_implicit; // use iterations for the variable density (implicit) YES/NO
    	double	dens_eps;      // stopping criterium for density iterations - pressure
    	bool		write_iterations; // write results during iterations to transport POS file YES/NO

    	// Transport
    	bool             transport_on;    // Compute transport YES/NO
    	int		 		 n_substances;    // # substances transported by water
    	char            *concentration_fname;// Name of file of concentration
        bool              sorption;     // Include sorption  YES/NO
        bool              dual_porosity;   // Include dual porosity YES/NO
        bool              reaction_on;     // Include reaction  YES/NO
    	char            *transport_bcd_fname;// Name of file of transport bcd


    	char            *transport_out_fname;// Name of file of trans. output
        char            *transport_out_im_fname;// Name of file of trans. immobile output
        char            *transport_out_sorp_fname;// Name of file of trans. output
        char            *transport_out_im_sorp_fname;// Name of file of trans. immobile output
    	char			**substance_name;	// Names of substances
        double          *substance_density_scale;
        int             transport_sub_problem;

        // Other
        struct Problem* problem;
        int sub_problem;	// 0-only transport,1-transport+dual porosity,
							// 2-transport+sorption
							// 3-transport+dual porosity+sorption
//REACTION
       struct Reaction *reaction;
       int n_reaction;
//BTC
        struct BTC		*btc;

//DECOVALEX
        struct FSection *fsec;

//PEPA
        int 	pepa; // It enables Pepa Chudoba's  crazy functions
        int 	type; // Type of crazy function


// NEW TRANSPORT

        VecScatter vconc_out_scatter;
        Mat tm; // PETSc transport matrix
        Mat bcm; // PETSc boundary condition matrix
        Vec *bcv; // boundary condition vector
        Vec *bcvcorr; // boundary condition correction vector

        Vec *vconc; // concentration vector
        Vec *vpconc; // previous concentration vector

        //Vec *vconc_im; // immobile concentration vector
        //Vec *vconc_so; // sorbed concentration vector
        //Vec *vconc_im_so; // immobile sorbed concentration vector

        Vec *vconc_out; // concentration vector output (gathered)


        int **d_row;  // diagonal row entries number in tm
        int **od_row; // off-diagonal row entries number in tm
        int **db_row; // diagonal column entries number in bcm
        int **odb_row; // off-diagonal column entries number in bcm
        int *l_row; // number of local rows in tm and bcm
        int *lb_col; // number of local columns in bcm
        
        int *row_4_el;
        int *el_4_loc;
        Distribution *el_ds;
};
//=============================================================================
// TRANSPORT MATRIX
//=============================================================================
struct TMatrix
{
        int* lenrow;
        int* irowst;
        int* jcn;
        double* val;
        int rows;
};

struct SVector{
        int* pos;
        double* val;
};


void transport_step_mpi(Mat *tm,Vec *conc,Vec *pconc,Vec *bc);

void alloc_transport(struct Problem *problem);
void transport_init(struct Problem *problem);
void make_transport(struct Transport *transport);
void convection(struct Transport *transport);
void transport_output(struct Transport *transport,double time,int frame);
void transport_output_init(struct Transport *transport);
void transport_output_finish(struct Transport *transport);
//void transport_partioning(struct Transport *transport,int np);

//DENSITY
int compare_dens_iter(struct Problem *problem);
void restart_iteration_C(struct Problem *problem);
void save_restart_iteration_H(struct Problem *problem);
void save_time_step_C(struct Problem *problem);


#endif /* TRANSPORT_H_ */
