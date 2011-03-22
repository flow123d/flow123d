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
#include "problem.h"
//#include "reaction.h"


struct BTC;

//=============================================================================
// TRANSPORT
//=============================================================================
#define MOBILE 0
#define IMMOBILE 1
#define MOBILE_SORB 2
#define IMMOBILE_SORB 3

#define MAX_PHASES 4




class ConvectionTransport {
public:
	ConvectionTransport(struct Problem *problem);
	~ConvectionTransport();
	void make_transport(); //
	void make_transport_partitioning(); //
//	void alloc_transport(struct Problem *problem);
	void transport_init(); //
	void read_initial_condition(); //
	void alloc_transport_vectors(); //
	void alloc_density_vectors(); //
	void alloc_transport_structs_mpi(); //
	void fill_transport_vectors_mpi(); //
	void create_transport_matrix_mpi(); //
	void transport_matrix_step_mpi(double time_step); //
	void calculate_bc_mpi(); //
	void transport_step_mpi(Mat *tm, Vec *conc, Vec *pconc, Vec *bc);
	void transport_dual_porosity( int elm_pos, MaterialDatabase::Iter material, int sbi); //
	void transport_sorption(int elm_pos, MaterialDatabase::Iter mtr, int sbi); //
	void compute_sorption(double conc_avg, vector<double> &sorp_coef, int sorp_type, double *concx, double *concx_sorb, double Nv,
	        double N); //
	void convection(); //
	void output_vector_gather(); //
//	void get_reaction(int i,oReaction *reaction); //
	//void transport_output(struct Transport *transport, double time, int frame);
//	void transport_output_init(struct Transport *transport);
//  void transport_output_finish(struct Transport *transport);
	//DENSITY
	int compare_dens_iter(); //
	void restart_iteration_C(); //
	void save_restart_iteration_H(); //
	void save_time_step_C(); //

    // global
    double ***out_conc;
	int		 		 n_substances;    // # substances transported by water
	char			**substance_name;	// Names of substances
	char            *transport_out_fname;// Name of file of trans. output
	int              dens_step;            //
	double 			update_dens_time;

protected:
    Mesh *mesh;
    MaterialDatabase *mat_base;


private:
	void subst_names(char *line); //
	void subst_scales(char *line); //

    //		double           stop_time;       // Number of time steps
    //		double           save_step;       // Step for outputing results
    		double           time_step;       // Time step for computation
    		double			 max_step;		// bounded by CFL

          //  struct TMatrix* tmatrix;

            // only local part
            double ***conc;
            double ***pconc;



            //double ****node_conc;	// zatim nepouzite


        	//Density
            bool density;			// Density Yes/NO
            double ***prev_conc;
            double *scalar_it;
            double          *substance_density_scale;

        	int		max_dens_it;	// maximum number of iterations in the variable-density iteration cycle
        	bool		dens_implicit; // use iterations for the variable density (implicit) YES/NO
        	double	dens_eps;      // stopping criterium for density iterations - pressure
        	bool		write_iterations; // write results during iterations to transport POS file YES/NO

        	// Transport
        	bool             transport_on;    // Compute transport YES/NO
        	//unsigned int 			n_elements; 		// number of elements

        	char            *concentration_fname;// Name of file of concentration
            bool              sorption;     // Include sorption  YES/NO
            bool              dual_porosity;   // Include dual porosity YES/NO
            bool              reaction_on;     // Include reaction  YES/NO
        	char            *transport_bcd_fname;// Name of file of transport bcd

            char            *transport_out_im_fname;// Name of file of trans. immobile output
            char            *transport_out_sorp_fname;// Name of file of trans. output
            char            *transport_out_im_sorp_fname;// Name of file of trans. immobile output


            int             transport_sub_problem;

            // Other
            struct Problem* problem;
            int sub_problem;	// 0-only transport,1-transport+dual porosity,
    							// 2-transport+sorption
    							// 3-transport+dual porosity+sorption
    //REACTION
            /*
           oReaction *react;
           int n_reaction;
           */
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


void alloc_transport(struct Problem *problem);
void transport_output(double ***out_conc,char **subst_names ,int n_subst,double time,int frame); //NO
void transport_output_init(char *transport_out_fname); // NO
void transport_output_finish(char *transport_out_fname); // NO
void transport_output(double ***out_conc,char **subst_names ,int n_subst,double time, int frame, char *transport_out_fname); //NO



#endif /* TRANSPORT_H_ */
