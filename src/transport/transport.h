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
#include "coupling/equation.hh"
#include "transport/sources.hh"
#include "materials.hh"
#include "input/accessors.hh"
#include "flow/mh_dofhandler.hh"


struct BTC;
class OutputTime;
class Mesh;
class Distribution;
class TimeMarks;
class MaterialDatabase;
class TransportSources;
class ConvectionTransport;

//=============================================================================
// TRANSPORT
//=============================================================================
#define MOBILE 0
#define IMMOBILE 1
#define MOBILE_SORB 2
#define IMMOBILE_SORB 3

#define MAX_PHASES 4

/**
 * TODO:
 * - doxy documentation
 * - remove boundary matrix, make method for assembly of boundary source vector directly (when bc or velocity has changed) or
 *   use rescaling for change of time step
 * - make separate method for changing time step (rescaling only)
 */




class ConvectionTransport : public EquationBase {
	friend class TransportSources;
public:
    /**
     * Constructor.
     */
	ConvectionTransport(TimeMarks &marks,  Mesh &init_mesh, MaterialDatabase &material_database, const Input::Record &in_rec);

	/**
	 * TODO: destructor
	 */
	virtual ~ConvectionTransport();

	/**
	 * Calculates one time step of explicit transport.
	 */
	void compute_one_step();

	//void transport_until_time(double time_interval);
	/**
	 * Use new flow field vector for construction of convection matrix.
	 * Updates CFL time step constrain.
	 */
	void set_flow_field_vector(const MH_DofHandler &dh);

    /**
     * Set time interval over which we should use fixed transport matrix. Rescale transport matrix.
     */
    void set_target_time(double target_time);

	/**
	 * Read time dependent boundary condition and update boundary source vectors.
	 */
	void read_bc_vector();
	/**
	 * Communicate parallel concentration vectors into sequential output vector.
	 */
	void output_vector_gather(); //
	/**
	 * Returns time step constrain given by CFL condition for the discretization of the
	 * convection term. The constrain depends on actual convection matrix assembled by
	 * create_transport_matrix_mpi()
	 */
	//double get_cfl_time_constrain();
	double ***get_concentration_matrix();
	double ***get_prev_concentration_matrix();
	void get_par_info(int * &el_4_loc, Distribution * &el_ds);
	bool get_dual_porosity();
	int get_n_substances();
	int *get_el_4_loc();
	virtual void get_parallel_solution_vector(Vec &vc);
	virtual void get_solution_vector(double* &vector, unsigned int &size);
	/**
	 * Return pointer to sequential arrays for output.
	 * TODO: Maybe this should be made by get_solution_vector, but here we have matrix of arrays.
	 */
	double ***get_out_conc();
    vector<string> &get_substance_names();
    TransportSources *transportsources;
    const MH_DofHandler *mh_dh;


private:


    /**
     * Assembly convection term part of the matrix and boundary matrix for application of boundary conditions.
     *
     * Discretization of the convection term use explicit time scheme and finite volumes with full upwinding.
     * We count on with exchange between dimensions and mixing on edges where more then two elements connect (can happen for 2D and 1D elements in
     * 3D embedding space)
     *
     * In order to get multiplication matrix for explicit transport one have to scale the convection part by the acctual time step and
     * add time term, i. e. unit matrix (see. transport_matrix_step_mpi)
     *
     * Updates CFL time step constrain.
     *
     * TODO: Remove assembling of boundary matrix and assembly directly boundary source vector.
     * TODO: Do not use velocty values stored in mesh. Use separate velocity vector.
     */
    void create_transport_matrix_mpi();
	void make_transport_partitioning(); //
//	void alloc_transport(struct Problem *problem);
	void read_initial_condition(string fname); //
	void read_concentration_sources();

	/**
	 * Compose file name for boundary condition at given level.
	 * For level -1 use original fixed file name.
	 */
	std::string make_bc_file_name(int level);
	/**
	 * Finish explicit transport matrix (time step scaling)
	 */
	void transport_matrix_step_mpi(double time_step); //

	void transport_dual_porosity( int elm_pos, MaterialDatabase::Iter material, int sbi); //
	void transport_sorption(int elm_pos, MaterialDatabase::Iter mtr, int sbi); //
	void compute_sorption(double conc_avg, vector<double> &sorp_coef, int sorp_type, double *concx, double *concx_sorb, double Nv,
	        double N); //
	//void compute_concentration_sources(int sbi);

//	void get_reaction(int i,oReaction *reaction); //
	//void transport_output(struct Transport *transport, double time, int frame);
//	void transport_output_init(struct Transport *transport);
//  void transport_output_finish(struct Transport *transport);
	//DENSITY
	int compare_dens_iter(); //
	void restart_iteration_C(); //
	void save_restart_iteration_H(); //
	void save_time_step_C(); //
	void alloc_transport_vectors(); //
	void alloc_density_vectors(); //
	void alloc_transport_structs_mpi(); //
	void subst_names(char *line); //
	void subst_scales(char *line); //

	//void get_names(string* array);

	bool is_convection_matrix_scaled;

    // TODO: Make simplified TimeGovernor and move following into it.
    //double time_step;                     ///< Time step for computation.
    //double max_step;		              ///< Time step constrain given by CFL condition.
    //unsigned int time_level;              ///< Number of computed time steps.

    //std::vector<double> bc_times;       ///< Times of reading time dependent boundary condition. Initial boundary condition is zero.
	TimeMark::Type bc_mark_type_;
    unsigned int bc_time_level;         ///< Index into bc_times vector.

    TimeMark::Type target_mark_type;    ///< TimeMark type for time marks denoting end of every time interval where transport matrix remains constant.
    double cfl_max_step;
            // only local part
            double ***conc;
            double ***pconc;
            double **cumulative_corr;

            // global

//        	std::string		transport_out_fname;// Name of file of trans. output
        	int              dens_step;            //
        	double 			update_dens_time;

            double ***out_conc;
            int              n_substances;    // # substances transported by water
            vector<string> substance_name;   // Names of substances
            string bc_fname; // name of input file with boundary conditions

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
            bool              sorption;     // Include sorption  YES/NO
            bool              dual_porosity;   // Include dual porosity YES/NO
            bool              reaction_on;     // Include reaction  YES/NO

            // Other
           // struct Problem* problem;
            int sub_problem;	// 0-only transport,1-transport+dual porosity,
    							// 2-transport+sorption
    							// 3-transport+dual porosity+sorption
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
            Vec *vcumulative_corr;


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
#endif /* TRANSPORT_H_ */
