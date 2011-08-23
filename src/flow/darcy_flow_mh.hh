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
 * @brief mixed-hybrid model of linear Darcy flow, possibly unsteady.
 *
 * Main object for mixed-hybrid discretization of the linear elliptic PDE (Laplace)
 * on  a multidimensional domain. Discretization of saturated Darcy flow using
 * RT0 approximation for the velocity
 *
 *
 *
 *  @author Jan Brezina
 *
 */

/*
 * list of files dependent on this one:
 *
 * posprocess.cc
 * problem.cc
 * main.hh
 * transport.cc
 */


#ifndef DARCY_FLOW_MH_HH
#define DARCY_FLOW_MH_HH

#include <petscmat.h>
#include "system/sys_vector.hh"
#include "time_governor.hh"
#include <field_p0.hh>
#include <materials.hh>
#include "equation.hh"

/// external types:
class LinSys;
struct Solver;
class Mesh;
class SchurComplement;
class Distribution;
class SparseGraph;


/**
 * @brief Mixed-hybrid model of linear Darcy flow, possibly unsteady.
 *
 * Abstract class for various implementations of Darcy flow. In future there should be
 * one further level of abstraction for general time dependent problem.
 *
 * maybe TODO:
 * split compute_one_step to :
 * 1) prepare_next_timestep
 * 2) actualize_solution - this is for iterative nonlinear solvers
 *
 */

class DarcyFlowMH : public EquationBase {
public:
    DarcyFlowMH(TimeMarks &marks, Mesh &mesh, MaterialDatabase &mat_base)
    : EquationBase(marks, mesh, mat_base)
    {}

    FieldP0<double>  * get_sources()
        { return sources; }

    void get_velocity_seq_vector(Vec &velocity_vec)
        { velocity_vec = velocity_vector; }

protected:
    void setup_velocity_vector() {
        double *velocity_array;
        unsigned int size;

        get_solution_vector(velocity_array, size);
        VecCreateSeqWithArray(PETSC_COMM_SELF, mesh_->n_sides, velocity_array, &velocity_vector);

    }

    virtual void postprocess() =0;

    //virtual void balance();
    //virtual void integrate_sources();

protected:
    FieldP0<double> *sources;
    Vec velocity_vector;
};


/**
 * @brief Mixed-hybrid of steady Darcy flow with sources and variable density.
 *
 * solve equations:
 * @f[
 *      q= -{\mathbf{K}} \nabla h -{\mathbf{K}} R \nabla z
 * @f]
 * @f[
 *      \mathrm{div} q = f
 * @f]
 *
 * where
 * - @f$ q @f$ is flux @f$[ms^{-1}]@f$ for 3d, @f$[m^2s^{-1}]@f$ for 2d and @f$[m^3s^{-1}]@f$ for 1d.
 * - @f$ \mathbf{K} @f$ is hydraulic tensor ( its orientation for 2d, 1d case is questionable )
 * - @f$ h = \frac{\pi}{\rho_0 g}+z @f$ is pressures head, @f$ \pi, \rho_0, g @f$ are the pressure, water density, and acceleration of gravity , respectively.
 *   Assumes gravity force acts counter to the direction of the @f$ z @f$ axis.
 * - @f$ R @f$ is destity or gravity variability coefficient. For density driven flow it should be
 * @f[
 *    R = \frac{\rho}{\rho_0} -1 = \rho_0^{-1}\sum_{i=1}^s c_i
 * @f]
 *   where @f$ c_i @f$ is concentration in @f$ kg m^{-3} @f$.
 */
class DarcyFlowMH_Steady : public DarcyFlowMH
{
public:
    DarcyFlowMH_Steady(TimeMarks &marks,Mesh &mesh, MaterialDatabase &mat_base_in);
    virtual void update_solution();
    virtual void get_solution_vector(double * &vec, unsigned int &vec_size);
    virtual void get_parallel_solution_vector(Vec &vector);
    virtual void postprocess() {};
    ~DarcyFlowMH_Steady();

protected:
    virtual void modify_system() {};
    void set_R() {};
    void prepare_parallel();
    void make_row_numberings();
    void preallocate_mh_matrix();
    void assembly_steady_mh_matrix();
    void make_schur0();
    void make_schur1();
    void make_schur2();


	int size;				// global size of MH matrix
	int  n_schur_compls;  	// number of shur complements to make
	double  *solution; 			// sequantial scattered solution vector

	struct Solver *solver;

	LinSys *schur0;  		// whole MH Linear System
	SchurComplement *schur1;  	// first schur compl.
	SchurComplement *schur2;	// second ..



	// parallel
	int np;  // number of procs
	int myp; // my proc number
	int	 lsize;				// local size of whole MH matrix
	Distribution *edge_ds; // optimal distribution of edges
	Distribution *el_ds; // optimal distribution of elements
	Distribution *side_ds; // optimal distribution of elements
	Distribution *rows_ds; // final distribution of rows of MH matrix

	int *el_4_loc;		// array of idexes of local elements (in ordering matching the optimal global)
	int *row_4_el;		// element index to matrix row
	int *side_id_4_loc;		// array of ids of local sides
	int	*side_row_4_id;		// side id to matrix row
	int *edge_4_loc;		// array of indexes of local edges
	int	*row_4_edge;		// edge index to matrix row
	//int *old_4_new;        // aux. array should be only part of parallel LinSys

	// MATIS related arrays
	int ndof_loc;                   // size of local block of MATIS matrix 
	int *global_row_4_sub_row;      // global dof index for subdomain index
	ISLocalToGlobalMapping map_side_local_to_global; ///< PETSC mapping form local SIDE indices of subdomain to global indices

	// gather of the solution
	Vec sol_vec;			// vector over solution array
	VecScatter par_to_all;
};


void make_element_connection_graph(Mesh *mesh, SparseGraph * &graph,bool neigh_on = false);
void id_maps(int n_ids, int *id_4_old, const Distribution &old_ds,
        int *loc_part, Distribution * &new_ds, int * &id_4_loc, int * &new_4_id);
void mat_count_off_proc_values(Mat m, Vec v);



/**
 * @brief Mixed-hybrid solution of unsteady Darcy flow.
 *
 * Standard discretization with time term and sources picewise constant
 * on the element. This leads to violation of the discrete maximum principle for
 * non-acute meshes or to too small timesteps. For simplicial meshes this can be solved by lumping to the edges. See DarcyFlowLMH_Unsteady.
 */

class DarcyFlowMH_Unsteady : public DarcyFlowMH_Steady
{
public:
    DarcyFlowMH_Unsteady(TimeMarks &marks,Mesh &mesh, MaterialDatabase &mat_base_in);
    DarcyFlowMH_Unsteady();
protected:
    virtual void modify_system();
    void setup_time_term();
private:
    Vec steady_diagonal;
    Vec steady_rhs;
    Vec new_diagonal;
    Vec previous_solution;

};

/**
 * @brief Edge lumped mixed-hybrid solution of unsteady Darcy flow.
 *
 * The time term and sources are evenly distributed form an element to its edges.
 * This applies directly to the second Schur complement. After this system for pressure traces is solved we reconstruct pressures and side flows as follows:
 *
 * -# Element pressure is  average of edge pressure. This is in fact same as the MH for steady case so we let SchurComplement class do its job.
 *
 * -# We let SchurComplement to reconstruct fluxes and then account time term and sources which are evenly distributed from an element to its sides.
 *    It can be proved, that this keeps continuity of the fluxes over the edges.
 *
 * This lumping technique preserves discrete maximum principle for any time step provided one use acute mesh. But in practice even worse meshes are tractable.
 */
class DarcyFlowLMH_Unsteady : public DarcyFlowMH_Steady
{
public:
    DarcyFlowLMH_Unsteady(TimeMarks &marks,Mesh &mesh, MaterialDatabase &mat_base_in);
    DarcyFlowLMH_Unsteady();
protected:
    virtual void modify_system();
    virtual void postprocess();
private:
    Vec steady_diagonal;
    Vec steady_rhs;
    Vec new_diagonal;
    Vec previous_solution;
    Vec time_term;
};

#endif  //DARCY_FLOW_MH_HH
//-----------------------------------------------------------------------------
// vim: set cindent:
