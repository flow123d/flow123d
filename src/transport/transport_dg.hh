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
 * @brief Discontinuous Galerkin method for equation of transport with dispersion.
 *  @author Jan Stebel
 */

#ifndef TRANSPORT_DG_HH_
#define TRANSPORT_DG_HH_

#include "transport_operator_splitting.hh"
#include "la/linsys.hh"
#include "fem/dofhandler.hh"
#include "fem/finite_element.hh"
#include "fem/fe_values.hh"




class TransportDG : public TransportBase
{
public:

    /**
     * Constructor.
     * @param marks
     * @param init_mesh
     * @param material_database
     */
    TransportDG(TimeMarks &marks,  Mesh &init_mesh, MaterialDatabase &material_database);

    /**
     * Computes solution in one time instant.
     */
	void update_solution();

	/**
	 * Returns the serialized solution array.
	 * @param vector the solution
	 * @param size   size of the array
	 */
	void get_solution_vector(double * &vector, unsigned int &size);

	/**
	 * Returns the (possibly) parallel solution vector.
	 * @param vector
	 */
	void get_parallel_solution_vector(Vec &vector);

	/**
	 * Updates the velocity field which determines some coefficients of the transport equation.
	 * (So far it does not work since the flow module returns a vector of zeros.)
	 * @param velocity_vector
	 */
	void set_velocity_field(Vec &velocity_vector);

	/**
	 * Postprocesses the solution and writes to output file.
	 */
	void output_data();

	/**
	 * Destructor.
	 */
	~TransportDG();

private:

	/**
	 * Assembles the mass matrix for the given dimension.
	 * The DOF handler and FiniteElement objects of specified dimension
     * must be passed as arguments.
	 * @param dh DOF handler.
	 * @param fe FiniteElement
	 */
	template<unsigned int dim>
	void assemble_mass_matrix(DOFHandler<dim,3> *dh, FiniteElement<dim,3> *fe);

	/**
	 * Assembles the stiffness matrix for elements of dimension @p dim.
	 * The DOF handler and FiniteElement objects of specified dimension
	 * must be passed as arguments.
	 * @param dh DOF handler.
	 * @param dh_sub DOF handler for sides.
	 * @param fe FiniteElement.
	 * @param fe_sub FiniteElement for sides.
	 */
	template<unsigned int dim>
	void assemble(DOFHandler<dim,3> *dh, DOFHandler<dim-1,3> *dh_sub, FiniteElement<dim,3> *fe, FiniteElement<dim-1,3> *fe_sub);


	/**
	 * Assembles the boundary fluxes into the stiffness matrix.
	 * The DOF handler and FiniteElement objects of specified dimension
	 * must be passed as arguments.
	 * @param dh DOF handler.
	 * @param dh_sub DOF handler for sides.
	 * @param fe FiniteElement.
	 * @param fe_sub FiniteElement for sides.
	 */
	template<unsigned int dim>
	void assemble_fluxes(DOFHandler<dim,3> *dh, DOFHandler<dim-1,3> *dh_sub, FiniteElement<dim,3> *fe, FiniteElement<dim-1,3> *fe_sub);

	/**
	 * Assembles the r.h.s. components corresponding to the Dirichlet
	 * boundary conditions.
	 * @param dh DOF handler.
     * @param fe FiniteElement.
	 */
	template<unsigned int dim>
	void set_boundary_conditions(DOFHandler<dim,3> *dh, FiniteElement<dim,3> *fe);

	/**
	 * Calculates the velocity field on a given @p dim dimensional cell.
	 * @param cell     The cell.
	 * @param velocity The computed velocity field (at quadrature points)
	 * @param fv       The FEValues class providing the quadrature points
	 *                 and the shape functions for velocity.
	 */
	template<unsigned int dim>
	void calculate_velocity(typename DOFHandler<dim,3>::CellIterator cell, std::vector<arma::vec3> &velocity, FEValuesBase<dim,3> &fv);

	/**
	 * Calculates the dispersivity (diffusivity) tensor from the velocity field.
	 * @param K        The computed dispersivity tensor.
	 * @param velocity The velocity field (at quadrature points).
	 */
	void calculate_dispersivity_tensor(std::vector<arma::mat33> &K, std::vector<arma::vec3> &velocity);

	/**
	 * Sets up some parameters of the DG method for two sides of a neighbour.
	 * @param n					The neighbour.
	 * @param s1				Side 1.
	 * @param s2				Side 2.
	 * @param n_points			Number of quadrature points.
	 * @param K					Dispersivity tensor.
	 * @param normal_vector		Normal vector to side 0 of the neighbour
	 * 							(assumed constant along the side).
	 * @param alpha				Penalty parameter that influences the continuity
	 * 							of the solution (large value=more continuity).
	 * @param advection			Coefficient of advection/transport (0=no advection).
	 * @param gamma				Computed penalty parameters.
	 * @param omega				Computed weights.
	 * @param transport_flux	Computed flux from side 1 to side 2.
	 */

	void set_DG_parameters(const Neighbour *n,
	        const int s1,
	        const int s2,
	        const unsigned int n_points,
	        const std::vector< std::vector<arma::mat33> > &K,
	        const arma::vec3 &normal_vector,
	        const double alpha,
	        const double advection,
	        double &gamma,
	        double *omega,
	        double &transport_flux);

	/**
	 * Sets up parameters of the DG method on a given boundary edge.
	 * Assupmtion is that the edge consists of only 1 side.
	 * @param edge				The edge.
	 * @param n_points			Number of quadrature points.
	 * @param K					Dispersivity tensor.
	 * @param normal_vector		Normal vector (assumed constant along the edge).
	 * @param alpha				Penalty parameter that influences the continuity
	 * 							of the solution (large value=more continuity).
	 * @param advection			Coefficient of advection/transport (0=no advection).
	 * @param gamma				Computed penalty parameters.
	 * @param omega				Computed weights.
	 */
	void set_DG_parameters_edge(const Edge *edge,
	            const unsigned int n_points,
	            const std::vector< vector<arma::mat33> > &K,
	            const arma::vec3 &normal_vector,
	            const double alpha,
	            const double advection,
	            double &gamma,
	            double *omega);

	/**
	 * Generates the file name for the time-dependent boundary condition.
	 * @param level
	 * @return
	 */
	string make_bc_file_name(int level);

	/**
	 * Reads the boundary condition.
	 */
	void read_bc_vector();

	/**
	 * Reads the initial condition.
	 */
	void read_initial_condition();

	/**
	 * Reads the names of transported substances.
	 */
	void read_subst_names();



	// physical parameters
	/**
	 * Longitudal dispersivity.
	 */
	double alphaL;

	/**
	 * Transversal dispersivity.
	 */
	double alphaT;

	/**
	 * Molecular diffusivity.
	 */
	double molecular_diffusivity;

	/**
	 * Tortuosity.
	 */
	double tortuosity;

	/**
	 * Number of transported substances.
	 */
	int n_substances;

	/**
	 * Names of transported substances.
	 */
	std::vector<string> substance_names;

	/**
	 * Time marks for boundary conditions.
	 */
	TimeMark::Type bc_mark_type_;

	/**
	 * Counter for the time level of the boundary condition.
	 * Value -1 means b.c. independent of time.
	 */
	unsigned int bc_time_level;

	/**
	 *  Boundary conditions for substance concentrations (parallel vector).
	 */
	Vec *bcv;

	/**
	 * Boundary conditions for substance concentrations (sequential array).
	 */
	double **bc_values;

	/**
	 * Vector of fluxes across element edges.
     */
	Vec flux_vector;

	Vec rhs;
	Mat stiffness_matrix;
	Mat mass_matrix;

	/**
	 * Linear algebra system for the transport equation.
	 */
	LinSys *ls;

	LinSys *ls_dt;

	/**
	 * Solver for the linear algebraic system.
	 */
	struct Solver *solver;

	/**
	 * Array for storing the output solution data.
	 */
	vector<double*> output_solution, output_cell_solution;

	/**
	 * Class for handling the solution output.
	 */
	OutputTime *transport_output;

	/**
	 * DOF handlers for 1D-3D.
	 */
	DOFHandler<1,3> *dof_handler1d;
	DOFHandler<2,3> *dof_handler2d;
	DOFHandler<3,3> *dof_handler3d;

	/**
	 * Finite elements for 1D-3D.
	 */
	FiniteElement<1,3> *fe1d;
	FiniteElement<2,3> *fe2d;
	FiniteElement<3,3> *fe3d;


	/**
	 * Penalty parameter.
	 */
	std::vector<double> gamma;

	// coefficient affecting inter-element continuity due to dispersion
	double alpha;


    // if flux through a boundary side is < -tol_flux_bc then the Dirichlet condition is applied
    const double tol_flux_bc;

    // coefficient of advection/transport (0=no advection)
    const double advection;

    /**
     * Indicates whether the fluxes have changed in the last time step.
     */
    bool flux_changed;
};





#endif /* TRANSPORT_DG_HH_ */
