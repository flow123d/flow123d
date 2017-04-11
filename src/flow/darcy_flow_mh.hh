/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    darcy_flow_mh.hh
 * @brief   mixed-hybrid model of linear Darcy flow, possibly unsteady.
 * @author  Jan Brezina
 *
 * Main object for mixed-hybrid discretization of the linear elliptic PDE (Laplace)
 * on  a multidimensional domain. Discretization of saturated Darcy flow using
 * RT0 approximation for the velocity
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

#include <memory>
#include "input/input_type_forward.hh"

#include <petscmat.h>
#include "system/sys_vector.hh"
#include "coupling/equation.hh"
#include "flow/mh_dofhandler.hh"

#include "fields/bc_field.hh"
#include "fields/field.hh"
#include "fields/field_set.hh"
#include "flow/darcy_flow_interface.hh"
#include "input/input_exception.hh"


/// external types:
class LinSys;
class LinSys_BDDC;
class Mesh;
class DarcyFlowMHOutput;
class Balance;
class AssemblyBase;

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
 * - @f$ h = \frac{\pi}{\rho_0 g}+z @f$ is pressure head, @f$ \pi, \rho_0, g @f$ are the pressure, water density, and acceleration of gravity , respectively.
 *   Assumes gravity force acts counter to the direction of the @f$ z @f$ axis.
 * - @f$ R @f$ is destity or gravity variability coefficient. For density driven flow it should be
 * @f[
 *    R = \frac{\rho}{\rho_0} -1 = \rho_0^{-1}\sum_{i=1}^s c_i
 * @f]
 *   where @f$ c_i @f$ is concentration in @f$ kg m^{-3} @f$.
 *
 * The time key is optional, when not specified the equation is forced to steady regime. Using Steady TimeGovernor which have no dt constraints.
 *
 *
 * TODO:
 * Make solution regular field (need FeSeystem and parallel DofHandler for edge pressures), then remove get_solution_vector from
 * Equation interface.
 */
/**
 * Model for transition coefficients due to Martin, Jaffre, Roberts (see manual for full reference)
 *
 * TODO:
 * - how we can reuse field values computed during assembly
 *
 */

class DarcyMH : public DarcyFlowInterface
{
public:
    TYPEDEF_ERR_INFO( EI_Reason, string);
    DECLARE_EXCEPTION(ExcSolverDiverge,
            << "Diverged nonlinear solver. Reason: " << EI_Reason::val
             );
    DECLARE_INPUT_EXCEPTION(ExcMissingTimeGovernor,
            << "Missing the key 'time', obligatory for the transient problems.");


    typedef std::vector<std::shared_ptr<AssemblyBase> > MultidimAssembly;

    /// Class with all fields used in the equation DarcyFlow.
    /// This is common to all implementations since this provides interface
    /// to this equation for possible coupling.
    class EqData : public FieldSet {
    public:

        /**
         * For compatibility with old BCD file we have to assign integer codes starting from 1.
         */
        enum BC_Type {
            none=0,
            dirichlet=1,
            total_flux=4,
            seepage=5,
            river=6
        };

        /// Return a Selection corresponding to enum BC_Type.
        static const Input::Type::Selection & get_bc_type_selection();

        /// Creation of all fields.
        EqData();


        Field<3, FieldValue<3>::TensorFixed > anisotropy;
        Field<3, FieldValue<3>::Scalar > conductivity;
        Field<3, FieldValue<3>::Scalar > cross_section;
        Field<3, FieldValue<3>::Scalar > water_source_density;
        Field<3, FieldValue<3>::Scalar > sigma;
        
        BCField<3, FieldValue<3>::Enum > bc_type; // Discrete need Selection for initialization
        BCField<3, FieldValue<3>::Scalar > bc_pressure; 
        BCField<3, FieldValue<3>::Scalar > bc_flux;
        BCField<3, FieldValue<3>::Scalar > bc_robin_sigma;
        BCField<3, FieldValue<3>::Scalar > bc_switch_pressure;
        
        Field<3, FieldValue<3>::Scalar > init_pressure;
        Field<3, FieldValue<3>::Scalar > storativity;

        /**
         * Gravity vector and constant shift of pressure potential. Used to convert piezometric head
         * to pressure head and vice versa.
         */
        arma::vec4 gravity_;
        arma::vec3 gravity_vec_;

        // Mirroring the following members of DarcyMH:
        Mesh *mesh;
        MH_DofHandler *mh_dh;


        uint water_balance_idx;
        unsigned int local_boundary_index;
        std::shared_ptr<Balance> balance;
        LinSys *lin_sys;
        
        unsigned int n_schur_compls;
        int is_linear;              ///< Hack fo BDDC solver.
        bool force_bc_switch;       ///< auxiliary flag for switchting Dirichlet like BC
        
        /// Idicator of dirichlet or neumann type of switch boundary conditions.
        std::vector<char> bc_switch_dirichlet;
    };

    /// Type of experimental Mortar-like method for non-compatible 1d-2d interaction.
    enum MortarMethod {
        NoMortar = 0,
        MortarP0 = 1,
        MortarP1 = 2
    };
    /// Selection for enum MortarMethod.
    static const Input::Type::Selection & get_mh_mortar_selection();





    DarcyMH(Mesh &mesh, const Input::Record in_rec);

    static const Input::Type::Record & type_field_descriptor();
    static const Input::Type::Record & get_input_type();

    const MH_DofHandler &get_mh_dofhandler()  override {
        double *array;
        unsigned int size;
        get_solution_vector(array, size);

        // here assume that velocity field is extended as constant
        // to the previous time, so here we set left bound of the interval where the velocity
        // has current value; this may not be good for every transport !!
        // we can resolve this when we use FieldFE to store computed velocities in few last steps and
        // let every equation set time according to nature of the time scheme

        // in particular this setting is necessary to prevent ConvectinTransport to recreate the transport matrix
        // every timestep ( this may happen for unsteady flow if we would use time->t() here since it returns infinity.
        mh_dh.set_solution(time_->last_t(), array, solution_precision());
       return mh_dh;
    }

    void init_eq_data();
    void initialize() override;
    virtual void initialize_specific();
    void zero_time_step() override;
    void update_solution() override;

    void get_solution_vector(double * &vec, unsigned int &vec_size) override;
    void get_parallel_solution_vector(Vec &vector) override;
    
    /// postprocess velocity field (add sources)
    virtual void prepare_new_time_step();
    virtual void postprocess();
    virtual void output_data() override;

    virtual ~DarcyMH() override;


protected:
    /**
     * Returns true is the fields involved in the time term have values that makes the time term zero.
     * For time_global==true, it returns true if there are no field descriptors in the input list, so the
     * fields )of the time ter) have their default values for whole simulation.
     * If time_global==false (default), only the actual values are considered.
     */
    virtual bool zero_time_term(bool time_global=false);

    /// Solve method common to zero_time_step and update solution.
    void solve_nonlinear();
    void make_serial_scatter();
    void modify_system();
    virtual void setup_time_term();


    //void prepare_parallel();
    //void make_row_numberings();
    /// Initialize global_row_4_sub_row.
    //void prepare_parallel_bddc();

    /**
     * Create and preallocate MH linear system (including matrix, rhs and solution vectors)
     */
    void create_linear_system(Input::AbstractRecord rec);

    /**
     * Read initial condition into solution vector.
     * Must be called after create_linear_system.
     *
     * For the LMH scheme we have to be able to save edge pressures in order to
     * restart simulation or use results of one simulation as initial condition for other one.
     */
    virtual void read_initial_condition();

    /**
     * Part of per element assembly that is specific for MH and LMH respectively.
     *
     * This implemnets MH case:
     * - compute conductivity scaling
     * - assembly source term
     * - no time term, managed by diagonal extraction etc.
     */
    //virtual void local_assembly_specific(AssemblyData &local_data);
   
    /**
     * Allocates linear system matrix for MH.
     * TODO:
     * - use general preallocation methods in DofHandler
     */
    void allocate_mh_matrix();

    /**
     * Assembles linear system matrix for MH.
     * Element by element assembly is done using dim-template assembly class.
     * Assembles only steady part of the equation.
     * TODO:
     * - include time term
     * - add support for Robin type sources
     * - support for nonlinear solvers - assembly either residual vector, matrix, or both (using FADBAD++)
     */
    void assembly_mh_matrix(MultidimAssembly& assembler);

    /// Source term is implemented differently in LMH version.
    virtual void assembly_source_term();

    /**
     * Assembly or update whole linear system.
     */
    virtual void assembly_linear_system();

    void set_mesh_data_for_bddc(LinSys_BDDC * bddc_ls);
    /**
     * Return a norm of residual vector.
     * TODO: Introduce Equation::compute_residual() updating
     * residual field, standard part of EqData.
     */
    virtual double solution_precision() const;
    
    /// Print darcy flow matrix in matlab format into a file.
    void print_matlab_matrix(string matlab_file);

    bool solution_changed_for_scatter;
    //Vec velocity_vector;
    MH_DofHandler mh_dh;    // provides access to seq. solution fluxes and pressures on sides

    MortarMethod mortar_method_;

    std::shared_ptr<Balance> balance_;
    /// index of water balance within the Balance object.
    unsigned int water_balance_idx_;

    DarcyFlowMHOutput *output_object;

	int size;				    // global size of MH matrix
	int  n_schur_compls;  	    // number of shur complements to make
	double  *solution; 			// sequantial scattered solution vector

	// Propagate test for the time term to the assembly.
	// This flag is necessary for switching BC to avoid setting zero neumann on the whole boundary in the steady case.
	bool use_steady_assembly_;
	bool data_changed_;

	// Setting of the nonlinear solver. TODO: Move to the solver class later on.
	double tolerance_;
	unsigned int min_n_it_;
	unsigned int max_n_it_;
	unsigned int nonlinear_iteration_; //< Actual number of completed nonlinear iterations, need to pass this information into assembly.


	LinSys *schur0;  		//< whole MH Linear System


	// gather of the solution
	Vec sol_vec;			                 //< vector over solution array
	VecScatter par_to_all;

	Vec steady_diagonal;
    Vec steady_rhs;
    Vec new_diagonal;
    Vec previous_solution;

	std::shared_ptr<EqData> data_;

    friend class DarcyFlowMHOutput;
    friend class P0_CouplingAssembler;
    friend class P1_CouplingAssembler;

private:
  /// Registrar of class to factory
  static const int registrar;
};


class P0_CouplingAssembler {
public:
	P0_CouplingAssembler(const DarcyMH &darcy)
	: darcy_(darcy),
	  master_list_(darcy.mesh_->master_elements),
	  intersections_(darcy.mesh_->intersections),
	  master_(nullptr),
	  tensor_average(2)
	{
		arma::mat master_map(1,2), slave_map(1,3);
		master_map.fill(1.0 / 2);
		slave_map.fill(1.0 / 3);

		tensor_average[0].push_back( master_map.t() * master_map );
		tensor_average[0].push_back( master_map.t() * slave_map );
		tensor_average[1].push_back( slave_map.t() * master_map );
		tensor_average[1].push_back( slave_map.t() * slave_map );
	}

	void assembly(LinSys &ls);
	void pressure_diff(int i_ele,
			vector<int> &dofs,
			unsigned int &ele_type,
			double &delta,
			arma::vec &dirichlet);
private:
	typedef vector<unsigned int> IsecList;

	const DarcyMH &darcy_;

	const vector<IsecList> &master_list_;
	const vector<Intersection> &intersections_;

	vector<IsecList>::const_iterator ml_it_;
	const Element *master_;

	/// Row matrices to compute element pressure as average of boundary pressures
	vector< vector< arma::mat > > tensor_average;
	/// measure of master element, should be sum of intersection measures
	double delta_0;
};



class P1_CouplingAssembler {
public:
	P1_CouplingAssembler(const DarcyMH &darcy)
	: darcy_(darcy),
	  intersections_(darcy.mesh_->intersections),
	  rhs(5),
	  dofs(5),
	  dirichlet(5)
	{
		rhs.zeros();
	}

	void assembly(LinSys &ls);
	void add_sides(const Element * ele, unsigned int shift, vector<int> &dofs, vector<double> &dirichlet);
private:

	const DarcyMH &darcy_;
	const vector<Intersection> &intersections_;

	arma::vec rhs;
	vector<int> dofs;
	vector<double> dirichlet;
};



void mat_count_off_proc_values(Mat m, Vec v);



/**
 * @brief Mixed-hybrid solution of unsteady Darcy flow.
 *
 * Standard discretization with time term and sources picewise constant
 * on the element. This leads to violation of the discrete maximum principle for
 * non-acute meshes or to too small timesteps. For simplicial meshes this can be solved by lumping to the edges. See DarcyFlowLMH_Unsteady.
 */
/*
class DarcyFlowMH_Unsteady : public DarcyFlowMH_Steady
{
public:

    DarcyFlowMH_Unsteady(Mesh &mesh, const Input::Record in_rec);
    //DarcyFlowMH_Unsteady();

    static const Input::Type::Record & get_input_type();
protected:
    void read_initial_condition() override;
    void modify_system() override;
    void setup_time_term();
    
private:
    /// Registrar of class to factory
    static const int registrar;


};
*/

#endif  //DARCY_FLOW_MH_HH
//-----------------------------------------------------------------------------
// vim: set cindent:

