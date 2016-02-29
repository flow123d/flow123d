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
#include "la/linsys_BDDC.hh"
#include "la/linsys_PETSC.hh"

#include "fields/bc_field.hh"
#include "fields/field.hh"
#include "fields/field_set.hh"
#include "flow/darcy_flow_interface.hh"


/// external types:
class LinSys;
class Mesh;
class SchurComplement;
class Distribution;
class SparseGraph;
class LocalToGlobalMap;
class DarcyFlowMHOutput;
class Balance;
class VectorSeqDouble;

template<unsigned int dim, unsigned int spacedim> class FE_RT0;
template<unsigned int degree, unsigned int dim, unsigned int spacedim> class FE_P_disc;
template<unsigned int dim, unsigned int spacedim> class MappingP1;
template<unsigned int dim, unsigned int spacedim> class FEValues;
template<unsigned int dim, unsigned int spacedim> class FESideValues;
template<unsigned int dim> class QGauss;

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

class DarcyFlowMH_Steady : public DarcyFlowInterface
{
public:

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
            //neumann=2,
            //robin=3,
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

        //FieldSet  time_term_fields;
        //FieldSet  main_matrix_fields;
        //FieldSet  rhs_fields;
    };

    /// Type of experimental Mortar-like method for non-compatible 1d-2d interaction.
    enum MortarMethod {
        NoMortar = 0,
        MortarP0 = 1,
        MortarP1 = 2
    };
    /// Selection for enum MortarMethod.
    static const Input::Type::Selection & get_mh_mortar_selection();








    DarcyFlowMH_Steady(Mesh &mesh, const Input::Record in_rec);

    static const Input::Type::Record & get_input_type();

    const MH_DofHandler &get_mh_dofhandler() {
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

    void update_solution() override;
    void zero_time_step() override;

    void get_solution_vector(double * &vec, unsigned int &vec_size) override;
    void get_parallel_solution_vector(Vec &vector) override;
    
    /// postprocess velocity field (add sources)
    virtual void postprocess();
    virtual void output_data() override;

    ~DarcyFlowMH_Steady();


protected:
    //class AssemblyBase;
    //template<unsigned int dim> class Assembly;
    
    struct AssemblyData
    {
        Mesh *mesh;
        EqData* data;
        MH_DofHandler *mh_dh;
    };
    
    class AssemblyBase
    {
    public:
        // assembly just A block of local matrix
        virtual void assembly_local_matrix(arma::mat &local_matrix, 
                                           ElementFullIter ele) = 0;

        // assembly compatible neighbourings
        virtual void assembly_local_vb(double *local_vb, 
                                       ElementFullIter ele,
                                       Neighbour *ngh) = 0;

        // compute velocity value in the barycenter
        // TOTO: implement and use general interpolations between discrete spaces
        virtual arma::vec3 make_element_vector(ElementFullIter ele) = 0;
    };
    
    template<unsigned int dim>
    class Assembly : public AssemblyBase
    {
    public:
        Assembly<dim>(AssemblyData ad);
        ~Assembly<dim>();
        void assembly_local_matrix(arma::mat &local_matrix, 
                                   ElementFullIter ele) override;
        void assembly_local_vb(double *local_vb, 
                               ElementFullIter ele,
                               Neighbour *ngh
                              ) override;
        arma::vec3 make_element_vector(ElementFullIter ele) override;

        // assembly volume integrals
        FE_RT0<dim,3> fe_rt_;
        MappingP1<dim,3> map_;
        QGauss<dim> quad_;
        FEValues<dim,3> fe_values_;
        
        // assembly face integrals (BC)
        QGauss<dim-1> side_quad_;
        FE_P_disc<0,dim,3> fe_p_disc_;
        FESideValues<dim,3> fe_side_values_;

        // Interpolation of velocity into barycenters
        QGauss<dim> velocity_interpolation_quad_;
        FEValues<dim,3> velocity_interpolation_fv_;

        // data shared by assemblers of different dimension
        AssemblyData ad_;

    };
    
    /*
    void setup_velocity_vector() {
        double *velocity_array;
        unsigned int size;

        get_solution_vector(velocity_array, size);
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, mesh_->n_sides(), velocity_array, &velocity_vector);

    }*/


    /// Solve method common to zero_time_step and update solution.
    void solve_nonlinear();
    void make_serial_scatter();
    virtual void modify_system()
    {  };
    virtual void setup_time_term()
    {  };


    void prepare_parallel();
    void make_row_numberings();
    /// Initialize global_row_4_sub_row.
    void prepare_parallel_bddc();

    /**
     * Create and preallocate MH linear system (including matrix, rhs and solution vectors)
     */
    void create_linear_system();

    /**
     * Read initial condition into solution vector.
     * Must be called after create_linear_system.
     *
     */
    virtual void read_initial_condition()
    {  };

    /**
     * Abstract assembly method used for both assembly and preallocation.
     * Assembly only steady part of the equation.
     * TODO:
     * - use general preallocation methods in DofHandler
     * - include alos time term
     * - add support for Robin type sources
     * - support for nonlinear solvers - assembly either residual vector, matrix, or both (using FADBAD++)
     *
     */
    void assembly_steady_mh_matrix();
    

    /// Source term is implemented differently in LMH version.
    virtual void assembly_source_term();

    /**
     * Assembly or update whole linear system.
     */
    void assembly_linear_system();

    void set_mesh_data_for_bddc(LinSys_BDDC * bddc_ls);
    /**
     * Return a norm of residual vector.
     * TODO: Introduce Equation::compute_residual() updating
     * residual field, standard part of EqData.
     */
    virtual double solution_precision() const;

    bool solution_changed_for_scatter;
    //Vec velocity_vector;
    MH_DofHandler mh_dh;    // provides access to seq. solution fluxes and pressures on sides

    MortarMethod mortar_method_;

    /// object for calculation and writing the water balance to file.
    boost::shared_ptr<Balance> balance_;
    /// index of water balance within the Balance object.
    unsigned int water_balance_idx_;

    DarcyFlowMHOutput *output_object;

	int size;				    // global size of MH matrix
	int  n_schur_compls;  	    // number of shur complements to make
	double  *solution; 			// sequantial scattered solution vector
	int is_linear_;             // Hack fo BDDC solver.
	bool use_steady_assembly_;   // Propagate test for the time term to the assembly.
	double tolerance_;
	unsigned int max_n_it_;

	unsigned int nonlinear_iteration_; //< Actual number of completed nonlinear iterations, need to pass this information into assembly.



	LinSys *schur0;  		//< whole MH Linear System

	//AssemblyData *assembly_data_;
	std::vector<AssemblyBase *> assembly_;
	
	// parallel
	Distribution *edge_ds;          //< optimal distribution of edges
	Distribution *el_ds;            //< optimal distribution of elements
	Distribution *side_ds;          //< optimal distribution of elements
	boost::shared_ptr<Distribution> rows_ds;          //< final distribution of rows of MH matrix

	int *el_4_loc;		        //< array of idexes of local elements (in ordering matching the optimal global)
	int *row_4_el;		        //< element index to matrix row
	int *side_id_4_loc;		//< array of ids of local sides
	int	*side_row_4_id;		//< side id to matrix row
	int *edge_4_loc;		//< array of indexes of local edges
	int	*row_4_edge;		//< edge index to matrix row

	/// Idicator of dirichlet or neumann type of switch boundary conditions.
	std::vector<char> bc_switch_dirichlet;

	/// Necessary only for BDDC solver.
    boost::shared_ptr<LocalToGlobalMap> global_row_4_sub_row;           //< global dof index for subdomain index

	// gather of the solution
	Vec sol_vec;			                 //< vector over solution array
	VecScatter par_to_all;
        
	EqData data_;

    friend class DarcyFlowMHOutput;
    friend class P0_CouplingAssembler;
    friend class P1_CouplingAssembler;

private:
  /// Registrar of class to factory
  static const int registrar;
};


class P0_CouplingAssembler {
public:
	P0_CouplingAssembler(const DarcyFlowMH_Steady &darcy)
	: darcy_(darcy),
	  master_list_(darcy.mesh_->master_elements),
	  intersections_(darcy.mesh_->intersections),
	  master_(nullptr),
	  tensor_average(2)
	{
		arma::mat master_map(1,2), slave_map(1,3);
		master_map.fill(1.0 / 2);
		slave_map.fill(1.0 / 3);

		tensor_average[0].push_back( trans( master_map ) * master_map );
		tensor_average[0].push_back( trans( master_map ) * slave_map );
		tensor_average[1].push_back( trans( slave_map ) * master_map );
		tensor_average[1].push_back( trans( slave_map ) * slave_map );
	}

	void assembly(LinSys &ls);
	void pressure_diff(int i_ele,
			vector<int> &dofs,
			unsigned int &ele_type,
			double &delta,
			arma::vec &dirichlet);
private:
	typedef vector<unsigned int> IsecList;

	const DarcyFlowMH_Steady &darcy_;

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
	P1_CouplingAssembler(const DarcyFlowMH_Steady &darcy)
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

	const DarcyFlowMH_Steady &darcy_;
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

    Vec steady_diagonal;
    Vec steady_rhs;
    Vec new_diagonal;
    Vec previous_solution;

};


#endif  //DARCY_FLOW_MH_HH
//-----------------------------------------------------------------------------
// vim: set cindent:

