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
 * @file    transport.h
 * @brief   
 * @todo
 * - remove transport_sources
 * - in create_transport_matric_mpi, there there is condition edge_flow > ZERO
 *   this makes matrix sparser, but can lead to elements without outflow and other problems
 *   when there are big differences in fluxes, more over it doesn't work if overall flow is very small
 */

#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include <petscmat.h>
#include "coupling/equation.hh"
#include "input/accessors.hh"
#include "flow/mh_dofhandler.hh"
#include "transport/transport_operator_splitting.hh"

#include "fields/field_algo_base.hh"
#include "fields/bc_field.hh"
#include "fields/field_values.hh"
#include "fields/multi_field.hh"
#include "fields/vec_seq_double.hh"


class SorptionImmob;
class OutputTime;
class Mesh;
class Distribution;
class ConvectionTransport;

//=============================================================================
// TRANSPORT
//=============================================================================

/**
 * TODO:
 * - doxy documentation
 * - make separate method for changing time step (rescaling only), reassembly matrix only when data are changed
 *
 * - needs support from EqData to determine next change of the data for 1) transport matrix 2) source vector
 *   this allows us precisely choose interval where to fix timestep
 *   : field->next_change_time() - returns time jump or actual time in case of time dep. field
 */



/**
 * Class that implements explicit finite volumes scheme with upwind. The timestep is given by CFL condition.
 *
 */
class ConvectionTransport : public TransportBase {
public:

    class EqData : public TransportBase::TransportEqData {
    public:

        EqData();
        virtual ~EqData() {};

        /// Override generic method in order to allow specification of the boundary conditions through the old bcd files.
        RegionSet read_boundary_list_item(Input::Record rec);

		/**
		 * Boundary conditions (Dirichlet) for concentrations.
		 * They are applied only on water inflow part of the boundary.
		 */
		BCField<3, FieldValue<3>::Vector> bc_conc;

		/// Initial concentrations.
		Field<3, FieldValue<3>::Vector> init_conc;

		Field<3, FieldValue<3>::Integer> region_id;
        MultiField<3, FieldValue<3>::Scalar>    conc_mobile;    ///< Calculated concentrations in the mobile zone.


        /// Fields indended for output, i.e. all input fields plus those representing solution.
        FieldSet output_fields;
    };

    /**
     * Constructor.
     */
        ConvectionTransport(Mesh &init_mesh, const Input::Record in_rec);
	/**
	 * TODO: destructor
	 */
	virtual ~ConvectionTransport();

	/**
	 * Initialize solution at zero time.
	 */
    void zero_time_step() override;
    /**
     * 
     */
    bool assess_time_constraint(double &time_constraint);
	/**
	 * Calculates one time step of explicit transport.
	 */
	void update_solution() override;

    /**
     * Set time interval which is considered as one time step by TransportOperatorSplitting.
     * In particular the velocity field dosn't change over this interval.
     *
     * Dependencies:
     *
     * velocity, porosity -> matrix, source_vector
     * matrix -> time_step
     *
     * data_read_times -> time_step (not necessary if we won't stick to jump times)
     * data -> source_vector
     * time_step -> scaling
     *
     *
     *
     */
    void set_target_time(double target_time);

    /**
     * Use Balance object from upstream equation (e.g. in various couplings) instead of own instance.
     */
    void set_balance_object(boost::shared_ptr<Balance> balance);

    const vector<unsigned int> &get_subst_idx()
	{ return subst_idx; }

    /**
     * Calculate quantities necessary for cumulative balance (over time).
     * This method is called at each (sub)iteration of the time loop.
     */
    void calculate_cumulative_balance();

    /**
     * Calculate instant quantities at output times.
     */
    void calculate_instant_balance();

	/**
	 * Communicate parallel concentration vectors into sequential output vector.
	 */
	void output_vector_gather();

    /**
     * @brief Write computed fields.
     */
    virtual void output_data() override;


	/**
	 * Getters.
	 */
	inline EqData *get_data() { return &data_; }

	inline std::shared_ptr<OutputTime> output_stream() { return output_stream_; }

	double **get_concentration_matrix();
	Vec *get_concentration_vector() { return vconc; }
	void get_par_info(int * &el_4_loc, Distribution * &el_ds);
	int *get_el_4_loc();
	int *get_row_4_el();

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
     */
    void create_transport_matrix_mpi();

    void make_transport_partitioning(); //
	void set_initial_condition();
	void read_concentration_sources();
	void set_boundary_conditions();
  
  //note: the source of concentration is multiplied by time interval (gives the mass, not the flow like before)
// 	void compute_concentration_sources(unsigned int sbi);

    //note: the source of concentration is multiplied by time interval (gives the mass, not the flow like before)
    void compute_concentration_sources();
    
	/**
	 * Finish explicit transport matrix (time step scaling)
	 */
	void transport_matrix_step_mpi(double time_step); //

    void alloc_transport_vectors();
    void alloc_transport_structs_mpi();



    /**
     *  Parameters of the equation, some are shared with other implementations since EqData is derived from TransportBase::TransportEqData
     */
    EqData data_;

    //@{
    /**
     * Flag indicates the state of object (transport matrix or source term).
     * If false, the object is freshly assembled and not rescaled.
     * If true, the object is scaled (not necessarily with the current time step).
     */
	bool is_convection_matrix_scaled, is_src_term_scaled;
    //@}
    
    double **sources_corr;
    Vec *v_sources_corr;
    

    TimeMark::Type target_mark_type;    ///< TimeMark type for time marks denoting end of every time interval where transport matrix remains constant.
    double cfl_max_step;    ///< Time step constraint coming from CFL condition.
    
    Vec vcfl_flow_,     ///< Parallel vector for flow contribution to CFL condition.
        vcfl_source_;   ///< Parallel vector for source term contribution to CFL condition.
    double *cfl_flow_, *cfl_source_;


    VecScatter vconc_out_scatter;
    Mat tm; // PETSc transport matrix
    Vec *v_tm_diag; // additions to PETSC transport matrix on the diagonal - from sources (for each substance)
    double **tm_diag;

    /// Time when the transport matrix was created.
    /// TODO: when we have our own classes for LA objects, we can use lazy dependence to check
    /// necessity for matrix update
    double transport_matrix_time;
    double transport_bc_time;   ///< Time of the last update of the boundary condition terms.

    /// Concentration vectors for mobile phase.
    Vec *vconc; // concentration vector
    /// Concentrations for phase, substance, element
    double **conc;

    ///
    Vec *vpconc; // previous concentration vector
    Vec *bcvcorr; // boundary condition correction vector
    Vec *vcumulative_corr;
    double **cumulative_corr;

    std::vector<VectorSeqDouble> out_conc;

	/// Record with output specification.
	Input::Record output_rec;

	std::shared_ptr<OutputTime> output_stream_;


	int *row_4_el;
	int *el_4_loc;
	Distribution *el_ds;

	/// List of indices used to call balance methods for a set of quantities.
	vector<unsigned int> subst_idx;

            friend class TransportOperatorSplitting;
};
#endif /* TRANSPORT_H_ */
