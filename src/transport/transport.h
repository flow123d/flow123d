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
#include "fields/bc_multi_field.hh"
#include "fields/field_values.hh"
#include "fields/multi_field.hh"
#include "fields/vec_seq_double.hh"
#include "io/equation_output.hh"

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
class ConvectionTransport : public ConcentrationTransportBase {
public:

    class EqData : public TransportEqData {
    public:

        EqData();
        virtual ~EqData() {};

        /// Override generic method in order to allow specification of the boundary conditions through the old bcd files.
        RegionSet read_boundary_list_item(Input::Record rec);

		/**
		 * Boundary conditions (Dirichlet) for concentrations.
		 * They are applied only on water inflow part of the boundary.
		 */
		BCMultiField<3, FieldValue<3>::Scalar> bc_conc;

		/// Initial concentrations.
		MultiField<3, FieldValue<3>::Scalar> init_conc;

		Field<3, FieldValue<3>::Integer> region_id;
        MultiField<3, FieldValue<3>::Scalar>    conc_mobile;    ///< Calculated concentrations in the mobile zone.


        /// Fields indended for output, i.e. all input fields plus those representing solution.
        EquationOutput output_fields;
    };


    typedef ConcentrationTransportBase FactoryBaseType;

    static const Input::Type::Record & get_input_type();

    static const IT::Selection & get_output_selection();

    /**
     * Constructor.
     */
    ConvectionTransport(Mesh &init_mesh, const Input::Record in_rec);
	/**
	 * TODO: destructor
	 */
	virtual ~ConvectionTransport();

	void initialize() override;

	/**
	 * Initialize solution at zero time.
	 */
    void zero_time_step() override;
    /**
      * Evaluates CFL condition.
      * Assembles the transport matrix and vector (including sources, bc terms).
      * @param time_constraint is the value CFL constraint (return parameter)
      * @return true if CFL is changed since previous step, false otherwise
      */
    bool evaluate_time_constraint(double &time_constraint) override;
	/**
	 * Calculates one time step of explicit transport.
	 */
	void update_solution() override;

	void calculate_concentration_matrix() override {};

	void update_after_reactions(bool solution_changed) override {};

    /**
     * Set time interval which is considered as one time step by TransportOperatorSplitting.
     * In particular the velocity field doesn't change over this interval.
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
    void set_target_time(double target_time) override;

    /**
     * Use Balance object from upstream equation (e.g. in various couplings) instead of own instance.
     */
    void set_balance_object(std::shared_ptr<Balance> balance) override;

    const vector<unsigned int> &get_subst_idx() override
	{ return subst_idx; }

    /**
     * Calculate quantities necessary for cumulative balance (over time).
     * This method is called at each (sub)iteration of the time loop.
     */
    void calculate_cumulative_balance() override;

    /**
     * Calculate instant quantities at output times.
     */
    void calculate_instant_balance() override;


    /**
     * @brief Write computed fields.
     */
    virtual void output_data() override;

    inline void set_velocity_field(const MH_DofHandler &dh) override
    { mh_dh=&dh; }

    void set_output_stream(std::shared_ptr<OutputTime> stream) override
    { output_stream_ = stream; }


	/**
	 * Getters.
	 */
	inline std::shared_ptr<OutputTime> output_stream() override
	{ return output_stream_; }

	double **get_concentration_matrix() override;

	const Vec &get_solution(unsigned int sbi) override
	{ return vconc[sbi]; }

	void get_par_info(int * &el_4_loc, Distribution * &el_ds) override;

	int *get_row_4_el() override;

    /// Returns number of transported substances.
    inline unsigned int n_substances() override
    { return substances_.size(); }

    /// Returns reference to the vector of substance names.
    inline SubstanceList &substances() override
    { return substances_; }

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
    void create_mass_matrix();

    void make_transport_partitioning(); //
	void set_initial_condition();
	void read_concentration_sources();
	void set_boundary_conditions();
  
    /** @brief Assembles concentration sources for each substance.
     * note: the source of concentration is multiplied by time interval (gives the mass, not the flow like before)
     */
    void compute_concentration_sources();
    
	/**
	 * Finish explicit transport matrix (time step scaling)
	 */
	void transport_matrix_step_mpi(double time_step); //

    void alloc_transport_vectors();
    void alloc_transport_structs_mpi();

	/**
	 * Communicate parallel concentration vectors into sequential output vector.
	 */
	void output_vector_gather();



    /// Registrar of class to factory
    static const int registrar;

    /**
     *  Parameters of the equation, some are shared with other implementations since EqData is derived from TransportBase::TransportEqData
     */
    EqData data_;

    //@{
    /**
     * Flag indicates the state of object (transport matrix or source or boundary term).
     * If false, the object is freshly assembled and not rescaled.
     * If true, the object is scaled (not necessarily with the current time step).
     */
	bool is_convection_matrix_scaled, is_src_term_scaled, is_bc_term_scaled;
	
	/// Flag indicates that porosity or cross_section changed during last time.
	bool is_mass_diag_changed;
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
    Vec mass_diag;  // diagonal entries in pass matrix (cross_section * porosity)
    Vec vpmass_diag;  // diagonal entries in mass matrix from last time (cross_section * porosity)
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

	/// Record with input specification.
	const Input::Record input_rec;

	std::shared_ptr<OutputTime> output_stream_;


	int *row_4_el;
	int *el_4_loc;
	Distribution *el_ds;

    /// Transported substances.
    SubstanceList substances_;

    /**
     * Temporary solution how to pass velocity field form the flow model.
     * TODO: introduce FieldDiscrete -containing true DOFHandler and data vector and pass such object together with other
     * data. Possibly make more general set_data method, allowing setting data given by name. needs support from EqDataBase.
     */
    const MH_DofHandler *mh_dh;

	/// List of indices used to call balance methods for a set of quantities.
	vector<unsigned int> subst_idx;

    friend class TransportOperatorSplitting;
};
#endif /* TRANSPORT_H_ */
