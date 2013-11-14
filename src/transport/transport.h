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
 *
 * TODO:
 * - remove transport_sources
 * - in create_transport_matric_mpi, there there is condition edge_flow > ZERO
 *   this makes matrix sparser, but can lead to elements without outflow and other problems
 *   when there are big differences in fluxes, more over it doesn't work if overall flow is very small
 *
 */

#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include <petscmat.h>
#include "coupling/equation.hh"
#include "input/accessors.hh"
#include "flow/mh_dofhandler.hh"
#include "transport/transport_operator_splitting.hh"

#include "fields/field_base.hh"
#include "fields/field_values.hh"


class OutputTime;
class Mesh;
class Distribution;
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
        static Input::Type::Selection sorption_type_selection;

        EqData();
        virtual ~EqData() {};

        /// Override generic method in order to allow specification of the boundary conditions through the old bcd files.
        RegionSet read_boundary_list_item(Input::Record rec);

        Field<3, FieldValue<3>::Scalar> por_imm;        ///< Immobile porosity
        Field<3, FieldValue<3>::Vector> alpha;          ///< Coefficients of non-equilibrium linear mobile-immobile exchange
        Field<3, FieldValue<3>::EnumVector> sorp_type;  ///< Type of sorption for each substance
        Field<3, FieldValue<3>::Vector> sorp_coef0;     ///< Coefficient of sorption for each substance
        Field<3, FieldValue<3>::Vector> sorp_coef1;     ///< Coefficient of sorption for each substance
        Field<3, FieldValue<3>::Scalar> phi;            ///< solid / solid mobile

        MultiField<3, FieldValue<3>::Scalar>    conc_mobile;    ///< Calculated concentrations in the mobile zone.
    };

    /**
     * Constructor.
     */
        ConvectionTransport(Mesh &init_mesh, const Input::Record &in_rec);
	/**
	 * TODO: destructor
	 */
	virtual ~ConvectionTransport();

	/**
	 * Calculates one time step of explicit transport.
	 */
	void compute_one_step();

	/**
	 * Use new flow field vector for construction of convection matrix.
	 * Updates CFL time step constrain.
	 * TODO: Just set the new velocity, postpone update till compute_one_step
	 */
	//void set_flow_field_vector(const MH_DofHandler &dh);
  
	/**
	 * Set cross section of fractures from the Flow equation.
	 *
	 * TODO: Make this and previous part of Transport interface in TransportBase.
	 */
	void set_cross_section_field(Field< 3, FieldValue<3>::Scalar >* cross_section);


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
	 * Communicate parallel concentration vectors into sequential output vector.
	 */
	void output_vector_gather(); //

    /**
     * @brief Write computed fields.
     */
    virtual void output_data();


	/**
	 * Getters.
	 */
	inline EqData *get_data() { return &data_; }

	double ***get_concentration_matrix();
	void get_par_info(int * &el_4_loc, Distribution * &el_ds);
	bool get_dual_porosity();
	int *get_el_4_loc();
	int *get_row_4_el();
	virtual void get_parallel_solution_vector(Vec &vc);
	virtual void get_solution_vector(double* &vector, unsigned int &size);

	TimeIntegrationScheme time_scheme() { return TimeIntegrationScheme::explicit_euler; }

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
	void compute_concentration_sources(unsigned int sbi);
	void compute_concentration_sources_for_mass_balance(unsigned int sbi);

	/**
	 * Finish explicit transport matrix (time step scaling)
	 */
	void transport_matrix_step_mpi(double time_step); //

	void transport_dual_porosity( int elm_pos, ElementFullIter elem, int sbi); //
	void transport_sorption(int elm_pos, ElementFullIter elem, int sbi); //
	void compute_sorption(double conc_avg, double sorp_coef0, double sorp_coef1, unsigned int sorp_type,
			double *concx, double *concx_sorb, double Nv, double N); //


    void alloc_transport_vectors();
    void alloc_transport_structs_mpi();

    /**
     * Implements the virtual method EquationForMassBalance::calc_fluxes().
     */
    void calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance);
    /**
     * Implements the virtual method EquationForMassBalance::calc_elem_sources().
     */
    void calc_elem_sources(vector<vector<double> > &mass, vector<vector<double> > &src_balance);


    /**
     *  Parameters of the equation, some are shared with other implementations since EqData is derived from TransportBase::TransportEqData
     */
    EqData data_;

    /**
     * Indicates if we finished the matrix and add vector by scaling with timestep factor.
     */
	bool is_convection_matrix_scaled, need_time_rescaling;

    bool              sorption;     // Include sorption  YES/NO
    bool              dual_porosity;   // Include dual porosity YES/NO
    int sub_problem;    // 0-only transport,1-transport+dual porosity,
                        // 2-transport+sorption
                        // 3-transport+dual porosity+sorption


    double *sources_corr;
    Vec v_sources_corr;
    
    ///temporary arrays to store constant values of fields over time interval
    //(avoiding calling "field.value()" to often)
    double **sources_density, 
           **sources_conc,
           **sources_sigma;

    TimeMark::Type target_mark_type;    ///< TimeMark type for time marks denoting end of every time interval where transport matrix remains constant.
    double cfl_max_step;
            // only local part



    VecScatter vconc_out_scatter;
    Mat tm; // PETSc transport matrix

    /// Time when the transport matrix was created.
    /// TODO: when we have our own classes for LA objects, we can use lazy dependence to check
    /// necessity for matrix update
    double transport_matrix_time;

    /// Concentration vectors for mobile phase.
    Vec *vconc; // concentration vector
    /// Concentrations for phase, substance, element
    double ***conc;

    ///
    Vec *vpconc; // previous concentration vector
    //double ***pconc;
    Vec *bcvcorr; // boundary condition correction vector
    Vec *vcumulative_corr;
    double **cumulative_corr;

    Vec *vconc_out; // concentration vector output (gathered)
    double ***out_conc;


            int *row_4_el;
            int *el_4_loc;
            Distribution *el_ds;

            friend class TransportOperatorSplitting;
};
#endif /* TRANSPORT_H_ */
