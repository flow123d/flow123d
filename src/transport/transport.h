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
 * - make separate method for changing time step (rescaling only)
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
	void set_flow_field_vector(const MH_DofHandler &dh);
  
	/**
	 * Set cross section of fractures from the Flow equation.
	 *
	 * TODO: Make this and previous part of Transport interface in TransportBase.
	 */
	void set_cross_section_field(Field< 3, FieldValue<3>::Scalar >* cross_section);


    /**
     * Set time interval over which we should use fixed transport matrix. Rescale transport matrix.
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
	//int get_n_substances();
	int *get_el_4_loc();
	int *get_row_4_el();
	virtual void get_parallel_solution_vector(Vec &vc);
	virtual void get_solution_vector(double* &vector, unsigned int &size);
	/**
	 * Return pointer to sequential arrays for output.
	 * TODO: Maybe this should be made by get_solution_vector, but here we have matrix of arrays.
	 */
	double ***get_out_conc();
	double ***get_conc();
    //vector<string> &get_substance_names();
    double *get_sources(int sbi);





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
	Vec compute_concentration_sources(unsigned int subst_i, double *conc);

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
     * Overriding the virtual method that is called by TransportBase::mass_balance() to get boundary balances over individual boundary regions.
     * TODO: more precise description
     */
    void calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance);
    /**
     * Overriding the virtual method that is called by TransportBase::mass_balance() to get source balances over individual boundary regions.
     * TODO: more precise description
     */
    void calc_elem_sources(vector<vector<double> > &mass, vector<vector<double> > &src_balance);


    /**
     *  Parameters of the equation, some are shared with other implementations since EqData is derived from TransportBase::TransportEqData
     */
    EqData data_;
    /**
     * Class for handling the solution output.
     */
    OutputTime *field_output;
    /**
     * Indicates if we finished the matrix and add vector by scaling with timestep factor.
     */
	bool is_convection_matrix_scaled, is_bc_vector_scaled;

    bool              sorption;     // Include sorption  YES/NO
    bool              dual_porosity;   // Include dual porosity YES/NO
    int sub_problem;    // 0-only transport,1-transport+dual porosity,
                        // 2-transport+sorption
                        // 3-transport+dual porosity+sorption


    double *sources_corr;
    Vec v_sources_corr;

    TimeMark::Type target_mark_type;    ///< TimeMark type for time marks denoting end of every time interval where transport matrix remains constant.
    double cfl_max_step;
            // only local part



            VecScatter vconc_out_scatter;
    Mat tm; // PETSc transport matrix

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
