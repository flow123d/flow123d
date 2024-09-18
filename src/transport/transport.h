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


#include <memory>                                     // for shared_ptr
#include <vector>                                     // for vector
#include <petscmat.h>
#include "fem/fe_values.hh"                           // for FEValues
#include "fields/field.hh"                            // for Field
#include "fields/bc_multi_field.hh"
#include "fields/field_values.hh"
#include "fields/multi_field.hh"
#include "la/vector_mpi.hh"
#include "fields/equation_output.hh"
#include "input/type_base.hh"                         // for Array
#include "input/type_generic.hh"                      // for Instance
#include "input/accessors.hh"
#include "system/index_types.hh"
#include "mesh/region.hh"                             // for RegionSet
#include "petscvec.h"                                 // for Vec, _p_Vec
#include "tools/time_marks.hh"                        // for TimeMark, TimeM...
#include "transport/substance.hh"                     // for SubstanceList
#include "transport/transport_operator_splitting.hh"
#include "quadrature/quadrature_lib.hh"
#include "la/vector_mpi.hh"

class OutputTime;
class Mesh;
class Distribution;
class Balance;
namespace Input {
	namespace Type {
		class Record;
		class Selection;
	}
}
template<unsigned int dim> class MassAssemblyConvection;
template<unsigned int dim> class InitCondAssemblyConvection;
template<unsigned int dim> class ConcSourcesBdrAssemblyConvection;
template<unsigned int dim> class MatrixMpiAssemblyConvection;
template< template<IntDim...> class DimAssembly> class GenericAssembly;


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
    class EqFields : public TransportEqFields {
    public:

    	EqFields();
        virtual ~EqFields() {};

        /// Calculate flux on given side point.
        inline double side_flux(SidePoint &side_p, FEValues<3> &fe_side_values) {
            return arma::dot(this->flow_flux(side_p), fe_side_values.normal_vector(0)) * fe_side_values.JxW(0);
        }

		/**
		 * Boundary conditions (Dirichlet) for concentrations.
		 * They are applied only on water inflow part of the boundary.
		 */
		BCMultiField<3, FieldValue<3>::Scalar> bc_conc;

		/// Initial concentrations.
		MultiField<3, FieldValue<3>::Scalar> init_conc;

		Field<3, FieldValue<3>::Scalar> region_id;
        Field<3, FieldValue<3>::Scalar> subdomain;

        MultiField<3, FieldValue<3>::Scalar>    conc_mobile;    ///< Calculated concentrations in the mobile zone.
        FieldFEScalarVec conc_mobile_fe;                        ///< Underlaying FieldFE for each substance of conc_mobile.

        /// Fields indended for output, i.e. all input fields plus those representing solution.
        EquationOutput output_fields;
    };


    class EqData : public EqDataBase {
    public:

        EqData() : EqDataBase(0), is_mass_diag_changed(false), cfl_source_(PETSC_COMM_WORLD), cfl_flow_(PETSC_COMM_WORLD) {}
        virtual ~EqData() {};

        /// Returns number of transported substances.
        inline unsigned int n_substances()
        { return substances_.size(); }

		inline void set_time_governor(TimeGovernor *time) {
		    ASSERT_PTR(time);
		    this->time_ = time;
		}

		void alloc_transport_structs_mpi(unsigned int lsize) {
		    this->cfl_flow_.resize(lsize);
		    this->cfl_source_.resize(lsize);
		}

        /**
         * Temporary solution how to pass velocity field form the flow model.
         * TODO: introduce FieldDiscrete -containing true DOFHandler and data vector and pass such object together with other
         * data. Possibly make more general set_data method, allowing setting data given by name. needs support from EqDataBase.
         */
        std::shared_ptr<DOFHandlerMultiDim> dh_;

        /// object for calculation and writing the mass balance to file, shared with EquationBase.
        std::shared_ptr<Balance> balance_;

        /// Flag indicates that porosity or cross_section changed during last time.
    	bool is_mass_diag_changed;

        /// Transported substances.
        SubstanceList substances_;

    	/// List of indices used to call balance methods for a set of quantities.
    	vector<unsigned int> subst_idx;

    	Vec mass_diag;  // diagonal entries in pass matrix (cross_section * porosity)

    	/// Flag indicates that sources part of equation was changed during last time.
    	bool sources_changed_;

        vector<VectorMPI> corr_vec;
        VectorMPI cfl_source_;      ///< Parallel vector for source term contribution to CFL condition.
        VectorMPI cfl_flow_;        ///< Parallel vector for flow contribution to CFL condition.
        Mat tm;                     ///< PETSc transport matrix
        vector<VectorMPI> tm_diag;  ///< additions to PETSC transport matrix on the diagonal - from sources (for each substance)
        Vec *bcvcorr;               ///< Boundary condition correction vector
        double transport_bc_time;   ///< Time of the last update of the boundary condition terms.
		TimeGovernor *time_;

	    /// Time when the transport matrix was created.
	    /// TODO: when we have our own classes for LA objects, we can use lazy dependence to check
	    /// necessity for matrix update
	    double transport_matrix_time;

		bool is_convection_matrix_scaled;   ///< Flag indicates the state of object

		/// Maximal number of edge sides (evaluate from dim 1,2,3)
		unsigned int max_edg_sides;

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

    /** Compute P0 interpolation of the solution (used reaction term).
     * Empty - solution is already P0 interpolation.
     */
    void compute_p0_interpolation() override {};

    /// Not used in this class.
	void update_after_reactions(bool) override {};

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
	{ return eq_data_->subst_idx; }

    /**
     * @brief Write computed fields.
     */
    virtual void output_data() override;

    void set_output_stream(std::shared_ptr<OutputTime> stream) override
    { output_stream_ = stream; }


	/**
	 * Getters.
	 */

    /// Getter for P0 interpolation by FieldFE.
	FieldFEScalarVec& get_p0_interpolation() override;

	Vec get_component_vec(unsigned int sbi) override;

    /// Returns number of transported substances.
    inline unsigned int n_substances() override
    { return eq_data_->substances_.size(); }

    /// Returns reference to the vector of substance names.
    inline SubstanceList &substances() override
    { return eq_data_->substances_; }

private:

    void alloc_transport_vectors();
    void alloc_transport_structs_mpi();



    /// Registrar of class to factory
    static const int registrar;

    /**
     *  Parameters of the equation, some are shared with other implementations since EqFields is derived from TransportBase::TransportEqFields
     */
    std::shared_ptr<EqFields> eq_fields_;
    std::shared_ptr<EqData> eq_data_;

    //@{
    /**
     * Flag indicates the state of object (transport matrix or source or boundary term).
     * If false, the object is freshly assembled and not rescaled.
     * If true, the object is scaled (not necessarily with the current time step).
     */
	bool is_src_term_scaled, is_bc_term_scaled;
	
    //@}

    TimeMark::Type target_mark_type;    ///< TimeMark type for time marks denoting end of every time interval where transport matrix remains constant.
    double cfl_max_step;    ///< Time step constraint coming from CFL condition.
    

    VecScatter vconc_out_scatter;
    Vec vpmass_diag;  // diagonal entries in mass matrix from last time (cross_section * porosity)

    ///
    Vec *vpconc; // previous concentration vector
    Vec *vcumulative_corr;

	/// Record with input specification.
	const Input::Record input_rec;

	std::shared_ptr<OutputTime> output_stream_;


    /// general assembly objects, hold assembly objects of appropriate dimension
    GenericAssembly< MassAssemblyConvection > * mass_assembly_;
    GenericAssembly< InitCondAssemblyConvection > * init_cond_assembly_;
    GenericAssembly< ConcSourcesBdrAssemblyConvection > * conc_sources_bdr_assembly_;
    GenericAssembly< MatrixMpiAssemblyConvection > * matrix_mpi_assembly_;

    friend class TransportOperatorSplitting;
};
#endif /* TRANSPORT_H_ */
