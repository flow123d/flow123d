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
 * @file    transport_operator_splitting.hh
 * @brief   
 */

#ifndef TRANSPORT_OPERATOR_SPLITTING_HH_
#define TRANSPORT_OPERATOR_SPLITTING_HH_

#include <boost/exception/info.hpp>             // for operator<<, error_inf...
#include <memory>                               // for shared_ptr
#include <vector>                               // for vector
#include "coupling/equation.hh"
#include "fields/field.hh"                      // for Field
#include "fields/field_values.hh"
#include "fields/field_set.hh"
#include "fields/multi_field.hh"
#include "transport/advection_process_base.hh"
#include "input/accessors.hh"                   // for Record
#include "input/type_base.hh"                   // for Array
#include "input/type_generic.hh"                // for Instance
#include "petscvec.h"                           // for Vec
#include "tools/time_governor.hh"               // for TimeGovernor, TimeGov...
#include "tools/time_marks.hh"                  // for TimeMarks
#include "system/index_types.hh"

/// external types:
class Mesh;
class ReactionTerm;
class Balance;
class Distribution;
class OutputTime;
class SubstanceList;
namespace Input {
	namespace Type {
		class Abstract;
		class Record;
	}
}






/**
 * Abstract interface class for implementations of transport equation within TransportOperatorSplitting.
 */
class ConcentrationTransportBase : public EquationBase {
public:

    /**
     * Constructor.
     */
    ConcentrationTransportBase(Mesh &init_mesh, const Input::Record in_rec)
	: EquationBase(init_mesh, in_rec) {};


	/// Common specification of the input record for secondary equations.
    static Input::Type::Abstract & get_input_type();


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
    virtual void set_target_time(double target_time) = 0;

    /**
     * Use Balance object from upstream equation (e.g. in various couplings) instead of own instance.
     */
    virtual void set_balance_object(std::shared_ptr<Balance> balance) = 0;
    
    /// Computes a constraint for time step.
    virtual bool evaluate_time_constraint(double &time_constraint) = 0;

    /// Return substance indices used in balance.
    virtual const vector<unsigned int> &get_subst_idx() = 0;

    /// Calculate the array of concentrations per element (for reactions).
    virtual void calculate_concentration_matrix() = 0;

    /// Perform changes to transport solution after reaction step.
    virtual void update_after_reactions(bool solution_changed) = 0;

    /// Setter for output stream.
    virtual void set_output_stream(std::shared_ptr<OutputTime> stream) = 0;

    /// Getter for output stream.
    virtual std::shared_ptr<OutputTime> output_stream() = 0;

    /// Getter for array of concentrations per element.
	virtual double **get_concentration_matrix() = 0;

	/// Return PETSc vector with solution for sbi-th substance.
	virtual const Vec &get_solution(unsigned int sbi) = 0;

	/// Return array of indices of local elements and parallel distribution of elements.
	virtual void get_par_info(LongIdx * &el_4_loc, Distribution * &el_ds) = 0;

	/// Return global array of order of elements within parallel vector.
	virtual LongIdx *get_row_4_el() = 0;

	/// Pass velocity from flow to transport.
    virtual void set_velocity_field(std::shared_ptr<FieldFE<3, FieldValue<3>::VectorFixed>> flux_field) = 0;

    /// Returns number of trnasported substances.
    virtual unsigned int n_substances() = 0;

    /// Returns reference to the vector of substnace names.
    virtual SubstanceList &substances() = 0;


};





/**
 * Class with fields that are common to all transport models.
 */
class TransportEqData : public FieldSet {
public:

	TransportEqData();
	inline virtual ~TransportEqData() {};

	/// Mobile porosity - usually saturated water content in the case of unsaturated flow model
	Field<3, FieldValue<3>::Scalar> porosity;

	/// Water content - result of unsaturated water flow model or porosity
	Field<3, FieldValue<3>::Scalar> water_content;

	/// Pointer to DarcyFlow field cross_section
	Field<3, FieldValue<3>::Scalar > cross_section;

	/// Concentration sources - density of substance source, only positive part is used.
	MultiField<3, FieldValue<3>::Scalar> sources_density;
	/// Concentration sources - Robin type, in_flux = sources_sigma * (sources_conc - mobile_conc)
	MultiField<3, FieldValue<3>::Scalar> sources_sigma;
	MultiField<3, FieldValue<3>::Scalar> sources_conc;

};




/**
 * @brief Empty transport class.
 */
class TransportNothing : public AdvectionProcessBase {
public:
    inline TransportNothing(Mesh &mesh_in)
    : AdvectionProcessBase(mesh_in, Input::Record() )

    {
        // make module solved for ever
        TimeGovernor::marks().new_mark_type();
        // auto eq_mark_type = TimeGovernor::marks().new_mark_type();
        time_= new TimeGovernor(TimeGovernor::inf_time, TimeGovernor::inf_time);
        time_->next_time();
    };

    inline virtual ~TransportNothing()
    {
        if(time_) delete time_;
    }

    inline void set_velocity_field(FMT_UNUSED std::shared_ptr<FieldFE<3, FieldValue<3>::VectorFixed>> flux_field) override {};

    inline virtual void output_data() override {};


};



/**
 * @brief Coupling of a transport model with a reaction model by operator splitting.
 *
 * Outline:
 * Transport model is any descendant of TransportBase (even TransportOperatorSplitting itself). This
 * should perform the transport possibly with diffusion and usually without coupling between substances and phases.
 *
 * Reaction is any descendant of the ReactionBase class. This represents reactions in general way of any coupling that
 * happens between substances and phases on one element or more generally on one DoF.
 */

class TransportOperatorSplitting : public AdvectionProcessBase {
public:
    typedef AdvectionProcessBase FactoryBaseType;

    /**
     * @brief Declare input record type for the equation TransportOperatorSplittiong.
     *
     * TODO: The question is if this should be a general coupling class
     * (e.g. allow coupling TranportDG with reactions even if it is not good idea for numerical reasons.)
     * To make this a coupling class we should modify all main input files for transport problems.
     */
    static const Input::Type::Record & get_input_type();

    /// Constructor.
    TransportOperatorSplitting(Mesh &init_mesh, const Input::Record in_rec);
    /// Destructor.
    virtual ~TransportOperatorSplitting();

    virtual void set_velocity_field(std::shared_ptr<FieldFE<3, FieldValue<3>::VectorFixed>> flux_field) override;

    void initialize() override;
    void zero_time_step() override;
    void update_solution() override;

    void compute_until_save_time();
    void compute_internal_step();
    void output_data() override;

   

private:
    /// Registrar of class to factory
    static const int registrar;

    std::shared_ptr<ConcentrationTransportBase> convection;
    std::shared_ptr<ReactionTerm> reaction;

    //double *** semchem_conc_ptr;   //dumb 3-dim array (for phases, which are not supported any more)
    //Semchem_interface *Semchem_reactions;
    
    double cfl_convection; ///< Time restriction due to transport
    double cfl_reaction;   ///< Time restriction due to reactions

};





#endif // TRANSPORT_OPERATOR_SPLITTING_HH_
