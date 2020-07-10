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
 * @file    heat_model.hh
 * @brief   Discontinuous Galerkin method for equation of transport with dispersion.
 * @author  Jan Stebel
 */

#ifndef HEAT_MODEL_HH_
#define HEAT_MODEL_HH_

#include <boost/exception/info.hpp>             // for operator<<, error_inf...
#include <memory>                               // for shared_ptr
#include <string>                               // for string
#include <vector>                               // for vector
#include <armadillo>
#include "advection_diffusion_model.hh"
#include "transport/advection_process_base.hh"
#include "transport/substance.hh"
#include "fields/field_values.hh"               // for FieldValue<>::Scalar
#include "fields/bc_field.hh"
#include "fields/bc_multi_field.hh"
#include "fields/field.hh"
#include "fields/multi_field.hh"
#include "fields/field_set.hh"
#include "input/type_base.hh"                   // for Array
#include "input/type_generic.hh"                // for Instance
#include "input/type_record.hh"                 // for Record
#include "input/type_selection.hh"              // for Selection

class Mesh;
class OutputTime;
namespace Input { class Record; }
template <int spacedim> class ElementAccessor;


// class HeatProcessBase : public EquationBase
// {
// public:
//     typedef HeatProcessBase FactoryBaseType;


//     HeatProcessBase(Mesh &mesh, const Input::Record in_rec)
//     : EquationBase(mesh, in_rec)
//     {};

//     /**
//      * This method takes sequential PETSc vector of side velocities and update
//      * transport matrix. The ordering is same as ordering of sides in the mesh.
//      * We just keep the pointer, but do not destroy the object.
//      *
//      * TODO: We should pass whole velocity field object (description of base functions and dof numbering) and vector.
//      */
//     virtual void set_velocity_changed() = 0;

//     /// Common specification of the input record for secondary equations.
//     static Input::Type::Abstract & get_input_type() {
//         return Input::Type::Abstract("Heat",
//                 "Equation for heat transfer.")
//                 .close();
//     }
// };

/*
class HeatNothing : public HeatProcessBase {
public:
    inline HeatNothing(Mesh &mesh_in)
    : HeatProcessBase(mesh_in, Input::Record() )

    {
        // make module solved for ever
        time_=new TimeGovernor();
        time_->next_time();
    };

    inline virtual ~HeatNothing()
    {}

    inline void set_velocity_changed() override {};

    inline virtual void output_data() override {};

};
*/

class HeatTransferModel : public AdvectionDiffusionModel, public AdvectionProcessBase {
public:

	class ModelEqData : public FieldSet {
	public:

		enum Heat_bc_types {
			bc_inflow,
			bc_dirichlet,
			bc_total_flux,
			bc_diffusive_flux
		};

		/// Type of boundary condition (see also BC_Type)
		BCMultiField<3, FieldValue<3>::Enum > bc_type;
		/// Dirichlet boundary condition for temperature.
		BCMultiField<3, FieldValue<3>::Scalar> bc_dirichlet_value;
		/// Flux value in total/diffusive flux b.c.
		BCMultiField<3, FieldValue<3>::Scalar > bc_flux;
		/// Transition coefficient in total/diffusive flux b.c.
		BCMultiField<3, FieldValue<3>::Scalar > bc_robin_sigma;
		/// Initial temperature.
		MultiField<3, FieldValue<3>::Scalar> init_condition;
		/// Porosity of solid.
		Field<3, FieldValue<3>::Scalar> porosity;
		/// Water content passed from Darcy flow model
		Field<3, FieldValue<3>::Scalar> water_content;
		/// Flow flux, can be result of water flow model.
    	Field<3, FieldValue<3>::VectorFixed> flow_flux;
		/// Density of fluid.
		Field<3, FieldValue<3>::Scalar> fluid_density;
		/// Heat capacity of fluid.
		Field<3, FieldValue<3>::Scalar> fluid_heat_capacity;
		/// Heat conductivity of fluid.
		Field<3, FieldValue<3>::Scalar> fluid_heat_conductivity;
		/// Density of solid.
		Field<3, FieldValue<3>::Scalar> solid_density;
		/// Heat capacity of solid.
		Field<3, FieldValue<3>::Scalar> solid_heat_capacity;
		/// Heat conductivity of solid.
		Field<3, FieldValue<3>::Scalar> solid_heat_conductivity;
		/// Longitudal heat dispersivity.
		Field<3, FieldValue<3>::Scalar> disp_l;
		/// Transversal heat dispersivity.
		Field<3, FieldValue<3>::Scalar> disp_t;
		/// Thermal source in fluid.
		Field<3, FieldValue<3>::Scalar> fluid_thermal_source;
		/// Thermal source in solid.
		Field<3, FieldValue<3>::Scalar> solid_thermal_source;
		/// Heat exchange rate in fluid.
		Field<3, FieldValue<3>::Scalar> fluid_heat_exchange_rate;
		/// Heat exchange rate in solid.
		Field<3, FieldValue<3>::Scalar> solid_heat_exchange_rate;
		/// Reference temperature in fluid.
		Field<3, FieldValue<3>::Scalar> fluid_ref_temperature;
		/// Reference temperature in solid.
		Field<3, FieldValue<3>::Scalar> solid_ref_temperature;

		/// Pointer to DarcyFlow field cross_section
		Field<3, FieldValue<3>::Scalar > cross_section;


		MultiField<3, FieldValue<3>::Scalar> output_field;


		/// @name Instances of FieldModel used in assembly methods
		// @{

		/// Velocity norm field.
        Field<3, FieldValue<3>::Scalar > v_norm;
		/// Field represents coefficients of mass matrix.
        Field<3, FieldValue<3>::Scalar > mass_matrix_coef;
		/// Field represents retardation coefficients due to sorption.
        MultiField<3, FieldValue<3>::Scalar> retardation_coef;
    	/// Concentration sources - density output
    	MultiField<3, FieldValue<3>::Scalar> sources_density_out;
    	/// Concentration sources - sigma output
    	MultiField<3, FieldValue<3>::Scalar> sources_sigma_out;
    	/// Concentration sources - concentration output
    	MultiField<3, FieldValue<3>::Scalar> sources_conc_out;
		/// Advection coefficients.
		MultiField<3, FieldValue<3>::VectorFixed> advection_coef;
		/// Diffusion coefficients.
		MultiField<3, FieldValue<3>::TensorFixed> diffusion_coef;

    	// @}



		ModelEqData();

		static  constexpr const char *  name() { return "Heat_AdvectionDiffusion"; }

		static string default_output_field() { return "\"temperature\""; }

        static const Input::Type::Selection & get_bc_type_selection();

		static IT::Selection get_output_selection();

		/**
		 * Initialize FieldModel instances.
		 */
		void initialize();

        /// Returns number of transported substances.
        inline unsigned int n_substances()
        { return 1; }


		/// @name Data of substances
		// @{

	    /// Transported substances.
	    SubstanceList substances_;

		/// List of indices used to call balance methods for a set of quantities.
		vector<unsigned int> subst_idx;

    	// @}
	};

	typedef AdvectionProcessBase FactoryBaseType;


	HeatTransferModel(Mesh &mesh, const Input::Record in_rec);

	void init_from_input(const Input::Record &) override {};

	~HeatTransferModel() override;

    /// Returns reference to the vector of substance names.
    inline SubstanceList &substances()
    { return data().substances_; }

    const vector<unsigned int> &get_subst_idx()
	{ return data().subst_idx; }


	/// Derived class should implement getter for ModelEqData instance.
	virtual ModelEqData &data() = 0;

protected:
	
	/**
	 * Create input type that can be passed to the derived class.
	 * @param implementation String characterizing the numerical method, e.g. DG, FEM, FVM.
	 * @param description    Comment used to describe the record key.
	 * @return
	 */
	static IT::Record get_input_type(const string &implementation, const string &description);

	void output_data() override;

	std::shared_ptr<OutputTime> &output_stream()
	{ return output_stream_; }

	virtual void calculate_cumulative_balance() = 0;

	/**
	 * Temporary method, sets balance object after construction of EqData object.
	 */
	void init_balance(const Input::Record &in_rec);

	std::shared_ptr<OutputTime> output_stream_;
};





#endif /* HEAT_MODEL_HH_ */
