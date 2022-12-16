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
 * @file    concentration_model.hh
 * @brief   Discontinuous Galerkin method for equation of transport with dispersion.
 * @author  Jan Stebel
 */

#ifndef CONC_DISP_MODEL_HH_
#define CONC_DISP_MODEL_HH_


#include <memory>                                     // for shared_ptr
#include <string>                                     // for string
#include <vector>                                     // for vector
#include <armadillo>
#include "advection_diffusion_model.hh"
#include "fields/field_values.hh"                     // for FieldValue<>::S...
#include "fields/bc_multi_field.hh"
#include "fields/field.hh"
#include "fields/multi_field.hh"

#include "input/accessors.hh"                  // for ExcAccessorForNullStorage
#include "input/accessors_impl.hh"             // for Record::val
#include "input/type_base.hh"                         // for Array
#include "input/type_generic.hh"                      // for Instance
#include "input/type_record.hh"                       // for Record
#include "input/type_selection.hh"                    // for Selection
#include "transport/substance.hh"                     // for SubstanceList
#include "transport/transport_operator_splitting.hh"  // for ConcentrationTr...

class Balance;
class Mesh;
class OutputTime;
namespace Input { class Record; }
template <int spacedim> class ElementAccessor;




class ConcDispersionModel : public AdvectionDiffusionModel, public ConcentrationTransportBase {
public:

	class ModelEqFields : public TransportEqFields {
	public:

		/// Type of boundary condition (see also BC_Type)
        BCMultiField<3, FieldValue<3>::Enum > bc_type;
		/// Prescribed concentration for Dirichlet/reference concentration for flux b.c.
		BCMultiField<3, FieldValue<3>::Scalar> bc_dirichlet_value;
		/// Flux value in total/diffusive flux b.c.
		BCMultiField<3, FieldValue<3>::Scalar > bc_flux;
		/// Transition coefficient in total/diffusive flux b.c.
		BCMultiField<3, FieldValue<3>::Scalar > bc_robin_sigma;
		/// Initial concentrations.
		MultiField<3, FieldValue<3>::Scalar> init_condition;

		/// Components of the full dispersion tensor
		Field<3, FieldValue<3>::Tensor > disp_uxx;
		Field<3, FieldValue<3>::Tensor > disp_uyy;
		Field<3, FieldValue<3>::Tensor > disp_uzz;
		Field<3, FieldValue<3>::Tensor > disp_uyz;
		Field<3, FieldValue<3>::Tensor > disp_uxz;
		Field<3, FieldValue<3>::Tensor > disp_uxy;
		/// Longitudal dispersivity (for each substance).

		//MultiField<3, FieldValue<3>::Scalar> disp_l;
		/// Transversal dispersivity (for each substance).
		//MultiField<3, FieldValue<3>::Scalar> disp_t;
		/// Molecular diffusivity (for each substance).
		//MultiField<3, FieldValue<3>::Tensor> diff_m;

	    Field<3, FieldValue<3>::Scalar > rock_density;      ///< Rock matrix density.
	    MultiField<3, FieldValue<3>::Scalar > sorption_coefficient;     ///< Coefficient of linear sorption.


		MultiField<3, FieldValue<3>::Scalar> output_field;


		/// @name Instances of FieldModel used in assembly methods
		// @{

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
		MultiField<3, FieldValue<3>::Vector> advection_coef;
		/// Diffusion coefficients.
		MultiField<3, FieldValue<3>::Tensor> diffusion_coef;
		/// Velocity norm field.
        Field<3, FieldValue<3>::Scalar > v_norm;

    	// @}

		enum Concentration_bc_types {
			bc_inflow,
			bc_dirichlet,
			bc_total_flux,
			bc_diffusive_flux
		};

		ModelEqFields();

        static const Input::Type::Selection & get_bc_type_selection();

		/**
		 * Initialize FieldModel instances.
		 */
		void initialize(Input::Record transport_rec);

	};

   	class ModelEqData {
   	public:

		ModelEqData() {}

		static constexpr const char * name() { return "Solute_ConcDispersion"; }

		static string default_output_field() { return "\"conc\""; }

		static IT::Selection get_output_selection();

        /// Returns number of transported substances.
        inline unsigned int n_substances()
        { return substances_.size(); }

        /// Returns reference to the vector of substance indices.
        const vector<unsigned int> &subst_idx()
    	{ return subst_idx_; }

        /// Returns reference to the vector of substance names.
        inline SubstanceList &substances()
        { return substances_; }


		/// @name Data of substances
		// @{

	    /// Transported substances.
	    SubstanceList substances_;

		/// List of indices used to call balance methods for a set of quantities.
		vector<unsigned int> subst_idx_;


    	// @}
	};


	typedef ConcentrationTransportBase FactoryBaseType;


	ConcDispersionModel(Mesh &mesh, const Input::Record &in_rec);

        void init_from_input(const Input::Record &in_rec) override;


	~ConcDispersionModel() override;

    /// Returns number of transported substances.
    inline unsigned int n_substances() override
    { return eq_data().n_substances(); }

    /// Returns reference to the vector of substance names.
    inline SubstanceList &substances() override
    { return eq_data().substances(); }


    // Methods inherited from ConcentrationTransportBase:

	// Must be implemented in descendants.
	void set_target_time(double) override {};

	void set_balance_object(std::shared_ptr<Balance> balance) override;

    const vector<unsigned int> &get_subst_idx() override
	{ return eq_data().subst_idx(); }

    void set_output_stream(std::shared_ptr<OutputTime> stream)
    { output_stream_ = stream; }


	/// Derived class should implement getter for ModelEqFields instance.
	virtual ModelEqFields &eq_fields() = 0;

	/// Derived class should implement getter for ModelEqData instance.
	virtual ModelEqData &eq_data() = 0;

protected:

	/**
	 * Create input type that can be passed to the derived class.
	 * @param implementation String characterizing the numerical method, e.g. DG, FEM, FVM.
	 * @param description    Comment used to describe the record key.
	 * @return
	 */
	static IT::Record get_input_type(const string &implementation, const string &description);

	/**
	 * Empty temporary method (must be implemented for continuity with HeatTransferModel)
	 */
	void init_balance(const Input::Record &in_rec);

	/// Density of liquid (a global constant).
	double solvent_density_;

	std::shared_ptr<OutputTime> output_stream_;
};








#endif /* CONC_DISP_MODEL_HH_ */
