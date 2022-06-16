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
 * @file    elasticity.hh
 * @brief   FEM for linear elasticity.
 * @author  Jan Stebel
 */

#ifndef ELASTICITY_HH_
#define ELASTICITY_HH_

#include "fields/bc_field.hh"
#include "fields/field.hh"
#include "fields/field_fe.hh"
#include "fields/multi_field.hh"
#include "la/linsys.hh"
#include "la/vector_mpi.hh"
#include "fields/equation_output.hh"
#include "coupling/equation.hh"
#include "fem/fe_values_views.hh"
#include "tools/mixed.hh"
#include "quadrature/quadrature_lib.hh"

class Distribution;
class OutputTime;
class DOFHandlerMultiDim;
template<unsigned int dim> class FiniteElement;
class Elasticity;
template<unsigned int dim> class StiffnessAssemblyElasticity;
template<unsigned int dim> class RhsAssemblyElasticity;
template<unsigned int dim> class ConstraintAssemblyElasticity;
template<unsigned int dim> class OutpuFieldsAssemblyElasticity;
template< template<IntDim...> class DimAssembly> class GenericAssembly;




class Elasticity : public EquationBase
{
public:

	class EqFields : public FieldSet {
	public:
      
        enum Bc_types {
          bc_type_displacement,
          bc_type_displacement_normal,
          bc_type_traction,
		  bc_type_stress,
        };

        EqFields();
        
        static const Input::Type::Selection & get_bc_type_selection();

        BCField<3, FieldValue<3>::Enum > bc_type;
        BCField<3, FieldValue<3>::VectorFixed> bc_displacement;
        BCField<3, FieldValue<3>::VectorFixed> bc_traction;
		BCField<3, FieldValue<3>::TensorFixed> bc_stress;
        Field<3, FieldValue<3>::VectorFixed> load;
        Field<3, FieldValue<3>::Scalar> young_modulus;
        Field<3, FieldValue<3>::Scalar> poisson_ratio;
		Field<3, FieldValue<3>::Scalar> fracture_sigma;    ///< Transition parameter for diffusive transfer on fractures.
        Field<3, FieldValue<3>::Scalar> roughness_angle;
        Field<3, FieldValue<3>::Scalar> roughness_height;
        Field<3, FieldValue<3>::TensorFixed> initial_stress;
		
		/// Pointer to DarcyFlow field cross_section
        Field<3, FieldValue<3>::Scalar > cross_section;
        Field<3, FieldValue<3>::Scalar > cross_section_min;
        Field<3, FieldValue<3>::Scalar > potential_load;   ///< Potential of an additional (external) load.
        Field<3, FieldValue<3>::Scalar > ref_potential_load; ///< Potential of reference external load on boundary. TODO: Switch to BCField when possible.
        Field<3, FieldValue<3>::Scalar> region_id;
        Field<3, FieldValue<3>::Scalar> subdomain;
        
        Field<3, FieldValue<3>::VectorFixed> output_field;
        Field<3, FieldValue<3>::TensorFixed> output_stress;
        Field<3, FieldValue<3>::Scalar> output_von_mises_stress;
        Field<3, FieldValue<3>::Scalar> output_mean_stress;
        Field<3, FieldValue<3>::Scalar> output_cross_section;
        Field<3, FieldValue<3>::Scalar> output_divergence;
        Field<3, FieldValue<3>::VectorFixed> output_displacement_jump;
        
		/// @name Instances of FieldModel used in assembly methods
		// @{

        Field<3, FieldValue<3>::Scalar > lame_mu;
        Field<3, FieldValue<3>::Scalar > lame_lambda;
        Field<3, FieldValue<3>::Scalar > dirichlet_penalty;

    	// @}

        std::shared_ptr<FieldFE<3, FieldValue<3>::VectorFixed> > output_field_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::TensorFixed> > output_stress_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > output_von_mises_stress_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > output_mean_stress_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > output_cross_section_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > output_div_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::VectorFixed> > output_displ_jump_ptr;

        EquationOutput output_fields;

	};

	class EqData {
	public:

		EqData()
        : ls(nullptr), constraint_matrix(nullptr), constraint_vec(nullptr) {}

		~EqData() {
		    if (ls!=nullptr) delete ls;
            if (constraint_matrix!=nullptr) MatDestroy(&constraint_matrix);
            if (constraint_vec!=nullptr) VecDestroy(&constraint_vec);
		}

		/// Create DOF handler objects
        void create_dh(Mesh * mesh, unsigned int fe_order);

        /// Objects for distribution of dofs.
        std::shared_ptr<DOFHandlerMultiDim> dh_;
        std::shared_ptr<DOFHandlerMultiDim> dh_scalar_;
        std::shared_ptr<DOFHandlerMultiDim> dh_vector_;
        std::shared_ptr<DOFHandlerMultiDim> dh_tensor_;

    	/// @name Solution of algebraic system
    	// @{

    	/// Linear algebraic system.
    	LinSys *ls;

        Mat constraint_matrix;
        Vec constraint_vec;

        // map local element -> constraint index
        std::map<LongIdx,LongIdx> constraint_idx;

    	// @}

    	/// Shared Balance object
    	std::shared_ptr<Balance> balance_;

	};


    /**
     * @brief Constructor.
     * @param init_mesh         computational mesh
     * @param in_rec            input record
     * @param tm                time governor (if nullptr then it is created from input record)
     */
    Elasticity(Mesh &init_mesh, const Input::Record in_rec, TimeGovernor *tm = nullptr);
    /**

     * @brief Declare input record type for the equation TransportDG.
     */
    static const Input::Type::Record & get_input_type();

    /**
     * @brief Initialize solution in the zero time.
     */
	void zero_time_step() override;

    /**
     * @brief Computes the solution in one time instant.
     */
	void update_solution() override;
    
    /// Pass to next time and update equation data.
    void next_time();
    
    /// Solve without updating time step and without output.
    void solve_linear_system();

	/**
	 * @brief Postprocesses the solution and writes to output file.
	 */
	void output_data();

	/**
	 * @brief Destructor.
	 */
	~Elasticity();

	void initialize() override;
    
	// Recompute fields for output (stress, divergence etc.)
	void update_output_fields();
    
    void set_potential_load(const Field<3, FieldValue<3>::Scalar> &potential,
                            const Field<3, FieldValue<3>::Scalar> &ref_potential)
    {
        eq_fields_->potential_load = potential;
        eq_fields_->ref_potential_load = ref_potential;
    }

    void calculate_cumulative_balance();

	const Vec &get_solution()
	{ return eq_data_->ls->get_solution(); }

	inline EqFields &eq_fields() { return *eq_fields_; }

	inline EqData &eq_data() { return *eq_data_; }

    
    
    typedef Elasticity FactoryBaseType;




private:
    /// Registrar of class to factory
    static const int registrar;

	void preallocate();


	void assemble_constraint_matrix();


	/// @name Physical parameters
	// @{

	/// Fields for model parameters.
	std::shared_ptr<EqFields> eq_fields_;

	/// Data for model parameters.
	std::shared_ptr<EqData> eq_data_;

    /// Indicator of contact conditions on fractures.
    bool has_contact_;

    
	// @}


	
    /// @name Output to file
	// @{

    std::shared_ptr<OutputTime> output_stream_;

	/// Record with input specification.
	Input::Record input_rec;


	// @}


    static constexpr const char *  name_ = "Mechanics_LinearElasticity";


    /// general assembly objects, hold assembly objects of appropriate dimension
    GenericAssembly< StiffnessAssemblyElasticity > * stiffness_assembly_;
    GenericAssembly< RhsAssemblyElasticity > * rhs_assembly_;
    GenericAssembly< ConstraintAssemblyElasticity > * constraint_assembly_;
    GenericAssembly< OutpuFieldsAssemblyElasticity > * output_fields_assembly_;

};


/*
 * TODO Remove these two methods after implementation new assembly algorithm in HM_Iterative class.
 */
double lame_mu(double young, double poisson);
double lame_lambda(double young, double poisson);






#endif /* ELASTICITY_HH_ */
