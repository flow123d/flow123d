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


namespace Mechanics {

/**
 * Auxiliary container class for Finite element and related objects of all dimensions.
 * Its purpose is to provide templated access to these objects, applicable in
 * the assembling methods.
 */
class FEObjects {
public:

	FEObjects(Mesh *mesh_, unsigned int fe_order);
	~FEObjects();

	template<unsigned int dim>
	inline std::shared_ptr<FiniteElement<dim>> fe();

	template<unsigned int dim>
	inline Quadrature *q() { return &(q_[dim]); }

	inline std::shared_ptr<DOFHandlerMultiDim> dh();
    inline std::shared_ptr<DOFHandlerMultiDim> dh_scalar();
    inline std::shared_ptr<DOFHandlerMultiDim> dh_tensor();
    
//     const FEValuesViews::Vector<dim,3> vec;

private:

        MixedPtr<FiniteElement> fe_;  ///< Finite elements for the solution of the mechanics equation.
	QGauss::array q_;


    std::shared_ptr<DiscreteSpace> ds_;
    
	/// Object for distribution of dofs.
	std::shared_ptr<DOFHandlerMultiDim> dh_;
    std::shared_ptr<DOFHandlerMultiDim> dh_scalar_;
    std::shared_ptr<DOFHandlerMultiDim> dh_tensor_;
};


} // namespace Mechanics





class Elasticity : public EquationBase
{
public:

	class EqData : public FieldSet {
	public:
      
        enum Bc_types {
          bc_type_displacement,
          bc_type_displacement_normal,
          bc_type_traction
        };

		EqData();
        
        static  constexpr const char *  name() { return "Mechanics_LinearElasticity"; }

        static string default_output_field() { return "\"displacement\""; }

        static const Input::Type::Selection & get_bc_type_selection();

        static IT::Selection get_output_selection();
        
        
        BCField<3, FieldValue<3>::Enum > bc_type;
        BCField<3, FieldValue<3>::VectorFixed> bc_displacement;
        BCField<3, FieldValue<3>::VectorFixed> bc_traction;
        Field<3, FieldValue<3>::VectorFixed> load;
        Field<3, FieldValue<3>::Scalar> young_modulus;
        Field<3, FieldValue<3>::Scalar> poisson_ratio;
		Field<3, FieldValue<3>::Scalar> fracture_sigma;    ///< Transition parameter for diffusive transfer on fractures.
		
		/// Pointer to DarcyFlow field cross_section
        Field<3, FieldValue<3>::Scalar > cross_section;
        Field<3, FieldValue<3>::Scalar > potential_load;   ///< Potential of an additional (external) load.
        Field<3, FieldValue<3>::Scalar> region_id;
        Field<3, FieldValue<3>::Scalar> subdomain;
        
        Field<3, FieldValue<3>::VectorFixed> output_field;
        Field<3, FieldValue<3>::TensorFixed> output_stress;
        Field<3, FieldValue<3>::Scalar> output_von_mises_stress;
        Field<3, FieldValue<3>::Scalar> output_cross_section;
        Field<3, FieldValue<3>::Scalar> output_divergence;
        
        std::shared_ptr<FieldFE<3, FieldValue<3>::VectorFixed> > output_field_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::TensorFixed> > output_stress_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > output_von_mises_stress_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > output_cross_section_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > output_div_ptr;

        EquationOutput output_fields;

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
	
    bool evaluate_time_constraint(double &)
    { return false; }

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
    
    void set_potential_load(const Field<3, FieldValue<3>::Scalar> &potential)
    { data_.potential_load = potential; }

    void calculate_cumulative_balance();

	const Vec &get_solution()
	{ return ls->get_solution(); }

    inline EqData &data() { return data_; }
    
    
    typedef Elasticity FactoryBaseType;




private:
    /// Registrar of class to factory
    static const int registrar;

    template<unsigned int dim>
    void compute_output_fields();

	void preallocate();

	/**
	 * @brief Assembles the stiffness matrix.
	 *
	 * This routine just calls assemble_volume_integrals(), assemble_fluxes_boundary(),
	 * assemble_fluxes_element_element() and assemble_fluxes_element_side() for each
	 * space dimension.
	 */
	void assemble_stiffness_matrix();

	/**
	 * @brief Assembles the volume integrals into the stiffness matrix.
	*/
	template<unsigned int dim>
	void assemble_volume_integrals();

	/**
	 * @brief Assembles the right hand side (forces, boundary conditions, tractions).
	 */
	void assemble_rhs();

	/**
	 * @brief Assembles the right hand side vector due to volume sources.
	 */
	template<unsigned int dim>
	void assemble_sources();
    
    /**
     * @brief Assembles the fluxes on the boundary.
     */
    template<unsigned int dim>
    void assemble_fluxes_boundary();

	/**
	 * @brief Assembles the fluxes between elements of different dimensions depending on displacement.
	 */
	template<unsigned int dim>
	void assemble_matrix_element_side();
    
    /** @brief Assemble fluxes between different dimensions that are independent of displacement. */
    template<unsigned int dim>
    void assemble_rhs_element_side();


	/**
	 * @brief Assembles the r.h.s. components corresponding to the Dirichlet boundary conditions
	 * for a given space dimension.
	 */
	template<unsigned int dim>
	void assemble_boundary_conditions();
    
    /** @brief Penalty to enforce boundary value in weak sense. */
    double dirichlet_penalty(SideIter side);
    
    


	/// @name Physical parameters
	// @{

	/// Field data for model parameters.
	EqData data_;
    
	// @}


	/// @name Parameters of the numerical method
	// @{

	/// Finite element objects
	Mechanics::FEObjects *feo;
    
	// @}



	/// @name Solution of algebraic system
	// @{

	/// Vector of right hand side.
	Vec rhs;

	/// The stiffness matrix.
	Mat stiffness_matrix;

	/// Linear algebra system for the transport equation.
	LinSys *ls;

	// @}


	/// @name Output to file
	// @{

    std::shared_ptr<OutputTime> output_stream_;

	/// Record with input specification.
	Input::Record input_rec;


	// @}




	/// @name Other
	// @{

    /// Indicates whether matrices have been preallocated.
    bool allocation_done;
    
    // @}
};



double lame_mu(double young, double poisson);
double lame_lambda(double young, double poisson);






#endif /* ELASTICITY_HH_ */
