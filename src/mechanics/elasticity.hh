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
#include "fields/multi_field.hh"
// #include "fields/vec_seq_double.hh"
#include "la/linsys.hh"
#include "la/vector_mpi.hh"
#include "flow/mh_dofhandler.hh"
#include "fields/equation_output.hh"
#include "fem/mapping_p1.hh"
#include "coupling/equation.hh"
#include "fem/fe_values_views.hh"

class Distribution;
class OutputTime;
class DOFHandlerMultiDim;
template<unsigned int dim, unsigned int spacedim> class FEValuesBase;
template<unsigned int dim> class FiniteElement;
template<unsigned int dim, unsigned int spacedim> class Mapping;
template<unsigned int dim> class Quadrature;
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
	inline FiniteElement<dim> *fe();

	template<unsigned int dim>
	inline FiniteElement<dim> *fe_rt();

	template<unsigned int dim>
	inline Quadrature<dim> *q();

	template<unsigned int dim>
	inline MappingP1<dim,3> *mapping();

	inline std::shared_ptr<DOFHandlerMultiDim> dh();
    
//     const FEValuesViews::Vector<dim,3> vec;

private:

	/// Finite elements for the solution of the advection-diffusion equation.
    FiniteElement<0> *fe0_;
	FiniteElement<1> *fe1_;
	FiniteElement<2> *fe2_;
	FiniteElement<3> *fe3_;

	/// Finite elements for the water velocity field.
	FiniteElement<1> *fe_rt1_;
	FiniteElement<2> *fe_rt2_;
	FiniteElement<3> *fe_rt3_;

	/// Quadratures used in assembling methods.
	Quadrature<0> *q0_;
	Quadrature<1> *q1_;
	Quadrature<2> *q2_;
	Quadrature<3> *q3_;

	/// Auxiliary mappings of reference elements.
	MappingP1<1,3> *map1_;
	MappingP1<2,3> *map2_;
	MappingP1<3,3> *map3_;

    std::shared_ptr<DiscreteSpace> ds_;
    
	/// Object for distribution of dofs.
	std::shared_ptr<DOFHandlerMultiDim> dh_;
};


class VolumeChange
{
public:
  
  void init(Mesh *mesh);
  
  double value(const ElementAccessor<3> &elem) const;

protected:
  
  void set_values(const vector<double> &new_values, double new_time);
  
  VectorMPI values_;
  VectorMPI old_values_;
  Mesh *mesh_;
  double time_;
  double old_time_;
  
  friend class ::Elasticity;
};



} // namespace Mechanics





class Elasticity : public EquationBase
{
public:

	class EqData : public FieldSet {
	public:
      
        enum Bc_types {
          bc_type_displacement,
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
        Field<3, FieldValue<3>::Scalar> biot_alpha;
		Field<3, FieldValue<3>::Scalar> fracture_sigma;    ///< Transition parameter for diffusive transfer on fractures.
		
		/// Pointer to DarcyFlow field cross_section
        Field<3, FieldValue<3>::Scalar > cross_section;
        Field<3, FieldValue<3>::Scalar> region_id;
        Field<3, FieldValue<3>::Scalar> subdomain;
        
        Field<3, FieldValue<3>::VectorFixed> output_field;

        EquationOutput output_fields;

	};



    /**
     * @brief Constructor.
     * @param init_mesh         computational mesh
     * @param in_rec            input record
     */
    Elasticity(Mesh &init_mesh, const Input::Record in_rec);
    /**

     * @brief Declare input record type for the equation TransportDG.
     */
    static const Input::Type::Record & get_input_type();

    /**
     * @brief Initialize solution in the zero time.
     */
	void zero_time_step() override;
	
    bool evaluate_time_constraint(double &time_constraint)
    { return false; }

    /**
     * @brief Computes the solution in one time instant.
     */
	void update_solution() override;

	/**
	 * @brief Postprocesses the solution and writes to output file.
	 */
	void output_data();

	/**
	 * @brief Destructor.
	 */
	~Elasticity();

	void initialize() override;

    void calculate_cumulative_balance();

	const Vec &get_solution()
	{ return ls->get_solution(); }

    void get_par_info(int * &el_4_loc, Distribution * &el_ds);

    int *get_row_4_el();
    
    inline EqData &data() { return data_; }
    
    inline void set_velocity_field(const MH_DofHandler &dh)
    {
        mh_dh = &dh;
        flux_changed = true;
    }
    
    const Mechanics::VolumeChange &get_volume_change() { return volume_change; }
    
    
    typedef Elasticity FactoryBaseType;




private:
    /// Registrar of class to factory
    static const int registrar;

	void output_vector_gather();

	void preallocate();

	/**
	 * @brief Assembles the mass matrix.
	 *
	 * The routine just calls templated method assemble_mass_matrix() for each
	 * space dimension.
	 */
// 	void assemble_mass_matrix();

	/**
	 * @brief Assembles the mass matrix for the given dimension.
	 */
// 	template<unsigned int dim>
// 	void assemble_mass_matrix();

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
	 * @brief Assembles the right hand side due to volume sources.
	 *
	 * This method just calls set_sources() for each space dimension.
	 */
	void set_sources();

	/**
	 * @brief Assembles the right hand side vector due to volume sources.
	 */
	template<unsigned int dim>
	void set_sources();
    
    /**
     * @brief Assembles the fluxes on the boundary.
     */
    template<unsigned int dim>
    void assemble_fluxes_boundary();

	/**
	 * @brief Assembles the fluxes between elements of different dimensions.
	 */
	template<unsigned int dim>
	void assemble_fluxes_element_side();


	/**
	 * @brief Assembles the r.h.s. components corresponding to the Dirichlet boundary conditions.
	 *
	 * The routine just calls templated method set_boundary_condition() for each space dimension.
	 */
	void set_boundary_conditions();

	/**
	 * @brief Assembles the r.h.s. components corresponding to the Dirichlet boundary conditions
	 * for a given space dimension.
	 */
	template<unsigned int dim>
	void set_boundary_conditions();

	/**
	 * @brief Calculates the velocity field on a given @p dim dimensional cell.
	 *
	 * @param cell     The cell.
	 * @param velocity The computed velocity field (at quadrature points).
	 * @param fv       The FEValues class providing the quadrature points
	 *                 and the shape functions for velocity.
	 */
	template<unsigned int dim>
	void calculate_velocity(const ElementAccessor<3> &cell, std::vector<arma::vec3> &velocity, FEValuesBase<dim,3> &fv);

	/**
	 * @brief Sets the initial condition.
	 */
	void set_initial_condition();

	/**
	 * @brief Assembles the auxiliary linear system to calculate the initial solution
	 * as L^2-projection of the prescribed initial condition.
	 */
	template<unsigned int dim>
	void prepare_initial_condition();
    
    
    void update_volume_change();
    
    template<unsigned int dim>
    void update_volume_change(vector<double> &divergence);



	/// @name Physical parameters
	// @{

	/// Field data for model parameters.
	EqData data_;
    
        /**
     * Temporary solution how to pass velocity field form the flow model.
     * TODO: introduce FieldDiscrete -containing true DOFHandler and data vector and pass such object together with other
     * data. Possibly make more general set_data method, allowing setting data given by name. needs support from EqDataBase.
     */
    const MH_DofHandler *mh_dh;
    
    Mechanics::VolumeChange volume_change;

	// @}


	/// @name Parameters of the numerical method
	// @{

	/// Finite element objects
	Mechanics::FEObjects *feo;
    
    const double gamma = 1e3;

	// @}



	/// @name Solution of algebraic system
	// @{

	/// Vector of right hand side.
	Vec rhs;

	/// The stiffness matrix.
	Mat stiffness_matrix;

	/// The mass matrix.
	Mat mass_matrix;
	
	/// Mass from previous time instant (necessary when coefficients of mass matrix change in time).
// 	Vec mass_vec;

	/// Linear algebra system for the transport equation.
	LinSys *ls;

	/// Linear algebra system for the time derivative (actually it is used only for handling the matrix structures).
	LinSys *ls_dt;

	// @}


	/// @name Output to file
	// @{

	/// Array for storing the output solution data.
	//vector<double*> output_solution;

	/// Vector of solution data.
	VectorMPI output_vec;
    
    std::shared_ptr<OutputTime> output_stream_;

	/// Record with input specification.
	Input::Record input_rec;


	// @}


	/// @name Auxiliary fields used during assembly
	// @{

	/// Mass matrix coefficients.
// 	vector<double> mm_coef;
	/// Advection coefficients.
	vector<arma::vec3> ad_coef;
	/// Diffusion coefficients.
	vector<arma::mat33> dif_coef;
	/// Advection coefficients on edges.
	vector<vector<arma::vec3> > ad_coef_edg;
	/// Diffusion coefficients on edges.
	vector<vector<arma::mat33> > dif_coef_edg;

	// @}




	/// @name Other
	// @{

    /// Indicates whether matrices have been preallocated.
    bool allocation_done;
    
    /// Index used to call balance methods for a set of quantities.
    unsigned int subst_idx;
    
    /// Indicator of change in advection vector field.
    bool flux_changed;

    // @}
};






#endif /* ELASTICITY_HH_ */
