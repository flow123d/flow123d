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
 * @brief Discontinuous Galerkin method for equation of transport with dispersion.
 *  @author Jan Stebel
 */

#ifndef TRANSPORT_DG_HH_
#define TRANSPORT_DG_HH_

#include "transport_operator_splitting.hh"
#include "fields/bc_field.hh"
#include "fields/field.hh"
#include "la/linsys.hh"
#include "flow/mh_dofhandler.hh"

class Distribution;
class OutputTime;
class DOFHandlerMultiDim;
template<unsigned int dim, unsigned int spacedim> class FEValuesBase;
template<unsigned int dim, unsigned int spacedim> class FiniteElement;
template<unsigned int dim, unsigned int spacedim> class Mapping;
template<unsigned int dim> class Quadrature;



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
	inline FiniteElement<dim,3> *fe();

	template<unsigned int dim>
	inline FiniteElement<dim,3> *fe_rt();

	template<unsigned int dim>
	inline Quadrature<dim> *q();

	template<unsigned int dim>
	inline Mapping<dim,3> *mapping();

	inline DOFHandlerMultiDim *dh();

private:

	/// Finite elements for the solution of the advection-diffusion equation.
	FiniteElement<1,3> *fe1_;
	FiniteElement<2,3> *fe2_;
	FiniteElement<3,3> *fe3_;

	/// Finite elements for the water velocity field.
	FiniteElement<1,3> *fe_rt1_;
	FiniteElement<2,3> *fe_rt2_;
	FiniteElement<3,3> *fe_rt3_;

	/// Quadratures used in assembling methods.
	Quadrature<0> *q0_;
	Quadrature<1> *q1_;
	Quadrature<2> *q2_;
	Quadrature<3> *q3_;

	/// Auxiliary mappings of reference elements.
	Mapping<0,3> *map0_;
	Mapping<1,3> *map1_;
	Mapping<2,3> *map2_;
	Mapping<3,3> *map3_;

	/// Object for distribution of dofs.
	DOFHandlerMultiDim *dh_;
};



/**
 * @brief Transport with dispersion implemented using discontinuous Galerkin method.
 *
 * TransportDG implements the discontinuous Galerkin method for the transport and diffusion of substances.
 * The concentration @f$ c_i ~[kg/m^3]@f$ of the i-th substance is governed by the advection-diffusion equation
 * @f[
 * 		\partial_t c_i + \mathbf v\cdot\nabla c_i - \mathrm{div}(D\nabla c_i) = F \mbox{ in }\Omega^d,
 * @f]
 * where @f$\mathbf v@f$ is the fluid velocity and @f$\Omega^d@f$ the @f$d@f$-dimensional domain, respectively.
 * The hydrodynamic dispersivity tensor @f$\mathbf D ~[m^2/s]@f$ is given by:
 * @f[
 *		\mathbf D = D_m\mathbf I + |\mathbf v|\left(\alpha_T\mathbf I + (\alpha_L-\alpha_T)\frac{\mathbf v\otimes\mathbf v}{|\mathbf v|^2}\right).
 * @f]
 * The molecular dispersivity @f$D_m~[m^2/s]@f$, as well as the longitudal and transversal dispersivity @f$\alpha_L,~\alpha_T~[m]@f$ are input parameters of the model.
 *
 * For lower dimensions @f$d=1,2@f$ the advection-diffusion equation is multiplied by the fracture cross-cut @f$\delta^d~[m^{3-d}]@f$.
 *
 * The boundary @f$\partial\Omega^d@f$ is divided into three disjoint parts @f$\Gamma^d_D\cup\Gamma^d_N\cup\Gamma^d_F@f$.
 * We prescribe the following boundary conditions:
 * @f{eqnarray*}{
 * 		c_i^d &= c_{iD}^d &\mbox{ on }\Gamma^d_D \mbox{ (Dirichlet)},\\
 * 		\mathbf D^d\nabla c_i^d\cdot\mathbf n &= 0 &\mbox{ on }\Gamma^d_N \mbox{ (Neumann)},
 * @f}
 * The transfer of mass through fractures is described by the transmission conditions on @f$\Gamma^d_F@f$:
 * @f[
 * 		-\mathbf D^d\nabla c_i^d\cdot\mathbf n = \sigma(c_i^d-c_i^{d-1}) + \left\{\begin{array}{cl}0 &\mbox{ if }\mathbf v^d\cdot\mathbf n\ge 0\\\mathbf v^d\cdot\mathbf n(c_i^{d-1}-c_i^d) & \mbox{ if }\mathbf v^d\cdot\mathbf n<0\end{array}\right.,\qquad
 * 		F^{d-1} = (\sigma + |\mathbf v^d\cdot\mathbf n|)(c_i^d-c_i^{d-1}).
 * @f]
 * Here @f$\mathbf n@f$ stands for the unit outward normal vector to @f$\partial\Omega^d@f$.
 * The coefficient @f$\sigma@f$ determines the transfer of mass through fractures due to diffusion.
 *
 * @ingroup transport_mod
 *
 */
template<class Model>
class TransportDG : public TransportBase, public Model
{
public:
	typedef AdvectionProcessBase FactoryBaseType;

	class EqData : public Model::ModelEqData {
	public:

        enum BC_Type {
            none,
            inflow,
            dirichlet,
            neumann,
            robin
        };
        static const Input::Type::Selection & get_bc_type_selection();

        static const Input::Type::Selection & get_output_selection();

		EqData();

		Field<3, FieldValue<3>::Vector> fracture_sigma;    ///< Transition parameter for diffusive transfer on fractures (for each substance).
		Field<3, FieldValue<3>::Vector> dg_penalty;        ///< Penalty enforcing inter-element continuity of solution (for each substance).

        BCField<3, FieldValue<3>::EnumVector > bc_type;    ///< Type of boundary condition (see also BC_Type)
        BCField<3, FieldValue<3>::Vector > bc_flux;        ///< Flux in Neumann or Robin b.c.
        BCField<3, FieldValue<3>::Vector > bc_robin_sigma; ///< Transition coefficient in Robin b.c.

        Field<3, FieldValue<3>::Integer> region_id;

	};



	enum DGVariant {
		// Non-symmetric weighted interior penalty DG
		non_symmetric = -1,

		// Incomplete weighted interior penalty DG
		incomplete = 0,

		// Symmetric weighted interior penalty DG
		symmetric = 1
	};

    /**
     * @brief Constructor.
     * @param init_mesh         computational mesh
     * @param in_rec            input record
     */
    TransportDG(Mesh &init_mesh, const Input::Record &in_rec);
    /**

     * @brief Declare input record type for the equation TransportDG.
     */
    static const Input::Type::Record & get_input_type();

    /**
     * @brief Input type for the DG variant selection.
     */
    static const Input::Type::Selection & get_dg_variant_selection_input_type();

    /**
     * @brief Initialize solution in the zero time.
     */
	void zero_time_step() override;

    /**
     * @brief Computes the solution in one time instant.
     */
	void update_solution() override;

	/**
	 * @brief Updates the velocity field which determines some coefficients of the transport equation.
	 * 
         * @param dh mixed hybrid dof handler
         * 
	 * (So far it does not work since the flow module returns a vector of zeros.)
	 * @param velocity_vector Input array of velocity values.
	 */
	virtual void set_velocity_field(const MH_DofHandler &dh);

	/**
	 * @brief Postprocesses the solution and writes to output file.
	 */
	void output_data();

	/**
	 * @brief Getter for field data.
	 */
	virtual EqData *get_data() { return &data_; }

	/**
	 * @brief Destructor.
	 */
	~TransportDG();

private:
    /// Registrar of class to factory
    static const int registrar;

	inline typename Model::ModelEqData &data() { return data_; }

	void output_vector_gather();

	void preallocate();

	/**
	 * @brief Assembles the mass matrix.
	 *
	 * The routine just calls templated method assemble_mass_matrix() for each
	 * space dimension.
	 */
	void assemble_mass_matrix();

	/**
	 * @brief Assembles the mass matrix for the given dimension.
	 */
	template<unsigned int dim>
	void assemble_mass_matrix();

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
	 * @brief Assembles the fluxes between elements of the same dimension.
	 */
	template<unsigned int dim>
	void assemble_fluxes_element_element();

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
	void calculate_velocity(const ElementFullIter &cell, std::vector<arma::vec3> &velocity, FEValuesBase<dim,3> &fv);

	/**
	 * @brief Calculates the dispersivity (diffusivity) tensor from the velocity field.
	 *
	 * @param K        The computed dispersivity tensor.
	 * @param velocity The velocity field (at quadrature points).
	 * @param Dm       Molecular diffusivities.
	 * @param alphaL   Longitudal dispersivities.
	 * @param alphaT   Transversal dispersivities.
	 * @param porosity  Porosities.
	 * @param cross_cut Cross-cuts of higher dimension.
	 */
	void calculate_dispersivity_tensor(arma::mat33 &K, const arma::vec3 &velocity,
			double Dm, double alphaL, double alphaT, double porosity,
			double cross_cut);

	/**
	 * @brief Sets up some parameters of the DG method for two sides of an edge.
	 *
	 * @param edg				The edge.
	 * @param s1				Side 1.
	 * @param s2				Side 2.
	 * @param K_size            Size of vector of tensors K.
	 * @param K1				Dispersivity tensors on side s1 (in quadrature points).
	 * @param K2				Dispersivity tensors on side s2 (in quadrature points).
	 * @param normal_vector		Normal vector to side 0 of the neighbour
	 * 							(assumed constant along the side).
	 * @param alpha1, alpha2	Penalty parameter that influences the continuity
	 * 							of the solution (large value=more continuity).
	 * @param gamma				Computed penalty parameters.
	 * @param omega				Computed weights.
	 * @param transport_flux	Computed flux from side s1 to side s2.
	 */
	void set_DG_parameters_edge(const Edge &edg,
	        const int s1,
	        const int s2,
	        const int K_size,
	        const std::vector<arma::mat33> &K1,
	        const std::vector<arma::mat33> &K2,
	        const std::vector<double> &fluxes,
	        const arma::vec3 &normal_vector,
	        const double alpha1,
	        const double alpha2,
	        double &gamma,
	        double *omega,
	        double &transport_flux);

	/**
	 * @brief Sets up parameters of the DG method on a given boundary edge.
	 *
	 * Assumption is that the edge consists of only 1 side.
	 * @param side       		The boundary side.
	 * @param K_size            Size of vector of tensors K.
	 * @param K					Dispersivity tensor.
	 * @param ad_vector         Advection vector.
	 * @param normal_vector		Normal vector (assumed constant along the edge).
	 * @param alpha				Penalty parameter that influences the continuity
	 * 							of the solution (large value=more continuity).
	 * @param gamma				Computed penalty parameters.
	 */
	void set_DG_parameters_boundary(const SideIter side,
			    const int K_size,
	            const std::vector<arma::mat33> &K,
	            const double flux,
	            const arma::vec3 &normal_vector,
	            const double alpha,
	            double &gamma);


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



	/// @name Physical parameters
	// @{

	/// Field data for model parameters.
	EqData data_;

	// @}


	/// @name Parameters of the numerical method
	// @{

	/// Finite element objects
	FEObjects *feo;

	/// Penalty parameters.
	std::vector<std::vector<double> > gamma;

	/// DG variant ((non-)symmetric/incomplete
	int dg_variant;

	/// Polynomial order of finite elements.
	unsigned int dg_order;

	// @}



	/// @name Solution of algebraic system
	// @{

	/// Vector of right hand side.
	Vec *rhs;

	/// The stiffness matrix.
	Mat *stiffness_matrix;

	/// The mass matrix.
	Mat mass_matrix;
	
	/// Mass from previous time instant (necessary when coefficients of mass matrix change in time).
	Vec *mass_vec;

	/// Linear algebra system for the transport equation.
	LinSys **ls;

	/// Linear algebra system for the time derivative (actually it is used only for handling the matrix structures).
	LinSys *ls_dt;

	// @}


	/// @name Output to file
	// @{

	/// Array for storing the output solution data.
	vector<double*> output_solution;

	/// Vector of solution data.
	vector<Vec> output_vec;

	/// Record with output specification.
	Input::Record output_rec;

	std::shared_ptr<OutputTime> output_stream;


	// @}


	/// @name Auxiliary fields used during assembly
	// @{

	/// Mass matrix coefficients.
	vector<double> mm_coef;
	/// Advection coefficients.
	vector<vector<arma::vec3> > ad_coef;
	/// Diffusion coefficients.
	vector<vector<arma::mat33> > dif_coef;
	/// Advection coefficients on edges.
	vector<vector<vector<arma::vec3> > > ad_coef_edg;
	/// Diffusion coefficients on edges.
	vector<vector<vector<arma::mat33> > > dif_coef_edg;
	/// List of indices used to call balance methods for a set of quantities.
	vector<unsigned int> subst_idx;

	// @}




	/// @name Other
	// @{

    /// Indicates whether matrices have been preallocated.
    bool allocation_done;

    // @}
};






#endif /* TRANSPORT_DG_HH_ */
