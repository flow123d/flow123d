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
 * @file    transport_dg.hh
 * @brief   Discontinuous Galerkin method for equation of transport with dispersion.
 * @author  Jan Stebel
 */

#ifndef TRANSPORT_DG_HH_
#define TRANSPORT_DG_HH_

#include <math.h>                              // for fabs
#include <string.h>                            // for memcpy
#include <algorithm>                           // for max
#include <boost/exception/info.hpp>            // for operator<<, error_info...
#include <string>                              // for operator<<
#include <vector>                              // for vector
#include <armadillo>
#include "fem/update_flags.hh"                 // for operator|
#include "fields/field_values.hh"              // for FieldValue<>::Scalar
#include "fields/field.hh"
#include "fields/multi_field.hh"
#include "la/vector_mpi.hh"
#include "fields/equation_output.hh"
#include "la/linsys.hh"
#include "input/accessors.hh"                  // for ExcAccessorForNullStorage
#include "input/accessors_impl.hh"             // for Record::val
#include "input/storage.hh"                    // for ExcStorageTypeMismatch
#include "input/type_base.hh"                  // for Array
#include "input/type_generic.hh"               // for Instance
#include "input/type_record.hh"                // for Record::ExcRecordKeyNo...
#include "system/index_types.hh"               // for LongIdx
#include "mesh/accessors.hh"                   // for ElementAccessor, SIdeIter
#include "mesh/elements.h"                     // for Element::dim, Element:...
#include "mesh/neighbours.h"                   // for Neighbour::element
#include "mpi.h"                               // for MPI_Comm_rank
#include "petscmat.h"                          // for Mat, MatDestroy
#include "petscvec.h"                          // for Vec, VecDestroy, VecSc...
#include "transport/concentration_model.hh"    // for ConcentrationTransport...
#include "transport/heat_model.hh"             // for HeatTransferModel, Hea...
#include "tools/mixed.hh"
class DiscreteSpace;
class Distribution;
class OutputTime;
class DOFHandlerMultiDim;
template<unsigned int dim, class Model> class AssemblyDG;
template<unsigned int dim, class Model> class MassAssemblyDG;
template<unsigned int dim, class Model> class StiffnessAssemblyDG;
template<unsigned int dim, class Model> class SourcesAssemblyDG;
template<unsigned int dim, class Model> class BdrConditionAssemblyDG;
template<unsigned int dim, class Model> class InitConditionAssemblyDG;
template< template<Dim...> class DimAssembly> class GenericAssembly;
template<unsigned int dim, unsigned int spacedim> class FEValuesBase;
template<unsigned int dim> class FiniteElement;
template<unsigned int dim, unsigned int spacedim> class Mapping;
class Quadrature;
namespace Input { namespace Type { class Selection; } }
class ElementCacheMap;

/*class FEObjects {
public:

	inline FEObjects(Mesh *mesh_, unsigned int fe_order)
	{
        fe = MixedPtr<FE_P_disc>(fe_order);
        fe_rt = MixedPtr<FE_RT0>();
        q = MixedPtr<QGauss>(2*fe_order);

        auto ds = std::make_shared<EqualOrderDiscreteSpace>(mesh_, fe);
        dh = std::make_shared<DOFHandlerMultiDim>(*mesh_);
        dh->distribute_dofs(ds_);
    }

	~FEObjects();

	MixedPtr<FiniteElement> fe;
	MixedPtr<FiniteElement> fe_rt;
	MixedPtr<Quadrature> q;

	/// Object for distribution of dofs.
	std::shared_ptr<DOFHandlerMultiDim> dh;
};*/



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
class TransportDG : public Model
{
public:

    template<unsigned int dim> using MassAssemblyDim = MassAssemblyDG<dim, Model>;
    template<unsigned int dim> using StiffnessAssemblyDim = StiffnessAssemblyDG<dim, Model>;
    template<unsigned int dim> using SourcesAssemblyDim = SourcesAssemblyDG<dim, Model>;
    template<unsigned int dim> using BdrConditionAssemblyDim = BdrConditionAssemblyDG<dim, Model>;
    template<unsigned int dim> using InitConditionAssemblyDim = InitConditionAssemblyDG<dim, Model>;

	class EqData : public Model::ModelEqData {
	public:

		EqData();

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
		void set_DG_parameters_boundary(Side side,
				    const int K_size,
		            const std::vector<arma::mat33> &K,
		            const double flux,
		            const arma::vec3 &normal_vector,
		            const double alpha,
		            double &gamma);


		/// Compute and return anisotropy of given element
		double elem_anisotropy(ElementAccessor<3> e) const;



		MultiField<3, FieldValue<3>::Scalar> fracture_sigma;    ///< Transition parameter for diffusive transfer on fractures (for each substance).
		MultiField<3, FieldValue<3>::Scalar> dg_penalty;        ///< Penalty enforcing inter-element continuity of solution (for each substance).
        Field<3, FieldValue<3>::Scalar> region_id;
        Field<3, FieldValue<3>::Scalar> subdomain;

        EquationOutput output_fields;


    	/// @name Parameters of the numerical method
    	// @{

    	/// Penalty parameters.
    	std::vector<std::vector<double> > gamma;

    	/// DG variant ((non-)symmetric/incomplete
    	int dg_variant;

    	/// Polynomial order of finite elements.
    	unsigned int dg_order;

    	// @}


        /// Auxiliary vectors for calculation of sources in balance due to retardation (e.g. sorption).
    	std::vector<Vec> ret_vec;

    	/// Linear algebra system for the transport equation.
    	LinSys **ls;

    	/// Linear algebra system for the time derivative (actually it is used only for handling the matrix structures).
    	LinSys **ls_dt;

    	/// @name Auxiliary fields used during assembly
    	// @{

    	/// Advection coefficients.
    	vector<vector<arma::vec3> > ad_coef;
    	/// Diffusion coefficients.
    	vector<vector<arma::mat33> > dif_coef;
    	/// Advection coefficients on edges.
    	vector<vector<vector<arma::vec3> > > ad_coef_edg;
    	/// Diffusion coefficients on edges.
    	vector<vector<vector<arma::mat33> > > dif_coef_edg;

    	// @}

        /// Object for distribution of dofs.
        std::shared_ptr<DOFHandlerMultiDim> dh_;

        /// general assembly objects, hold assembly objects of appropriate dimension
        GenericAssembly< MassAssemblyDim > * mass_assembly_;
        GenericAssembly< StiffnessAssemblyDim > * stiffness_assembly_;
        GenericAssembly< SourcesAssemblyDim > * sources_assembly_;
        GenericAssembly< BdrConditionAssemblyDim > * bdr_cond_assembly_;
        GenericAssembly< InitConditionAssemblyDim > * init_cond_assembly_;
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
    TransportDG(Mesh &init_mesh, const Input::Record in_rec);
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
	
    bool evaluate_time_constraint(double &)
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
	~TransportDG();

	void initialize() override;

    void calculate_cumulative_balance();

	const Vec &get_solution(unsigned int sbi)
	{ return data_->ls[sbi]->get_solution(); }

	double **get_concentration_matrix()
	{ return solution_elem_; }

	void calculate_concentration_matrix();

	void update_after_reactions(bool solution_changed);

    void get_par_info(LongIdx * &el_4_loc, Distribution * &el_ds);

    LongIdx *get_row_4_el();

    /// Access to balance object of Model
    inline std::shared_ptr<Balance> balance() const {
        return Model::balance_;
    }

    /// Return vector of substances indices
    inline const vector<unsigned int> subst_idx() const {
        return Model::subst_idx;
    }




private:
    /// Registrar of class to factory
    static const int registrar;

	inline typename Model::ModelEqData &data() { return *data_; }

	void preallocate();

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
// 	void calculate_dispersivity_tensor(arma::mat33 &K, const arma::vec3 &velocity,
// 			double Dm, double alphaL, double alphaT, double porosity,
// 			double cross_cut);

	/**
	 * @brief Sets the initial condition.
	 *
	 * This method just calls AssemblyDG::prepare_initial_condition() for each elements.
	 */
	void set_initial_condition();


	/// @name Physical parameters
	// @{

	/// Field data for model parameters.
	std::shared_ptr<EqData> data_;

	// @}



	/// @name Solution of algebraic system
	// @{

	/// Vector of right hand side.
	std::vector<Vec> rhs;

	/// The stiffness matrix.
	std::vector<Mat> stiffness_matrix;

	/// The mass matrix.
	std::vector<Mat> mass_matrix;
	
	/// Mass from previous time instant (necessary when coefficients of mass matrix change in time).
	std::vector<Vec> mass_vec;
    
	/// Element averages of solution (the array is passed to reactions in operator splitting).
	double **solution_elem_;

	// @}


	/// @name Output to file
	// @{

	/// Array for storing the output solution data.
	//vector<double*> output_solution;

	/// Vector of solution data.
	vector<VectorMPI> output_vec;

	/// Record with input specification.
	Input::Record input_rec;


	// @}


	/// @name Auxiliary fields used during assembly
	// @{

	/// Temporary values of increments due to retardation (e.g. sorption)
    vector<double> ret_sources, ret_sources_prev;

	// @}




	/// @name Other
	// @{

    /// Indicates whether matrices have been preallocated.
    bool allocation_done;

    // @}

};






#endif /* TRANSPORT_DG_HH_ */
