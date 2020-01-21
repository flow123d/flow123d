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
 * @file    assembly_dg_new.hh
 * @brief
 */

#ifndef ASSEMBLY_DG_NEW_HH_
#define ASSEMBLY_DG_NEW_HH_

#include "transport/transport_dg.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_values.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fields/eval_subset.hh"
#include "fields/eval_points.hh"



template <class MultidimAssembly>
class GenericAssembly
{
private:
	enum IntegralType {
        none, bulk, edge, ngh_lower_dim, ngh_higher_dim, boundary
    };

    struct IntegralData {
	    IntegralData() : elm_idx(0), side_idx(0), integral(IntegralType::none), data_size(0) {}

	    void set(unsigned int elm, unsigned int side, IntegralType itg, unsigned int size)
	    { ASSERT_DBG(itg!=IntegralType::none); elm_idx=elm; side_idx=side; integral=itg; data_size=size; }

	    void reset()
	    { integral = IntegralType::none; }

	    unsigned int elm_idx;
	    unsigned int side_idx;
	    IntegralType integral;
	    unsigned int data_size;
	};
public:
    /// Constructor
	GenericAssembly(MultidimAssembly &multidim_assembly)
    : multidim_assembly_(multidim_assembly), integral_size_(0) {
	    eval_points_ = std::make_shared<EvalPoints>();
	    std::vector<const Quadrature *> quads = { std::get<0>(multidim_assembly_)->quad_, std::get<1>(multidim_assembly_)->quad_,
	                                              std::get<2>(multidim_assembly_)->quad_, std::get<0>(multidim_assembly_)->quad_low_,
                                                  std::get<1>(multidim_assembly_)->quad_low_, std::get<2>(multidim_assembly_)->quad_low_ };
        bulk_integral_[0] = eval_points_->add_bulk<1>(*quads[0]);
        bulk_integral_[1] = eval_points_->add_bulk<2>(*quads[1]);
        bulk_integral_[2] = eval_points_->add_bulk<3>(*quads[2]);
        edge_integral_[0] = eval_points_->add_edge<1>(*quads[3]);
        edge_integral_[1] = eval_points_->add_edge<2>(*quads[4]);
        edge_integral_[2] = eval_points_->add_edge<3>(*quads[5]);
        coupling_integral_[0] = eval_points_->add_coupling<2>(*quads[4]);
        coupling_integral_[1] = eval_points_->add_coupling<3>(*quads[5]);
        //boundary_integral_
	}

	void add_compute_volume_integrals(DHCellAccessor cell) {
		unsigned int data_size = eval_points_->subset_size( cell.dim(), bulk_integral_[cell.dim()-1]->get_subset_idx() );
		integral_data_[integral_size_].set(cell.elm_idx(), 0, IntegralType::bulk, data_size);
		integral_size_++;
	}

	void add_compute_fluxes_element_element(DHCellSide edge_side) {
		unsigned int data_size = eval_points_->subset_size( edge_side.dim(), edge_integral_[edge_side.dim()-1]->get_subset_idx() ) / (edge_side.dim()+1);
		integral_data_[integral_size_].set(edge_side.elem_idx(), edge_side.side_idx(), IntegralType::edge, data_size);
		integral_size_++;
	}

	void add_compute_fluxes_element_side(DHCellAccessor cell) {
		unsigned int data_size = eval_points_->subset_size( cell.dim(), coupling_integral_[cell.dim()-1]->get_subset_low_idx() );
		integral_data_[integral_size_].set(cell.elm_idx(), 0, IntegralType::ngh_lower_dim, data_size);
		integral_size_++;
	}

	void add_compute_fluxes_element_side(DHCellSide ngh_side) {
		unsigned int data_size = eval_points_->subset_size( ngh_side.dim(), coupling_integral_[ngh_side.dim()-1]->get_subset_high_idx() ) / (ngh_side.dim()+1);
		integral_data_[integral_size_].set(ngh_side.elem_idx(), ngh_side.side_idx(), IntegralType::ngh_higher_dim, data_size);
		integral_size_++;
	}

	void insert_eval_points_from_integral_data() {
	    for (unsigned int i=0; i<integral_size_; ++i) {
	        // add data to cache if there is free space, else return
	        integral_data_[integral_size_].reset();
	    }
	    integral_size_ = 0;
	}

private:
    /// Assembly object
    MultidimAssembly &multidim_assembly_;

    std::array<std::shared_ptr<BulkIntegral>, 3> bulk_integral_;          ///< Bulk integrals of elements of dimensions 1, 2, 3
    std::array<std::shared_ptr<EdgeIntegral>, 3> edge_integral_;          ///< Edge integrals between elements of dimensions 1, 2, 3
    std::array<std::shared_ptr<CouplingIntegral>, 2> coupling_integral_;  ///< Coupling integrals between elements of dimensions 1-2, 2-3
    std::array<std::shared_ptr<BoundaryIntegral>, 3> boundary_integral_;  ///< Boundary integrals betwwen elements of dimensions 1, 2, 3 and boundaries
    std::shared_ptr<EvalPoints> eval_points_;                             ///< EvalPoints object shared by all integrals

    std::array<IntegralData, 22> integral_data_;
    unsigned int integral_size_;
};


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim, class Model>
class AssemblyDGNew
{
public:
    typedef typename TransportDG<Model>::EqData EqDataDG;

    /// Constructor.
    AssemblyDGNew(std::shared_ptr<EqDataDG> data, TransportDG<Model> &model)
    : fe_(make_shared< FE_P_disc<dim> >(data->dg_order)), fe_low_(make_shared< FE_P_disc<dim-1> >(data->dg_order)),
      fe_rt_(new FE_RT0<dim>), fe_rt_low_(new FE_RT0<dim-1>),
      quad_(new QGauss(dim, 2*data->dg_order)),
	  quad_low_(new QGauss(dim-1, 2*data->dg_order)),
      model_(model), data_(data), fv_rt_(*quad_, *fe_rt_, update_values | update_gradients | update_quadrature_points),
      fe_values_(*quad_, *fe_, update_values | update_gradients | update_JxW_values | update_quadrature_points),
      fv_rt_vb_(nullptr), fe_values_vb_(nullptr),
      fe_values_side_(*quad_low_, *fe_, update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points),
      fsv_rt_(*quad_low_, *fe_rt_, update_values | update_quadrature_points) {

        if (dim>1) {
            fv_rt_vb_ = new FEValues<dim-1,3>(*quad_low_, *fe_rt_low_, update_values | update_quadrature_points);
            fe_values_vb_ = new FEValues<dim-1,3>(*quad_low_, *fe_low_,
                    update_values | update_gradients | update_JxW_values | update_quadrature_points);
        }
        ndofs_ = fe_->n_dofs();
        qsize_ = quad_->size();
        qsize_lower_dim_ = quad_low_->size();
        dof_indices_.resize(ndofs_);
        loc_dof_indices_.resize(ndofs_);
        side_dof_indices_vb_.resize(2*ndofs_);
    }

    /// Destructor.
    ~AssemblyDGNew() {
        delete fe_rt_;
        delete fe_rt_low_;
        delete quad_;
        delete quad_low_;
        if (fv_rt_vb_!=nullptr) delete fv_rt_vb_;
        if (fe_values_vb_!=nullptr) delete fe_values_vb_;

        for (unsigned int i=0; i<data_->ad_coef_edg.size(); i++)
        {
            delete fe_values_vec_[i];
        }
    }

    /// Initialize auxiliary vectors and other data members
    void initialize() {
        local_matrix_.resize(4*ndofs_*ndofs_);
        local_retardation_balance_vector_.resize(ndofs_);
        local_mass_balance_vector_.resize(ndofs_);
        local_rhs_.resize(ndofs_);
        local_source_balance_vector_.resize(ndofs_);
        local_source_balance_rhs_.resize(ndofs_);
        local_flux_balance_vector_.resize(ndofs_);
        velocity_.resize(qsize_);
        side_velocity_vec_.resize(data_->ad_coef_edg.size());
        sources_conc_.resize(model_.n_substances(), std::vector<double>(qsize_));
        sources_density_.resize(model_.n_substances(), std::vector<double>(qsize_));
        sources_sigma_.resize(model_.n_substances(), std::vector<double>(qsize_));
        sigma_.resize(qsize_lower_dim_);
        csection_.resize(qsize_lower_dim_);
        csection_higher_.resize(qsize_lower_dim_);
        dg_penalty_.resize(data_->ad_coef_edg.size());
        bc_values_.resize(qsize_lower_dim_);
        bc_fluxes_.resize(qsize_lower_dim_);
        bc_ref_values_.resize(qsize_lower_dim_);
        init_values_.resize(model_.n_substances(), std::vector<double>(qsize_));

        mm_coef_.resize(qsize_);
        ret_coef_.resize(model_.n_substances());
        for (unsigned int sbi=0; sbi<model_.n_substances(); sbi++)
        {
            ret_coef_[sbi].resize(qsize_);
        }

        for (unsigned int sid=0; sid<data_->ad_coef_edg.size(); sid++)
        {
            side_dof_indices_.push_back( vector<LongIdx>(ndofs_) );
            fe_values_vec_.push_back(new FESideValues<dim,3>(*quad_low_, *fe_,
                    update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points));
        }

        // index 0 = element with lower dimension,
        // index 1 = side of element with higher dimension
        fv_sb_.resize(2);
        fv_sb_[0] = fe_values_vb_;
        fv_sb_[1] = &fe_values_side_;
    }

//private:
	/**
	 * @brief Calculates the velocity field on a given cell.
	 *
	 * @param cell       The cell.
	 * @param velocity   The computed velocity field (at quadrature points).
	 * @param point_list The quadrature points.
	 */
    void calculate_velocity(const ElementAccessor<3> &cell, vector<arma::vec3> &velocity,
                            const std::vector<arma::vec::fixed<3>> &point_list)
    {
        velocity.resize(point_list.size());
        model_.velocity_field_ptr()->value_list(point_list, cell, velocity);
    }


    shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.
    shared_ptr<FiniteElement<dim-1>> fe_low_;   ///< Finite element for the solution of the advection-diffusion equation (dim-1).
    FiniteElement<dim> *fe_rt_;                 ///< Finite element for the water velocity field.
    FiniteElement<dim-1> *fe_rt_low_;           ///< Finite element for the water velocity field (dim-1).
    Quadrature *quad_;                     ///< Quadrature used in assembling methods.
    Quadrature *quad_low_;               ///< Quadrature used in assembling methods (dim-1).

    /// Reference to model (we must use common ancestor of concentration and heat model)
    TransportDG<Model> &model_;

    /// Data object shared with TransportDG
    std::shared_ptr<EqDataDG> data_;

    unsigned int ndofs_;                                      ///< Number of dofs
    unsigned int qsize_;                                      ///< Size of quadrature of actual dim
    unsigned int qsize_lower_dim_;                            ///< Size of quadrature of dim-1
    FEValues<dim,3> fv_rt_;                                   ///< FEValues of object (of RT0 finite element type)
    FEValues<dim,3> fe_values_;                               ///< FEValues of object (of P disc finite element type)
    FEValues<dim-1,3> *fv_rt_vb_;                             ///< FEValues of dim-1 object (of RT0 finite element type)
    FEValues<dim-1,3> *fe_values_vb_;                         ///< FEValues of dim-1 object (of P disc finite element type)
    FESideValues<dim,3> fe_values_side_;                      ///< FESideValues of object (of P disc finite element type)
    FESideValues<dim,3> fsv_rt_;                              ///< FESideValues of object (of RT0 finite element type)
    vector<FESideValues<dim,3>*> fe_values_vec_;              ///< Vector of FESideValues of object (of P disc finite element types)
    vector<FEValuesSpaceBase<3>*> fv_sb_;                     ///< Auxiliary vector, holds FEValues objects for assemble element-side

    vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
    vector<LongIdx> loc_dof_indices_;                         ///< Vector of local DOF indices
    vector< vector<LongIdx> > side_dof_indices_;              ///< Vector of vectors of side DOF indices
    vector<LongIdx> side_dof_indices_vb_;                     ///< Vector of side DOF indices (assemble element-side fluxex)
    vector<PetscScalar> local_matrix_;                        ///< Auxiliary vector for assemble methods
    vector<PetscScalar> local_retardation_balance_vector_;    ///< Auxiliary vector for assemble mass matrix.
    vector<PetscScalar> local_mass_balance_vector_;           ///< Same as previous.
    vector<PetscScalar> local_rhs_;                           ///< Auxiliary vector for set_sources method.
    vector<PetscScalar> local_source_balance_vector_;         ///< Auxiliary vector for set_sources method.
    vector<PetscScalar> local_source_balance_rhs_;            ///< Auxiliary vector for set_sources method.
    vector<PetscScalar> local_flux_balance_vector_;           ///< Auxiliary vector for set_boundary_conditions method.
    PetscScalar local_flux_balance_rhs_;                      ///< Auxiliary variable for set_boundary_conditions method.
    vector<arma::vec3> velocity_;                             ///< Auxiliary results.
    vector<arma::vec3> velocity_higher_;                      ///< Velocity results of higher dim element (element-side computation).
    vector<vector<arma::vec3> > side_velocity_vec_;           ///< Vector of velocities results.
    vector<vector<double> > sources_conc_;                    ///< Auxiliary vectors for set_sources method.
    vector<vector<double> > sources_density_;                 ///< Auxiliary vectors for set_sources method.
    vector<vector<double> > sources_sigma_;                   ///< Auxiliary vectors for assemble volume integrals and set_sources method.
    vector<double> sigma_;                                    ///< Auxiliary vector for assemble boundary fluxes (robin sigma), element-side fluxes (frac sigma) and set boundary conditions method
    vector<double> csection_;                                 ///< Auxiliary vector for assemble boundary fluxes, element-side fluxes and set boundary conditions
    vector<double> csection_higher_;                          ///< Auxiliary vector for assemble element-side fluxes
    vector<vector<double> > dg_penalty_;                      ///< Auxiliary vectors for assemble element-element fluxes
    vector<double> bc_values_;                                ///< Auxiliary vector for set boundary conditions method
    vector<double> bc_fluxes_;                                ///< Same as previous
    vector<double> bc_ref_values_;                            ///< Same as previous
    std::vector<std::vector<double> > init_values_;           ///< Auxiliary vectors for prepare initial condition

	/// Mass matrix coefficients.
	vector<double> mm_coef_;
	/// Retardation coefficient due to sorption.
	vector<vector<double> > ret_coef_;

	/// @name Auxiliary variables used during element-element assembly
	// @{

	double gamma_l, omega[2], transport_flux, delta[2], delta_sum;
    double aniso1, aniso2;
    int sid, s1, s2;

	// @}

	/// @name Auxiliary variables used during element-side assembly
	// @{

    unsigned int n_dofs[2], n_indices;
    double comm_flux[2][2];

	// @}

	/// @name Auxiliary variables used during set sources
	// @{

    double source;

	// @}

    friend class TransportDG<Model>;

};



#endif /* ASSEMBLY_DG_NEW_HH_ */
