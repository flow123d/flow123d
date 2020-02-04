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
#include "fields/field_value_cache.hh"



/// Allow set mask of active integrals.
enum ActiveIntegrals {
    none     =      0,
    bulk     = 0x0001,
    edge     = 0x0002,
    coupling = 0x0004,
    boundary = 0x0008
};


template < template<Dim...> class DimAssembly>
class GenericAssembly
{
private:
    struct BulkIntegralData {
    	BulkIntegralData() : data_size(0) {}

	    void reset()
	    { data_size = 0; }

	    DHCellAccessor cell;
	    unsigned int data_size;
	    unsigned int subset_index;
	};

    struct EdgeIntegralData {
    	EdgeIntegralData() : data_size(0) {}

	    void reset()
	    { data_size = 0; }

	    DHCellSide side;
	    unsigned int data_size;
	    unsigned int start_point;
	    unsigned int subset_index;
	};

    struct CouplingIntegralData {
       	CouplingIntegralData() : bulk_data_size(0), side_data_size(0) {}

        void reset()
        { bulk_data_size = 0; side_data_size=0; }

        DHCellAccessor cell;
        unsigned int bulk_data_size;
	    unsigned int bulk_subset_index;

        DHCellSide side;
        unsigned int side_data_size;
	    unsigned int side_start_point;
	    unsigned int side_subset_index;
    };

public:

    /// Constructor
    GenericAssembly(std::shared_ptr<DimAssembly<0>> assembly0, std::shared_ptr<DimAssembly<1>> assembly1,
                    std::shared_ptr<DimAssembly<2>> assembly2, std::shared_ptr<DimAssembly<3>> assembly3 )
    : multidim_assembly_(assembly0, assembly1, assembly2, assembly3),
      active_integrals_(ActiveIntegrals::none), integrals_size_({0, 0, 0, 0})
    {
        eval_points_ = std::make_shared<EvalPoints>();
        // first step - create integrals, then - initialize cache
        multidim_assembly_.get<1>()->create_integrals(eval_points_);
        multidim_assembly_.get<2>()->create_integrals(eval_points_);
        multidim_assembly_.get<3>()->create_integrals(eval_points_);
        multidim_assembly_.get<1>()->data_->element_cache_map_.init(eval_points_);
    }

	inline void set_active(int active) {
	    active_integrals_ = active;
	}

    /// Call initialize method of inner AssemblyDGNew objects.
    void initialize() {
        multidim_assembly_.get<1>()->initialize();
        multidim_assembly_.get<2>()->initialize();
        multidim_assembly_.get<3>()->initialize();
    }

    void assemble_stiffness_matrix(std::shared_ptr<DOFHandlerMultiDim> dh) {
        START_TIMER("assemble_stiffness");
        ElementCacheMap &el_cache_map = multidim_assembly_.get<1>()->data_->element_cache_map_;
        for (auto cell : dh->local_range() )
        {
            // generic_assembly.check_integral_data();
            if (active_integrals_ & ActiveIntegrals::bulk)
        	    if (cell.is_own()) { // Not ghost
                    this->add_volume_integral(cell);
                    el_cache_map.add(cell);
        	    }

            for( DHCellSide cell_side : cell.side_range() ) {
                if (active_integrals_ & ActiveIntegrals::boundary)
                    if (cell.is_own()) // Not ghost
                        if ( (cell_side.side().edge()->n_sides == 1) && (cell_side.side().cond() != NULL) ) {
                            this->add_boundary_integral(cell_side);
                            el_cache_map.add(cell_side);
                            continue;
                        }
                if (active_integrals_ & ActiveIntegrals::edge)
                    if ( (cell_side.n_edge_sides() >= 2) && (cell_side.edge_sides().begin()->element().idx() == cell.elm_idx()))
                        for( DHCellSide edge_side : cell_side.edge_sides() ) {
                            this->add_edge_integral(edge_side);
                            el_cache_map.add(edge_side);
                        }
            }

	        if (active_integrals_ & ActiveIntegrals::coupling)
                for( DHCellSide neighb_side : cell.neighb_sides() ) { // cell -> elm lower dim, neighb_side -> elm higher dim
                    if (cell.dim() != neighb_side.dim()-1) continue;
                    this->add_compute_fluxes_element_side(cell, neighb_side);
                    el_cache_map.add(cell);
                    el_cache_map.add(neighb_side);
                }

            this->insert_eval_points_from_integral_data(el_cache_map);
            multidim_assembly_.get<1>()->data_->cache_update(el_cache_map, dh->mesh());
        }
        END_TIMER("assemble_stiffness");
    }

    void add_volume_integral(const DHCellAccessor &cell) {
        bulk_integral_data_[ integrals_size_[0] ].cell = cell;
        switch (cell.dim()) {
        case 1:
        	bulk_integral_data_[ integrals_size_[0] ].subset_index = std::get<1>(multidim_assembly_)->bulk_integral_->get_subset_idx();
            break;
        case 2:
            bulk_integral_data_[ integrals_size_[0] ].subset_index = std::get<2>(multidim_assembly_)->bulk_integral_->get_subset_idx();
            break;
        case 3:
            bulk_integral_data_[ integrals_size_[0] ].subset_index = std::get<3>(multidim_assembly_)->bulk_integral_->get_subset_idx();
            break;
    	}
        bulk_integral_data_[ integrals_size_[0] ].data_size = eval_points_->subset_size( cell.dim(), bulk_integral_data_[ integrals_size_[0] ].subset_index );
        integrals_size_[0]++;
    }

    void add_edge_integral(const DHCellSide &edge_side) {
    	edge_integral_data_[ integrals_size_[1] ].side = edge_side;
        switch (edge_side.dim()) {
        case 1:
            edge_integral_data_[ integrals_size_[1] ].subset_index = std::get<1>(multidim_assembly_)->edge_integral_->get_subset_idx();
            break;
        case 2:
            edge_integral_data_[ integrals_size_[1] ].subset_index = std::get<2>(multidim_assembly_)->edge_integral_->get_subset_idx();
            break;
        case 3:
            edge_integral_data_[ integrals_size_[1] ].subset_index = std::get<3>(multidim_assembly_)->edge_integral_->get_subset_idx();
            break;
    	}
        edge_integral_data_[ integrals_size_[1] ].data_size =
                eval_points_->subset_size( edge_side.dim(), edge_integral_data_[ integrals_size_[1] ].subset_index ) / (edge_side.dim() +1);
        edge_integral_data_[ integrals_size_[1] ].start_point = edge_integral_data_[ integrals_size_[1] ].data_size * edge_side.side_idx();
        integrals_size_[1]++;
    }

    void add_compute_fluxes_element_side(const DHCellAccessor &cell, const DHCellSide &ngh_side) {
    	coupling_integral_data_[ integrals_size_[2] ].cell = cell;
    	coupling_integral_data_[ integrals_size_[2] ].side = ngh_side;
        switch (cell.dim()) {
        case 1:
        	coupling_integral_data_[ integrals_size_[2] ].bulk_subset_index = std::get<2>(multidim_assembly_)->coupling_integral_->get_subset_low_idx();
        	coupling_integral_data_[ integrals_size_[2] ].side_subset_index = std::get<2>(multidim_assembly_)->coupling_integral_->get_subset_high_idx();
            break;
        case 2:
            coupling_integral_data_[ integrals_size_[2] ].bulk_subset_index = std::get<3>(multidim_assembly_)->coupling_integral_->get_subset_low_idx();
        	coupling_integral_data_[ integrals_size_[2] ].side_subset_index = std::get<3>(multidim_assembly_)->coupling_integral_->get_subset_high_idx();
            break;
    	}
        coupling_integral_data_[ integrals_size_[2] ].bulk_data_size =
        		eval_points_->subset_size( cell.dim(), coupling_integral_data_[ integrals_size_[2] ].bulk_subset_index );
        coupling_integral_data_[ integrals_size_[2] ].side_data_size =
        		eval_points_->subset_size( ngh_side.dim(), coupling_integral_data_[ integrals_size_[2] ].side_subset_index );
        coupling_integral_data_[ integrals_size_[2] ].side_start_point = coupling_integral_data_[ integrals_size_[2] ].side_data_size * ngh_side.side_idx();
        integrals_size_[2]++;
    }

    void add_boundary_integral(const DHCellSide &bdr_side) {
    	boundary_integral_data_[ integrals_size_[3] ].side = bdr_side;
        switch (bdr_side.dim()) {
        case 1:
        	boundary_integral_data_[ integrals_size_[3] ].subset_index = std::get<1>(multidim_assembly_)->edge_integral_->get_subset_idx();
            break;
        case 2:
        	boundary_integral_data_[ integrals_size_[3] ].subset_index = std::get<2>(multidim_assembly_)->edge_integral_->get_subset_idx();
            break;
        case 3:
        	boundary_integral_data_[ integrals_size_[3] ].subset_index = std::get<3>(multidim_assembly_)->edge_integral_->get_subset_idx();
            break;
    	}
        boundary_integral_data_[ integrals_size_[3] ].data_size =
                eval_points_->subset_size( bdr_side.dim(), boundary_integral_data_[ integrals_size_[3] ].subset_index ) / (bdr_side.dim() +1);
        boundary_integral_data_[ integrals_size_[3] ].start_point = boundary_integral_data_[ integrals_size_[3] ].data_size * bdr_side.side_idx();
        integrals_size_[3]++;
    }

    void insert_eval_points_from_integral_data(ElementCacheMap &el_cache_map) {
        for (unsigned int i=0; i<integrals_size_[0]; ++i) {
            // add data to cache if there is free space, else return
            el_cache_map.mark_used_eval_points(bulk_integral_data_[i].cell, bulk_integral_data_[i].subset_index, bulk_integral_data_[i].data_size);
            bulk_integral_data_[i].reset();
        }
        integrals_size_[0] = 0;
        for (unsigned int i=0; i<integrals_size_[1]; ++i) {
            // add data to cache if there is free space, else return
            el_cache_map.mark_used_eval_points(edge_integral_data_[i].side.cell(), edge_integral_data_[i].subset_index, edge_integral_data_[i].data_size, edge_integral_data_[i].start_point);
            edge_integral_data_[i].reset();
        }
        integrals_size_[1] = 0;
        for (unsigned int i=0; i<integrals_size_[2]; ++i) {
            // add data to cache if there is free space, else return
            el_cache_map.mark_used_eval_points(coupling_integral_data_[i].cell, coupling_integral_data_[i].bulk_subset_index, coupling_integral_data_[i].bulk_data_size);
            el_cache_map.mark_used_eval_points(coupling_integral_data_[i].side.cell(), coupling_integral_data_[i].side_subset_index, coupling_integral_data_[i].side_data_size, coupling_integral_data_[i].side_start_point);
            coupling_integral_data_[i].reset();
        }
        integrals_size_[2] = 0;
        for (unsigned int i=0; i<integrals_size_[3]; ++i) {
            // add data to cache if there is free space, else return
            el_cache_map.mark_used_eval_points(boundary_integral_data_[i].side.cell(), boundary_integral_data_[i].subset_index, boundary_integral_data_[i].data_size, boundary_integral_data_[i].start_point);
            boundary_integral_data_[i].reset();
        }
        integrals_size_[3] = 0;
    }

private:
    /// Assembly object
    MixedPtr<DimAssembly> multidim_assembly_;

    /// Holds mask of active integrals.
    int active_integrals_;

    std::shared_ptr<EvalPoints> eval_points_;                     ///< EvalPoints object shared by all integrals

    // Following variables hold data of all integrals depending of actual computed element.
    std::array<BulkIntegralData, 1>     bulk_integral_data_;      ///< Holds data for computing bulk integrals.
    std::array<EdgeIntegralData, 18>    edge_integral_data_;      ///< Holds data for computing edge integrals.
    std::array<CouplingIntegralData, 6> coupling_integral_data_;  ///< Holds data for computing couplings integrals.
    std::array<EdgeIntegralData, 4>     boundary_integral_data_;  ///< Holds data for computing boundary integrals.
    std::array<unsigned int, 4>         integrals_size_;          ///< Holds used sizes of previous integral data types
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

    void create_integrals(std::shared_ptr<EvalPoints> eval_points) {
        bulk_integral_ = eval_points->add_bulk<dim>(*quad_);
        edge_integral_ = eval_points->add_edge<dim>(*quad_low_);
        if (dim>1) coupling_integral_ = eval_points->add_coupling<dim>(*quad_low_);
        boundary_integral_ = eval_points->add_boundary<dim>(*quad_low_);
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
    /*void calculate_velocity(const ElementAccessor<3> &cell, vector<arma::vec3> &velocity,
                            const std::vector<arma::vec::fixed<3>> &point_list)
    {
        velocity.resize(point_list.size());
        model_.velocity_field_ptr()->value_list(point_list, cell, velocity);
    }*/

    std::shared_ptr<BulkIntegral> bulk_integral_;          ///< Bulk integrals of elements of given dimension
    std::shared_ptr<EdgeIntegral> edge_integral_;          ///< Edge integrals between sides of elements of given dimension
    std::shared_ptr<CouplingIntegral> coupling_integral_;  ///< Coupling integrals between elements of given dimension and sides of elements of dim+1 dimension
    std::shared_ptr<BoundaryIntegral> boundary_integral_;  ///< Boundary integrals betwwen sides of elements of given dimension and mesh boundary

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

/// Template specialization of dim=0
template <class Model>
class AssemblyDGNew<0, Model>
{
public:
    typedef typename TransportDG<Model>::EqData EqDataDG;

    /// Constructor.
    AssemblyDGNew(std::shared_ptr<EqDataDG> data, TransportDG<Model> &model) {}

    /// Destructor.
    ~AssemblyDGNew() {}

    void initialize() {}

    friend class TransportDG<Model>;

};



#endif /* ASSEMBLY_DG_NEW_HH_ */
