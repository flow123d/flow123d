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
 * @file    assembly_dg.hh
 * @brief
 */

#ifndef ASSEMBLY_DG_HH_
#define ASSEMBLY_DG_HH_

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


/// Set of all used integral necessary in assemblation
struct AssemblyIntegrals {
    std::array<std::shared_ptr<BulkIntegral>, 3> bulk_;          ///< Bulk integrals of elements of dimensions 1, 2, 3
    std::array<std::shared_ptr<EdgeIntegral>, 3> edge_;          ///< Edge integrals between elements of dimensions 1, 2, 3
    std::array<std::shared_ptr<CouplingIntegral>, 2> coupling_;  ///< Coupling integrals between elements of dimensions 1-2, 2-3
    std::array<std::shared_ptr<BoundaryIntegral>, 3> boundary_;  ///< Boundary integrals betwwen elements of dimensions 1, 2, 3 and boundaries
};



/**
 * @brief Generic class of assemblation.
 *
 * Class
 *  - holds assemblation structures (EvalPoints, Integral objects, Integral data table).
 *  - associates assemblation objects specified by dimension
 *  - provides general assemble method
 *  - provides methods that allow construction of element patches
 */
template < template<IntDim...> class DimAssembly>
class GenericAssembly
{
private:
    struct BulkIntegralData {
        BulkIntegralData() {}

        DHCellAccessor cell;
        unsigned int subset_index;
    };

    struct EdgeIntegralData {
    	EdgeIntegralData()
    	: edge_side_range(make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide() ), make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide() )) {}

    	RangeConvert<DHEdgeSide, DHCellSide> edge_side_range;
        unsigned int subset_index;
	};

    struct CouplingIntegralData {
       	CouplingIntegralData() {}

        DHCellAccessor cell;
	    unsigned int bulk_subset_index;
        DHCellSide side;
	    unsigned int side_subset_index;
    };

    struct BoundaryIntegralData {
    	BoundaryIntegralData() {}

	    DHCellSide side;
	    unsigned int subset_index;
	};

public:

    /// Constructor
    GenericAssembly( typename DimAssembly<1>::EqDataDG *eq_data, int active_integrals )
    : multidim_assembly_(eq_data),
      active_integrals_(active_integrals), integrals_size_({0, 0, 0, 0})
    {
        eval_points_ = std::make_shared<EvalPoints>();
        // first step - create integrals, then - initialize cache
        multidim_assembly_[1_d]->create_integrals(eval_points_, integrals_, active_integrals_);
        multidim_assembly_[2_d]->create_integrals(eval_points_, integrals_, active_integrals_);
        multidim_assembly_[3_d]->create_integrals(eval_points_, integrals_, active_integrals_);
        element_cache_map_.init(eval_points_);
    }

    inline MixedPtr<DimAssembly, 1> multidim_assembly() const {
        return multidim_assembly_;
    }

    inline std::shared_ptr<EvalPoints> eval_points() const {
        return eval_points_;
    }

	/**
	 * @brief General assemble methods.
	 *
	 * Loops through local cells and calls assemble methods of assembly
	 * object of each cells over space dimension.
	 */
    void assemble(std::shared_ptr<DOFHandlerMultiDim> dh) {
        unsigned int i;
        multidim_assembly_[1_d]->begin();
        for (auto cell : dh->local_range() )
        {
            this->add_integrals_of_computing_step(cell);

            if ( cell.is_own() && (active_integrals_ & ActiveIntegrals::bulk) ) {
                START_TIMER("assemble_volume_integrals");
                for (i=0; i<integrals_size_[0]; ++i) { // volume integral
                    switch (bulk_integral_data_[i].cell.dim()) {
                    case 1:
                        multidim_assembly_[1_d]->assemble_volume_integrals(bulk_integral_data_[i].cell);
                        break;
                    case 2:
                        multidim_assembly_[2_d]->assemble_volume_integrals(bulk_integral_data_[i].cell);
                        break;
                    case 3:
                        multidim_assembly_[3_d]->assemble_volume_integrals(bulk_integral_data_[i].cell);
                        break;
                    }
                }
                END_TIMER("assemble_volume_integrals");
            }

            if ( cell.is_own() && (active_integrals_ & ActiveIntegrals::boundary) ) {
                START_TIMER("assemble_fluxes_boundary");
                for (i=0; i<integrals_size_[3]; ++i) { // boundary integral
                    switch (boundary_integral_data_[i].side.dim()) {
                    case 1:
                        multidim_assembly_[1_d]->assemble_fluxes_boundary(boundary_integral_data_[i].side);
                        break;
                    case 2:
                        multidim_assembly_[2_d]->assemble_fluxes_boundary(boundary_integral_data_[i].side);
                        break;
                    case 3:
                        multidim_assembly_[3_d]->assemble_fluxes_boundary(boundary_integral_data_[i].side);
                        break;
                    }
                }
                END_TIMER("assemble_fluxes_boundary");
            }

            if (active_integrals_ & ActiveIntegrals::edge) {
                START_TIMER("assemble_fluxes_elem_elem");
                for (i=0; i<integrals_size_[1]; ++i) { // edge integral
                    switch (edge_integral_data_[i].edge_side_range.begin()->dim()) {
                    case 1:
                        multidim_assembly_[1_d]->assemble_fluxes_element_element(edge_integral_data_[i].edge_side_range);
                        break;
                    case 2:
                        multidim_assembly_[2_d]->assemble_fluxes_element_element(edge_integral_data_[i].edge_side_range);
                        break;
                    case 3:
                        multidim_assembly_[3_d]->assemble_fluxes_element_element(edge_integral_data_[i].edge_side_range);
                        break;
                    }
                }
                END_TIMER("assemble_fluxes_elem_elem");
            }

            if (active_integrals_ & ActiveIntegrals::coupling) {
                START_TIMER("assemble_fluxes_elem_side");
                for (i=0; i<integrals_size_[2]; ++i) { // coupling integral
                    switch (coupling_integral_data_[i].side.dim()) {
                    case 2:
                        multidim_assembly_[2_d]->assemble_fluxes_element_side(coupling_integral_data_[i].cell, coupling_integral_data_[i].side);
                        break;
                    case 3:
                        multidim_assembly_[3_d]->assemble_fluxes_element_side(coupling_integral_data_[i].cell, coupling_integral_data_[i].side);
                        break;
                    }
                }
                END_TIMER("assemble_fluxes_elem_side");
            }
        }
        multidim_assembly_[1_d]->end();
    }

private:
    /// Mark eval points in table of Element cache map.
    void insert_eval_points_from_integral_data() {
        for (unsigned int i=0; i<integrals_size_[0]; ++i) {
            // add data to cache if there is free space, else return
        	unsigned int data_size = eval_points_->subset_size( bulk_integral_data_[i].cell.dim(), bulk_integral_data_[i].subset_index );
        	element_cache_map_.mark_used_eval_points(bulk_integral_data_[i].cell, bulk_integral_data_[i].subset_index, data_size);
        }

        for (unsigned int i=0; i<integrals_size_[1]; ++i) {
            // add data to cache if there is free space, else return
        	for (DHCellSide edge_side : edge_integral_data_[i].edge_side_range) {
        	    unsigned int side_dim = edge_side.dim();
                unsigned int data_size = eval_points_->subset_size( side_dim, edge_integral_data_[i].subset_index ) / (side_dim+1);
                unsigned int start_point = data_size * edge_side.side_idx();
                element_cache_map_.mark_used_eval_points(edge_side.cell(), edge_integral_data_[i].subset_index, data_size, start_point);
        	}
        }

        for (unsigned int i=0; i<integrals_size_[2]; ++i) {
            // add data to cache if there is free space, else return
            unsigned int bulk_data_size = eval_points_->subset_size( coupling_integral_data_[i].cell.dim(), coupling_integral_data_[i].bulk_subset_index );
            element_cache_map_.mark_used_eval_points(coupling_integral_data_[i].cell, coupling_integral_data_[i].bulk_subset_index, bulk_data_size);

            unsigned int side_dim = coupling_integral_data_[i].side.dim();
            unsigned int side_data_size = eval_points_->subset_size( side_dim, coupling_integral_data_[i].side_subset_index ) / (side_dim+1);
            unsigned int start_point = side_data_size * coupling_integral_data_[i].side.side_idx();
            element_cache_map_.mark_used_eval_points(coupling_integral_data_[i].side.cell(), coupling_integral_data_[i].side_subset_index, side_data_size, start_point);
        }

        for (unsigned int i=0; i<integrals_size_[3]; ++i) {
            // add data to cache if there is free space, else return
        	unsigned int side_dim = boundary_integral_data_[i].side.dim();
            unsigned int data_size = eval_points_->subset_size( side_dim, boundary_integral_data_[i].subset_index ) / (side_dim+1);
            unsigned int start_point = data_size * boundary_integral_data_[i].side.side_idx();
            element_cache_map_.mark_used_eval_points(boundary_integral_data_[i].side.cell(), boundary_integral_data_[i].subset_index, data_size, start_point);
        }
    }

    /**
     * Add data of integrals to appropriate structure and register elements to ElementCacheMap.
     *
     * Types of used integrals must be set in data member \p active_integrals_.
     */
    void add_integrals_of_computing_step(DHCellAccessor cell) {
        for (unsigned int i=0; i<4; i++) integrals_size_[i] = 0; // clean integral data from previous step
        element_cache_map_.start_elements_update();

        // generic_assembly.check_integral_data();
        if (active_integrals_ & ActiveIntegrals::bulk)
    	    if (cell.is_own()) { // Not ghost
                this->add_volume_integral(cell);
                element_cache_map_.add(cell);
    	    }

        for( DHCellSide cell_side : cell.side_range() ) {
            if (active_integrals_ & ActiveIntegrals::boundary)
                if (cell.is_own()) // Not ghost
                    if ( (cell_side.side().edge().n_sides() == 1) && (cell_side.side().is_boundary()) ) {
                        this->add_boundary_integral(cell_side);
                        element_cache_map_.add(cell_side);
                        continue;
                    }
            if (active_integrals_ & ActiveIntegrals::edge)
                if ( (cell_side.n_edge_sides() >= 2) && (cell_side.edge_sides().begin()->element().idx() == cell.elm_idx())) {
                    this->add_edge_integral(cell_side.edge_sides());
                	for( DHCellSide edge_side : cell_side.edge_sides() ) {
                		element_cache_map_.add(edge_side);
                    }
                }
        }

        if (active_integrals_ & ActiveIntegrals::coupling)
            for( DHCellSide neighb_side : cell.neighb_sides() ) { // cell -> elm lower dim, neighb_side -> elm higher dim
                if (cell.dim() != neighb_side.dim()-1) continue;
                this->add_coupling_integral(cell, neighb_side);
                element_cache_map_.add(cell);
                element_cache_map_.add(neighb_side);
            }

        element_cache_map_.prepare_elements_to_update();
        this->insert_eval_points_from_integral_data();
        element_cache_map_.create_elements_points_map();
        // not used yet: TODO need fix in MultiField, HeatModel ...; need better access to EqData
        //multidim_assembly_[1]->data_->cache_update(element_cache_map_);
        element_cache_map_.finish_elements_update();
    }

    /// Add data of volume integral to appropriate data structure.
    void add_volume_integral(const DHCellAccessor &cell) {
        bulk_integral_data_[ integrals_size_[0] ].cell = cell;
        bulk_integral_data_[ integrals_size_[0] ].subset_index = integrals_.bulk_[cell.dim()-1]->get_subset_idx();
        integrals_size_[0]++;
    }

    /// Add data of edge integral to appropriate data structure.
    void add_edge_integral(RangeConvert<DHEdgeSide, DHCellSide> edge_side_range) {
    	edge_integral_data_[ integrals_size_[1] ].edge_side_range = edge_side_range;
    	edge_integral_data_[ integrals_size_[1] ].subset_index = integrals_.edge_[edge_side_range.begin()->dim()-1]->get_subset_idx();
        integrals_size_[1]++;
    }

    /// Add data of coupling integral to appropriate data structure.
    void add_coupling_integral(const DHCellAccessor &cell, const DHCellSide &ngh_side) {
    	coupling_integral_data_[ integrals_size_[2] ].cell = cell;
    	coupling_integral_data_[ integrals_size_[2] ].side = ngh_side;
    	coupling_integral_data_[ integrals_size_[2] ].bulk_subset_index = integrals_.coupling_[cell.dim()-1]->get_subset_low_idx();
    	coupling_integral_data_[ integrals_size_[2] ].side_subset_index = integrals_.coupling_[cell.dim()-1]->get_subset_high_idx();
        integrals_size_[2]++;
    }

    /// Add data of boundary integral to appropriate data structure.
    void add_boundary_integral(const DHCellSide &bdr_side) {
    	boundary_integral_data_[ integrals_size_[3] ].side = bdr_side;
    	boundary_integral_data_[ integrals_size_[3] ].subset_index = integrals_.boundary_[bdr_side.dim()-1]->get_subset_idx();
        integrals_size_[3]++;
    }


    /// Assembly object
    MixedPtr<DimAssembly, 1> multidim_assembly_;

    /// Holds mask of active integrals.
    int active_integrals_;

    AssemblyIntegrals integrals_;                                 ///< Holds integral objects.
    std::shared_ptr<EvalPoints> eval_points_;                     ///< EvalPoints object shared by all integrals
    ElementCacheMap element_cache_map_;                           ///< ElementCacheMap according to EvalPoints

    // Following variables hold data of all integrals depending of actual computed element.
    std::array<BulkIntegralData, 1>     bulk_integral_data_;      ///< Holds data for computing bulk integrals.
    std::array<EdgeIntegralData, 4>     edge_integral_data_;      ///< Holds data for computing edge integrals.
    std::array<CouplingIntegralData, 6> coupling_integral_data_;  ///< Holds data for computing couplings integrals.
    std::array<BoundaryIntegralData, 4> boundary_integral_data_;  ///< Holds data for computing boundary integrals.
    std::array<unsigned int, 4>         integrals_size_;          ///< Holds used sizes of previous integral data types
};


/**
 * Base class define empty methods, these methods can be overwite in descendants.
 */
template <unsigned int dim>
class AssemblyBase
{
public:
	/// Constructor
	AssemblyBase(unsigned int quad_order) {
        quad_ = new QGauss(dim, 2*quad_order);
        quad_low_ = new QGauss(dim-1, 2*quad_order);
	}

	// Destructor
    virtual ~AssemblyBase() {
        delete quad_;
        delete quad_low_;
    }

    /// Assembles the volume integrals on cell.
    virtual void assemble_volume_integrals(FMT_UNUSED DHCellAccessor cell) {}

    /// Assembles the fluxes on the boundary.
    virtual void assemble_fluxes_boundary(FMT_UNUSED DHCellSide cell_side) {}

    /// Assembles the fluxes between sides on the edge.
    virtual void assemble_fluxes_element_element(FMT_UNUSED RangeConvert<DHEdgeSide, DHCellSide> edge_side_range) {}

    /// Assembles the fluxes between elements of different dimensions.
    virtual void assemble_fluxes_element_side(FMT_UNUSED DHCellAccessor cell_lower_dim, FMT_UNUSED DHCellSide neighb_side) {}

    /// Method prepares object before assemblation (e.g. balance, ...).
    virtual void begin() {}

    /// Method finishes object after assemblation (e.g. balance, ...).
    virtual void end() {}

    /// Create integrals according to dim of assembly object
    void create_integrals(std::shared_ptr<EvalPoints> eval_points, AssemblyIntegrals &integrals, int active_integrals) {
    	if (active_integrals & ActiveIntegrals::bulk)
    	    integrals.bulk_[dim-1] = eval_points->add_bulk<dim>(*quad_);
    	if (active_integrals & ActiveIntegrals::edge)
    	    integrals.edge_[dim-1] = eval_points->add_edge<dim>(*quad_low_);
       	if ((dim>1) && (active_integrals & ActiveIntegrals::coupling))
       	    integrals.coupling_[dim-2] = eval_points->add_coupling<dim>(*quad_low_);
       	if (active_integrals & ActiveIntegrals::boundary)
       	    integrals.boundary_[dim-1] = eval_points->add_boundary<dim>(*quad_low_);
    }

protected:
    Quadrature *quad_;                                     ///< Quadrature used in assembling methods.
    Quadrature *quad_low_;                                 ///< Quadrature used in assembling methods (dim-1).
};


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim, class Model>
class MassAssemblyDG : public AssemblyBase<dim>
{
public:
    typedef typename TransportDG<Model>::EqData EqDataDG;

    /// Constructor.
    MassAssemblyDG(EqDataDG *data)
    : AssemblyBase<dim>(data->dg_order), model_(nullptr), data_(data) {}

    /// Destructor.
    ~MassAssemblyDG() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(TransportDG<Model> &model) {
        this->model_ = &model;

        fe_ = std::make_shared< FE_P_disc<dim> >(data_->dg_order);
        fe_values_.initialize(*this->quad_, *fe_, update_values | update_gradients | update_JxW_values | update_quadrature_points);
        ndofs_ = fe_->n_dofs();
        qsize_ = this->quad_->size();
        dof_indices_.resize(ndofs_);
        local_matrix_.resize(4*ndofs_*ndofs_);
        local_retardation_balance_vector_.resize(ndofs_);
        local_mass_balance_vector_.resize(ndofs_);

        mm_coef_.resize(qsize_);
        ret_coef_.resize(model_->n_substances());
        for (unsigned int sbi=0; sbi<model_->n_substances(); sbi++)
        {
            ret_coef_[sbi].resize(qsize_);
        }
    }


    /// Assemble integral over element
    void assemble_volume_integrals(DHCellAccessor cell) override
    {
        ASSERT_EQ_DBG(cell.dim(), dim).error("Dimension of element mismatch!");
        ElementAccessor<3> elm = cell.elm();

        fe_values_.reinit(elm);
        cell.get_dof_indices(dof_indices_);

        model_->compute_mass_matrix_coefficient(fe_values_.point_list(), elm, mm_coef_);
        model_->compute_retardation_coefficient(fe_values_.point_list(), elm, ret_coef_);

        for (unsigned int sbi=0; sbi<model_->n_substances(); ++sbi)
        {
            // assemble the local mass matrix
            for (unsigned int i=0; i<ndofs_; i++)
            {
                for (unsigned int j=0; j<ndofs_; j++)
                {
                    local_matrix_[i*ndofs_+j] = 0;
                    for (unsigned int k=0; k<qsize_; k++)
                        local_matrix_[i*ndofs_+j] += (mm_coef_[k]+ret_coef_[sbi][k])*fe_values_.shape_value(j,k)*fe_values_.shape_value(i,k)*fe_values_.JxW(k);
                }
            }

            for (unsigned int i=0; i<ndofs_; i++)
            {
                local_mass_balance_vector_[i] = 0;
                local_retardation_balance_vector_[i] = 0;
                for (unsigned int k=0; k<qsize_; k++)
                {
                    local_mass_balance_vector_[i] += mm_coef_[k]*fe_values_.shape_value(i,k)*fe_values_.JxW(k);
                    local_retardation_balance_vector_[i] -= ret_coef_[sbi][k]*fe_values_.shape_value(i,k)*fe_values_.JxW(k);
                }
            }

            model_->balance()->add_mass_values(model_->get_subst_idx()[sbi], cell, cell.get_loc_dof_indices(),
                                               local_mass_balance_vector_, 0);

            data_->ls_dt[sbi]->mat_set_values(ndofs_, &(dof_indices_[0]), ndofs_, &(dof_indices_[0]), &(local_matrix_[0]));
            VecSetValues(data_->ret_vec[sbi], ndofs_, &(dof_indices_[0]), &(local_retardation_balance_vector_[0]), ADD_VALUES);
        }
    }

    /// Implements @p AssemblyBase::begin.
    void begin() override
    {
        model_->balance()->start_mass_assembly( model_->subst_idx() );
    }

    /// Implements @p AssemblyBase::end.
    void end() override
    {
        model_->balance()->finish_mass_assembly( model_->subst_idx() );
    }


    private:
    	/**
    	 * @brief Calculates the velocity field on a given cell.
    	 *
    	 * @param cell       The cell.
    	 * @param velocity   The computed velocity field (at quadrature points).
    	 * @param point_list The quadrature points.
    	 */
        void calculate_velocity(const ElementAccessor<3> &cell, vector<arma::vec3> &velocity,
                                const Armor::array &point_list)
        {
            velocity.resize(point_list.size());
            model_->velocity_field_ptr()->value_list(point_list, cell, velocity);
        }

        shared_ptr<FiniteElement<dim>> fe_;                    ///< Finite element for the solution of the advection-diffusion equation.

        /// Pointer to model (we must use common ancestor of concentration and heat model)
        TransportDG<Model> *model_;

        /// Data object shared with TransportDG
        EqDataDG *data_;

        unsigned int ndofs_;                                      ///< Number of dofs
        unsigned int qsize_;                                      ///< Size of quadrature of actual dim
        FEValues<3> fe_values_;                                   ///< FEValues of object (of P disc finite element type)

        vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
        vector<PetscScalar> local_matrix_;                        ///< Auxiliary vector for assemble methods
        vector<PetscScalar> local_retardation_balance_vector_;    ///< Auxiliary vector for assemble mass matrix.
        vector<PetscScalar> local_mass_balance_vector_;           ///< Same as previous.

    	/// Mass matrix coefficients.
    	vector<double> mm_coef_;
    	/// Retardation coefficient due to sorption.
    	vector<vector<double> > ret_coef_;

        friend class TransportDG<Model>;
        template < template<IntDim...> class DimAssembly>
        friend class GenericAssembly;

};


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim, class Model>
class StiffnessAssemblyDG : public AssemblyBase<dim>
{
public:
    typedef typename TransportDG<Model>::EqData EqDataDG;

    /// Constructor.
    StiffnessAssemblyDG(EqDataDG *data)
    : AssemblyBase<dim>(data->dg_order), fe_rt_(nullptr), model_(nullptr), data_(data) {}

    /// Destructor.
    ~StiffnessAssemblyDG() {
        if (fe_rt_==nullptr) return; // uninitialized object

    	delete fe_rt_;
        delete fe_rt_low_;
    }

    /// Initialize auxiliary vectors and other data members
    void initialize(TransportDG<Model> &model) {
        this->model_ = &model;

        fe_ = std::make_shared< FE_P_disc<dim> >(data_->dg_order);
        fe_low_ = std::make_shared< FE_P_disc<dim-1> >(data_->dg_order);
        fe_rt_ = new FE_RT0<dim>();
        fe_rt_low_ = new FE_RT0<dim-1>();
        fv_rt_.initialize(*this->quad_, *fe_rt_, update_values | update_gradients | update_quadrature_points);
        fe_values_.initialize(*this->quad_, *fe_, update_values | update_gradients | update_JxW_values | update_quadrature_points);
        if (dim>1) {
            fv_rt_vb_.initialize(*this->quad_low_, *fe_rt_low_, update_values | update_quadrature_points);
            fe_values_vb_.initialize(*this->quad_low_, *fe_low_, update_values | update_gradients | update_JxW_values | update_quadrature_points);
        }
        fe_values_side_.initialize(*this->quad_low_, *fe_, update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
        fsv_rt_.initialize(*this->quad_low_, *fe_rt_, update_values | update_quadrature_points);
        ndofs_ = fe_->n_dofs();
        qsize_ = this->quad_->size();
        qsize_lower_dim_ = this->quad_low_->size();
        dof_indices_.resize(ndofs_);
        side_dof_indices_vb_.resize(2*ndofs_);
        local_matrix_.resize(4*ndofs_*ndofs_);
        local_retardation_balance_vector_.resize(ndofs_);
        local_mass_balance_vector_.resize(ndofs_);
        velocity_.resize(qsize_);
        side_velocity_vec_.resize(data_->ad_coef_edg.size());
        sources_sigma_.resize(model_->n_substances(), std::vector<double>(qsize_));
        sigma_.resize(qsize_lower_dim_);
        csection_.resize(qsize_lower_dim_);
        csection_higher_.resize(qsize_lower_dim_);
        dg_penalty_.resize(data_->ad_coef_edg.size());

        mm_coef_.resize(qsize_);
        ret_coef_.resize(model_->n_substances());
        for (unsigned int sbi=0; sbi<model_->n_substances(); sbi++)
        {
            ret_coef_[sbi].resize(qsize_);
        }

        fe_values_vec_.resize(data_->ad_coef_edg.size());
        for (unsigned int sid=0; sid<data_->ad_coef_edg.size(); sid++)
        {
            side_dof_indices_.push_back( vector<LongIdx>(ndofs_) );
            fe_values_vec_[sid].initialize(*this->quad_low_, *fe_,
                    update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
        }

        // index 0 = element with lower dimension,
        // index 1 = side of element with higher dimension
        fv_sb_.resize(2);
        fv_sb_[0] = &fe_values_vb_;
        fv_sb_[1] = &fe_values_side_;
    }


    /// Assembles the volume integrals into the stiffness matrix.
    void assemble_volume_integrals(DHCellAccessor cell) override
    {
        ASSERT_EQ_DBG(cell.dim(), dim).error("Dimension of element mismatch!");
        if (!cell.is_own()) return;

        ElementAccessor<3> elm = cell.elm();

        fe_values_.reinit(elm);
        fv_rt_.reinit(elm);
        cell.get_dof_indices(dof_indices_);

        calculate_velocity(elm, velocity_, fv_rt_.point_list());
        model_->compute_advection_diffusion_coefficients(fe_values_.point_list(), velocity_, elm, data_->ad_coef, data_->dif_coef);
        model_->compute_sources_sigma(fe_values_.point_list(), elm, sources_sigma_);

        // assemble the local stiffness matrix
        for (unsigned int sbi=0; sbi<model_->n_substances(); sbi++)
        {
            for (unsigned int i=0; i<ndofs_; i++)
                for (unsigned int j=0; j<ndofs_; j++)
                    local_matrix_[i*ndofs_+j] = 0;

            for (unsigned int k=0; k<qsize_; k++)
            {
                for (unsigned int i=0; i<ndofs_; i++)
                {
                    arma::vec3 Kt_grad_i = data_->dif_coef[sbi][k].t()*fe_values_.shape_grad(i,k);
                    double ad_dot_grad_i = arma::dot(data_->ad_coef[sbi][k], fe_values_.shape_grad(i,k));

                    for (unsigned int j=0; j<ndofs_; j++)
                        local_matrix_[i*ndofs_+j] += (arma::dot(Kt_grad_i, fe_values_.shape_grad(j,k))
                                                  -fe_values_.shape_value(j,k)*ad_dot_grad_i
                                                  +sources_sigma_[sbi][k]*fe_values_.shape_value(j,k)*fe_values_.shape_value(i,k))*fe_values_.JxW(k);
                }
            }
            data_->ls[sbi]->mat_set_values(ndofs_, &(dof_indices_[0]), ndofs_, &(dof_indices_[0]), &(local_matrix_[0]));
        }
    }


    /// Assembles the fluxes on the boundary.
    void assemble_fluxes_boundary(DHCellSide cell_side) override
    {
        ASSERT_EQ_DBG(cell_side.dim(), dim).error("Dimension of element mismatch!");
        if (!cell_side.cell().is_own()) return;

        Side side = cell_side.side();
        const DHCellAccessor &cell = cell_side.cell();

        ElementAccessor<3> elm_acc = cell.elm();
        cell.get_dof_indices(dof_indices_);
        fe_values_side_.reinit(side);
        fsv_rt_.reinit(side);

        calculate_velocity(elm_acc, velocity_, fsv_rt_.point_list());
        model_->compute_advection_diffusion_coefficients(fe_values_side_.point_list(), velocity_, elm_acc, data_->ad_coef, data_->dif_coef);
        arma::uvec bc_type;
        model_->get_bc_type(side.cond().element_accessor(), bc_type);
        data_->cross_section.value_list(fe_values_side_.point_list(), elm_acc, csection_);

        for (unsigned int sbi=0; sbi<model_->n_substances(); sbi++)
        {
            std::fill(local_matrix_.begin(), local_matrix_.end(), 0);

            // On Neumann boundaries we have only term from integrating by parts the advective term,
            // on Dirichlet boundaries we additionally apply the penalty which enforces the prescribed value.
            double side_flux = 0;
            for (unsigned int k=0; k<qsize_lower_dim_; k++)
                side_flux += arma::dot(data_->ad_coef[sbi][k], fe_values_side_.normal_vector(k))*fe_values_side_.JxW(k);
            double transport_flux = side_flux/side.measure();

            if (bc_type[sbi] == AdvectionDiffusionModel::abc_dirichlet)
            {
                // set up the parameters for DG method
                double gamma_l;
                data_->set_DG_parameters_boundary(side, qsize_lower_dim_, data_->dif_coef[sbi], transport_flux, fe_values_side_.normal_vector(0), data_->dg_penalty[sbi].value(elm_acc.centre(), elm_acc), gamma_l);
                data_->gamma[sbi][side.cond_idx()] = gamma_l;
                transport_flux += gamma_l;
            }

            // fluxes and penalty
            for (unsigned int k=0; k<qsize_lower_dim_; k++)
            {
                double flux_times_JxW;
                if (bc_type[sbi] == AdvectionDiffusionModel::abc_total_flux)
                {
                    //sigma_ corresponds to robin_sigma
                    model_->get_flux_bc_sigma(sbi, fe_values_side_.point_list(), side.cond().element_accessor(), sigma_);
                    flux_times_JxW = csection_[k]*sigma_[k]*fe_values_side_.JxW(k);
                }
                else if (bc_type[sbi] == AdvectionDiffusionModel::abc_diffusive_flux)
                {
                    model_->get_flux_bc_sigma(sbi, fe_values_side_.point_list(), side.cond().element_accessor(), sigma_);
                    flux_times_JxW = (transport_flux + csection_[k]*sigma_[k])*fe_values_side_.JxW(k);
                }
                else if (bc_type[sbi] == AdvectionDiffusionModel::abc_inflow && side_flux < 0)
                    flux_times_JxW = 0;
                else
                    flux_times_JxW = transport_flux*fe_values_side_.JxW(k);

                for (unsigned int i=0; i<ndofs_; i++)
                {
                    for (unsigned int j=0; j<ndofs_; j++)
                    {
                        // flux due to advection and penalty
                        local_matrix_[i*ndofs_+j] += flux_times_JxW*fe_values_side_.shape_value(i,k)*fe_values_side_.shape_value(j,k);

                        // flux due to diffusion (only on dirichlet and inflow boundary)
                        if (bc_type[sbi] == AdvectionDiffusionModel::abc_dirichlet)
                            local_matrix_[i*ndofs_+j] -= (arma::dot(data_->dif_coef[sbi][k]*fe_values_side_.shape_grad(j,k),fe_values_side_.normal_vector(k))*fe_values_side_.shape_value(i,k)
                                    + arma::dot(data_->dif_coef[sbi][k]*fe_values_side_.shape_grad(i,k),fe_values_side_.normal_vector(k))*fe_values_side_.shape_value(j,k)*data_->dg_variant
                                    )*fe_values_side_.JxW(k);
                    }
                }
            }

            data_->ls[sbi]->mat_set_values(ndofs_, &(dof_indices_[0]), ndofs_, &(dof_indices_[0]), &(local_matrix_[0]));
        }
    }


    /// Assembles the fluxes between elements of the same dimension.
    void assemble_fluxes_element_element(RangeConvert<DHEdgeSide, DHCellSide> edge_side_range) override {
        ASSERT_EQ_DBG(edge_side_range.begin()->element().dim(), dim).error("Dimension of element mismatch!");

   	    sid=0;
        for( DHCellSide edge_side : edge_side_range )
        {
            auto dh_edge_cell = data_->dh_->cell_accessor_from_element( edge_side.elem_idx() );
            ElementAccessor<3> edg_elm = dh_edge_cell.elm();
            dh_edge_cell.get_dof_indices(side_dof_indices_[sid]);
            fe_values_vec_[sid].reinit(edge_side.side());
            fsv_rt_.reinit(edge_side.side());
            calculate_velocity(edg_elm, side_velocity_vec_[sid], fsv_rt_.point_list());
            model_->compute_advection_diffusion_coefficients(fe_values_vec_[sid].point_list(), side_velocity_vec_[sid], edg_elm, data_->ad_coef_edg[sid], data_->dif_coef_edg[sid]);
            dg_penalty_[sid].resize(model_->n_substances());
            for (unsigned int sbi=0; sbi<model_->n_substances(); sbi++)
                dg_penalty_[sid][sbi] = data_->dg_penalty[sbi].value(edg_elm.centre(), edg_elm);
            ++sid;
        }
        arma::vec3 normal_vector = fe_values_vec_[0].normal_vector(0);

        // fluxes and penalty
        for (unsigned int sbi=0; sbi<model_->n_substances(); sbi++)
        {
            vector<double> fluxes(edge_side_range.begin()->n_edge_sides());
            double pflux = 0, nflux = 0; // calculate the total in- and out-flux through the edge
            sid=0;
            for( DHCellSide edge_side : edge_side_range )
            {
                fluxes[sid] = 0;
                for (unsigned int k=0; k<qsize_lower_dim_; k++)
                    fluxes[sid] += arma::dot(data_->ad_coef_edg[sid][sbi][k], fe_values_vec_[sid].normal_vector(k))*fe_values_vec_[sid].JxW(k);
                fluxes[sid] /= edge_side.measure();
                if (fluxes[sid] > 0)
                    pflux += fluxes[sid];
                else
                    nflux += fluxes[sid];
                ++sid;
            }

            s1=0;
            for( DHCellSide edge_side1 : edge_side_range )
            {
                s2=-1; // need increment at begin of loop (see conditionally 'continue' directions)
                for( DHCellSide edge_side2 : edge_side_range )
                {
                    s2++;
                    if (s2<=s1) continue;
                    ASSERT(edge_side1.is_valid()).error("Invalid side of edge.");

                    arma::vec3 nv = fe_values_vec_[s1].normal_vector(0);

                    // set up the parameters for DG method
                    // calculate the flux from edge_side1 to edge_side2
                    if (fluxes[s2] > 0 && fluxes[s1] < 0)
                        transport_flux = fluxes[s1]*fabs(fluxes[s2]/pflux);
                    else if (fluxes[s2] < 0 && fluxes[s1] > 0)
                        transport_flux = fluxes[s1]*fabs(fluxes[s2]/nflux);
                    else
                        transport_flux = 0;

                    gamma_l = 0.5*fabs(transport_flux);

                    delta[0] = 0;
                    delta[1] = 0;
                    for (unsigned int k=0; k<qsize_lower_dim_; k++)
                    {
                        delta[0] += dot(data_->dif_coef_edg[s1][sbi][k]*normal_vector,normal_vector);
                        delta[1] += dot(data_->dif_coef_edg[s2][sbi][k]*normal_vector,normal_vector);
                    }
                    delta[0] /= qsize_lower_dim_;
                    delta[1] /= qsize_lower_dim_;

                    delta_sum = delta[0] + delta[1];

//                        if (delta_sum > numeric_limits<double>::epsilon())
                    if (fabs(delta_sum) > 0)
                    {
                        omega[0] = delta[1]/delta_sum;
                        omega[1] = delta[0]/delta_sum;
                        double local_alpha = max(dg_penalty_[s1][sbi], dg_penalty_[s2][sbi]);
                        double h = edge_side1.diameter();
                        aniso1 = data_->elem_anisotropy(edge_side1.element());
                        aniso2 = data_->elem_anisotropy(edge_side2.element());
                        gamma_l += local_alpha/h*aniso1*aniso2*(delta[0]*delta[1]/delta_sum);
                    }
                    else
                        for (int i=0; i<2; i++) omega[i] = 0;
                    // end of set up the parameters for DG method

                    int sd[2]; bool is_side_own[2];
                    sd[0] = s1; is_side_own[0] = edge_side1.cell().is_own();
                    sd[1] = s2; is_side_own[1] = edge_side2.cell().is_own();

#define AVERAGE(i,k,side_id)  (fe_values_vec_[sd[side_id]].shape_value(i,k)*0.5)
#define WAVERAGE(i,k,side_id) (arma::dot(data_->dif_coef_edg[sd[side_id]][sbi][k]*fe_values_vec_[sd[side_id]].shape_grad(i,k),nv)*omega[side_id])
#define JUMP(i,k,side_id)     ((side_id==0?1:-1)*fe_values_vec_[sd[side_id]].shape_value(i,k))

                    // For selected pair of elements:
                    for (int n=0; n<2; n++)
                    {
                        if (!is_side_own[n]) continue;

                        for (int m=0; m<2; m++)
                        {
                            for (unsigned int i=0; i<fe_values_vec_[sd[n]].n_dofs(); i++)
                                for (unsigned int j=0; j<fe_values_vec_[sd[m]].n_dofs(); j++)
                                    local_matrix_[i*fe_values_vec_[sd[m]].n_dofs()+j] = 0;

                            for (unsigned int k=0; k<qsize_lower_dim_; k++)
                            {
                                double flux_times_JxW = transport_flux*fe_values_vec_[0].JxW(k);
                                double gamma_times_JxW = gamma_l*fe_values_vec_[0].JxW(k);

                                for (unsigned int i=0; i<fe_values_vec_[sd[n]].n_dofs(); i++)
                                {
                                    double flux_JxW_jump_i = flux_times_JxW*JUMP(i,k,n);
                                    double gamma_JxW_jump_i = gamma_times_JxW*JUMP(i,k,n);
                                    double JxW_jump_i = fe_values_vec_[0].JxW(k)*JUMP(i,k,n);
                                    double JxW_var_wavg_i = fe_values_vec_[0].JxW(k)*WAVERAGE(i,k,n)*data_->dg_variant;

                                    for (unsigned int j=0; j<fe_values_vec_[sd[m]].n_dofs(); j++)
                                    {
                                        int index = i*fe_values_vec_[sd[m]].n_dofs()+j;

                                        // flux due to transport (applied on interior edges) (average times jump)
                                        local_matrix_[index] += flux_JxW_jump_i*AVERAGE(j,k,m);

                                        // penalty enforcing continuity across edges (applied on interior and Dirichlet edges) (jump times jump)
                                        local_matrix_[index] += gamma_JxW_jump_i*JUMP(j,k,m);

                                        // terms due to diffusion
                                        local_matrix_[index] -= WAVERAGE(j,k,m)*JxW_jump_i;
                                        local_matrix_[index] -= JUMP(j,k,m)*JxW_var_wavg_i;
                                    }
                                }
                            }
                            data_->ls[sbi]->mat_set_values(fe_values_vec_[sd[n]].n_dofs(), &(side_dof_indices_[sd[n]][0]), fe_values_vec_[sd[m]].n_dofs(), &(side_dof_indices_[sd[m]][0]), &(local_matrix_[0]));
                        }
                    }
#undef AVERAGE
#undef WAVERAGE
#undef JUMP
                }
            s1++;
            }
        }
    }


    /// Assembles the fluxes between elements of different dimensions.
    void assemble_fluxes_element_side(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) override {
        if (dim == 1) return;
        ASSERT_EQ_DBG(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

        // Note: use data members csection_ and velocity_ for appropriate quantities of lower dim element

        ElementAccessor<3> elm_lower_dim = cell_lower_dim.elm();
        n_indices = cell_lower_dim.get_dof_indices(dof_indices_);
        for(unsigned int i=0; i<n_indices; ++i) {
            side_dof_indices_vb_[i] = dof_indices_[i];
        }
        fe_values_vb_.reinit(elm_lower_dim);
        n_dofs[0] = fv_sb_[0]->n_dofs();

        DHCellAccessor cell_higher_dim = data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
        ElementAccessor<3> elm_higher_dim = cell_higher_dim.elm();
        n_indices = cell_higher_dim.get_dof_indices(dof_indices_);
        for(unsigned int i=0; i<n_indices; ++i) {
            side_dof_indices_vb_[i+n_dofs[0]] = dof_indices_[i];
        }
        fe_values_side_.reinit(neighb_side.side());
        n_dofs[1] = fv_sb_[1]->n_dofs();

        // Testing element if they belong to local partition.
        bool own_element_id[2];
        own_element_id[0] = cell_lower_dim.is_own();
        own_element_id[1] = cell_higher_dim.is_own();

        fsv_rt_.reinit(neighb_side.side());
        fv_rt_vb_.reinit(elm_lower_dim);
        calculate_velocity(elm_higher_dim, velocity_higher_, fsv_rt_.point_list());
        calculate_velocity(elm_lower_dim, velocity_, fv_rt_vb_.point_list());
        model_->compute_advection_diffusion_coefficients(fe_values_vb_.point_list(), velocity_, elm_lower_dim, data_->ad_coef_edg[0], data_->dif_coef_edg[0]);
        model_->compute_advection_diffusion_coefficients(fe_values_vb_.point_list(), velocity_higher_, elm_higher_dim, data_->ad_coef_edg[1], data_->dif_coef_edg[1]);
        data_->cross_section.value_list(fe_values_vb_.point_list(), elm_lower_dim, csection_);
        data_->cross_section.value_list(fe_values_vb_.point_list(), elm_higher_dim, csection_higher_);

        for (unsigned int sbi=0; sbi<model_->n_substances(); sbi++) // Optimize: SWAP LOOPS
        {
            for (unsigned int i=0; i<n_dofs[0]+n_dofs[1]; i++)
                for (unsigned int j=0; j<n_dofs[0]+n_dofs[1]; j++)
                    local_matrix_[i*(n_dofs[0]+n_dofs[1])+j] = 0;

            // sigma_ corresponds to frac_sigma
            data_->fracture_sigma[sbi].value_list(fe_values_vb_.point_list(), elm_lower_dim, sigma_);

            // set transmission conditions
            for (unsigned int k=0; k<qsize_lower_dim_; k++)
            {
                // The communication flux has two parts:
                // - "diffusive" term containing sigma
                // - "advective" term representing usual upwind
                //
                // The calculation differs from the reference manual, since ad_coef and dif_coef have different meaning
                // than b and A in the manual.
                // In calculation of sigma there appears one more csection_lower in the denominator.
                double sigma = sigma_[k]*arma::dot(data_->dif_coef_edg[0][sbi][k]*fe_values_side_.normal_vector(k),fe_values_side_.normal_vector(k))*
                        2*csection_higher_[k]*csection_higher_[k]/(csection_[k]*csection_[k]);

                double transport_flux = arma::dot(data_->ad_coef_edg[1][sbi][k], fe_values_side_.normal_vector(k));

                comm_flux[0][0] =  (sigma-min(0.,transport_flux))*fv_sb_[0]->JxW(k);
                comm_flux[0][1] = -(sigma-min(0.,transport_flux))*fv_sb_[0]->JxW(k);
                comm_flux[1][0] = -(sigma+max(0.,transport_flux))*fv_sb_[0]->JxW(k);
                comm_flux[1][1] =  (sigma+max(0.,transport_flux))*fv_sb_[0]->JxW(k);

                for (int n=0; n<2; n++)
                {
                    if (!own_element_id[n]) continue;

                    for (unsigned int i=0; i<n_dofs[n]; i++)
                        for (int m=0; m<2; m++)
                            for (unsigned int j=0; j<n_dofs[m]; j++)
                                local_matrix_[(i+n*n_dofs[0])*(n_dofs[0]+n_dofs[1]) + m*n_dofs[0] + j] +=
                                        comm_flux[m][n]*fv_sb_[m]->shape_value(j,k)*fv_sb_[n]->shape_value(i,k);
                }
            }
            data_->ls[sbi]->mat_set_values(n_dofs[0]+n_dofs[1], &(side_dof_indices_vb_[0]), n_dofs[0]+n_dofs[1], &(side_dof_indices_vb_[0]), &(local_matrix_[0]));
        }
    }


private:
	/**
	 * @brief Calculates the velocity field on a given cell.
	 *
	 * @param cell       The cell.
	 * @param velocity   The computed velocity field (at quadrature points).
	 * @param point_list The quadrature points.
	 */
    void calculate_velocity(const ElementAccessor<3> &cell, vector<arma::vec3> &velocity,
                            const Armor::array &point_list)
    {
        velocity.resize(point_list.size());
        model_->velocity_field_ptr()->value_list(point_list, cell, velocity);
    }

    shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.
    shared_ptr<FiniteElement<dim-1>> fe_low_;   ///< Finite element for the solution of the advection-diffusion equation (dim-1).
    FiniteElement<dim> *fe_rt_;                 ///< Finite element for the water velocity field.
    FiniteElement<dim-1> *fe_rt_low_;           ///< Finite element for the water velocity field (dim-1).

    /// Pointer to model (we must use common ancestor of concentration and heat model)
    TransportDG<Model> *model_;

    /// Data object shared with TransportDG
    EqDataDG *data_;

    unsigned int ndofs_;                                      ///< Number of dofs
    unsigned int qsize_;                                      ///< Size of quadrature of actual dim
    unsigned int qsize_lower_dim_;                            ///< Size of quadrature of dim-1
    FEValues<3> fv_rt_;                                       ///< FEValues of object (of RT0 finite element type)
    FEValues<3> fe_values_;                                   ///< FEValues of object (of P disc finite element type)
    FEValues<3> fv_rt_vb_;                                    ///< FEValues of dim-1 object (of RT0 finite element type)
    FEValues<3> fe_values_vb_;                                ///< FEValues of dim-1 object (of P disc finite element type)
    FEValues<3> fe_values_side_;                              ///< FEValues of object (of P disc finite element type)
    FEValues<3> fsv_rt_;                                      ///< FEValues of object (of RT0 finite element type)
    vector<FEValues<3>> fe_values_vec_;                       ///< Vector of FEValues of object (of P disc finite element types)
    vector<FEValues<3>*> fv_sb_;                              ///< Auxiliary vector, holds FEValues objects for assemble element-side

    vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
    vector< vector<LongIdx> > side_dof_indices_;              ///< Vector of vectors of side DOF indices
    vector<LongIdx> side_dof_indices_vb_;                     ///< Vector of side DOF indices (assemble element-side fluxex)
    vector<PetscScalar> local_matrix_;                        ///< Auxiliary vector for assemble methods
    vector<PetscScalar> local_retardation_balance_vector_;    ///< Auxiliary vector for assemble mass matrix.
    vector<PetscScalar> local_mass_balance_vector_;           ///< Same as previous.
    vector<arma::vec3> velocity_;                             ///< Auxiliary results.
    vector<arma::vec3> velocity_higher_;                      ///< Velocity results of higher dim element (element-side computation).
    vector<vector<arma::vec3> > side_velocity_vec_;           ///< Vector of velocities results.
    vector<vector<double> > sources_sigma_;                   ///< Auxiliary vectors for assemble volume integrals and set_sources method.
    vector<double> sigma_;                                    ///< Auxiliary vector for assemble boundary fluxes (robin sigma), element-side fluxes (frac sigma) and set boundary conditions method
    vector<double> csection_;                                 ///< Auxiliary vector for assemble boundary fluxes, element-side fluxes and set boundary conditions
    vector<double> csection_higher_;                          ///< Auxiliary vector for assemble element-side fluxes
    vector<vector<double> > dg_penalty_;                      ///< Auxiliary vectors for assemble element-element fluxes

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

    friend class TransportDG<Model>;
    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim, class Model>
class SourcesAssemblyDG : public AssemblyBase<dim>
{
public:
    typedef typename TransportDG<Model>::EqData EqDataDG;

    /// Constructor.
    SourcesAssemblyDG(EqDataDG *data)
    : AssemblyBase<dim>(data->dg_order), model_(nullptr), data_(data) {}

    /// Destructor.
    ~SourcesAssemblyDG() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(TransportDG<Model> &model) {
        this->model_ = &model;

        fe_ = std::make_shared< FE_P_disc<dim> >(data_->dg_order);
        fe_values_.initialize(*this->quad_, *fe_, update_values | update_gradients | update_JxW_values | update_quadrature_points);
        ndofs_ = fe_->n_dofs();
        qsize_ = this->quad_->size();
        dof_indices_.resize(ndofs_);
        loc_dof_indices_.resize(ndofs_);
        local_rhs_.resize(ndofs_);
        local_source_balance_vector_.resize(ndofs_);
        local_source_balance_rhs_.resize(ndofs_);
        sources_conc_.resize(model_->n_substances(), std::vector<double>(qsize_));
        sources_density_.resize(model_->n_substances(), std::vector<double>(qsize_));
        sources_sigma_.resize(model_->n_substances(), std::vector<double>(qsize_));
    }


    /// Assemble integral over element
    void assemble_volume_integrals(DHCellAccessor cell) override
    {
    	ASSERT_EQ_DBG(cell.dim(), dim).error("Dimension of element mismatch!");

        ElementAccessor<3> elm = cell.elm();

        fe_values_.reinit(elm);
        cell.get_dof_indices(dof_indices_);

        model_->compute_source_coefficients(fe_values_.point_list(), elm, sources_conc_, sources_density_, sources_sigma_);

        // assemble the local stiffness matrix
        for (unsigned int sbi=0; sbi<model_->n_substances(); sbi++)
        {
            fill_n( &(local_rhs_[0]), ndofs_, 0 );
            local_source_balance_vector_.assign(ndofs_, 0);
            local_source_balance_rhs_.assign(ndofs_, 0);

            // compute sources
            for (unsigned int k=0; k<qsize_; k++)
            {
                source = (sources_density_[sbi][k] + sources_conc_[sbi][k]*sources_sigma_[sbi][k])*fe_values_.JxW(k);

                for (unsigned int i=0; i<ndofs_; i++)
                    local_rhs_[i] += source*fe_values_.shape_value(i,k);
            }
            data_->ls[sbi]->rhs_set_values(ndofs_, &(dof_indices_[0]), &(local_rhs_[0]));

            for (unsigned int i=0; i<ndofs_; i++)
            {
                for (unsigned int k=0; k<qsize_; k++)
                    local_source_balance_vector_[i] -= sources_sigma_[sbi][k]*fe_values_.shape_value(i,k)*fe_values_.JxW(k);

                local_source_balance_rhs_[i] += local_rhs_[i];
            }
            model_->balance()->add_source_values(model_->get_subst_idx()[sbi], elm.region().bulk_idx(),
                                                 cell.get_loc_dof_indices(),
                                                 local_source_balance_vector_, local_source_balance_rhs_);
        }
    }

    /// Implements @p AssemblyBase::begin.
    void begin() override
    {
        model_->balance()->start_source_assembly( model_->subst_idx() );
    }

    /// Implements @p AssemblyBase::end.
    void end() override
    {
        model_->balance()->finish_source_assembly( model_->subst_idx() );
    }


    private:
    	/**
    	 * @brief Calculates the velocity field on a given cell.
    	 *
    	 * @param cell       The cell.
    	 * @param velocity   The computed velocity field (at quadrature points).
    	 * @param point_list The quadrature points.
    	 */
        void calculate_velocity(const ElementAccessor<3> &cell, vector<arma::vec3> &velocity,
                                const Armor::array &point_list)
        {
            velocity.resize(point_list.size());
            model_->velocity_field_ptr()->value_list(point_list, cell, velocity);
        }

        shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.

        /// Pointer to model (we must use common ancestor of concentration and heat model)
        TransportDG<Model> *model_;

        /// Data object shared with TransportDG
        EqDataDG *data_;

        unsigned int ndofs_;                                      ///< Number of dofs
        unsigned int qsize_;                                      ///< Size of quadrature of actual dim
        FEValues<3> fe_values_;                                   ///< FEValues of object (of P disc finite element type)

        vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
        vector<LongIdx> loc_dof_indices_;                         ///< Vector of local DOF indices
        vector<PetscScalar> local_rhs_;                           ///< Auxiliary vector for set_sources method.
        vector<PetscScalar> local_source_balance_vector_;         ///< Auxiliary vector for set_sources method.
        vector<PetscScalar> local_source_balance_rhs_;            ///< Auxiliary vector for set_sources method.
        vector<vector<double> > sources_conc_;                    ///< Auxiliary vectors for set_sources method.
        vector<vector<double> > sources_density_;                 ///< Auxiliary vectors for set_sources method.
        vector<vector<double> > sources_sigma_;                   ///< Auxiliary vectors for assemble volume integrals and set_sources method.

        /// @name Auxiliary variables used during set sources
    	// @{

        double source;

    	// @}

        friend class TransportDG<Model>;
        template < template<IntDim...> class DimAssembly>
        friend class GenericAssembly;

};


/**
 * Assembles the r.h.s. components corresponding to the Dirichlet boundary conditions..
 */
template <unsigned int dim, class Model>
class BdrConditionAssemblyDG : public AssemblyBase<dim>
{
public:
    typedef typename TransportDG<Model>::EqData EqDataDG;

    /// Constructor.
    BdrConditionAssemblyDG(EqDataDG *data)
    : AssemblyBase<dim>(data->dg_order), fe_rt_(nullptr), model_(nullptr), data_(data) {}

    /// Destructor.
    ~BdrConditionAssemblyDG() {
        if (fe_rt_==nullptr) return; // uninitialized object

        delete fe_rt_;
    }

    /// Initialize auxiliary vectors and other data members
    void initialize(TransportDG<Model> &model) {
        this->model_ = &model;

        fe_ = std::make_shared< FE_P_disc<dim> >(data_->dg_order);
        fe_rt_ = new FE_RT0<dim>();
        fe_values_side_.initialize(*this->quad_low_, *fe_, update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
        fsv_rt_.initialize(*this->quad_low_, *fe_rt_, update_values | update_quadrature_points);
        ndofs_ = fe_->n_dofs();
        qsize_ = this->quad_->size();
        qsize_lower_dim_ = this->quad_low_->size();
        dof_indices_.resize(ndofs_);
        local_rhs_.resize(ndofs_);
        local_flux_balance_vector_.resize(ndofs_);
        velocity_.resize(qsize_);
        sigma_.resize(qsize_lower_dim_);
        csection_.resize(qsize_lower_dim_);
        bc_values_.resize(qsize_lower_dim_);
        bc_fluxes_.resize(qsize_lower_dim_);
        bc_ref_values_.resize(qsize_lower_dim_);
    }


    /// Assemble integral over element
    void assemble_volume_integrals(DHCellAccessor cell) override
    {
        ElementAccessor<3> elm = cell.elm();
        if (elm->boundary_idx_ == nullptr) return;

        for (DHCellSide dh_side : cell.side_range())
        {
            if (dh_side.n_edge_sides() > 1) continue;
            // skip edges lying not on the boundary
            if (! dh_side.side().is_boundary()) continue;

            const unsigned int cond_idx = dh_side.side().cond_idx();

            ElementAccessor<3> bc_elm = dh_side.cond().element_accessor();

            arma::uvec bc_type;
            model_->get_bc_type(bc_elm, bc_type);

            fe_values_side_.reinit(dh_side.side());
            fsv_rt_.reinit(dh_side.side());
            calculate_velocity(elm, velocity_, fsv_rt_.point_list());

            cell.get_dof_indices(dof_indices_);

            model_->compute_advection_diffusion_coefficients(fe_values_side_.point_list(), velocity_, elm, data_->ad_coef, data_->dif_coef);
            data_->cross_section.value_list(fe_values_side_.point_list(), elm, csection_);

            for (unsigned int sbi=0; sbi<model_->n_substances(); sbi++)
            {
                fill_n(&(local_rhs_[0]), ndofs_, 0);
                local_flux_balance_vector_.assign(ndofs_, 0);
                local_flux_balance_rhs_ = 0;

                // The b.c. data are fetched for all possible b.c. types since we allow
                // different bc_type for each substance.
                data_->bc_dirichlet_value[sbi].value_list(fe_values_side_.point_list(), bc_elm, bc_values_);

                double side_flux = 0;
                for (unsigned int k=0; k<qsize_lower_dim_; k++)
                    side_flux += arma::dot(data_->ad_coef[sbi][k], fe_values_side_.normal_vector(k))*fe_values_side_.JxW(k);
                double transport_flux = side_flux/dh_side.measure();

                if (bc_type[sbi] == AdvectionDiffusionModel::abc_inflow && side_flux < 0)
                {
                    for (unsigned int k=0; k<qsize_lower_dim_; k++)
                    {
                        double bc_term = -transport_flux*bc_values_[k]*fe_values_side_.JxW(k);
                        for (unsigned int i=0; i<ndofs_; i++)
                            local_rhs_[i] += bc_term*fe_values_side_.shape_value(i,k);
                    }
                    for (unsigned int i=0; i<ndofs_; i++)
                        local_flux_balance_rhs_ -= local_rhs_[i];
                }
                else if (bc_type[sbi] == AdvectionDiffusionModel::abc_dirichlet)
                {
                    for (unsigned int k=0; k<qsize_lower_dim_; k++)
                    {
                        double bc_term = data_->gamma[sbi][cond_idx]*bc_values_[k]*fe_values_side_.JxW(k);
                        arma::vec3 bc_grad = -bc_values_[k]*fe_values_side_.JxW(k)*data_->dg_variant*(arma::trans(data_->dif_coef[sbi][k])*fe_values_side_.normal_vector(k));
                        for (unsigned int i=0; i<ndofs_; i++)
                            local_rhs_[i] += bc_term*fe_values_side_.shape_value(i,k)
                                    + arma::dot(bc_grad,fe_values_side_.shape_grad(i,k));
                    }
                    for (unsigned int k=0; k<qsize_lower_dim_; k++)
                    {
                        for (unsigned int i=0; i<ndofs_; i++)
                        {
                            local_flux_balance_vector_[i] += (arma::dot(data_->ad_coef[sbi][k], fe_values_side_.normal_vector(k))*fe_values_side_.shape_value(i,k)
                                    - arma::dot(data_->dif_coef[sbi][k]*fe_values_side_.shape_grad(i,k),fe_values_side_.normal_vector(k))
                                    + data_->gamma[sbi][cond_idx]*fe_values_side_.shape_value(i,k))*fe_values_side_.JxW(k);
                        }
                    }
                    if (model_->time().tlevel() > 0)
                        for (unsigned int i=0; i<ndofs_; i++)
                            local_flux_balance_rhs_ -= local_rhs_[i];
                }
                else if (bc_type[sbi] == AdvectionDiffusionModel::abc_total_flux)
                {
                	model_->get_flux_bc_data(sbi, fe_values_side_.point_list(), bc_elm, bc_fluxes_, sigma_, bc_ref_values_);
                    for (unsigned int k=0; k<qsize_lower_dim_; k++)
                    {
                        double bc_term = csection_[k]*(sigma_[k]*bc_ref_values_[k]+bc_fluxes_[k])*fe_values_side_.JxW(k);
                        for (unsigned int i=0; i<ndofs_; i++)
                            local_rhs_[i] += bc_term*fe_values_side_.shape_value(i,k);
                    }

                    for (unsigned int i=0; i<ndofs_; i++)
                    {
                        for (unsigned int k=0; k<qsize_lower_dim_; k++)
                            local_flux_balance_vector_[i] += csection_[k]*sigma_[k]*fe_values_side_.JxW(k)*fe_values_side_.shape_value(i,k);
                        local_flux_balance_rhs_ -= local_rhs_[i];
                    }
                }
                else if (bc_type[sbi] == AdvectionDiffusionModel::abc_diffusive_flux)
                {
                	model_->get_flux_bc_data(sbi, fe_values_side_.point_list(), bc_elm, bc_fluxes_, sigma_, bc_ref_values_);
                    for (unsigned int k=0; k<qsize_lower_dim_; k++)
                    {
                        double bc_term = csection_[k]*(sigma_[k]*bc_ref_values_[k]+bc_fluxes_[k])*fe_values_side_.JxW(k);
                        for (unsigned int i=0; i<ndofs_; i++)
                            local_rhs_[i] += bc_term*fe_values_side_.shape_value(i,k);
                    }

                    for (unsigned int i=0; i<ndofs_; i++)
                    {
                        for (unsigned int k=0; k<qsize_lower_dim_; k++)
                            local_flux_balance_vector_[i] += csection_[k]*(arma::dot(data_->ad_coef[sbi][k], fe_values_side_.normal_vector(k)) + sigma_[k])*fe_values_side_.JxW(k)*fe_values_side_.shape_value(i,k);
                        local_flux_balance_rhs_ -= local_rhs_[i];
                    }
                }
                else if (bc_type[sbi] == AdvectionDiffusionModel::abc_inflow && side_flux >= 0)
                {
                    for (unsigned int k=0; k<qsize_lower_dim_; k++)
                    {
                        for (unsigned int i=0; i<ndofs_; i++)
                            local_flux_balance_vector_[i] += arma::dot(data_->ad_coef[sbi][k], fe_values_side_.normal_vector(k))*fe_values_side_.JxW(k)*fe_values_side_.shape_value(i,k);
                    }
                }
                data_->ls[sbi]->rhs_set_values(ndofs_, &(dof_indices_[0]), &(local_rhs_[0]));

                model_->balance()->add_flux_values(model_->get_subst_idx()[sbi], dh_side,
                                              cell.get_loc_dof_indices(),
                                              local_flux_balance_vector_, local_flux_balance_rhs_);
            }
        }
    }

    /// Implements @p AssemblyBase::begin.
    void begin() override
    {
        model_->balance()->start_flux_assembly( model_->subst_idx() );
    }

    /// Implements @p AssemblyBase::end.
    void end() override
    {
        model_->balance()->finish_flux_assembly( model_->subst_idx() );
    }


    private:
    	/**
    	 * @brief Calculates the velocity field on a given cell.
    	 *
    	 * @param cell       The cell.
    	 * @param velocity   The computed velocity field (at quadrature points).
    	 * @param point_list The quadrature points.
    	 */
        void calculate_velocity(const ElementAccessor<3> &cell, vector<arma::vec3> &velocity,
                                const Armor::array &point_list)
        {
            velocity.resize(point_list.size());
            model_->velocity_field_ptr()->value_list(point_list, cell, velocity);
        }

        shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.
        FiniteElement<dim> *fe_rt_;                 ///< Finite element for the water velocity field.

        /// Pointer to model (we must use common ancestor of concentration and heat model)
        TransportDG<Model> *model_;

        /// Data object shared with TransportDG
        EqDataDG *data_;

        unsigned int ndofs_;                                      ///< Number of dofs
        unsigned int qsize_;                                      ///< Size of quadrature of actual dim
        unsigned int qsize_lower_dim_;                            ///< Size of quadrature of dim-1
        FEValues<3> fe_values_side_;                              ///< FEValues of object (of P disc finite element type)
        FEValues<3> fsv_rt_;                                      ///< FEValues of object (of RT0 finite element type)

        vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
        vector<PetscScalar> local_rhs_;                           ///< Auxiliary vector for set_sources method.
        vector<PetscScalar> local_flux_balance_vector_;           ///< Auxiliary vector for set_boundary_conditions method.
        PetscScalar local_flux_balance_rhs_;                      ///< Auxiliary variable for set_boundary_conditions method.
        vector<arma::vec3> velocity_;                             ///< Auxiliary results.
        vector<double> sigma_;                                    ///< Auxiliary vector for assemble boundary fluxes (robin sigma), element-side fluxes (frac sigma) and set boundary conditions method
        vector<double> csection_;                                 ///< Auxiliary vector for assemble boundary fluxes, element-side fluxes and set boundary conditions
        vector<double> bc_values_;                                ///< Auxiliary vector for set boundary conditions method
        vector<double> bc_fluxes_;                                ///< Same as previous
        vector<double> bc_ref_values_;                            ///< Same as previous

        friend class TransportDG<Model>;
        template < template<IntDim...> class DimAssembly>
        friend class GenericAssembly;

};


/**
 * Auxiliary container class sets the initial condition.
 */
template <unsigned int dim, class Model>
class InitConditionAssemblyDG : public AssemblyBase<dim>
{
public:
    typedef typename TransportDG<Model>::EqData EqDataDG;

    /// Constructor.
    InitConditionAssemblyDG(EqDataDG *data)
    : AssemblyBase<dim>(data->dg_order), model_(nullptr), data_(data) {}

    /// Destructor.
    ~InitConditionAssemblyDG() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(TransportDG<Model> &model) {
        this->model_ = &model;

        fe_ = std::make_shared< FE_P_disc<dim> >(data_->dg_order);
        fe_values_.initialize(*this->quad_, *fe_, update_values | update_gradients | update_JxW_values | update_quadrature_points);
        ndofs_ = fe_->n_dofs();
        qsize_ = this->quad_->size();
        dof_indices_.resize(ndofs_);
        local_matrix_.resize(4*ndofs_*ndofs_);
        local_rhs_.resize(ndofs_);
        init_values_.resize(model_->n_substances(), std::vector<double>(qsize_));
    }


    /// Assemble integral over element
    void assemble_volume_integrals(DHCellAccessor cell) override
    {
        ASSERT_EQ_DBG(cell.dim(), dim).error("Dimension of element mismatch!");

        ElementAccessor<3> elem = cell.elm();
        cell.get_dof_indices(dof_indices_);
        fe_values_.reinit(elem);
        model_->compute_init_cond(fe_values_.point_list(), elem, init_values_);

        for (unsigned int sbi=0; sbi<model_->n_substances(); sbi++)
        {
            for (unsigned int i=0; i<ndofs_; i++)
            {
                local_rhs_[i] = 0;
                for (unsigned int j=0; j<ndofs_; j++)
                    local_matrix_[i*ndofs_+j] = 0;
            }

            for (unsigned int k=0; k<qsize_; k++)
            {
                double rhs_term = init_values_[sbi][k]*fe_values_.JxW(k);

                for (unsigned int i=0; i<ndofs_; i++)
                {
                    for (unsigned int j=0; j<ndofs_; j++)
                        local_matrix_[i*ndofs_+j] += fe_values_.shape_value(i,k)*fe_values_.shape_value(j,k)*fe_values_.JxW(k);

                    local_rhs_[i] += fe_values_.shape_value(i,k)*rhs_term;
                }
            }
            data_->ls[sbi]->set_values(ndofs_, &(dof_indices_[0]), ndofs_, &(dof_indices_[0]), &(local_matrix_[0]), &(local_rhs_[0]));
        }
    }


    private:
    	/**
    	 * @brief Calculates the velocity field on a given cell.
    	 *
    	 * @param cell       The cell.
    	 * @param velocity   The computed velocity field (at quadrature points).
    	 * @param point_list The quadrature points.
    	 */
        void calculate_velocity(const ElementAccessor<3> &cell, vector<arma::vec3> &velocity,
                                const Armor::array &point_list)
        {
            velocity.resize(point_list.size());
            model_->velocity_field_ptr()->value_list(point_list, cell, velocity);
        }

        shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.

        /// Pointer to model (we must use common ancestor of concentration and heat model)
        TransportDG<Model> *model_;

        /// Data object shared with TransportDG
        EqDataDG *data_;

        unsigned int ndofs_;                                      ///< Number of dofs
        unsigned int qsize_;                                      ///< Size of quadrature of actual dim
        FEValues<3> fe_values_;                                   ///< FEValues of object (of P disc finite element type)

        vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
        vector<PetscScalar> local_matrix_;                        ///< Auxiliary vector for assemble methods
        vector<PetscScalar> local_rhs_;                           ///< Auxiliary vector for set_sources method.
        std::vector<std::vector<double> > init_values_;           ///< Auxiliary vectors for prepare initial condition

        friend class TransportDG<Model>;
        template < template<IntDim...> class DimAssembly>
        friend class GenericAssembly;

};



#endif /* ASSEMBLY_DG_HH_ */

