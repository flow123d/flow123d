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


/**
 * Base (non-templated) class of container for Finite element and related objects.
 */
class AssemblyDGBase
{
public:
    virtual ~AssemblyDGBase() {}

    virtual void initialize() = 0;

    virtual void assemble_mass_matrix(DHCellAccessor cell) = 0;

    virtual void assemble_volume_integrals(DHCellAccessor cell) = 0;

    virtual void assemble_fluxes_boundary(DHCellAccessor cell) = 0;

    virtual void assemble_fluxes_element_element(DHCellAccessor cell) = 0;

    virtual void assemble_fluxes_element_side(DHCellAccessor cell_lower_dim) = 0;

    virtual void set_sources(DHCellAccessor cell) = 0;

    virtual void set_boundary_conditions(DHCellAccessor cell) = 0;

    virtual void prepare_initial_condition(DHCellAccessor cell) = 0;
};


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim, class Model>
class AssemblyDG : public AssemblyDGBase
{
public:
    typedef typename TransportDG<Model>::EqData EqDataDG;

    /// Constructor.
    AssemblyDG(std::shared_ptr<EqDataDG> data, TransportDG<Model> &model)
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
    ~AssemblyDG() {
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
    void initialize() override {
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

    /// Assemble integral over element
    void assemble_mass_matrix(DHCellAccessor cell) override
    {
        ASSERT_EQ_DBG(cell.dim(), dim).error("Dimension of element mismatch!");
        ElementAccessor<3> elm = cell.elm();

        fe_values_.reinit(elm);
        cell.get_dof_indices(dof_indices_);

        model_.compute_mass_matrix_coefficient(fe_values_.point_list(), elm, mm_coef_);
        model_.compute_retardation_coefficient(fe_values_.point_list(), elm, ret_coef_);

        for (unsigned int sbi=0; sbi<model_.n_substances(); ++sbi)
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

            model_.balance()->add_mass_matrix_values(model_.get_subst_idx()[sbi], elm.region().bulk_idx(), dof_indices_, local_mass_balance_vector_);
            data_->ls_dt[sbi]->mat_set_values(ndofs_, &(dof_indices_[0]), ndofs_, &(dof_indices_[0]), &(local_matrix_[0]));
            VecSetValues(data_->ret_vec[sbi], ndofs_, &(dof_indices_[0]), &(local_retardation_balance_vector_[0]), ADD_VALUES);
        }
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
        model_.compute_advection_diffusion_coefficients(fe_values_.point_list(), velocity_, elm, data_->ad_coef, data_->dif_coef);
        model_.compute_sources_sigma(fe_values_.point_list(), elm, sources_sigma_);

        // assemble the local stiffness matrix
        for (unsigned int sbi=0; sbi<model_.n_substances(); sbi++)
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
    void assemble_fluxes_boundary(DHCellAccessor cell) override
    {
        ASSERT_EQ_DBG(cell.dim(), dim).error("Dimension of element mismatch!");
        if (!cell.is_own()) return;

        for( DHCellSide cell_side : cell.side_range() )
        {
            Side side = cell_side.side();
            if (side.edge()->n_sides > 1) continue;
            // check spatial dimension
            if (side.dim() != dim-1) continue;
            // skip edges lying not on the boundary
            if (side.cond() == NULL) continue;

            ElementAccessor<3> elm_acc = cell.elm();
            cell.get_dof_indices(dof_indices_);
            fe_values_side_.reinit(elm_acc, side.side_idx());
            fsv_rt_.reinit(elm_acc, side.side_idx());

            calculate_velocity(elm_acc, velocity_, fsv_rt_.point_list());
            model_.compute_advection_diffusion_coefficients(fe_values_side_.point_list(), velocity_, elm_acc, data_->ad_coef, data_->dif_coef);
            arma::uvec bc_type;
            model_.get_bc_type(side.cond()->element_accessor(), bc_type);
            data_->cross_section.value_list(fe_values_side_.point_list(), elm_acc, csection_);

            for (unsigned int sbi=0; sbi<model_.n_substances(); sbi++)
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
                        model_.get_flux_bc_sigma(sbi, fe_values_side_.point_list(), side.cond()->element_accessor(), sigma_);
                        flux_times_JxW = csection_[k]*sigma_[k]*fe_values_side_.JxW(k);
                    }
                    else if (bc_type[sbi] == AdvectionDiffusionModel::abc_diffusive_flux)
                    {
                        model_.get_flux_bc_sigma(sbi, fe_values_side_.point_list(), side.cond()->element_accessor(), sigma_);
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
    }


    /// Assembles the fluxes between elements of the same dimension.
    void assemble_fluxes_element_element(DHCellAccessor cell) override
    {
        ASSERT_EQ_DBG(cell.dim(), dim).error("Dimension of element mismatch!");

        // assemble integral over sides
        for( DHCellSide cell_side : cell.side_range() )
        {
            if (cell_side.n_edge_sides() < 2) continue;
            bool unique_edge = (cell_side.edge_sides().begin()->element().idx() != cell.elm_idx());
      	    if ( unique_edge ) continue;
       	    sid=0;
            for( DHCellSide edge_side : cell_side.edge_sides() )
            {
                auto dh_edge_cell = data_->dh_->cell_accessor_from_element( edge_side.elem_idx() );
                ElementAccessor<3> edg_elm = dh_edge_cell.elm();
                dh_edge_cell.get_dof_indices(side_dof_indices_[sid]);
                fe_values_vec_[sid]->reinit(edg_elm, edge_side.side_idx());
                fsv_rt_.reinit(edg_elm, edge_side.side_idx());
                calculate_velocity(edg_elm, side_velocity_vec_[sid], fsv_rt_.point_list());
                model_.compute_advection_diffusion_coefficients(fe_values_vec_[sid]->point_list(), side_velocity_vec_[sid], edg_elm, data_->ad_coef_edg[sid], data_->dif_coef_edg[sid]);
                dg_penalty_[sid].resize(model_.n_substances());
                for (unsigned int sbi=0; sbi<model_.n_substances(); sbi++)
                    dg_penalty_[sid][sbi] = data_->dg_penalty[sbi].value(edg_elm.centre(), edg_elm);
                ++sid;
            }
            arma::vec3 normal_vector = fe_values_vec_[0]->normal_vector(0);

            // fluxes and penalty
            for (unsigned int sbi=0; sbi<model_.n_substances(); sbi++)
            {
                vector<double> fluxes(cell_side.n_edge_sides());
                double pflux = 0, nflux = 0; // calculate the total in- and out-flux through the edge
                sid=0;
                for( DHCellSide edge_side : cell_side.edge_sides() )
                {
                    fluxes[sid] = 0;
                    for (unsigned int k=0; k<qsize_lower_dim_; k++)
                        fluxes[sid] += arma::dot(data_->ad_coef_edg[sid][sbi][k], fe_values_vec_[sid]->normal_vector(k))*fe_values_vec_[sid]->JxW(k);
                    fluxes[sid] /= edge_side.measure();
                    if (fluxes[sid] > 0)
                        pflux += fluxes[sid];
                    else
                        nflux += fluxes[sid];
                    ++sid;
                }

                s1=0;
                for( DHCellSide edge_side1 : cell_side.edge_sides() )
                {
                    s2=-1; // need increment at begin of loop (see conditionally 'continue' directions)
                    for( DHCellSide edge_side2 : cell_side.edge_sides() )
                    {
                        s2++;
                        if (s2<=s1) continue;
                        ASSERT(edge_side1.is_valid()).error("Invalid side of edge.");

                        arma::vec3 nv = fe_values_vec_[s1]->normal_vector(0);

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

#define AVERAGE(i,k,side_id)  (fe_values_vec_[sd[side_id]]->shape_value(i,k)*0.5)
#define WAVERAGE(i,k,side_id) (arma::dot(data_->dif_coef_edg[sd[side_id]][sbi][k]*fe_values_vec_[sd[side_id]]->shape_grad(i,k),nv)*omega[side_id])
#define JUMP(i,k,side_id)     ((side_id==0?1:-1)*fe_values_vec_[sd[side_id]]->shape_value(i,k))

                        // For selected pair of elements:
                        for (int n=0; n<2; n++)
                        {
                            if (!is_side_own[n]) continue;

                            for (int m=0; m<2; m++)
                            {
                                for (unsigned int i=0; i<fe_values_vec_[sd[n]]->n_dofs(); i++)
                                    for (unsigned int j=0; j<fe_values_vec_[sd[m]]->n_dofs(); j++)
                                        local_matrix_[i*fe_values_vec_[sd[m]]->n_dofs()+j] = 0;

                                for (unsigned int k=0; k<qsize_lower_dim_; k++)
                                {
                                    double flux_times_JxW = transport_flux*fe_values_vec_[0]->JxW(k);
                                    double gamma_times_JxW = gamma_l*fe_values_vec_[0]->JxW(k);

                                    for (unsigned int i=0; i<fe_values_vec_[sd[n]]->n_dofs(); i++)
                                    {
                                        double flux_JxW_jump_i = flux_times_JxW*JUMP(i,k,n);
                                        double gamma_JxW_jump_i = gamma_times_JxW*JUMP(i,k,n);
                                        double JxW_jump_i = fe_values_vec_[0]->JxW(k)*JUMP(i,k,n);
                                        double JxW_var_wavg_i = fe_values_vec_[0]->JxW(k)*WAVERAGE(i,k,n)*data_->dg_variant;

                                        for (unsigned int j=0; j<fe_values_vec_[sd[m]]->n_dofs(); j++)
                                        {
                                            int index = i*fe_values_vec_[sd[m]]->n_dofs()+j;

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
                                data_->ls[sbi]->mat_set_values(fe_values_vec_[sd[n]]->n_dofs(), &(side_dof_indices_[sd[n]][0]), fe_values_vec_[sd[m]]->n_dofs(), &(side_dof_indices_[sd[m]][0]), &(local_matrix_[0]));
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

    }


    /// Assembles the fluxes between elements of different dimensions.
    void assemble_fluxes_element_side(DHCellAccessor cell_lower_dim) override
    {
        if (dim == 1) return;
        ASSERT_EQ_DBG(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

        // Note: use data members csection_ and velocity_ for appropriate quantities of lower dim element

        for( DHCellSide neighb_side : cell_lower_dim.neighb_sides() )
        {
            // skip neighbours of different dimension
            if (cell_lower_dim.dim() != dim-1) continue;

            ElementAccessor<3> elm_lower_dim = cell_lower_dim.elm();
            n_indices = cell_lower_dim.get_dof_indices(dof_indices_);
            for(unsigned int i=0; i<n_indices; ++i) {
                side_dof_indices_vb_[i] = dof_indices_[i];
            }
            fe_values_vb_->reinit(elm_lower_dim);
            n_dofs[0] = fv_sb_[0]->n_dofs();

            DHCellAccessor cell_higher_dim = data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
            ElementAccessor<3> elm_higher_dim = cell_higher_dim.elm();
            n_indices = cell_higher_dim.get_dof_indices(dof_indices_);
            for(unsigned int i=0; i<n_indices; ++i) {
                side_dof_indices_vb_[i+n_dofs[0]] = dof_indices_[i];
            }
            fe_values_side_.reinit(elm_higher_dim, neighb_side.side_idx());
            n_dofs[1] = fv_sb_[1]->n_dofs();

            // Testing element if they belong to local partition.
            bool own_element_id[2];
            own_element_id[0] = cell_lower_dim.is_own();
            own_element_id[1] = cell_higher_dim.is_own();

            fsv_rt_.reinit(elm_higher_dim, neighb_side.side_idx());
            fv_rt_vb_->reinit(elm_lower_dim);
            calculate_velocity(elm_higher_dim, velocity_higher_, fsv_rt_.point_list());
            calculate_velocity(elm_lower_dim, velocity_, fv_rt_vb_->point_list());
            model_.compute_advection_diffusion_coefficients(fe_values_vb_->point_list(), velocity_, elm_lower_dim, data_->ad_coef_edg[0], data_->dif_coef_edg[0]);
            model_.compute_advection_diffusion_coefficients(fe_values_vb_->point_list(), velocity_higher_, elm_higher_dim, data_->ad_coef_edg[1], data_->dif_coef_edg[1]);
            data_->cross_section.value_list(fe_values_vb_->point_list(), elm_lower_dim, csection_);
            data_->cross_section.value_list(fe_values_vb_->point_list(), elm_higher_dim, csection_higher_);

            for (unsigned int sbi=0; sbi<model_.n_substances(); sbi++) // Optimize: SWAP LOOPS
            {
                for (unsigned int i=0; i<n_dofs[0]+n_dofs[1]; i++)
                    for (unsigned int j=0; j<n_dofs[0]+n_dofs[1]; j++)
                        local_matrix_[i*(n_dofs[0]+n_dofs[1])+j] = 0;

                // sigma_ corresponds to frac_sigma
                data_->fracture_sigma[sbi].value_list(fe_values_vb_->point_list(), elm_lower_dim, sigma_);

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

    }


    /// Assembles the right hand side vector due to volume sources.
    void set_sources(DHCellAccessor cell) override
    {
    	ASSERT_EQ_DBG(cell.dim(), dim).error("Dimension of element mismatch!");

        ElementAccessor<3> elm = cell.elm();

        fe_values_.reinit(elm);
        cell.get_dof_indices(dof_indices_);

        model_.compute_source_coefficients(fe_values_.point_list(), elm, sources_conc_, sources_density_, sources_sigma_);

        // assemble the local stiffness matrix
        for (unsigned int sbi=0; sbi<model_.n_substances(); sbi++)
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
                // temporary, TODO: replace with LocDofVec in balance
                loc_dof_indices_[i] = cell.get_loc_dof_indices()[i];
            }
            model_.balance()->add_source_values(model_.get_subst_idx()[sbi], elm.region().bulk_idx(), loc_dof_indices_,
                                               local_source_balance_vector_, local_source_balance_rhs_);
        }
    }


    void set_boundary_conditions(DHCellAccessor cell) override
    {
        if (cell.elm()->boundary_idx_ == nullptr) return;

        for (unsigned int si=0; si<cell.elm()->n_sides(); si++)
        {
            const Edge *edg = cell.elm().side(si)->edge();
            if (edg->n_sides > 1) continue;
            // skip edges lying not on the boundary
            if (edg->side(0)->cond() == NULL) continue;


            Side side = *(edg->side(0));
            ElementAccessor<3> elm = model_.mesh().element_accessor( side.element().idx() );
            ElementAccessor<3> ele_acc = side.cond()->element_accessor();

            arma::uvec bc_type;
            model_.get_bc_type(ele_acc, bc_type);

            fe_values_side_.reinit(elm, side.side_idx());
            fsv_rt_.reinit(elm, side.side_idx());
            calculate_velocity(elm, velocity_, fsv_rt_.point_list());

            model_.compute_advection_diffusion_coefficients(fe_values_side_.point_list(), velocity_, side.element(), data_->ad_coef, data_->dif_coef);
            data_->cross_section.value_list(fe_values_side_.point_list(), side.element(), csection_);

            DHCellAccessor dh_side_cell = data_->dh_->cell_accessor_from_element( side.element().idx() );
            dh_side_cell.get_dof_indices(dof_indices_);

            for (unsigned int sbi=0; sbi<model_.n_substances(); sbi++)
            {
                fill_n(&(local_rhs_[0]), ndofs_, 0);
                local_flux_balance_vector_.assign(ndofs_, 0);
                local_flux_balance_rhs_ = 0;

                // The b.c. data are fetched for all possible b.c. types since we allow
                // different bc_type for each substance.
                data_->bc_dirichlet_value[sbi].value_list(fe_values_side_.point_list(), ele_acc, bc_values_);

                double side_flux = 0;
                for (unsigned int k=0; k<qsize_lower_dim_; k++)
                    side_flux += arma::dot(data_->ad_coef[sbi][k], fe_values_side_.normal_vector(k))*fe_values_side_.JxW(k);
                double transport_flux = side_flux/side.measure();

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
                        double bc_term = data_->gamma[sbi][side.cond_idx()]*bc_values_[k]*fe_values_side_.JxW(k);
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
                                    + data_->gamma[sbi][side.cond_idx()]*fe_values_side_.shape_value(i,k))*fe_values_side_.JxW(k);
                        }
                    }
                    if (model_.time().tlevel() > 0)
                        for (unsigned int i=0; i<ndofs_; i++)
                            local_flux_balance_rhs_ -= local_rhs_[i];
                }
                else if (bc_type[sbi] == AdvectionDiffusionModel::abc_total_flux)
                {
                	model_.get_flux_bc_data(sbi, fe_values_side_.point_list(), ele_acc, bc_fluxes_, sigma_, bc_ref_values_);
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
                	model_.get_flux_bc_data(sbi, fe_values_side_.point_list(), ele_acc, bc_fluxes_, sigma_, bc_ref_values_);
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

                model_.balance()->add_flux_matrix_values(model_.get_subst_idx()[sbi], side, dof_indices_, local_flux_balance_vector_);
                model_.balance()->add_flux_vec_value(model_.get_subst_idx()[sbi], side, local_flux_balance_rhs_);
            }
        }
    }


	/**
	 * @brief Assembles the auxiliary linear system to calculate the initial solution
	 * as L^2-projection of the prescribed initial condition.
	 */
    void prepare_initial_condition(DHCellAccessor cell) override
    {
        ASSERT_EQ_DBG(cell.dim(), dim).error("Dimension of element mismatch!");

        ElementAccessor<3> elem = cell.elm();
        cell.get_dof_indices(dof_indices_);
        fe_values_.reinit(elem);
        model_.compute_init_cond(fe_values_.point_list(), elem, init_values_);

        for (unsigned int sbi=0; sbi<model_.n_substances(); sbi++)
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



#endif /* FE_VALUE_HANDLER_HH_ */
