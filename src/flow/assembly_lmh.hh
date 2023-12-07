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
 * @file    assembly_lmh.hh
 * @brief
 */

#ifndef ASSEMBLY_LMH_HH_
#define ASSEMBLY_LMH_HH_

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "flow/darcy_flow_lmh.hh"
#include "fields/field_value_cache.hh"

#include "system/index_types.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "mesh/neighbours.h"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_values_views.hh"
#include "fem/fe_system.hh"
#include "quadrature/quadrature_lib.hh"

#include "la/linsys_PETSC.hh"
#include "la/schur.hh"
#include "la/local_system.hh"
#include "la/local_constraint.hh"

#include "coupling/balance.hh"


template <unsigned int dim>
class ReadInitCondAssemblyLMH : public AssemblyBase<dim>
{
public:
    typedef typename DarcyLMH::EqFields EqFields;
    typedef typename DarcyLMH::EqData EqData;

    static constexpr const char * name() { return "ReadInitCondAssemblyLMH"; }

    /// Constructor.
    ReadInitCondAssemblyLMH(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : AssemblyBase<dim>(0, fe_values), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = ActiveIntegrals::bulk;
        this->used_fields_ += eq_fields_->init_pressure;
    }

    /// Destructor.
    virtual ~ReadInitCondAssemblyLMH() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        ASSERT_EQ(cell.dim(), dim).error("Dimension of element mismatch!");

        l_indices_ = cell.get_loc_dof_indices();
        ASSERT(l_indices_.n_elem == cell.elm().element()->n_sides());

        // set initial condition
        auto p = *( this->bulk_points(element_patch_idx).begin() );
        double init_value = eq_fields_->init_pressure(p);

        for (unsigned int i=0; i<cell.elm()->n_sides(); i++) {
             double init_value_on_edge = init_value / cell.elm().side(i)->edge().n_sides();
             eq_data_->p_edge_solution.add(l_indices_[i], init_value_on_edge);
        }
    }

    /// Implements @p AssemblyBase::end.
    void end() override
    {
        eq_data_->p_edge_solution.ghost_to_local_begin();
        eq_data_->p_edge_solution.ghost_to_local_end();
        eq_data_->p_edge_solution.local_to_ghost_begin();
        eq_data_->p_edge_solution.local_to_ghost_end();
    }



protected:
    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    /// Data objects shared with Flow equation
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Vector of pre-computed local DOF indices
    LocDofVec l_indices_;

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};

template <unsigned int dim>
class MHMatrixAssemblyLMH : public AssemblyBase<dim>
{
public:
    typedef typename DarcyLMH::EqFields EqFields;
    typedef typename DarcyLMH::EqData EqData;

    DECLARE_EXCEPTION( ExcBCNotSupported, << "BC type not supported.\n" );

    static constexpr const char * name() { return "MHMatrixAssemblyLMH"; }

    /// Constructor.
    MHMatrixAssemblyLMH(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : AssemblyBase<dim>(0, fe_values), eq_fields_(eq_fields), eq_data_(eq_data), quad_rt_(dim, 2) {
        this->active_integrals_ = (ActiveIntegrals::bulk | ActiveIntegrals::coupling | ActiveIntegrals::boundary);
        this->used_fields_ += eq_fields_->cross_section;
        this->used_fields_ += eq_fields_->conductivity;
        this->used_fields_ += eq_fields_->anisotropy;
        this->used_fields_ += eq_fields_->sigma;
        this->used_fields_ += eq_fields_->water_source_density;
        this->used_fields_ += eq_fields_->extra_source;
        this->used_fields_ += eq_fields_->storativity;
        this->used_fields_ += eq_fields_->extra_storativity;
        this->used_fields_ += eq_fields_->bc_type;
        this->used_fields_ += eq_fields_->bc_pressure;
        this->used_fields_ += eq_fields_->bc_flux;
        this->used_fields_ += eq_fields_->bc_robin_sigma;
        this->used_fields_ += eq_fields_->bc_switch_pressure;
    }

    /// Destructor.
    virtual ~MHMatrixAssemblyLMH() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        //this->balance_ = eq_data_->balance_;
        this->element_cache_map_ = element_cache_map;

        fe_ = std::make_shared< FE_P_disc<dim> >(0);
        fe_values_side_.initialize(*this->quad_low_, *fe_, update_normal_vectors);

        fe_values_.initialize(quad_rt_, fe_rt_, update_values | update_JxW_values | update_quadrature_points);

        // local numbering of dofs for MH system
        // note: this shortcut supposes that the fe_system is the same on all elements
        // TODO the function should be DiscreteSpace.fe(ElementAccessor)
        // and correct form fe(ad_->dh_->own_range().begin()->elm()) (see documentation of DiscreteSpace::fe)
        auto fe = eq_data_->dh_->ds()->fe()[Dim<dim>{}];
        FESystem<dim>* fe_system = dynamic_cast<FESystem<dim>*>(fe.get());
        eq_data_->loc_side_dofs[dim-1] = fe_system->fe_dofs(0);
        eq_data_->loc_ele_dof[dim-1] = fe_system->fe_dofs(1)[0];
        eq_data_->loc_edge_dofs[dim-1] = fe_system->fe_dofs(2);

        // reconstructed vector for the velocity and pressure
        // length = local Schur complement offset in LocalSystem
        eq_data_->schur_offset_[dim-1] = eq_data_->loc_edge_dofs[dim-1][0];
        reconstructed_solution_.zeros(eq_data_->schur_offset_[dim-1]);
    }


    /**
     * Integral over element.
     *
     * Override in descendants.
     */
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        ASSERT_EQ(cell.dim(), dim).error("Dimension of element mismatch!");

        // evaluation point
        auto p = *( this->bulk_points(element_patch_idx).begin() );
        bulk_local_idx_ = cell.local_idx();

        this->asm_sides(cell, p, eq_fields_->conductivity(p));
        this->asm_element();
        this->asm_source_term_darcy(cell, p);
    }


    /**
     * Assembles between boundary element and corresponding side on bulk element.
     *
     * Override in descendants.
     */
    inline void boundary_side_integral(DHCellSide cell_side)
    {
        ASSERT_EQ(cell_side.dim(), dim).error("Dimension of element mismatch!");
        if (!cell_side.cell().is_own()) return;

        auto p_side = *( this->boundary_points(cell_side).begin() );
        auto p_bdr = p_side.point_bdr(cell_side.cond().element_accessor() );
        ElementAccessor<3> b_ele = cell_side.side().cond().element_accessor(); // ??

        precompute_boundary_side(cell_side, p_side, p_bdr);

    	if (type_==DarcyLMH::EqFields::seepage) {
    	    this->use_dirichlet_switch(cell_side, b_ele, p_bdr);
    	}

        this->boundary_side_integral_in(cell_side, b_ele, p_bdr);
    }


    /**
     * Assembles the fluxes between elements of different dimensions.
     *
     * Common in all descendants.
     */
    inline void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
        if (dim == 1) return;
        ASSERT_EQ(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

        unsigned int neigh_idx = ngh_idx(cell_lower_dim, neighb_side); // TODO use better evaluation of neighbour_idx
        unsigned int loc_dof_higher = (2*(cell_lower_dim.dim()+1) + 1) + neigh_idx; // loc dof of higher ele edge
        bulk_local_idx_ = cell_lower_dim.local_idx();

        // Evaluation points
        auto p_high = *( this->coupling_points(neighb_side).begin() );
        auto p_low = p_high.lower_dim(cell_lower_dim);

        fe_values_side_.reinit(neighb_side.side());
        nv_ = fe_values_side_.normal_vector(0);

        ngh_value_ = eq_fields_->sigma(p_low) *
                        2*eq_fields_->conductivity(p_low) *
                        arma::dot(eq_fields_->anisotropy(p_low)*nv_, nv_) *
						eq_fields_->cross_section(p_high) * // cross-section of higher dim. (2d)
						eq_fields_->cross_section(p_high) /
						eq_fields_->cross_section(p_low) *      // crossection of lower dim.
                        neighb_side.measure();

        eq_data_->loc_system_[bulk_local_idx_].add_value(eq_data_->loc_ele_dof[dim-2], eq_data_->loc_ele_dof[dim-2], -ngh_value_);
        eq_data_->loc_system_[bulk_local_idx_].add_value(eq_data_->loc_ele_dof[dim-2], loc_dof_higher,                ngh_value_);
        eq_data_->loc_system_[bulk_local_idx_].add_value(loc_dof_higher,               eq_data_->loc_ele_dof[dim-2],  ngh_value_);
        eq_data_->loc_system_[bulk_local_idx_].add_value(loc_dof_higher,               loc_dof_higher,               -ngh_value_);

//             // update matrix for weights in BDDCML
//             if ( typeid(*eq_data_->lin_sys) == typeid(LinSys_BDDC) ) {
//                 int ind = eq_data_->loc_system_[bulk_local_idx_].row_dofs[p];
//                // there is -value on diagonal in block C!
//                static_cast<LinSys_BDDC*>(eq_data_->lin_sys)->diagonal_weights_set_value( ind, -value );
//             }
    }


    /// Implements @p AssemblyBase::begin.
    void begin() override
    {
        this->begin_mh_matrix();
    }


    /// Implements @p AssemblyBase::end.
    void end() override
    {
        this->end_mh_matrix();
    }

protected:
    /// Common code of begin method of MH matrix assembly (Darcy and Richards)
    void begin_mh_matrix()
    {
        // DebugOut() << "assembly_mh_matrix \n";
        // set auxiliary flag for switchting Dirichlet like BC
        eq_data_->force_no_neumann_bc = eq_data_->use_steady_assembly_ && (eq_data_->nonlinear_iteration_ == 0);

        eq_data_->reset();

        eq_data_->balance_->start_flux_assembly(eq_data_->water_balance_idx);
        eq_data_->balance_->start_source_assembly(eq_data_->water_balance_idx);
        eq_data_->balance_->start_mass_assembly(eq_data_->water_balance_idx);

        this->set_dofs();
    }


    /// Common code of end method of MH matrix assembly (Darcy and Richards)
    void end_mh_matrix()
    {
	    std::sort(eq_data_->loc_constraint_.begin(), eq_data_->loc_constraint_.end());
        eq_data_->loc_constraint_.emplace_back(uint(-1), 0, 0.0); // add constraint of invalid element to end of vector
        unsigned int i_constr = 0;
        for ( DHCellAccessor dh_cell : eq_data_->dh_cr_->own_range() ) {
            this->set_loc_schur(dh_cell);
            while (eq_data_->loc_constraint_[i_constr].i_element() == dh_cell.local_idx()) {
                this->loc_schur_.set_solution(eq_data_->loc_constraint_[i_constr]);
            	i_constr++;
            }
            eq_data_->loc_system_[dh_cell.local_idx()].compute_schur_complement(
                    eq_data_->schur_offset_[dh_cell.dim()-1], this->loc_schur_, true);

            // for seepage BC, save local system
            if (eq_data_->save_local_system_[dh_cell.local_idx()])
            	eq_data_->seepage_bc_systems[dh_cell.elm_idx()] = eq_data_->loc_system_[dh_cell.local_idx()];

            this->loc_schur_.eliminate_solution();
            eq_data_->lin_sys_schur->set_local_system(this->loc_schur_, eq_data_->dh_cr_->get_local_to_global_map());

            // TODO:
            // if (mortar_assembly)
            //     mortar_assembly->assembly(dh_cell);
        }

        eq_data_->balance_->finish_mass_assembly(eq_data_->water_balance_idx);
        eq_data_->balance_->finish_source_assembly(eq_data_->water_balance_idx);
        eq_data_->balance_->finish_flux_assembly(eq_data_->water_balance_idx);

        eq_data_->loc_constraint_.clear();
    }


    /// Common code of begin method of Reconstruct Schur assembly (Darcy and Richards)
    void begin_reconstruct_schur()
    {
        this->eq_data_->full_solution.zero_entries();
        this->eq_data_->p_edge_solution.local_to_ghost_begin();
        this->eq_data_->p_edge_solution.local_to_ghost_end();
    }


    /// Common code of end method of Reconstruct Schur assembly (Darcy and Richards)
    void end_reconstruct_schur()
    {
        for ( DHCellAccessor dh_cell : this->eq_data_->dh_cr_->own_range() ) {
            this->bulk_local_idx_ = dh_cell.local_idx();
            this->set_loc_schur(dh_cell);
            arma::vec schur_solution = this->eq_data_->p_edge_solution.get_subvec(this->loc_schur_.row_dofs);
            // reconstruct the velocity and pressure
            this->eq_data_->loc_system_[this->bulk_local_idx_].reconstruct_solution_schur(this->eq_data_->schur_offset_[dh_cell.dim()-1],
                                                                                  schur_solution, this->reconstructed_solution_);

            this->reconstructed_solution_ += this->eq_data_->postprocess_solution_[this->bulk_local_idx_];

            this->eq_data_->full_solution.set_subvec(this->eq_data_->loc_system_[this->bulk_local_idx_].row_dofs.head(
                    this->eq_data_->schur_offset_[dh_cell.dim()-1]), this->reconstructed_solution_);
            this->eq_data_->full_solution.set_subvec(this->eq_data_->loc_system_[this->bulk_local_idx_].row_dofs.tail(
                    this->loc_schur_.row_dofs.n_elem), schur_solution);
        }

        this->eq_data_->full_solution.local_to_ghost_begin();
        this->eq_data_->full_solution.local_to_ghost_end();

        eq_data_->loc_constraint_.clear();
    }

    /// Part of cell_integral method, common in all descendants
    inline void asm_sides(const DHCellAccessor& cell, BulkPoint &p, double conductivity)
    {
        auto ele = cell.elm();
        double scale_sides = 1 / eq_fields_->cross_section(p) / conductivity;

        fe_values_.reinit(ele);
        auto velocity = fe_values_.vector_view(0);

        for (unsigned int k=0; k<fe_values_.n_points(); k++)
            for (unsigned int i=0; i<fe_values_.n_dofs(); i++){
                double rhs_val = arma::dot(eq_data_->gravity_vec_, velocity.value(i,k))
                           * fe_values_.JxW(k);
                eq_data_->loc_system_[bulk_local_idx_].add_value(i, rhs_val);

                for (unsigned int j=0; j<fe_values_.n_dofs(); j++){
                    double mat_val =
                        arma::dot( velocity.value(i,k), //TODO: compute anisotropy before
                                   (eq_fields_->anisotropy(p)).i() * velocity.value(j,k)
                                 )
                        * scale_sides * fe_values_.JxW(k);

                    eq_data_->loc_system_[bulk_local_idx_].add_value(i, j, mat_val);
                }
            }

    // assemble matrix for weights in BDDCML
    // approximation to diagonal of
    // S = -C - B*inv(A)*B'
    // as
    // diag(S) ~ - diag(C) - 1./diag(A)
    // the weights form a partition of unity to average a discontinuous solution from neighbouring subdomains
    // to a continuous one
    // it is important to scale the effect - if conductivity is low for one subdomain and high for the other,
    // trust more the one with low conductivity - it will be closer to the truth than an arithmetic average
//            if ( typeid(*ad_->lin_sys) == typeid(LinSys_BDDC) ) {
//                const arma::mat& local_matrix = loc_system_.get_matrix();
//                for(unsigned int i=0; i < ndofs; i++) {
//                    double val_side = local_matrix(i,i);
//                    double val_edge = -1./local_matrix(i,i);
//
//                    unsigned int side_row = loc_system_.row_dofs[loc_side_dofs[i]];
//                    unsigned int edge_row = loc_system_.row_dofs[loc_edge_dofs[i]];
//                    static_cast<LinSys_BDDC*>(ad_->lin_sys)->diagonal_weights_set_value( side_row, val_side );
//                    static_cast<LinSys_BDDC*>(ad_->lin_sys)->diagonal_weights_set_value( edge_row, val_edge );
//                }
//            }
    }

    /// Part of cell_integral method, common in all descendants
    inline void asm_element()
    {
        // set block B, B': element-side, side-element

        for(unsigned int side = 0; side < eq_data_->loc_side_dofs[dim-1].size(); side++){
            eq_data_->loc_system_[bulk_local_idx_].add_value(eq_data_->loc_ele_dof[dim-1], eq_data_->loc_side_dofs[dim-1][side], -1.0);
            eq_data_->loc_system_[bulk_local_idx_].add_value(eq_data_->loc_side_dofs[dim-1][side], eq_data_->loc_ele_dof[dim-1], -1.0);
        }

//            if ( typeid(*ad_->lin_sys) == typeid(LinSys_BDDC) ) {
//                double val_ele =  1.;
//                static_cast<LinSys_BDDC*>(ad_->lin_sys)->
//                                diagonal_weights_set_value( loc_system_.row_dofs[loc_ele_dof], val_ele );
//            }
    }

    /// Part of cell_integral method, specialized in Darcy equation
    inline void asm_source_term_darcy(const DHCellAccessor& cell, BulkPoint &p)
    {
        // compute lumped source
        uint n_sides = cell.elm()->n_sides();
        coef_ = (1.0 / n_sides) * cell.elm().measure() * eq_fields_->cross_section(p);
        source_term_ = coef_ * (eq_fields_->water_source_density(p) + eq_fields_->extra_source(p));

        // in unsteady, compute time term
        time_term_diag_ = 0.0;
        time_term_ = 0.0;
        time_term_rhs_ = 0.0;

        if(! eq_data_->use_steady_assembly_)
        {
            time_term_ = coef_ * (eq_fields_->storativity(p) + eq_fields_->extra_storativity(p));
        }

        const DHCellAccessor cr_cell = cell.cell_with_other_dh(eq_data_->dh_cr_.get());
        auto loc_dof_vec = cr_cell.get_loc_dof_indices();

        for (unsigned int i=0; i<n_sides; i++)
        {
            if(! eq_data_->use_steady_assembly_)
            {
                time_term_diag_ = time_term_ / eq_data_->time_step_;
                time_term_rhs_ = time_term_diag_ * eq_data_->p_edge_solution_previous_time.get(loc_dof_vec[i]);

                eq_data_->balance_->add_mass_values(eq_data_->water_balance_idx, cell,
                                  {eq_data_->loc_system_[bulk_local_idx_].row_dofs[eq_data_->loc_edge_dofs[dim-1][i]]}, {time_term_}, 0);
            }
            else
            {
                // Add zeros explicitely to keep the sparsity pattern.
                // Otherwise Petsc would compress out the zeros in FinishAssembly.
                eq_data_->balance_->add_mass_values(eq_data_->water_balance_idx, cell,
                                  {eq_data_->loc_system_[bulk_local_idx_].row_dofs[eq_data_->loc_edge_dofs[dim-1][i]]}, {0}, 0);
            }

            eq_data_->loc_system_[bulk_local_idx_].add_value(eq_data_->loc_edge_dofs[dim-1][i], eq_data_->loc_edge_dofs[dim-1][i],
                            -time_term_diag_,
                            -source_term_ - time_term_rhs_);

            eq_data_->balance_->add_source_values(eq_data_->water_balance_idx, cell.elm().region().bulk_idx(),
                            {eq_data_->loc_system_[bulk_local_idx_].row_dofs[eq_data_->loc_edge_dofs[dim-1][i]]}, {0}, {source_term_});
        }
    }

    /// Precompute values used in boundary side integral on given DHCellSide
    inline void precompute_boundary_side(DHCellSide &cell_side, BoundaryPoint &p_side, BulkPoint &p_bdr)
    {
        bulk_local_idx_ = cell_side.cell().local_idx();

        cross_section_ = eq_fields_->cross_section(p_side);
        type_ = (DarcyLMH::EqFields::BC_Type)eq_fields_->bc_type(p_bdr);

        sidx_ = cell_side.side_idx();
        side_row_ = eq_data_->loc_side_dofs[dim-1][sidx_];    //local
        edge_row_ = eq_data_->loc_edge_dofs[dim-1][sidx_];    //local
    }

    /// Update BC switch dirichlet in MH matrix assembly if BC type is seepage
    inline void use_dirichlet_switch(DHCellSide &cell_side, const ElementAccessor<3> &b_ele, BulkPoint &p_bdr)
    {
        char & switch_dirichlet = eq_data_->bc_switch_dirichlet[b_ele.idx()];
        if (switch_dirichlet) {
            // check and possibly switch to flux BC
            // The switch raise error on the corresponding edge row.
            // Magnitude of the error is abs(solution_flux - side_flux).

            // try reconstructing local system for seepage BC
            auto reconstr_solution = this->load_local_system(cell_side.cell());
            double side_flux = -eq_fields_->bc_flux(p_bdr) * b_ele.measure() * cross_section_;

            if ( reconstr_solution[side_row_] < side_flux) {
                //DebugOut().fmt("x: {}, to neum, p: {} f: {} -> f: {}\n", b_ele.centre()[0], bc_pressure, reconstructed_solution_[side_row_], side_flux);
                switch_dirichlet = 0;
            }
        } else {
            // check and possibly switch to  pressure BC
            // TODO: What is the appropriate DOF in not local?
            // The switch raise error on the corresponding side row.
            // Magnitude of the error is abs(solution_head - bc_pressure)
            // Since usually K is very large, this error would be much
            // higher then error caused by the inverse switch, this
            // cause that a solution  with the flux violating the
            // flux inequality leading may be accepted, while the error
            // in pressure inequality is always satisfied.

            const DHCellAccessor cr_cell = cell_side.cell().cell_with_other_dh(eq_data_->dh_cr_.get());
            auto loc_dof_vec = cr_cell.get_loc_dof_indices();
            double solution_head = eq_data_->p_edge_solution.get(loc_dof_vec[sidx_]);
            double bc_pressure = eq_fields_->bc_switch_pressure(p_bdr);

            if ( solution_head > bc_pressure) {
                //DebugOut().fmt("x: {}, to dirich, p: {} -> p: {} f: {}\n",b_ele.centre()[0], solution_head, bc_pressure, bc_flux);
                switch_dirichlet=1;
            }
        }
    }

    /**
     * Common code of boundary_side_integral.
     */
    inline void boundary_side_integral_in(DHCellSide cell_side, const ElementAccessor<3> &b_ele, BulkPoint &p_bdr)
    {
        if ( type_ == DarcyLMH::EqFields::none) {
            // homogeneous neumann
        } else if ( type_ == DarcyLMH::EqFields::dirichlet ) {
            eq_data_->loc_constraint_.emplace_back( bulk_local_idx_, sidx_, eq_fields_->bc_pressure(p_bdr) );

        } else if ( type_ == DarcyLMH::EqFields::total_flux) {
            // internally we work with outward flux
            eq_data_->loc_system_[bulk_local_idx_].add_value(edge_row_, edge_row_,
                    -b_ele.measure() * eq_fields_->bc_robin_sigma(p_bdr) * cross_section_,
                    (-eq_fields_->bc_flux(p_bdr) - eq_fields_->bc_robin_sigma(p_bdr) * eq_fields_->bc_pressure(p_bdr)) * b_ele.measure() * cross_section_);
        } else if (type_==DarcyLMH::EqFields::seepage) {
            eq_data_->is_linear=false;

            double bc_pressure = eq_fields_->bc_switch_pressure(p_bdr);
            double side_flux = -eq_fields_->bc_flux(p_bdr) * b_ele.measure() * cross_section_;

            eq_data_->save_local_system_[bulk_local_idx_] = (bool) eq_data_->bc_switch_dirichlet[b_ele.idx()];

            // ** Apply BCUpdate BC type. **
            // Force Dirichlet type during the first iteration of the unsteady case.
            if (eq_data_->save_local_system_[bulk_local_idx_] || eq_data_->force_no_neumann_bc ) {
                //DebugOut().fmt("x: {}, dirich, bcp: {}\n", b_ele.centre()[0], bc_pressure);
                eq_data_->loc_constraint_.emplace_back(bulk_local_idx_, sidx_, bc_pressure);
            } else {
                //DebugOut()("x: {}, neuman, q: {}  bcq: {}\n", b_ele.centre()[0], side_flux, bc_flux);
                eq_data_->loc_system_[bulk_local_idx_].add_value(edge_row_, side_flux);
            }

        } else if (type_==DarcyLMH::EqFields::river) {
            eq_data_->is_linear=false;

            double bc_pressure = eq_fields_->bc_pressure(p_bdr);
            double bc_switch_pressure = eq_fields_->bc_switch_pressure(p_bdr);
            double bc_flux = -eq_fields_->bc_flux(p_bdr);
            double bc_sigma = eq_fields_->bc_robin_sigma(p_bdr);

            const DHCellAccessor cr_cell = cell_side.cell().cell_with_other_dh(eq_data_->dh_cr_.get());
            auto loc_dof_vec = cr_cell.get_loc_dof_indices();
            double solution_head = eq_data_->p_edge_solution.get(loc_dof_vec[sidx_]);

            // Force Robin type during the first iteration of the unsteady case.
            if (solution_head > bc_switch_pressure  || eq_data_->force_no_neumann_bc) {
                // Robin BC
                //DebugOut().fmt("x: {}, robin, bcp: {}\n", b_ele.centre()[0], bc_pressure);
                eq_data_->loc_system_[bulk_local_idx_].add_value(edge_row_, edge_row_,
                                        -b_ele.measure() * bc_sigma * cross_section_,
                                        b_ele.measure() * cross_section_ * (bc_flux - bc_sigma * bc_pressure)  );
            } else {
                // Neumann BC
                //DebugOut().fmt("x: {}, neuman, q: {}  bcq: {}\n", b_ele.centre()[0], bc_switch_pressure, bc_pressure);
                double bc_total_flux = bc_flux + bc_sigma*(bc_switch_pressure - bc_pressure);

                eq_data_->loc_system_[bulk_local_idx_].add_value(edge_row_, bc_total_flux * b_ele.measure() * cross_section_);
            }
        }
        else {
            THROW( ExcBCNotSupported() );
        }

        eq_data_->balance_->add_flux_values(eq_data_->water_balance_idx, cell_side,
                                      {eq_data_->loc_system_[bulk_local_idx_].row_dofs[eq_data_->loc_side_dofs[dim-1][sidx_]]},
                                      {1}, 0);
    }

    /// Precompute loc_system and loc_schur data members.
    void set_dofs() {
        unsigned int size, loc_size, elm_dim;;
        for ( DHCellAccessor dh_cell : eq_data_->dh_->local_range() ) {
            const ElementAccessor<3> ele = dh_cell.elm();

            elm_dim = ele.dim();
            size = (elm_dim+1) + 1 + (elm_dim+1); // = n_sides + 1 + n_sides
            loc_size = size + ele->n_neighs_vb();
            LocDofVec dofs(loc_size);

            // add full vec indices
            dofs.head(dh_cell.n_dofs()) = dh_cell.get_loc_dof_indices();

            if(ele->n_neighs_vb() != 0)
            {
                //D, E',E block: compatible connections: element-edge
                unsigned int i = 0;

                for ( DHCellSide neighb_side : dh_cell.neighb_sides() ) {

                    // read neighbor dofs (dh dofhandler)
                    // neighbor cell owning neighb_side
                    DHCellAccessor dh_neighb_cell = neighb_side.cell();

                    // read neighbor dofs (dh_cr dofhandler)
                    // neighbor cell owning neighb_side
                    DHCellAccessor dh_neighb_cr_cell = dh_neighb_cell.cell_with_other_dh(eq_data_->dh_cr_.get());

                    // local index of pedge dof in local system
                    const unsigned int p = size+i;
                    // local index of pedge dof on neighboring cell
                    const unsigned int t = dh_neighb_cell.n_dofs() - dh_neighb_cr_cell.n_dofs() + neighb_side.side().side_idx();
                    dofs[p] = dh_neighb_cell.get_loc_dof_indices()[t];
                    i++;
                }
            }
            eq_data_->loc_system_[dh_cell.local_idx()].reset(dofs, dofs);

            // Set side-edge (flux-lambda) terms
            for (DHCellSide dh_side : dh_cell.side_range()) {
                unsigned int sidx = dh_side.side_idx();
                // side-edge (flux-lambda) terms
                eq_data_->loc_system_[dh_cell.local_idx()].add_value(eq_data_->loc_side_dofs[elm_dim-1][sidx], eq_data_->loc_edge_dofs[elm_dim-1][sidx], 1.0);
                eq_data_->loc_system_[dh_cell.local_idx()].add_value(eq_data_->loc_edge_dofs[elm_dim-1][sidx], eq_data_->loc_side_dofs[elm_dim-1][sidx], 1.0);
            }
        }
    }


    /// Precompute loc_schur data member of given cell.
    void set_loc_schur(const DHCellAccessor dh_cr_cell) {
        const ElementAccessor<3> ele = dh_cr_cell.elm();

        unsigned int loc_size_schur = ele->n_sides() + ele->n_neighs_vb();
        LocDofVec dofs_schur(loc_size_schur);

        // add schur vec indices
        dofs_schur.head(dh_cr_cell.n_dofs()) = dh_cr_cell.get_loc_dof_indices();

        if(ele->n_neighs_vb() != 0)
        {
            //D, E',E block: compatible connections: element-edge
            unsigned int i = 0;

            for ( DHCellSide neighb_side : dh_cr_cell.neighb_sides() ) {

                // read neighbor dofs (dh dofhandler)
                // neighbor cell owning neighb_side
                DHCellAccessor dh_neighb_cr_cell = neighb_side.cell();

                // local index of pedge dof in local schur system
                const unsigned int tt = dh_cr_cell.n_dofs()+i;
                dofs_schur[tt] = dh_neighb_cr_cell.get_loc_dof_indices()[neighb_side.side().side_idx()];
                i++;
            }
        }


        loc_schur_.reset(dofs_schur, dofs_schur);
    }


    /// Temporary method find neighbour index in higher-dim cell
    inline unsigned int ngh_idx(DHCellAccessor &dh_cell, DHCellSide &neighb_side) {
        for (uint n_i=0; n_i<dh_cell.elm()->n_neighs_vb(); ++n_i) {
            auto side = dh_cell.elm()->neigh_vb[n_i]->side();
            if ( (side->elem_idx() == neighb_side.elem_idx()) && (side->side_idx() == neighb_side.side_idx()) ) return n_i;
        }
        ASSERT_PERMANENT(false)(dh_cell.elm_idx())(neighb_side.side_idx()).error("Side is not a part of neighbour!\n");
        return 0;
    }


    /** Loads the local system from a map: element index -> LocalSystem,
     * if it exits, or if the full solution is not yet reconstructed,
     * and reconstructs the full solution on the element.
     * Currently used only for seepage BC.
     */
    arma::vec load_local_system(const DHCellAccessor& dh_cell)
    {
        // possibly find the corresponding local system
        auto ls = eq_data_->seepage_bc_systems.find(dh_cell.elm_idx());
        if (ls != eq_data_->seepage_bc_systems.end())
        {
            const DHCellAccessor cr_cell = dh_cell.cell_with_other_dh(eq_data_->dh_cr_.get());
            auto loc_dof_vec = cr_cell.get_loc_dof_indices();
            arma::vec schur_solution = eq_data_->p_edge_solution.get_subvec(loc_dof_vec);
            // reconstruct the velocity and pressure
            ls->second.reconstruct_solution_schur(eq_data_->schur_offset_[dim-1], schur_solution, reconstructed_solution_);

        	unsigned int pos_in_cache = this->element_cache_map_->position_in_cache(dh_cell.elm_idx());
        	auto p = *( this->bulk_points(pos_in_cache).begin() );
            postprocess_velocity_darcy(dh_cell, p, reconstructed_solution_);

            eq_data_->bc_fluxes_reconstruted[bulk_local_idx_] = true;
        }

        return reconstructed_solution_;
    }


    /**
     * Precompute edge_scale and edge_source_term.
     *
     * This method must be calls in methods postprocess_velocity_darcy and postprocess_velocity_richards
     */
    void postprocess_velocity(const DHCellAccessor& dh_cell, BulkPoint &p)
    {
        edge_scale_ = dh_cell.elm().measure()
                          * eq_fields_->cross_section(p)
                          / dh_cell.elm()->n_sides();

        edge_source_term_ = edge_scale_ *
                ( eq_fields_->water_source_density(p)
                + eq_fields_->extra_source(p));
    }


    /// Postprocess velocity during loading of local system and after calculating of cell integral.
    void postprocess_velocity_darcy(const DHCellAccessor& dh_cell, BulkPoint &p, arma::vec& solution)
    {
    	postprocess_velocity(dh_cell, p);

        time_term_ = 0.0;
        const DHCellAccessor cr_cell = dh_cell.cell_with_other_dh(eq_data_->dh_cr_.get());
        auto loc_dof_vec = cr_cell.get_loc_dof_indices();

        for (unsigned int i=0; i<dh_cell.elm()->n_sides(); i++) {

            if( ! eq_data_->use_steady_assembly_)
            {
                double new_pressure = eq_data_->p_edge_solution.get(loc_dof_vec[i]);
                double old_pressure = eq_data_->p_edge_solution_previous_time.get(loc_dof_vec[i]);
                time_term_ = edge_scale_ * (eq_fields_->storativity(p) + eq_fields_->extra_storativity(p))
                             / eq_data_->time_step_ * (new_pressure - old_pressure);
            }
            solution[eq_data_->loc_side_dofs[dim-1][i]] += edge_source_term_ - time_term_;
        }
    }


    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    /// Data objects shared with DarcyFlow
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Assembly volume integrals
    FE_RT0<dim> fe_rt_;
    QGauss quad_rt_;
    FEValues<3> fe_values_;

    shared_ptr<FiniteElement<dim>> fe_;                    ///< Finite element for the solution of the advection-diffusion equation.
    FEValues<3> fe_values_side_;                           ///< FEValues of object (of P disc finite element type)

    /// Vector for reconstructed solution (velocity and pressure on element) from Schur complement.
    arma::vec reconstructed_solution_;

    double coef_, source_term_;                            ///< Variables used in compute lumped source.
    double time_term_diag_, time_term_, time_term_rhs_;    ///< Variables used in compute time term (unsteady)
    double cross_section_;                                 ///< Precomputed cross-section value
    unsigned int bulk_local_idx_;                          ///< Local idx of bulk element
    unsigned int sidx_, side_row_, edge_row_;              ///< Helper indices in boundary assembly
    DarcyLMH::EqFields::BC_Type type_;                     ///< Type of boundary condition
    arma::vec3 nv_;                                        ///< Normal vector
    double ngh_value_;                                     ///< Precomputed ngh value
    double edge_scale_, edge_source_term_;                 ///< Precomputed values in postprocess_velocity

    LocalSystem loc_schur_;

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};

template <unsigned int dim>
class ReconstructSchurAssemblyLMH : public MHMatrixAssemblyLMH<dim>
{
public:
    typedef typename DarcyLMH::EqFields EqFields;
    typedef typename DarcyLMH::EqData EqData;

    static constexpr const char * name() { return "ReconstructSchurAssemblyLMH"; }

    /// Constructor.
    ReconstructSchurAssemblyLMH(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : MHMatrixAssemblyLMH<dim>(eq_fields, eq_data, fe_values) {}

    /// Integral over element.
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        ASSERT_EQ(cell.dim(), dim).error("Dimension of element mismatch!");

        // evaluation point
        auto p = *( this->bulk_points(element_patch_idx).begin() );
        this->bulk_local_idx_ = cell.local_idx();

        { // postprocess the velocity
            this->eq_data_->postprocess_solution_[this->bulk_local_idx_].zeros(this->eq_data_->schur_offset_[dim-1]);
            this->postprocess_velocity_darcy(cell, p, this->eq_data_->postprocess_solution_[this->bulk_local_idx_]);
        }
    }


    /// Assembles between boundary element and corresponding side on bulk element.
    inline void boundary_side_integral(FMT_UNUSED DHCellSide cell_side)
    {}

    inline void dimjoin_intergral(FMT_UNUSED DHCellAccessor cell_lower_dim, FMT_UNUSED DHCellSide neighb_side)
    {}


    /// Implements @p AssemblyBase::begin.
    void begin() override
    {
        this->begin_reconstruct_schur();
    }


    /// Implements @p AssemblyBase::end.
    void end() override
    {
        this->end_reconstruct_schur();
    }

};

#endif /* ASSEMBLY_LMH_HH_ */

