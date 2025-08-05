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
 * @file    assembly_richards.hh
 * @brief
 */

#ifndef ASSEMBLY_RICHARDS_HH_
#define ASSEMBLY_RICHARDS_HH_

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "flow/assembly_lmh.hh"
#include "flow/soil_models.hh"
#include "fem/element_cache_map.hh"


template <unsigned int dim>
class InitCondPostprocessAssembly : public AssemblyBase<dim>
{
public:
    typedef typename RichardsLMH::EqFields EqFields;
    typedef typename RichardsLMH::EqData EqData;

    static constexpr const char * name() { return "InitCondPostprocessAssembly"; }

    /// Constructor.
    InitCondPostprocessAssembly(EqFields *eq_fields, EqData *eq_data, AssemblyInternals *asm_internals)
    : AssemblyBase<dim>(0, asm_internals), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = ActiveIntegrals::bulk;
        this->used_fields_ += this->eq_fields_->storativity;
        this->used_fields_ += this->eq_fields_->extra_storativity;
        this->used_fields_ += this->eq_fields_->genuchten_n_exponent;
        this->used_fields_ += this->eq_fields_->genuchten_p_head_scale;
        this->used_fields_ += this->eq_fields_->water_content_residual;
        this->used_fields_ += this->eq_fields_->water_content_saturated;
        this->used_fields_ += this->eq_fields_->conductivity;
        this->used_fields_ += this->eq_fields_->cross_section;
    }

    /// Destructor.
    ~InitCondPostprocessAssembly() {}

    /// Initialize auxiliary vectors and other data members
    void initialize() {}

    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx) {
        ASSERT_EQ(cell.dim(), dim).error("Dimension of element mismatch!");

        edge_indices_ = cell.get_loc_dof_indices();
        cr_disc_dofs_ = cell.cell_with_other_dh(this->eq_data_->dh_cr_disc_.get()).get_loc_dof_indices();
        const DHCellAccessor dh_cell = cell.cell_with_other_dh(this->eq_data_->dh_.get());

        auto p = *( this->bulk_points(element_patch_idx).begin() );
        bool genuchten_on = reset_soil_model(cell, p);
        storativity_ = this->eq_fields_->storativity(p)
                         + this->eq_fields_->extra_storativity(p);
        VectorMPI water_content_vec = this->eq_fields_->water_content_ptr->vec();
        double diagonal_coef = cell.elm().measure() * eq_fields_->cross_section(p) / cell.elm()->n_sides();

        for (unsigned int i=0; i<cell.elm()->n_sides(); i++) {
            const int local_side = cr_disc_dofs_[i];
            capacity = 0;
            water_content = 0;
            phead = this->eq_data_->p_edge_solution.get( edge_indices_[i] );

            if (genuchten_on) {
                fadbad::B<double> x_phead(phead);
                fadbad::B<double> evaluated( this->eq_data_->soil_model_->water_content_diff(x_phead) );
                evaluated.diff(0,1);
                water_content = evaluated.val();
                capacity = x_phead.d(0);
            }
            this->eq_data_->capacity.set( cr_disc_dofs_[i], capacity + storativity_ );
            water_content_vec.set( cr_disc_dofs_[i], water_content + storativity_ * phead);

            this->eq_data_->balance_->add_mass_values(eq_data_->water_balance_idx, dh_cell, {local_side},
                                               {0.0}, diagonal_coef*(water_content + storativity_ * phead) );
        }
    }

    /// Implements @p AssemblyBase::begin.
    void begin() override
    {
        this->eq_data_->balance_->start_mass_assembly(this->eq_data_->water_balance_idx);
    }


    /// Implements @p AssemblyBase::end.
    void end() override
    {
        this->eq_data_->balance_->finish_mass_assembly(this->eq_data_->water_balance_idx);
    }


private:
    bool reset_soil_model(const DHCellAccessor& cell, BulkPoint &p) {
        bool genuchten_on = (this->eq_fields_->genuchten_p_head_scale.field_result({cell.elm().region()}) != result_zeros);
        if (genuchten_on) {
            SoilData soil_data;
            soil_data.n     = this->eq_fields_->genuchten_n_exponent(p);
            soil_data.alpha = this->eq_fields_->genuchten_p_head_scale(p);
            soil_data.Qr    = this->eq_fields_->water_content_residual(p);
            soil_data.Qs    = this->eq_fields_->water_content_saturated(p);
            soil_data.Ks    = this->eq_fields_->conductivity(p);
            //soil_data.cut_fraction = 0.99; // set by model

            this->eq_data_->soil_model_->reset(soil_data);
        }
        return genuchten_on;
    }


    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    /// Data objects shared with ConvectionTransport
    EqFields *eq_fields_;
    EqData *eq_data_;

    LocDofVec cr_disc_dofs_;                 ///< Vector of local DOF indices pre-computed on different DOF handlers
    LocDofVec edge_indices_;                 ///< Dofs of discontinuous fields on element edges.
    double storativity_;
    double capacity, water_content, phead;

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;
};


template <unsigned int dim>
class MHMatrixAssemblyRichards : public MHMatrixAssemblyLMH<dim>
{
public:
    typedef typename RichardsLMH::EqFields EqFields;
    typedef typename RichardsLMH::EqData EqData;

    static constexpr const char * name() { return "MHMatrixAssemblyRichards"; }

    MHMatrixAssemblyRichards(EqFields *eq_fields, EqData *eq_data, AssemblyInternals *asm_internals)
    : MHMatrixAssemblyLMH<dim>(eq_fields, eq_data, asm_internals), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = (ActiveIntegrals::bulk | ActiveIntegrals::coupling | ActiveIntegrals::boundary);
        this->used_fields_ += eq_fields_->cross_section;
        this->used_fields_ += eq_fields_->conductivity;
        this->used_fields_ += eq_fields_->anisotropy;
        this->used_fields_ += eq_fields_->sigma;
        this->used_fields_ += eq_fields_->water_source_density;
        this->used_fields_ += eq_fields_->extra_source;
        this->used_fields_ += eq_fields_->storativity;
        this->used_fields_ += eq_fields_->extra_storativity;
        this->used_fields_ += eq_fields_->genuchten_n_exponent;
        this->used_fields_ += eq_fields_->genuchten_p_head_scale;
        this->used_fields_ += eq_fields_->water_content_residual;
        this->used_fields_ += eq_fields_->water_content_saturated;
        this->used_fields_ += eq_fields_->bc_type;
        this->used_fields_ += eq_fields_->bc_pressure;
        this->used_fields_ += eq_fields_->bc_flux;
        this->used_fields_ += eq_fields_->bc_pressure;
        this->used_fields_ += eq_fields_->bc_robin_sigma;
        this->used_fields_ += eq_fields_->bc_switch_pressure;
    }

    /// Destructor.
    ~MHMatrixAssemblyRichards() {}

    /// Initialize auxiliary vectors and other data members
    void initialize() {
        //this->balance_ = eq_data_->balance_;
    	MHMatrixAssemblyLMH<dim>::initialize();
    }


    /// Integral over element.
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        ASSERT_EQ(cell.dim(), dim).error("Dimension of element mismatch!");

        // evaluation point
        auto p = *( this->bulk_points(element_patch_idx).begin() );
        this->bulk_local_idx_ = cell.local_idx();

        this->asm_sides(cell, p, this->compute_conductivity(cell, p));
        this->asm_element();
        this->asm_source_term_richards(cell, p);
    }


    /// Assembles between boundary element and corresponding side on bulk element.
    inline void boundary_side_integral(DHCellSide cell_side)
    {
        ASSERT_EQ(cell_side.dim(), dim).error("Dimension of element mismatch!");
        if (!cell_side.cell().is_own()) return;

        auto p_side = *( this->boundary_points(cell_side).begin() );
        auto p_bdr = p_side.point_bdr(cell_side.cond().element_accessor() );
        ElementAccessor<3> b_ele = cell_side.side().cond().element_accessor(); // ??

        this->precompute_boundary_side(cell_side, p_side, p_bdr);

    	if (this->type_==DarcyLMH::EqFields::seepage) {
    	    this->use_dirichlet_switch(cell_side, b_ele, p_bdr);
    	}

        this->boundary_side_integral_in(cell_side, b_ele, p_bdr);
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
    /// Part of cell_integral method, specialized in Richards equation
    inline void asm_source_term_richards(const DHCellAccessor& cell, BulkPoint &p)
    {
        update_water_content(cell, p);
        const ElementAccessor<3> ele = cell.elm();

        // set lumped source
        diagonal_coef_ = ele.measure() * eq_fields_->cross_section(p) / ele->n_sides();
        source_diagonal_ = diagonal_coef_ * ( eq_fields_->water_source_density(p) + eq_fields_->extra_source(p));

        VectorMPI water_content_vec = eq_fields_->water_content_ptr->vec();

        const DHCellAccessor cr_cell = cell.cell_with_other_dh(eq_data_->dh_cr_.get());
        auto loc_dof_vec = cr_cell.get_loc_dof_indices();

        for (unsigned int i=0; i<ele->n_sides(); i++)
        {

            const int local_side = cr_disc_dofs_[i];
            /*if (this->dirichlet_edge[i] == 0)*/ { //TODO Fix condition evaluating dirichlet_edge

                water_content_diff_ = -water_content_vec.get(local_side) + eq_data_->water_content_previous_time.get(local_side);
                mass_diagonal_ = diagonal_coef_ * eq_data_->capacity.get(local_side);

                /*
                DebugOut().fmt("w diff: {:g}  mass: {:g} w prev: {:g} w time: {:g} c: {:g} p: {:g} z: {:g}",
                      water_content_diff,
                      mass_diagonal * eq_data_->p_edge_solution[local_edge],
                     -eq_data_->water_content_previous_it[local_side],
                      eq_data_->water_content_previous_time[local_side],
                      capacity,
                      eq_data_->p_edge_solution[local_edge],
                      ele.centre()[0] );
                */


                mass_rhs_ = mass_diagonal_ * eq_data_->p_edge_solution.get( loc_dof_vec[i] ) / eq_data_->time_step_
                                  + diagonal_coef_ * water_content_diff_ / eq_data_->time_step_;

                /*
                DBGCOUT(<< "source [" << loc_system_.row_dofs[this->loc_edge_dofs[i]] << ", " << loc_system_.row_dofs[this->loc_edge_dofs[i]]
                             << "] mat: " << -mass_diagonal/this->eq_data_->time_step_
                             << " rhs: " << -source_diagonal_ - mass_rhs
                             << "\n");
                */
                eq_data_->loc_system_[cell.local_idx()].add_value(eq_data_->loc_edge_dofs[dim-1][i], eq_data_->loc_edge_dofs[dim-1][i],
                                            -mass_diagonal_/eq_data_->time_step_,
                                            -source_diagonal_ - mass_rhs_);
            }

            eq_data_->balance_->add_mass_values(eq_data_->water_balance_idx, cell, {local_side},
                                               {0.0}, diagonal_coef_*water_content_vec.get(local_side));
            eq_data_->balance_->add_source_values(eq_data_->water_balance_idx, ele.region().bulk_idx(),
                                                {this->eq_data_->loc_system_[cell.local_idx()].row_dofs[eq_data_->loc_edge_dofs[dim-1][i]]},
                                                {0},{source_diagonal_});
        }
    }

    bool reset_soil_model(const DHCellAccessor& cell, BulkPoint &p) {
        bool genuchten_on = (this->eq_fields_->genuchten_p_head_scale.field_result({cell.elm().region()}) != result_zeros);
        if (genuchten_on) {
            SoilData soil_data;
            soil_data.n     = this->eq_fields_->genuchten_n_exponent(p);
            soil_data.alpha = this->eq_fields_->genuchten_p_head_scale(p);
            soil_data.Qr    = this->eq_fields_->water_content_residual(p);
            soil_data.Qs    = this->eq_fields_->water_content_saturated(p);
            soil_data.Ks    = this->eq_fields_->conductivity(p);
            //soil_data.cut_fraction = 0.99; // set by model

            this->eq_data_->soil_model_->reset(soil_data);
        }
        return genuchten_on;
    }


    void update_water_content(const DHCellAccessor& cell, BulkPoint &p) {
        edge_indices_ = cell.cell_with_other_dh(this->eq_data_->dh_cr_.get()).get_loc_dof_indices();
        cr_disc_dofs_ = cell.cell_with_other_dh(this->eq_data_->dh_cr_disc_.get()).get_loc_dof_indices();

        bool genuchten_on = reset_soil_model(cell, p);
        storativity_ = this->eq_fields_->storativity(p)
                         + this->eq_fields_->extra_storativity(p);
        VectorMPI water_content_vec = this->eq_fields_->water_content_ptr->vec();

        for (unsigned int i=0; i<cell.elm()->n_sides(); i++) {
            capacity = 0;
            water_content = 0;
            phead = this->eq_data_->p_edge_solution.get( edge_indices_[i] );

            if (genuchten_on) {
                fadbad::B<double> x_phead(phead);
                fadbad::B<double> evaluated( this->eq_data_->soil_model_->water_content_diff(x_phead) );
                evaluated.diff(0,1);
                water_content = evaluated.val();
                capacity = x_phead.d(0);
            }
            this->eq_data_->capacity.set( cr_disc_dofs_[i], capacity + storativity_ );
            water_content_vec.set( cr_disc_dofs_[i], water_content + storativity_ * phead);
        }
    }

    /// Precompute conductivity on bulk point.
    double compute_conductivity(const DHCellAccessor& cell, BulkPoint &p)
    {
        bool genuchten_on = reset_soil_model(cell, p);

        double conductivity = 0;
        if (genuchten_on) {
            const DHCellAccessor cr_cell = cell.cell_with_other_dh(eq_data_->dh_cr_.get());
            auto loc_dof_vec = cr_cell.get_loc_dof_indices();

            for (unsigned int i=0; i<cell.elm()->n_sides(); i++)
            {
                double phead = eq_data_->p_edge_solution.get( loc_dof_vec[i] );
                conductivity += eq_data_->soil_model_->conductivity(phead);
            }
            conductivity /= cell.elm()->n_sides();
        }
        else {
            conductivity = eq_fields_->conductivity(p);
        }
        return conductivity;
    }


    /// Postprocess velocity after calculating of cell integral.
    void postprocess_velocity_richards(const DHCellAccessor& dh_cell, BulkPoint &p, arma::vec& solution)
    {
        this->postprocess_velocity(dh_cell, p);

        this->update_water_content(dh_cell, p);

        VectorMPI water_content_vec = eq_fields_->water_content_ptr->vec();

        for (unsigned int i=0; i<dh_cell.elm()->n_sides(); i++) {
            water_content = water_content_vec.get( this->cr_disc_dofs_[i] );
            water_content_previous_time = eq_data_->water_content_previous_time.get( this->cr_disc_dofs_[i] );

            solution[eq_data_->loc_side_dofs[dim-1][i]]
                += this->edge_source_term_ - this->edge_scale_ * (water_content - water_content_previous_time) / eq_data_->time_step_;
        }

        IntIdx p_dof = dh_cell.cell_with_other_dh(eq_data_->dh_p_.get()).get_loc_dof_indices()(0);
        eq_fields_->conductivity_ptr->vec().set( p_dof, compute_conductivity(dh_cell, p) );
    }


    /// Data objects shared with ConvectionTransport
    EqFields *eq_fields_;
    EqData *eq_data_;

    LocDofVec cr_disc_dofs_;                       ///< Vector of local DOF indices pre-computed on different DOF handlers
    LocDofVec edge_indices_;                       ///< Dofs of discontinuous fields on element edges.
    double storativity_;
    double capacity, phead;
    double water_content, water_content_previous_time;
    double diagonal_coef_, source_diagonal_;
    double water_content_diff_, mass_diagonal_, mass_rhs_;

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;
};


template <unsigned int dim>
class ReconstructSchurAssemblyRichards : public MHMatrixAssemblyRichards<dim>
{
public:
    typedef typename RichardsLMH::EqFields EqFields;
    typedef typename RichardsLMH::EqData EqData;

    static constexpr const char * name() { return "ReconstructSchurAssemblyRichards"; }

    ReconstructSchurAssemblyRichards(EqFields *eq_fields, EqData *eq_data, AssemblyInternals *asm_internals)
    : MHMatrixAssemblyRichards<dim>(eq_fields, eq_data, asm_internals) {
    }

    /// Integral over element.
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        ASSERT_EQ(cell.dim(), dim).error("Dimension of element mismatch!");

        // evaluation point
        auto p = *( this->bulk_points(element_patch_idx).begin() );
        this->bulk_local_idx_ = cell.local_idx();

        { // postprocess the velocity
            this->eq_data_->postprocess_solution_[this->bulk_local_idx_].zeros(this->eq_data_->schur_offset_[dim-1]);
            this->postprocess_velocity_richards(cell, p, this->eq_data_->postprocess_solution_[this->bulk_local_idx_]);
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
protected:
    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;
};



#endif /* ASSEMBLY_RICHARDS_HH_ */

