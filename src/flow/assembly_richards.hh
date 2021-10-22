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
#include "flow/richards_lmh.hh"
#include "flow/assembly_lmh.hh"
#include "flow/soil_models.hh"
//#include "fem/fe_p.hh"
//#include "fem/fe_values.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fields/field_value_cache.hh"


template <unsigned int dim>
class ReadInitCondAssemblyRichards : public ReadInitCondAssemblyLMH<dim>
{
public:
    typedef typename RichardsLMH::EqFields EqFields;
    typedef typename RichardsLMH::EqData EqData;

    static constexpr const char * name() { return "ReadInitCondAssemblyLMH"; }

    /// Constructor.
    ReadInitCondAssemblyRichards(EqFields *eq_fields, EqData *eq_data)
    : ReadInitCondAssemblyLMH<dim>(eq_fields, eq_data), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = ActiveIntegrals::bulk;
        this->used_fields_ += this->eq_fields_->storativity;
        this->used_fields_ += this->eq_fields_->extra_storativity;
        this->used_fields_ += this->eq_fields_->genuchten_n_exponent;
        this->used_fields_ += this->eq_fields_->genuchten_p_head_scale;
        this->used_fields_ += this->eq_fields_->water_content_residual;
        this->used_fields_ += this->eq_fields_->water_content_saturated;
        this->used_fields_ += this->eq_fields_->conductivity;
    }

    /// Destructor.
    ~ReadInitCondAssemblyRichards() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        //this->balance_ = eq_data_->balance_;
        this->element_cache_map_ = element_cache_map;
        genuchten_on = false;
    }


private:
    void reset_soil_model(const DHCellAccessor& cell, BulkPoint &p) {
        genuchten_on = (this->eq_fields_->genuchten_p_head_scale.field_result({cell.elm().region()}) != result_zeros);
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
    }


    void update_water_content(const DHCellAccessor& cell, BulkPoint &p) override {
        edge_indices_ = cell.cell_with_other_dh(this->eq_data_->dh_cr_.get()).get_loc_dof_indices();
        cr_disc_dofs_ = cell.cell_with_other_dh(this->eq_data_->dh_cr_disc_.get()).get_loc_dof_indices();

        reset_soil_model(cell, p);
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

    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    LocDofVec cr_disc_dofs_;                 ///< Vector of local DOF indices pre-computed on different DOF handlers
    LocDofVec edge_indices_;                 ///< Dofs of discontinuous fields on element edges.
    bool genuchten_on;
    double storativity_;
    double capacity, water_content, phead;

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;
};


#endif /* ASSEMBLY_RICHARDS_HH_ */

