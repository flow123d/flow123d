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
 * @file    assembly_reaction.hh
 * @brief
 */

#ifndef ASSEMBLY_REACTION_HH_
#define ASSEMBLY_REACTION_HH_

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "reaction/dual_porosity.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "quadrature/quadrature_lib.hh"
//#include "coupling/balance.hh"
#include "fields/field_value_cache.hh"



template <unsigned int dim>
class InitConditionAssemblyDp : public AssemblyBase<dim>
{
public:
    typedef typename DualPorosity::EqFields EqFields;
    typedef typename DualPorosity::EqData EqData;

    static constexpr const char * name() { return "InitConditionAssemblyDp"; }

    /// Constructor.
    InitConditionAssemblyDp(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(0), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = ActiveIntegrals::bulk;
        this->used_fields_ += eq_fields_->init_conc_immobile;
    }

    /// Destructor.
    ~InitConditionAssemblyDp() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        //this->balance_ = eq_data_->balance_;
        this->element_cache_map_ = element_cache_map;
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        ASSERT_EQ_DBG(cell.dim(), dim).error("Dimension of element mismatch!");

        dof_p0_ = cell.get_loc_dof_indices()[0];
        auto p = *( this->bulk_points(element_patch_idx).begin() );

        //setting initial solid concentration for substances involved in adsorption
        for (unsigned int sbi = 0; sbi < eq_data_->substances_.size(); sbi++)
        {
            eq_fields_->conc_immobile_fe[sbi]->vec().set( dof_p0_, eq_fields_->init_conc_immobile[sbi](p) );
        }
    }


private:
    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    IntIdx dof_p0_;                                     ///< Index of local DOF

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;
};

template <unsigned int dim>
class ReactionAssemblyDp : public AssemblyBase<dim>
{
public:
    typedef typename DualPorosity::EqFields EqFields;
    typedef typename DualPorosity::EqData EqData;

    static constexpr const char * name() { return "InitConditionAssemblyDp"; }

    /// Constructor.
    ReactionAssemblyDp(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(0), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = ActiveIntegrals::bulk;
        this->used_fields_ += eq_fields_->porosity;
        this->used_fields_ += eq_fields_->porosity_immobile;
        this->used_fields_ += eq_fields_->diffusion_rate_immobile;
    }

    /// Destructor.
    ~ReactionAssemblyDp() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        //this->balance_ = eq_data_->balance_;
        this->element_cache_map_ = element_cache_map;
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        ASSERT_EQ_DBG(cell.dim(), dim).error("Dimension of element mismatch!");

        auto p = *( this->bulk_points(element_patch_idx).begin() );

        // if porosity_immobile == 0 then mobile concentration stays the same
        // and immobile concentration cannot change
        por_immob_ = eq_fields_->porosity_immobile(p);
        if (por_immob_ == 0.0) return;

        // get data from fields
        dof_p0_ = cell.get_loc_dof_indices()[0];
        por_mob_ = eq_fields_->porosity(p);
        arma::Col<double> diff_vec(eq_data_->substances_.size());
        for (sbi_=0; sbi_<eq_data_->substances_.size(); sbi_++) // Optimize: SWAP LOOPS
            diff_vec[sbi_] = eq_fields_->diffusion_rate_immobile[sbi_](p);

        temp_exponent_ = (por_mob_ + por_immob_) / (por_mob_ * por_immob_) * eq_data_->time_->dt();

        for (sbi_ = 0; sbi_ < eq_data_->substances_.size(); sbi_++) //over all substances
        {
        	exponent_ = diff_vec[sbi_] * temp_exponent_;
            //previous values
            previous_conc_mob_ = eq_fields_->conc_mobile_fe[sbi_]->vec().get(dof_p0_);
            previous_conc_immob_ = eq_fields_->conc_immobile_fe[sbi_]->vec().get(dof_p0_);

            // ---compute average concentration------------------------------------------
            conc_average_ = ((por_mob_ * previous_conc_mob_) + (por_immob_ * previous_conc_immob_))
                           / (por_mob_ + por_immob_);

            conc_max_ = std::max(previous_conc_mob_-conc_average_, previous_conc_immob_-conc_average_);

            // the following 2 conditions guarantee:
            // 1) stability of forward Euler's method
            // 2) the error of forward Euler's method will not be large
            if(eq_data_->time_->dt() <= por_mob_*por_immob_/(max(diff_vec)*(por_mob_+por_immob_)) &&
                conc_max_ <= (2*eq_data_->scheme_tolerance_/(exponent_*exponent_)*conc_average_))               // forward euler
            {
            	temp_ = diff_vec[sbi_]*(previous_conc_immob_ - previous_conc_mob_) * eq_data_->time_->dt();
                // ---compute concentration in mobile area
                conc_mob_ = temp_ / por_mob_ + previous_conc_mob_;

                // ---compute concentration in immobile area
                conc_immob_ = -temp_ / por_immob_ + previous_conc_immob_;
            }
            else                                                        //analytic solution
            {
                temp_ = exp(-exponent_);
                // ---compute concentration in mobile area
                conc_mob_ = (previous_conc_mob_ - conc_average_) * temp_ + conc_average_;

                // ---compute concentration in immobile area
                conc_immob_ = (previous_conc_immob_ - conc_average_) * temp_ + conc_average_;
            }

            eq_fields_->conc_mobile_fe[sbi_]->vec().set(dof_p0_, conc_mob_);
            eq_fields_->conc_immobile_fe[sbi_]->vec().set(dof_p0_, conc_immob_);
        }
    }


private:
    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int sbi_;                                ///< Index of substance
    IntIdx dof_p0_;                                   ///< Index of local DOF
    double conc_average_;                             ///< weighted (by porosity) average of concentration
    double conc_mob_, conc_immob_;                    ///< new mobile and immobile concentration
    double previous_conc_mob_, previous_conc_immob_;  ///< mobile and immobile concentration in previous time step
    double conc_max_;                                 ///< difference between concentration and average concentration
    double por_mob_, por_immob_;                      ///< mobile and immobile porosity
    double exponent_, temp_exponent_, temp_;          ///< Precomputed values

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;
};

#endif /* ASSEMBLY_REACTION_HH_ */

