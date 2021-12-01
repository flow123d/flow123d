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

//        shared_ptr<FE_P<dim>> fe_p = std::make_shared< FE_P<dim> >(1);
//        fe_ = std::make_shared<FESystem<dim>>(fe_p, FEVector, 3);
//        fv_.initialize(*this->quad_, *fe_,
//        		update_values | update_gradients | update_quadrature_points);
//        fsv_.initialize(*this->quad_low_, *fe_,
//        		update_values | update_normal_vectors | update_quadrature_points);
//        n_dofs_ = fe_->n_dofs();
//        vec_view_ = &fv_.vector_view(0);
//        //        if (dim>1) ??
//        vec_view_side_ = &fsv_.vector_view(0);
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

#endif /* ASSEMBLY_REACTION_HH_ */

