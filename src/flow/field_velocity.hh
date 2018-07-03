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
 * @file    field_velocity.hh
 * @brief   
 */

#ifndef FIELD_VELOCITY_HH_
#define FIELD_VELOCITY_HH_

#include "mh_dofhandler.hh"

#include "fem/mapping_p1.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/fe_rt.hh"

#include "fem/singularity.hh"
#include "fem/xfe_values.hh"

#include "quadrature/quadrature_lib.hh"
#include "quadrature/qxfem.hh"
#include "quadrature/qxfem_factory.hh"

#include "mesh/ref_element.hh"

/// Internal dim dependent part of FieldVelocity class.
template<int dim>
class FieldVelocityInternal
{
public:
    static const unsigned int spacedim = 3;
    typedef typename Space<spacedim>::Point Point;
    typedef typename arma::vec3 Value;
    
    FieldVelocityInternal(MH_DofHandler* mh_dh, bool enr=false, bool reg=true)
    : mh_dh_(mh_dh), fe_p0_(0), enr_(enr), reg_(reg)
    {
        fv_rt_xfem_ = std::make_shared<XFEValues<dim,3>>(map_, fe_rt_, fe_p0_, update_values | update_jacobians |
                                                           update_inverse_jacobians | update_quadrature_points);
    }
    
    Value value_vector(const ElementAccessor<spacedim> &elm, const Point &p){
        //HACK to get accessor; only for SINGLE processor
        auto ele_ac = mh_dh_->accessor(elm.idx());
        update_quad(elm, p);
        Value flux; flux.zeros();
        
        if(enr_ && ele_ac.is_enriched()){
            flux =  value_vector_xfem(ele_ac);
        }
        else if(reg_){
            flux = value_vector_regular(elm);
        }
        
        return flux;
    }
    
    void update_quad(const ElementAccessor<spacedim> &ele, const Point &p){
        arma::vec unit_p = map_.project_real_to_unit(p,map_.element_map(ele));
        quad_.resize(1);
        quad_.set_point(0, RefElement<dim>::bary_to_local(unit_p));
        quad_.set_real_point(0,p);
        quad_.set_weight(0,1.0);
//         DBGCOUT(<< "p: [" << p(0) << " " << p(1) << " " << p(2) << "]\n");
//         if(dim == 2) DBGCOUT(<< "p: [" << unit_p(0) << " " << unit_p(1) << " " << unit_p(2) << "]\n");
    }
    
    Value value_vector_regular(const ElementAccessor<spacedim> &ele){
        Value flux;
        flux.zeros();
        
        fv_rt_ = std::make_shared<FEValues<dim,3>>(map_, quad_, fe_rt_, update_values);
        ElementAccessor<3> ele_nonconst = ele;
        fv_rt_->reinit(ele_nonconst);
        auto velocity = fv_rt_->vector_view(0);
        
        for (unsigned int li = 0; li < ele->n_sides(); li++) {
//             DBGCOUT(<< "ele " << ele->index() << " flux " << mh_dh_->side_flux( *(ele->side( li ) ) ) 
//                 << " [" << fv_rt_->shape_vector(li,0)(0) << " " << fv_rt_->shape_vector(li,0)(1) << " " << fv_rt_->shape_vector(li,0)(2) << "]\n");
            flux += mh_dh_->side_flux( *(ele.side( li ) ) ) * velocity.value(li,0);
        }
        
//         DBGCOUT(<< "ele " << ele->index() << " flux [" << flux(0) << " " << flux(1) << " " << flux(2) << "]\n");
        return flux;
    }
    
    Value value_vector_xfem(LocalElementAccessorBase<3> ele_ac){
        Value flux;
        flux.zeros();
        ElementAccessor<3> ele = ele_ac.element_accessor();
        
        XFEMElementSingularData<dim> * xdata = ele_ac.xfem_data_sing<dim>();

        //TODO select enr type
        if(mh_dh_->single_enr)
            fv_rt_xfem_->reinit(ele, *xdata, quad_);
        auto velocity = fv_rt_xfem_->vector_view(0);
        
        int dofs[fv_rt_xfem_->n_dofs()];
        int ndofs = ele_ac.get_dofs_vel(dofs);
        ASSERT_DBG(fv_rt_xfem_->n_dofs() == ndofs);
        
    //     DBGCOUT("####################### DOFS\n");
    //     for(int i =0; i < ndofs; i++) cout << dofs[i] << " ";
    //     cout << "\n";
        
        int li = 0;
        if(!reg_) li = fv_rt_xfem_->n_regular_dofs();  //skip regular part
        for (; li < ndofs; li++) {
            flux += mh_dh_->mh_solution[dofs[li]]
                        * velocity.value(li,0);
        }
        
        return flux;
    }
    
private:
    MH_DofHandler* mh_dh_;
    
    MappingP1<dim,3> map_;
    QXFEM<dim,3> quad_;
    
    // assembly volume integrals
    FE_RT0<dim> fe_rt_;
    FE_P_disc<dim> fe_p0_;
    std::shared_ptr<FEValues<dim,3>> fv_rt_;
    
    std::shared_ptr<XFEValues<dim,3>> fv_rt_xfem_;
    
    bool enr_, reg_;
};


/**
 * Temporary class creating "FieldFE" for velocity with old MH Dofhandler.
 *
 */
class FieldVelocity : public FieldAlgorithmBase<3, FieldValue<3>::VectorFixed>
{
public:
    typedef typename FieldAlgorithmBase<3, FieldValue<3>::VectorFixed>::Point Point;
    typedef FieldAlgorithmBase<3, FieldValue<3>::VectorFixed> FactoryBaseType;

    FieldVelocity(MH_DofHandler* mh_dh, Field<3, FieldValue<3>::Scalar>* cross_section,
                  bool enriched_part=false, bool regular_part=true)
    : FieldAlgorithmBase<3, FieldValue<3>::VectorFixed>(0),
      fe_val_1d_(mh_dh),
      fe_val_2d_(mh_dh, enriched_part, regular_part),
      fe_val_3d_(mh_dh, enriched_part, regular_part),
      cross_section_(cross_section)
    {}

    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    typename FieldValue<3>::VectorFixed::return_type const &value(const Point &p,
                                                                  const ElementAccessor<3> &elm)
    {
        arma::vec3 val; val.zeros();
        switch(elm.dim()){
            case 1: val = fe_val_1d_.value_vector(elm,p); break;
            case 2: val = fe_val_2d_.value_vector(elm,p); break;
            case 3: val = fe_val_3d_.value_vector(elm,p);
        }
        
        val = val / cross_section_->value(p, elm);
        
//         val.print(cout,"val");
        for (unsigned int i=0; i<3; i++)
            this->value_(i,0) = val(i);
        
        return this->r_value_;
    }

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    void value_list (const std::vector< Point >  &point_list,
                     const ElementAccessor<3> &elm,
                     std::vector<typename FieldValue<3>::VectorFixed::return_type>  &value_list)
    {
        ASSERT_DBG(0).error("Not implemented!");
    }


    virtual ~FieldVelocity()
    {}

private:
    
    FieldVelocityInternal<1> fe_val_1d_;
    FieldVelocityInternal<2> fe_val_2d_;
    FieldVelocityInternal<3> fe_val_3d_;
    
    Field<3, FieldValue<3>::Scalar>* cross_section_;
    
    /// Registrar of class to factory
    static const int registrar;
};

#endif /* FIELD_VELOCITY_HH_ */
