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

#include "fields/fe_value_handler.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "quadrature/quadrature_lib.hh"


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class AssemblyDG
{
public:

    /// Constructor.
    AssemblyDG(unsigned int fe_order) {
        fe_ = new FE_P_disc<dim>(fe_order);
        fe_rt_ = new FE_RT0<dim>;
        quad_ = new QGauss<dim>(2*fe_order);
       	mapping_ = new MappingP1<dim,3>;
    }

    /// Set FieldEvaluate object
    void set_field_evaluation( std::shared_ptr<FieldEvaluate<dim, 3, FieldValue<3>::VectorFixed>> field_eval) {
        this->field_eval_ = field_eval;
    }

    /// Destructor.
    ~AssemblyDG() {
        delete fe_;
        delete fe_rt_;
        delete quad_;
        delete mapping_;
    }

private:

    FiniteElement<dim> *fe_;     ///< Finite element for the solution of the advection-diffusion equation.
    FiniteElement<dim> *fe_rt_;  ///< Finite element for the water velocity field.
    Quadrature<dim> *quad_;      ///< Quadrature used in assembling methods.
    MappingP1<dim,3> *mapping_;  ///< Auxiliary mapping of reference elements.

    /// Evaluate values of velocity field.
    std::shared_ptr<FieldEvaluate<dim, 3, FieldValue<3>::VectorFixed>> field_eval_;

    friend class FEObjects;
};



#endif /* FE_VALUE_HANDLER_HH_ */
