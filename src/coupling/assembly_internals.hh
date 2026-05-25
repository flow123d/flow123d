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
 * @file    assembly_internals.hh
 * @brief
 */

#ifndef ASSEMBLY_INTERNALS_HH_
#define ASSEMBLY_INTERNALS_HH_

#include "fem/eval_points.hh"
#include "fem/element_cache_map.hh"
#include "fem/patch_fe_values.hh"



/// Holds common data shared between GenericAssemblz and Assembly<dim> classes.
struct AssemblyInternals {
public:
    AssemblyInternals()
    : eval_points_(std::make_shared<EvalPoints>()) {}

    AssemblyInternals(MixedPtr<FiniteElement> fe)
    : eval_points_(std::make_shared<EvalPoints>()), fe_values_(fe) {}

    std::shared_ptr<EvalPoints> eval_points_;                     ///< EvalPoints object shared by all integrals
    ElementCacheMap element_cache_map_;                           ///< ElementCacheMap according to EvalPoints
    PatchFEValues<3> fe_values_;                                  ///< Common FEValues object over all dimensions
};


#endif /* ASSEMBLY_INTERNALS_HH_ */
