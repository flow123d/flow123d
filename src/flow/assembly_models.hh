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
 * @file    assembly_models.hh
 * @brief   Functors of FieldModels used in Darcy flow module.
 * @author  David Flanderka
 *
 */

#ifndef ASSEMBLY_MODELS_HH
#define ASSEMBLY_MODELS_HH

#include <armadillo>

using Sclr = double;
using Vect = arma::vec3;

// Functor computing velocity (flux / cross_section)
struct fn_mh_velocity {
	inline Vect operator() (Vect flux, Sclr csec) {
        return flux / csec;
    }
};


// Functor computing piezo_head_p0
struct fn_mh_piezohead {
	inline Sclr operator() (Vect gravity, Vect coords, Sclr pressure) {
        return arma::dot((-1*gravity), coords) + pressure;
    }
};

#endif  //ASSEMBLY_MODELS_HH
