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
 * @file    index_types.hh
 * @brief   
 */

#ifndef INDEX_TYPES_HH
#define INDEX_TYPES_HH

#include <armadillo>

/// Define type that represents indices of large arrays (elements, nodes, dofs etc.)
typedef int LongIdx;
typedef int IntIdx;

typedef arma::Col<LongIdx> GlobalDofVec;
typedef arma::Col<IntIdx> LocDofVec;

#define MPI_LONG_IDX MPI_INT

#endif // INDEX_TYPES_HH
