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
 * @file    eval_subset.cc
 * @brief
 * @author  David Flanderka
 */

#include "fields/eval_subset.hh"
#include "fields/eval_points.hh"


/******************************************************************************
 * Implementation of EvalSubset methods.
 */
void EvalSubset::print_bulk_points() {
    std::cout << "Print bulk points:" << std::endl;
    for (unsigned int i=0; i<this->point_indices_[0].size(); ++i)
        std::cout << "--- bulk point:" << std::endl << this->eval_points().local_point( this->point_indices_[0][i] );
}

void EvalSubset::print_side_points(unsigned int permutation) {
    std::cout << "Print side points with permutation: " << permutation << std::endl;
    unsigned int point_size = this->point_indices_[permutation].size();
    unsigned int points_per_side = point_size / (dim+1);
    for (unsigned int i=0; i<point_size; ++i)
        std::cout << "--- side point (side " << (i / points_per_side) << ")" << std::endl
	        << this->eval_points().local_point( this->point_indices_[permutation][i] );
}



/******************************************************************************
 * Explicit instantiation of templates
 */
