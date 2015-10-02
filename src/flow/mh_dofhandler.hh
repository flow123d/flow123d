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
 * @file    mh_dofhandler.hh
 * @brief   
 */

#ifndef MH_DOFHANDLER_HH_
#define MH_DOFHANDLER_HH_

#include <vector>
#include "mesh/mesh_types.hh"

using namespace std;

class Mesh;
class Side;
class SideIter;


/// temporary solution to provide access to results
/// from DarcyFlowMH independent of mesh
class MH_DofHandler {
public:
    void reinit(Mesh *mesh);

    void set_solution( double time, double * solution, double precision);

    inline double time_changed() const
        { return time_; }

    unsigned int side_dof(const SideIter side) const;

    /// temporary replacement for DofHandler accessor, flux through given side
    double side_flux(const Side &side) const;

    /// temporary replacement for DofHandler accessor, scalar (pressure) on edge of the side
    double side_scalar(const Side &side) const;

    /// temporary replacement for DofHandler accessor, scalar (pressure) on element
    double element_scalar( ElementFullIter &ele ) const;

    inline double precision() const { return solution_precision; };


protected:
    vector< vector<unsigned int> > elem_side_to_global;
    double * mh_solution;
    double solution_precision;
    double time_;
};

#endif /* MH_DOFHANDLER_HH_ */
