/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief ???
 *
 */

#include "system.hh"
#include "mesh.h"
#include "problem.h"
#include "boundaries.h"
#include "calculation.h"
#include "local_matrix.h"
#include "constantdb.h"
#include "mesh/ini_constants_mesh.hh"

/**
 * calculation_mh(struct Problem *problem)
 */
void calculation_mh(struct Problem *problem) {
    struct Side *sde;

    ASSERT(!(problem == NULL), "NULL as argument of function calculation_mh()\n");

    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    FOR_SIDES(sde) {
        calc_side_metrics(sde);
    }

    edge_calculation_mh(mesh);
    element_calculation_mh(mesh);
    side_calculation_mh(mesh, problem);
    boundary_calculation_mh(mesh);
    local_matrices_mh(mesh);
}

/**
 * calculation_unsteady(struct Problem *problem)
 */
void calculation_unsteady(struct Problem *problem) {
    element_calculation_unsteady(problem);
}
//-----------------------------------------------------------------------------
// vim: set cindent:
