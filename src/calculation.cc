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
