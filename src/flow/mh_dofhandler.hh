/*
 * mh_dofhandler.hh
 *
 *  Created on: Jun 18, 2012
 *      Author: jb
 */

#ifndef MH_DOFHANDLER_HH_
#define MH_DOFHANDLER_HH_

#include <vector>
using namespace std;

class Mesh;
class Side;
class SideIter;

/// temporary solution to provide access to results
/// from DarcyFlowMH independent of mesh
class MH_DofHandler {
public:
    void reinit(Mesh *mesh);

    void set_solution( double * solution);

    unsigned int side_dof(const SideIter side) const;

    /// temporary replacement for DofHandler accessor, flux through given side
    double side_flux(const Side &side) const;

    /// temporary replacement for DofHandler accessor, scalar (pressure) on edge of the side
    double side_scalar(const Side &side) const;

protected:
    vector< vector<unsigned int> > elem_side_to_global;
    double * mh_solution;
};

#endif /* MH_DOFHANDLER_HH_ */
