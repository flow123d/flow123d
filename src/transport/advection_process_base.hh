/*
 * advection_process_base.hh
 *
 *  Created on: Sep 17, 2015
 *      Author: jb
 */

#ifndef SRC_TRANSPORT_ADVECTION_PROCESS_BASE_HH_
#define SRC_TRANSPORT_ADVECTION_PROCESS_BASE_HH_

#include "coupling/equation.hh"
#include "input/input_type_forward.hh"

class Balance;
class Mesh;
class MH_DofHandler;
class SubstanceList;


/**
 * Abstract interface class for secondary equations in HC_ExplicitCoupling.
 */
class AdvectionProcessBase : public EquationBase {

public:
    AdvectionProcessBase(Mesh &mesh, const Input::Record in_rec)
    : EquationBase(mesh, in_rec)
    {};

    /**
     * This method takes sequential PETSc vector of side velocities and update
     * transport matrix. The ordering is same as ordering of sides in the mesh.
     * We just keep the pointer, but do not destroy the object.
     *
     * TODO: We should pass whole velocity field object (description of base functions and dof numbering) and vector.
     */
    virtual void set_velocity_field(const MH_DofHandler &dh) = 0;

    /// Common specification of the input record for secondary equations.
    static Input::Type::Abstract & get_input_type() {
        return Input::Type::Abstract("AdvectionProcess",
                "Abstract advection process. In particular: transport of substances or heat transfer.")
                .close();
    }


};




#endif /* SRC_TRANSPORT_ADVECTION_PROCESS_BASE_HH_ */
