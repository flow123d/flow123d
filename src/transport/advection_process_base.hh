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

class Mesh;

/**
 * Abstract interface class for secondary equations in HC_ExplicitCoupling.
 */
class AdvectionProcessBase : public EquationBase {

public:
    AdvectionProcessBase(Mesh &mesh, const Input::Record in_rec)
    : EquationBase(mesh, in_rec)
    {};

    /// Common specification of the input record for secondary equations.
    static Input::Type::Abstract & get_input_type() {
        return Input::Type::Abstract("AdvectionProcess",
                "Abstract advection process. In particular: transport of substances or heat transfer.")
                .close();
    }


};




#endif /* SRC_TRANSPORT_ADVECTION_PROCESS_BASE_HH_ */
