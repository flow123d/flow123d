/*
 * darcy_flow_interface.hh
 *
 *  Created on: Sep 17, 2015
 *      Author: jb
 */

#ifndef SRC_FLOW_DARCY_FLOW_INTERFACE_HH_
#define SRC_FLOW_DARCY_FLOW_INTERFACE_HH_

#include "input/input_type_forward.hh"

class MH_DofHandler;

class DarcyFlowInterface : public EquationBase {
public:
    /// Typedef for usage of Input::Factory in child classes.
    typedef DarcyFlowInterface FactoryBaseType;

    static Input::Type::Abstract & get_input_type() {
        return Input::Type::Abstract("DarcyFlow",
                "Darcy flow model. Abstraction of various porous media flow models.")
                .close();
    }

    DarcyFlowInterface(Mesh &mesh, const Input::Record in_rec)
    : EquationBase(mesh, in_rec)
    {}

    virtual const MH_DofHandler &get_mh_dofhandler() =0;

    virtual ~DarcyFlowInterface()
    {}
};




#endif /* SRC_FLOW_DARCY_FLOW_INTERFACE_HH_ */
