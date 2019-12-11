/*
 * darcy_flow_interface.hh
 *
 *  Created on: Sep 17, 2015
 *      Author: jb
 */

#ifndef SRC_FLOW_DARCY_FLOW_INTERFACE_HH_
#define SRC_FLOW_DARCY_FLOW_INTERFACE_HH_

#include "input/input_type_forward.hh"
#include "coupling/equation.hh"
#include "fields/field_values.hh"

template <int spacedim, class Value> class FieldFE;

class DarcyFlowInterface : public EquationBase {
public:
    /// Typedef for usage of Input::Factory in child classes.
    typedef DarcyFlowInterface FactoryBaseType;

    static Input::Type::Abstract & get_input_type() {
        return Input::Type::Abstract("DarcyFlow",
                "Darcy flow model. Abstraction of various porous media flow models.")
                .close();
    }
    
    /// Type of experimental Mortar-like method for non-compatible 1d-2d interaction.
    enum MortarMethod {
        NoMortar = 0,
        MortarP0 = 1,
        MortarP1 = 2
    };

    DarcyFlowInterface(Mesh &mesh, const Input::Record in_rec)
    : EquationBase(mesh, in_rec)
    {}

    /// Return last time of TimeGovernor.
    virtual double last_t() =0;

    // TODO: remove! Due to MH and LMH Darcy flow versions.
    virtual std::shared_ptr< FieldFE<3, FieldValue<3>::VectorFixed> > get_velocity_field()
    { return nullptr; }
    
    virtual ~DarcyFlowInterface()
    {}
};




#endif /* SRC_FLOW_DARCY_FLOW_INTERFACE_HH_ */
