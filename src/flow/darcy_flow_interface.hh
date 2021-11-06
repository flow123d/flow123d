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

class DarcyFlowInterface : public EquationBase {
public:
    /// Typedef for usage of Input::Factory in child classes.
    typedef DarcyFlowInterface FactoryBaseType;

    DECLARE_EXCEPTION( ExcBddcmlNotSupported, << "Flow123d was not build with BDDCML support.\n" );
    DECLARE_EXCEPTION( ExcUnknownSolver, << "Unknown solver type. Internal error.\n" );

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
    
    virtual ~DarcyFlowInterface()
    {}
};




#endif /* SRC_FLOW_DARCY_FLOW_INTERFACE_HH_ */
