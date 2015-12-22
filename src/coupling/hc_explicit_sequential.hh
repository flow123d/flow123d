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
 * @file    hc_explicit_sequential.hh
 * @brief   
 * @author  Jan Brezina
 */

#ifndef HC_EXPLICIT_SEQUENTIAL_HH_
#define HC_EXPLICIT_SEQUENTIAL_HH_

#include "input/input_type_forward.hh"
#include "input/accessors_forward.hh"

#include "coupling/equation.hh"

class DarcyFlowInterface;
class Mesh;
class EquationBase;
class AdvectionProcessBase;


/**
 * TODO: should be derived from EquationBase in order to chain couplings
 */
class CouplingBase {
public:
    static Input::Type::Abstract & get_input_type();

};


/**
 * @brief Class for solution of steady or unsteady flow with sequentially coupled explicit transport.
 *
 */
class HC_ExplicitSequential : public CouplingBase {
public:
    static const Input::Type::Record & get_input_type();

    HC_ExplicitSequential(Input::Record in_record);
    void run_simulation();
    ~HC_ExplicitSequential();

private:

    static const int registrar;

    /// mesh common to darcy flow and transport
    Mesh *mesh;

    /// steady or unsteady water flow simulator based on MH scheme
    std::shared_ptr<DarcyFlowInterface> water;

    /// explicit transport with chemistry through operator splitting
    std::shared_ptr<AdvectionProcessBase> secondary_eq;

};

#endif /* HC_EXPLICIT_SEQUENTIAL_HH_ */
