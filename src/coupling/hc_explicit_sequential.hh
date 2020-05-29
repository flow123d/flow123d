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

#include <memory>
#include <string>
#include <vector>
#include "input/input_type_forward.hh"
#include "input/accessors_forward.hh"
#include "coupling/equation.hh"

class DarcyFlowInterface;
class Mesh;
class AdvectionProcessBase;
class FieldCommon;


/**
 * TODO: should be derived from EquationBase in order to chain couplings
 */
class CouplingBase : public EquationBase {
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
    typedef std::shared_ptr<AdvectionProcessBase> AdvectionPtr;

    struct AdvectionData {
        AdvectionData(AdvectionPtr p)
        : process(p), velocity_changed(false), velocity_time(0.0)
        {}

        AdvectionPtr process;
        bool velocity_changed;
        double velocity_time;
    };

    /**
     * Create an advection process for given input key.
     */
    AdvectionPtr make_advection_process(std::string process_key);

    /**
     * Perform a single time step of given advection process.
     */
    void advection_process_step(AdvectionData &pdata);

    static const int registrar;

    ///
    Input::Record in_record_;

    /// mesh common to darcy flow and transport
    Mesh *mesh;

    /// steady or unsteady water flow simulator based on MH scheme
    std::shared_ptr<DarcyFlowInterface> water;

    /// solute transport with chemistry through operator splitting
    std::vector<AdvectionData> processes_;

    ///
    double min_velocity_time;

    bool is_end_all_;

    FieldCommon *water_content_saturated_;
    FieldCommon *water_content_p0_;
};

#endif /* HC_EXPLICIT_SEQUENTIAL_HH_ */
