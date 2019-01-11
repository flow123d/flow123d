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
 * @file    hm_iterative.hh
 * @brief   
 * @author  Jan Stebel
 */

#ifndef HM_ITERATIVE_HH_
#define HM_ITERATIVE_HH_

#include <memory>
#include <string>
#include <vector>
#include "input/input_type_forward.hh"
#include "input/accessors_forward.hh"
#include "coupling/equation.hh"
#include "flow/darcy_flow_interface.hh"

class Mesh;
class FieldCommon;
class RichardsLMH;


/**
 * @brief Class for solution of fully coupled flow and mechanics, solution by fixed-stress iterative splitting.
 *
 */
class HM_Iterative : public DarcyFlowInterface {
public:
    static const Input::Type::Record & get_input_type();

    HM_Iterative(Mesh &mesh, Input::Record in_record);
    void zero_time_step() override;
    void update_solution() override;
    const MH_DofHandler & get_mh_dofhandler() override {};
    ~HM_Iterative();

private:


    static const int registrar;

    /// steady or unsteady water flow simulator based on MH scheme
    std::shared_ptr<RichardsLMH> flow_;

    /// solute transport with chemistry through operator splitting
    std::shared_ptr<EquationBase> mechanics_;
    
    double beta_;
    
    unsigned int min_it_;
    
    unsigned int max_it_;
    
    double a_tol_;
    
    double r_tol_;

};

#endif /* HC_EXPLICIT_SEQUENTIAL_HH_ */
