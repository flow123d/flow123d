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
 * @file    darcy_flow_mh_output.hh
 * @brief   Output class for darcy_flow_mh model.
 * @author  Jan Brezina
 */

#ifndef DARCY_FLOW_MH_OUTPUT_XFEM_HH_
#define DARCY_FLOW_MH_OUTPUT_XFEM_HH_

#include <stdio.h>                       // for sprintf
#include <string.h>                      // for memcpy
#include <boost/exception/info.hpp>      // for operator<<, error_info::erro...
#include <cmath>                         // for pow
#include <iosfwd>                        // for ofstream
#include <memory>                        // for shared_ptr
#include <vector>                        // for vector
#include <armadillo>
#include "fem/fe_p.hh"                   // for FE_P_disc
#include "fem/mapping_p1.hh"             // for MappingP1
#include "fields/equation_output.hh"     // for EquationOutput
#include "fields/field.hh"               // for Field
#include "fields/field_set.hh"           // for FieldSet
#include "fields/field_values.hh"        // for FieldValue<>::Scalar, FieldV...
#include "fields/vec_seq_double.hh"      // for VectorSeqDouble
#include "input/type_base.hh"            // for Array
#include "input/type_generic.hh"         // for Instance
#include "petscvec.h"                    // for Vec, _p_Vec
#include "system/exceptions.hh"          // for ExcAssertMsg::~ExcAssertMsg

#include "flow/darcy_flow_mh_output.hh"

class DOFHandlerMultiDim;
class DarcyMH;
class Mesh;
class OutputTime;
class FieldVelocity;
namespace Input {
	class Record;
	namespace Type {
		class Record;
	}
}


/**
 * Actually this class only collect former code from postprocess.*
 * This code depends on values stored in mesh and has to be changed to use fields or other data provided by
 * interface to darcy_flow. We have to relay on the interface in order to allow different implementation of darcy_flow.
 *
 * Principal functionalities of current postprocess are:
 * - move computed values into mesh - this can be removed after we get new output classes,
 *   other dependent code can use directly field interface
 * - compute interpolation to nodes - this can be external and operate directly on fields
 * - compute water balances over materials and boundary parts - is feasible through field interface
 *
 * Further functionality of this class should be:
 * - prepare data for general output classes
 *
 */
class DarcyFlowMHOutputXFEM : public DarcyFlowMHOutput {
public:
    DarcyFlowMHOutputXFEM(DarcyMH *flow, Input::Record in_rec) ;
    ~DarcyFlowMHOutputXFEM();


private:
    void make_element_vector(ElementSetRef element_indices) override;
  
    /**
     * Temporary hack.
     * Calculate approximation of L2 norm for:
     * 1) difference between regularized pressure and analytical solution (using FunctionPython)
     * 2) difference between RT velocities and analytical solution
     * 3) difference of divergence
     *
     * TODO:
     * 1) implement field objects
     * 2) implement DG_P2 finite elements
     * 3) implement pressure postprocessing (result is DG_P2 field)
     * 4) implement calculation of L2 norm for two field (compute the norm and values on individual elements as P0 field)
     */
    void compute_l2_difference() override;

    std::shared_ptr<FieldVelocity> field_velocity;
    
    std::shared_ptr<FieldVelocity> field_velocity_enr_part;
    std::shared_ptr<FieldVelocity> field_velocity_reg_part;

    void prepare_specific_output() override;
};


#endif /* DARCY_FLOW_MH_OUTPUT_XFEM_HH_ */

