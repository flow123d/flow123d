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

#ifndef DARCY_FLOW_MH_OUTPUT_HH_
#define DARCY_FLOW_MH_OUTPUT_HH_

#include <stdio.h>                       // for sprintf
#include <string.h>                      // for memcpy

#include <cmath>                         // for pow
#include <iosfwd>                        // for ofstream
#include <memory>                        // for shared_ptr
#include <vector>                        // for vector
#include <armadillo>
#include "fields/equation_output.hh"     // for EquationOutput
#include "fields/field.hh"               // for Field
#include "fields/field_set.hh"           // for FieldSet
#include "fields/field_values.hh"        // for FieldValue<>::Scalar, FieldV...
#include "fields/field_python.hh"
#include "la/vector_mpi.hh"              // for VectorMPI
#include "input/type_base.hh"            // for Array
#include "input/type_generic.hh"         // for Instance
#include "petscvec.h"                    // for Vec, _p_Vec
#include "system/exceptions.hh"          // for ExcAssertMsg::~ExcAssertMsg

class DOFHandlerMultiDim;
class DarcyLMH;
class Mesh;
class OutputTime;
namespace Input {
	class Record;
	namespace Type {
		class Record;
	}
}

template<unsigned int spacedim> class FEValues;
template <int spacedim, class Value> class FieldFE;
template<unsigned int dim> class L2DifferenceAssembly;
template<unsigned int dim> class OutputInternalFlowAssembly;
template< template<IntDim...> class DimAssembly> class GenericAssembly;

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
class DarcyFlowMHOutput {
public:

    /// Standard quantities for output in DarcyFlowMH.
	class OutputFields : public EquationOutput {
	public:

		OutputFields();

	    Field<3, FieldValue<3>::Scalar> subdomain;
	    Field<3, FieldValue<3>::Scalar> region_id;
	};

    /// Specific quantities for output in DarcyFlowMH - error estimates etc.
    class OutputSpecificFields : public EquationOutput {
        public:
            OutputSpecificFields();
            
            Field<3, FieldValue<3>::Scalar> velocity_diff;
            Field<3, FieldValue<3>::Scalar> pressure_diff;
            Field<3, FieldValue<3>::Scalar> div_diff;
    };

    /// Output specific field stuff
    class DiffEqData {
    public:
        DiffEqData() {}

        /// Returns pair { quad_order_asm, quad_order_fields}
        inline std::vector<unsigned int> quad_order() const {
            return {2, 2};
        }

        double pressure_error[3], velocity_error[3], div_error[3];
        double mask_vel_error;

        // Temporary objects holding pointers to appropriate FieldFE
        // TODO remove after final fix of equations
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> pressure_diff_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> vel_diff_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> div_diff_ptr;

        std::shared_ptr<SubDOFHandlerMultiDim> dh_;

        std::vector<int> velocity_mask;
        std::shared_ptr<DarcyLMH::EqData> flow_data_;
    };

    class RawOutputEqData {
    public:
    	RawOutputEqData() {}

        /// Returns pair { quad_order_asm, quad_order_fields}
        inline std::vector<unsigned int> quad_order() const {
            return {0, 0};
        }

        ofstream raw_output_file;                        ///< Raw data output file.
        std::vector< std::string > raw_output_strings_;  ///< Output lines of cells.
        TimeGovernor *time_;                             ///< Time is shared with flow equation.

        std::shared_ptr<DarcyLMH::EqData> flow_data_;
    };

    DarcyFlowMHOutput(DarcyLMH *flow, Input::Record in_rec) ;
    virtual ~DarcyFlowMHOutput();

    static const Input::Type::Instance & get_input_type(FieldSet& eq_data, const std::string &equation_name);
    static const Input::Type::Instance & get_input_type_specific();

    /** \brief Calculate values for output.  **/
    void output();

    //const OutputFields &get_output_fields() { return output_fields; }



protected:
    typedef const vector<unsigned int> & ElementSetRef;

    virtual void prepare_output(Input::Record in_rec);
    void prepare_specific_output(Input::Record in_rec);
    void set_specific_output_python_fields();
    
    template <class FieldType>
    void set_ref_solution(const Input::Record &rec, Field<3, FieldType> &output_field, std::vector<std::string> reg);


    DarcyLMH *darcy_flow;
    Mesh *mesh_;

    /// Specific experimental error computing.
    bool compute_errors_;
    

    /** Pressure head (in [m]) interpolated into nodes. Provides P1 approximation. Indexed by element-node numbering.*/
    VectorMPI corner_pressure;
    /** Pressure head (in [m]) in barycenters of elements (or equivalently mean pressure over every element). Indexed by element indexes in the mesh.*/
    VectorMPI ele_pressure;
    /** Piezo-metric head (in [m]) in barycenter of elements (or equivalently mean pressure over every element). Indexed by element indexes in the mesh.*/
    VectorMPI ele_piezo_head;

    // integrals of squared differences on individual elements - error indicators, can be written out into VTK files
    std::vector<double>     l2_diff_pressure, l2_diff_velocity, l2_diff_divergence;

    std::shared_ptr<DOFHandlerMultiDim> dh_;
    std::shared_ptr<DiscreteSpace> ds;

    OutputFields output_fields;
    OutputSpecificFields output_specific_fields;
    
    std::shared_ptr<OutputTime> output_stream;

    /// Output specific field stuff
    bool is_output_specific_fields;
    std::shared_ptr<DarcyLMH::EqFields> flow_eq_fields_;
    std::shared_ptr<DiffEqData> diff_eq_data_;
    std::shared_ptr<RawOutputEqData> raw_eq_data_;
    
    //MixedPtr<FE_P_disc> fe_p0;

    /// general assembly objects, hold assembly objects of appropriate dimension
    GenericAssembly< L2DifferenceAssembly > * l2_difference_assembly_;
    GenericAssembly< OutputInternalFlowAssembly > * output_internal_assembly_;

};


#endif /* DARCY_FLOW_MH_OUTPUT_HH_ */

