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
#include <boost/exception/info.hpp>      // for operator<<, error_info::erro...
#include <cmath>                         // for pow
#include <iosfwd>                        // for ofstream
#include <memory>                        // for shared_ptr
#include <vector>                        // for vector
#include <armadillo>
#include "fem/fe_p.hh"                   // for FE_P_disc
#include "fem/fe_rt.hh"                  // for FE_RT0
#include "fem/fe_values.hh"              // for FEValues
#include "quadrature/quadrature_lib.hh"  // for QGauss
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
class DarcyFlowInterface;
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

    DarcyFlowMHOutput(DarcyFlowInterface *flow, Input::Record in_rec) ;
    virtual ~DarcyFlowMHOutput();

    static const Input::Type::Instance & get_input_type(FieldSet& eq_data, const std::string &equation_name);
    static const Input::Type::Instance & get_input_type_specific();

    /** \brief Calculate values for output.  **/
    void output();

    //const OutputFields &get_output_fields() { return output_fields; }



protected:
    typedef const vector<unsigned int> & ElementSetRef;

    virtual void prepare_output(Input::Record in_rec);
    virtual void prepare_specific_output(Input::Record in_rec);
    
    void output_internal_flow_data();

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
    void compute_l2_difference();


    DarcyFlowInterface *darcy_flow;
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

    /// Raw data output file.
    ofstream raw_output_file;

    /// Output specific field stuff
    bool is_output_specific_fields;
    struct DiffData {
        double pressure_error[3], velocity_error[3], div_error[3];
        double mask_vel_error;

        // Temporary objects holding pointers to appropriate FieldFE
        // TODO remove after final fix of equations
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> pressure_diff_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> vel_diff_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> div_diff_ptr;

        std::shared_ptr<SubDOFHandlerMultiDim> dh_;

        std::vector<int> velocity_mask;
        DarcyMH::EqData* data_;
    } diff_data;
    
    //MixedPtr<FE_P_disc> fe_p0;
    
    /// Struct containing all dim dependent FE classes needed for output
    /// (and for computing solution error).
    struct FEData{
        FEData();
        
        const unsigned int order; // order of Gauss quadrature
        QGauss::array quad;
        MixedPtr<FE_P_disc> fe_p1;

        // following is used for calculation of postprocessed pressure difference
        // and comparison to analytical solution
        MixedPtr<FE_P_disc> fe_p0;
        std::vector<FEValues<3>> fe_values;
        
        // FEValues for velocity.
        MixedPtr<FE_RT0> fe_rt;
        std::vector<FEValues<3>> fv_rt;
    };
    
    FEData fe_data;
    
    /// Computes L2 error on an element.
    void l2_diff_local(DHCellAccessor dh_cell,
                      FEValues<3> &fe_values, FEValues<3> &fv_rt,
                      FieldPython<3, FieldValue<3>::Vector > &anal_sol,  DiffData &result);
};


#endif /* DARCY_FLOW_MH_OUTPUT_HH_ */

