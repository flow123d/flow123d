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
#include "fem/mapping_p1.hh"             // for MappingP1
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
class DarcyMH;
class Mesh;
class OutputTime;
namespace Input {
	class Record;
	namespace Type {
		class Record;
	}
}

template<unsigned int dim, unsigned int spacedim> class FEValues;
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

	    Field<3, FieldValue<3>::Scalar> field_ele_pressure;
	    Field<3, FieldValue<3>::Scalar> field_node_pressure;
	    Field<3, FieldValue<3>::Scalar> field_ele_piezo_head;
	    Field<3, FieldValue<3>::VectorFixed> field_ele_flux;
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

    DarcyFlowMHOutput(DarcyMH *flow, Input::Record in_rec) ;
    virtual ~DarcyFlowMHOutput();

    static const Input::Type::Instance & get_input_type();
    static const Input::Type::Instance & get_input_type_specific();

    /** \brief Calculate values for output.  **/
    void output();

    //const OutputFields &get_output_fields() { return output_fields; }



protected:
    typedef const vector<unsigned int> & ElementSetRef;

    virtual void prepare_output(Input::Record in_rec);
    virtual void prepare_specific_output(Input::Record in_rec);
    
    void make_side_flux();
    void make_element_scalar(ElementSetRef element_indices);
    
    /** Computes fluxes at the barycenters of elements.
     *  TODO:
     *  We use FEValues to get fluxes at the barycenters of elements,
     *  but we still use MHDofHandler. Once we are able to make output routines
     *  parallel, we can use simply FieldFE for velocity here.
     */
    void make_element_vector(ElementSetRef element_indices);

    //void make_sides_scalar();
    /**
     * \brief Calculate nodes scalar,
     * store it in double* node_scalars instead of node->scalar
     *  */
    void make_node_scalar_param();
    //void make_node_scalar();
    void make_corner_scalar(vector<double> &node_scalar);
    //void make_neighbour_flux();
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


    DarcyMH *darcy_flow;
    Mesh *mesh_;

    /// Specific experimental error computing.
    bool compute_errors_;
    

    /** Pressure head (in [m]) interpolated into nodes. Provides P1 approximation. Indexed by element-node numbering.*/
    VectorMPI corner_pressure;
    /** Pressure head (in [m]) in barycenters of elements (or equivalently mean pressure over every element). Indexed by element indexes in the mesh.*/
    VectorMPI ele_pressure;
    /** Piezo-metric head (in [m]) in barycenter of elements (or equivalently mean pressure over every element). Indexed by element indexes in the mesh.*/
    VectorMPI ele_piezo_head;

    /** Average flux in barycenter of every element. Indexed as elements in the mesh. */
    // TODO: Definitely we need more general (templated) implementation of Output that accept arbitrary containers. So
    // that we can pass there directly vector< arma:: vec3 >
    VectorMPI ele_flux;

    // Temporary objects holding pointers to appropriate FieldFE
    // TODO remove after final fix of equations
    std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> ele_pressure_ptr;   ///< Field of pressure head in barycenters of elements.
    std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> ele_piezo_head_ptr; ///< Field of piezo-metric head in barycenters of elements.
    std::shared_ptr<FieldFE<3, FieldValue<3>::VectorFixed>> ele_flux_ptr;  ///< Field of flux in barycenter of every element.

    // A vector of all element indexes
    std::vector<unsigned int> all_element_idx_;

    // integrals of squared differences on individual elements - error indicators, can be written out into VTK files
    std::vector<double>     l2_diff_pressure, l2_diff_velocity, l2_diff_divergence;

    std::shared_ptr<DOFHandlerMultiDim> dh_;
    FE_P_disc<0> fe0; //TODO temporary solution - add support of FEData<0>
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
        VectorMPI pressure_diff;
        VectorMPI velocity_diff;
        VectorMPI div_diff;

        // Temporary objects holding pointers to appropriate FieldFE
        // TODO remove after final fix of equations
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> pressure_diff_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> vel_diff_ptr;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> div_diff_ptr;

        double * solution;
        const MH_DofHandler * dh;

        std::vector<int> velocity_mask;
        DarcyMH *darcy;
        DarcyMH::EqData *data_;
    } diff_data;
    
    
    /// Struct containing all dim dependent FE classes needed for output
    /// (and for computing solution error).
    template<int dim> struct FEData{
        FEData();
        
        // we create trivial Dofhandler , for P0 elements, to get access to, FEValues on individual elements
        // this we use to integrate our own functions - difference of postprocessed pressure and analytical solution
        FE_P_disc<dim> fe_p0;
        FE_P_disc<dim> fe_p1;

        const unsigned int order; // order of Gauss quadrature
        QGauss quad;

        MappingP1<dim,3> mapp;

        FEValues<dim,3> fe_values;
        
        // FEValues for velocity.
        FE_RT0<dim> fe_rt;
        FEValues<dim, 3> fv_rt;
    };
    
    FEData<1> fe_data_1d;
    FEData<2> fe_data_2d;
    FEData<3> fe_data_3d;
    
    /// Computes L2 error on an element.
    template <int dim>
    void l2_diff_local(ElementAccessor<3> &ele,
                      FEValues<dim,3> &fe_values, FEValues<dim,3> &fv_rt,
                      FieldPython<3, FieldValue<3>::Vector > &anal_sol,  DiffData &result);
};


#endif /* DARCY_FLOW_MH_OUTPUT_HH_ */

