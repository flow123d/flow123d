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

#include "mesh/mesh.h"
#include <string>
#include <vector>

#include "io/output_time.hh"
#include "input/input_type_forward.hh"
#include "input/accessors.hh"

#include "fem/mapping_p1.hh"
#include "fem/fe_p.hh"

#include "fields/vec_seq_double.hh"
#include "io/equation_output.hh"


class DarcyMH;
class OutputTime;
class DOFHandlerMultiDim;


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

	class OutputFields : public EquationOutput {
	public:

		OutputFields();

	    Field<3, FieldValue<3>::Scalar> field_ele_pressure;
	    Field<3, FieldValue<3>::Scalar> field_node_pressure;
	    Field<3, FieldValue<3>::Scalar> field_ele_piezo_head;
	    Field<3, FieldValue<3>::VectorFixed> field_ele_flux;
	    Field<3, FieldValue<3>::Integer> subdomain;
	    Field<3, FieldValue<3>::Integer> region_id;

	    Field<3, FieldValue<3>::Scalar> velocity_diff;
	    Field<3, FieldValue<3>::Scalar> pressure_diff;
	    Field<3, FieldValue<3>::Scalar> div_diff;

	    FieldSet error_fields_for_output;
	};

    DarcyFlowMHOutput(DarcyMH *flow, Input::Record in_rec) ;
    ~DarcyFlowMHOutput();

    static const Input::Type::Instance & get_input_type();
    static const Input::Type::Record & get_input_type_specific();

    /** \brief Calculate values for output.  **/
    void output();

    //const OutputFields &get_output_fields() { return output_fields; }



private:
    typedef const vector<unsigned int> & ElementSetRef;

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
    void make_node_scalar_param(ElementSetRef element_indices);
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
    vector<double> corner_pressure;
    /** Pressure head (in [m]) in barycenters of elements (or equivalently mean pressure over every element). Indexed by element indexes in the mesh.*/
    VectorSeqDouble ele_pressure;
    /** Piezo-metric head (in [m]) in barycenter of elements (or equivalently mean pressure over every element). Indexed by element indexes in the mesh.*/
    VectorSeqDouble ele_piezo_head;

    /** Average flux in barycenter of every element. Indexed as elements in the mesh. */
    // TODO: Definitely we need more general (templated) implementation of Output that accept arbitrary containers. So
    // that we can pass there directly vector< arma:: vec3 >
    VectorSeqDouble ele_flux;

    // A vector of all element indexes
    std::vector<unsigned int> all_element_idx_;

    // integrals of squared differences on individual elements - error indicators, can be written out into VTK files
    std::vector<double>     l2_diff_pressure, l2_diff_velocity, l2_diff_divergence;

    Vec vec_corner_pressure;
    DOFHandlerMultiDim *dh;
    MappingP1<1,3> map1;
    MappingP1<2,3> map2;
    MappingP1<3,3> map3;
    FE_P_disc<1,1,3> fe1;
    FE_P_disc<1,2,3> fe2;
    FE_P_disc<1,3,3> fe3;

    OutputFields output_fields;

    std::shared_ptr<OutputTime> output_stream;

    /// Raw data output file.
    ofstream raw_output_file;
};


#endif /* DARCY_FLOW_MH_OUTPUT_HH_ */

