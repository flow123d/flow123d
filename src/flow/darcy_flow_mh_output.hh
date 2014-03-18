/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id: darcy_flow_mh.hh 877 2011-02-04 13:13:25Z jakub.sistek $
 * $Revision: 877 $
 * $LastChangedBy: jakub.sistek $
 * $LastChangedDate: 2011-02-04 14:13:25 +0100 (Fri, 04 Feb 2011) $
 *
 * @file
 * @brief Output class for darcy_flow_mh model.
 *
 *  @author Jan Brezina
 *
 */

#ifndef DARCY_FLOW_MH_OUTPUT_HH_
#define DARCY_FLOW_MH_OUTPUT_HH_

#include "mesh/mesh.h"
#include <string>
#include <vector>
#include "io/output.h"

#include "input/input_type.hh"
#include "input/accessors.hh"

#include "fem/mapping_p1.hh"
#include "fem/fe_p.hh"


class DarcyFlowMH;
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

	class OutputFields : public FieldSet {
	public:

		OutputFields();

	    Field<3, FieldValue<3>::Scalar> field_ele_pressure;
	    Field<3, FieldValue<3>::Scalar> field_node_pressure;
	    Field<3, FieldValue<3>::Scalar> field_ele_piezo_head;
	    Field<3, FieldValue<3>::VectorFixed> field_ele_flux;

	};

    DarcyFlowMHOutput(DarcyFlowMH *flow, Input::Record in_rec) ;
    ~DarcyFlowMHOutput();

    static Input::Type::Record input_type;

    void postprocess();

    /** \brief Calculate values for output.  **/
    void output();


private:
    void make_side_flux();
    void make_element_scalar();
    void make_element_vector();

    //void make_element_vector_line(ElementFullIter, arma::vec3 &vec);
    //void make_element_vector_triangle(ElementFullIter, arma::vec3 &vec);
    //void make_element_vector_tetrahedron(ElementFullIter, arma::vec3 &vec);

    void make_sides_scalar();
    /**
     * \brief Calculate nodes scalar,
     * store it in double* node_scalars instead of node->scalar
     *  */
    void make_node_scalar_param(double scalars[]);
    void make_node_scalar();
    void make_corner_scalar(double *node_scalar, double *corner_scalar);
    void make_neighbour_flux();
    //void make_previous_scalar();
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

    /**
     * Calculate and output water balance over material subdomains and boudary fluxes.
     * Works only for steady flow.
     *
     * TODO:
     * - fix it also for unsteady flow
     * - create separate class for this caculations and output
     * - create class for output of tables with support to output into various file formats
     *   like GNUplot of excel/open calc
     **/
    void water_balance();
    double calc_water_balance();

    DarcyFlowMH *darcy_flow;
    Mesh *mesh_;

    TimeMark::Type output_mark_type;

    /// Accessor to the input record for the DarcyFlow output.
    Input::Record   in_rec_;


    /** This we need to allow piezo output and nead not to modify all test outputs. It should be replaced by
     *  more general scheme, where you can switch every output field on or off.
     */
    bool output_piezo_head;

    /** Pressure head (in [m]) interpolated into nodes. Provides P1 approximation. Indexed by node indexes in mesh.*/
    double *node_pressure;
    /** Pressure head (in [m]) interpolated into nodes. Provides P1 approximation. Indexed by element-node numbering.*/
    double *corner_pressure;
    /** Pressure head (in [m]) in barycenters of elements (or equivalently mean pressure over every element). Indexed by element indexes in the mesh.*/
    double *ele_pressure;
    /** Piezo-metric head (in [m]) in barycenter of elements (or equivalently mean pressure over every element). Indexed by element indexes in the mesh.*/
    double *ele_piezo_head;

    /** Average flux in barycenter of every element. Indexed as elements in the mesh. */
    // TODO: Definitely we need more general (templated) implementation of Output that accept arbitrary containers. So
    // that we can pass there directly vector< arma:: vec3 >
    double *ele_flux;

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


    Input::Record output_rec;

    /// Temporary solution for writing balance into separate file.
    FILE *balance_output_file;
    /// Raw data output file.
    FILE *raw_output_file;
};


#endif /* DARCY_FLOW_MH_OUTPUT_HH_ */

