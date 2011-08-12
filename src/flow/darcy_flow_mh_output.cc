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
 * @ingroup flow
 *
 *  @author Jan Brezina
 *
 */

#include "flow/darcy_flow_mh.hh"
#include "flow/darcy_flow_mh_output.hh"
#include "field_p0.hh"

#include "io/output.h"


DarcyFlowMHOutput::DarcyFlowMHOutput(DarcyFlowMH *flow)
: darcy_flow(flow), mesh_(&darcy_flow->mesh())
{
    // setup output
    string output_file = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Output", "Output_file", "\\"));
    DBGMSG("create output\n");
    output_writer = new OutputTime(mesh_, output_file);

    // set output time marks
    TimeMarks &marks=darcy_flow->time().marks();
    output_mark_type = marks.new_strict_mark_type();
    marks.add_time_marks(0.0, OptGetDbl("Global", "Save_step", "1.0"), darcy_flow->time().end_time(), output_mark_type );
    DBGMSG("end create output\n");

    ele_scalars = new double[mesh_->n_elements()];
    node_scalars = new double[mesh_->node_vector.size()];
    element_vectors = new OutVector();
}

DarcyFlowMHOutput::~DarcyFlowMHOutput(){
    delete [] node_scalars;
    delete [] ele_scalars;
    delete element_vectors;
    delete output_writer;
};

//=============================================================================
// CONVERT SOLUTION, CALCULATE BALANCES, ETC...
//=============================================================================

void DarcyFlowMHOutput::postprocess() {

    make_side_flux();

    make_element_scalar();
    make_element_scalar(ele_scalars);

    make_element_vector();
    make_sides_scalar();

    /** new version of make_node_scalar */
    make_node_scalar_param(node_scalars);
    make_node_scalar();


    //make_node_vector( mesh );
    make_neighbour_flux();
    //make_previous_scalar();
    xprintf(Msg, "O.K.\n")/*orig verb 2*/;
    water_balance();
    //  if (problem->transport_on == true)
    //         transport( problem );
    xprintf(Msg, "Postprocessing phase O.K.\n")/*orig verb 2*/;
}

void DarcyFlowMHOutput::output()
{
    std::string nodeName = "pressure_nodes";
    std::string nodeUnit = "L";
    std::string eleScalarName = "pressure_elements";
    std::string eleScalarUnit = "L";
    std::string eleVectorName = "velocity_elements";
    std::string eleVectorUnit = "L/T";

    unsigned int result = 0;

    if (darcy_flow->time().is_current(output_mark_type)) {
        result = output_writer->register_node_data(nodeName, nodeUnit, node_scalars, mesh_->node_vector.size());
        xprintf(Msg, "Register_node_data - result: %i, node size: %i\n", result,  mesh_->node_vector.size());

        result = output_writer->register_elem_data(eleScalarName, eleScalarUnit, ele_scalars, mesh_->n_elements());
        xprintf(Msg, "Register_elem_data scalars - result: %i\n", result);

        element_vectors->vectors = new VectorFloatVector;

        element_vectors->vectors->reserve(mesh_->n_elements());
        FOR_ELEMENTS(mesh_, ele) {
            /* Add vector */
            vector<double> vec;
            vec.reserve(3);
            vec.push_back(ele->vector[0]);
            vec.push_back(ele->vector[1]);
            vec.push_back(ele->vector[2]);
            element_vectors->vectors->push_back(vec);
        }
        result = output_writer->register_elem_data(eleVectorName, eleVectorUnit, *element_vectors->vectors);
        xprintf(Msg, "Register_elem_data vectors - result: %i\n", result);

        output_writer->write_data(darcy_flow->solved_time());

        if(element_vectors->vectors != NULL) {
            delete element_vectors->vectors;
        }
    }
}

//=============================================================================
// FILL TH "FLUX" FIELD FOR ALL SIDES IN THE MESH
//=============================================================================

void DarcyFlowMHOutput::make_side_flux() {
    int li, soi;
    unsigned int sol_size;
    double *sol;
    struct Side *sde;

    soi = 0;
    darcy_flow->get_solution_vector(sol, sol_size);
    FOR_ELEMENTS(mesh_, ele)
        for (li = 0; li < ele->n_sides; li++) {
            sde = ele->side[ li ];
            sde->flux = sol[ soi++ ];
            //if( fabs( sde->flux ) < ZERO )
            //  sde->flux = 0.0;
        }
}
//=============================================================================
// FILL TH "SCALAR" FIELD FOR ALL ELEMENTS IN THE MESH
//=============================================================================

void DarcyFlowMHOutput::make_element_scalar() {
    int soi;
    unsigned int sol_size;
    double *sol;

    soi = mesh_->n_sides;
    darcy_flow->get_solution_vector(sol, sol_size);
    FOR_ELEMENTS(mesh_,ele) ele->scalar = sol[ soi++ ];
}

void DarcyFlowMHOutput::make_element_scalar(double* scalars) {
    int soi;
    unsigned int sol_size;
    double *sol;
    int ele_index = 0; //!< index of each element */

    for (int i = 0; i < mesh_->n_elements(); i++){
        scalars[i] = 0.0;
    };
//    ele_index = mesh->element->ele_index;

    soi = mesh_->n_sides;
    darcy_flow->get_solution_vector(sol, sol_size);
    FOR_ELEMENTS(mesh_, ele){
        ele_index = mesh_->element.index(ele);
//        xprintf(Msg, "ele_index: %i\n", ele_index);
        scalars[ele_index] = sol[ soi++ ];
    }
}

/****
 * compute Darcian velocity in centre of elements
 *
 */
void DarcyFlowMHOutput::make_element_vector() {
    //FILE *out;

    // upravy Ji -- pomocny tisk
    //out = xfopen( "pomout2.txt", "wt" );
    //xfprintf( out, "Pomocny tisk bazovych funkci po vypoctu\n\n");

    FOR_ELEMENTS(mesh_, ele) {
        switch (ele->type) {
        case LINE:
            make_element_vector_line(ele);
            break;
        case TRIANGLE:
            make_element_vector_triangle(ele);
            //xfprintf( out, "%d \n", ele->id);
            //xfprintf( out, "Plocha %12.8f\n", ele->metrics);
            //xfprintf( out, "Teziste\n");
            //xfprintf( out, "%12.8f %12.8f %12.8f \n", ele->centre[0], ele->centre[1], ele->centre[2]);
            //xfprintf( out, "Bazove funkce\n");
            //xfprintf( out, "%12.8f %12.8f %12.8f \n", ele->bas_alfa[0], ele->bas_beta[0], ele->bas_gama[0] );
            //xfprintf( out, "%12.8f %12.8f %12.8f \n", ele->bas_alfa[1], ele->bas_beta[1], ele->bas_gama[1] );
            //xfprintf( out, "%12.8f %12.8f %12.8f \n", ele->bas_alfa[2], ele->bas_beta[2], ele->bas_gama[2] );
            //xfprintf( out, " \n");
            //xfprintf( out, "Pretoky\n");
            //xfprintf( out, "%12.8f %12.8f %12.8f \n", ele->side[0]->flux, ele->side[1]->flux, ele->side[2]->flux);
            //xfprintf( out, "Vektor v tezisti\n");
            //xfprintf( out, "%12.8f %12.8f %12.8f \n",  ele->vector[0], ele->vector[1], ele->vector[2]);
            //xfprintf( out, " \n");
            //xfprintf( out, " \n\n");
            break;
        case TETRAHEDRON:
            make_element_vector_tetrahedron(ele);
            break;
        }
        ele->v_length = vector_length(ele->vector);
    }
    //xfclose( out );
}
//=============================================================================
//
//=============================================================================

void DarcyFlowMHOutput::make_element_vector_line(ElementFullIter ele) {
    double darcy_vel = (ele->side[1]->flux - ele->side[0]->flux) / 2.0 / ele->material->size;

    // normalize element vector [node 0, node 1]
    ele->vector[ 0 ] = (ele->node[1]->getX() - ele->node[0]->getX());
    ele->vector[ 1 ] = (ele->node[1]->getY() - ele->node[0]->getY());
    ele->vector[ 2 ] = (ele->node[1]->getZ() - ele->node[0]->getZ());
    normalize_vector(ele->vector);

    ele->vector[ 0 ] *= darcy_vel;
    ele->vector[ 1 ] *= darcy_vel;
    ele->vector[ 2 ] *= darcy_vel;

}
//=============================================================================
//
//=============================================================================

void DarcyFlowMHOutput::make_element_vector_triangle(ElementFullIter ele) {
    double bas[ 3 ][ 3 ];
    double ex[ 3 ];
    double ey[ 3 ];
    double ez[ 3 ];
    double tmp[ 3 ];
    double X[ 3 ];
    double mid[3];
    double u[3];
    int i, li;

    // begin -- upravy Ji. -- prepocet vektoru do teziste elementu
    // 23.2.2007
    ex[ 0 ] = ele->node[1]->getX() - ele->node[0]->getX();
    ex[ 1 ] = ele->node[1]->getY() - ele->node[0]->getY();
    ex[ 2 ] = ele->node[1]->getZ() - ele->node[0]->getZ();

    tmp[ 0 ] = ele->node[2]->getX() - ele->node[0]->getX();
    tmp[ 1 ] = ele->node[2]->getY() - ele->node[0]->getY();
    tmp[ 2 ] = ele->node[2]->getZ() - ele->node[0]->getZ();
    vector_product(ex, tmp, ez);
    normalize_vector(ex);
    normalize_vector(ez);
    vector_product(ez, ex, ey);
    normalize_vector(ey);

    u[0] = ele->centre[0] - ele->node[0]->getX();
    u[1] = ele->centre[1] - ele->node[0]->getY();
    u[2] = ele->centre[2] - ele->node[0]->getZ();
    mid[0] = scalar_product(ex, u);
    mid[1] = scalar_product(ey, u);
    mid[2] = scalar_product(ez, u);


    for (li = 0; li < 3; li++)
        ele->vector[ li ] = 0;
    for (i = 0; i < 3; i++) {
        bas[ i ][ 0 ] = ele->bas_gama[ i ] * (mid[0] - ele->bas_alfa[ i ]);
        bas[ i ][ 1 ] = ele->bas_gama[ i ] * (mid[1] - ele->bas_beta[ i ]);
        bas[ i ][ 2 ] = 0;
        for (li = 0; li < 2; li++)
            ele->vector[ li ] += ele->side[ i ]->flux * bas[ i ][ li ] / ele->material->size;
    }

    /*for (li = 0; li < 3; li++){
            A[ 0 ][ li ] = ex[ li ];
            A[ 1 ][ li ] = ey[ li ];
            A[ 2 ][ li ] = ez[ li ];
    }
    matrix_x_matrix(ele->vector, 1, 3, A, 3, 3, X);
    for (li = 0; li < 3; li++)
            ele->vector[ li ] = X[ li ];

     */
    X[0] = ele->vector[0] * ex[0] +
            ele->vector[1] * ey[0] +
            ele->vector[2] * ez[0];
    X[1] = ele->vector[0] * ex[1] +
            ele->vector[1] * ey[1] +
            ele->vector[2] * ez[1];
    X[2] = ele->vector[0] * ex[2] +
            ele->vector[1] * ey[2] +
            ele->vector[2] * ez[2];
    for (li = 0; li < 3; li++)
        ele->vector[ li ] = X[ li ];
    // end -- upravy Ji. -- prepocet vektoru do teziste elementu

}
//=============================================================================
//
//=============================================================================

void DarcyFlowMHOutput::make_element_vector_tetrahedron(ElementFullIter ele) {
    double bas[ 4 ][ 3 ];
    int i, li;

    for (li = 0; li < 3; li++)
        ele->vector[ li ] = 0;
    for (i = 0; i < 4; i++) {
        bas[ i ][ 0 ] = ele->bas_delta[ i ] * (ele->centre[ 0 ]
                - ele->bas_alfa[ i ]);
        bas[ i ][ 1 ] = ele->bas_delta[ i ] * (ele->centre[ 1 ]
                - ele->bas_beta[ i ]);
        bas[ i ][ 2 ] = ele->bas_delta[ i ] * (ele->centre[ 2 ]
                - ele->bas_gama[ i ]);
        for (li = 0; li < 3; li++)
            ele->vector[ li ] += ele->side[ i ]->flux * bas[ i ][ li ];
    }
}
//=============================================================================
// FILL TH "SCALAR" FIELD FOR ALL INTERNAL SIDES IN THE MESH
//=============================================================================

void DarcyFlowMHOutput::make_sides_scalar() {
    double *sol;
    int soi, si;
    unsigned int sol_size;
    struct Side *sde;

    soi = mesh_->n_sides + mesh_->n_elements();
    darcy_flow->get_solution_vector(sol, sol_size);

    FOR_EDGES(mesh_, edg) {
        for (si = 0; si < edg->n_sides; si++) {
            sde = edg->side[ si ];
            sde->scalar = sol[ soi ];
        }
        soi++;
    }
}
//=============================================================================
// pressure interpolation
//
// some sort of weighted average over elements and sides/edges of an node
//
// TODO:
// questions:
// - explain details of averaging and motivation for this type
// - edge scalars are stronger since thay are computed twice by each side
// - why division by (nod.aux-1)
// - ?implement without TED, TSD global arrays with single loop over elements, sides and one loop over nodes for normalization
//
//=============================================================================

void DarcyFlowMHOutput::make_node_scalar_param(double* scalars) {
    F_ENTRY_P("nodes"+mesh_->node_vector.size());

    double dist; //!< tmp variable for storing particular distance node --> element, node --> side*/

    /** Iterators */
    NodeIter node;
    ElementIter ele;
    struct Side* side;

    int n_nodes = mesh_->node_vector.size(); //!< number of nodes in the mesh */
    xprintf(Msg,"n_nodes: %i\n", n_nodes);
    int node_index = 0; //!< index of each node */

    int* sum_elements = new int [n_nodes]; //!< sum elements joined to node */
    int* sum_sides = new int [n_nodes]; //!< sum sides joined to node */
    double* sum_ele_dist = new double [n_nodes]; //!< sum distances to all joined elements */
    double* sum_side_dist = new double [n_nodes]; //!<  Sum distances to all joined sides */

    /** tmp variables, will be replaced by ini keys
     * TODO include them into ini file*/
    bool count_elements = true; //!< scalar is counted via elements*/
    bool count_sides = true; //!< scalar is counted via sides */


    /** init arrays */
    for (int i = 0; i < n_nodes; i++){
        sum_elements[i] = 0;
        sum_sides[i] = 0;
        sum_ele_dist[i] = 0.0;
        sum_side_dist[i] = 0.0;
        scalars[i] = 0.0;
    };

    /**first pass - calculate sums (weights)*/
    if (count_elements){
        FOR_ELEMENTS(mesh_, ele)
            for (int li = 0; li < ele->n_nodes; li++) {
                node = ele->node[li]; //!< get Node pointer from element */
                node_index = mesh_->node_vector.index(node); //!< get nod index from mesh */

                dist = sqrt(
                        ((node->getX() - ele->centre[ 0 ])*(node->getX() - ele->centre[ 0 ])) +
                        ((node->getY() - ele->centre[ 1 ])*(node->getY() - ele->centre[ 1 ])) +
                        ((node->getZ() - ele->centre[ 2 ])*(node->getZ() - ele->centre[ 2 ]))
                );
                sum_ele_dist[node_index] += dist;
                sum_elements[node_index]++;
            }
    }
    if (count_sides){
        FOR_SIDES(mesh_, side) {
            for (int li = 0; li < side->n_nodes; li++) {
                node = side->node[li];//!< get Node pointer from element */
                node_index = mesh_->node_vector.index(node); //!< get nod index from mesh */
                dist = sqrt(
                        ((node->getX() - side->centre[ 0 ])*(node->getX() - side->centre[ 0 ])) +
                        ((node->getY() - side->centre[ 1 ])*(node->getY() - side->centre[ 1 ])) +
                        ((node->getZ() - side->centre[ 2 ])*(node->getZ() - side->centre[ 2 ]))
                );

                sum_side_dist[node_index] += dist;
                sum_sides[node_index]++;
            }
        }
    }

    /**second pass - calculate scalar  */
    if (count_elements){
        FOR_ELEMENTS(mesh_, ele)
            for (int li = 0; li < ele->n_nodes; li++) {
                node = ele->node[li];//!< get Node pointer from element */
                node_index = mesh_->node_vector.index(node); //!< get nod index from mesh */

                /**TODO - calculate it again or store it in prior pass*/
                dist = sqrt(
                        ((node->getX() - ele->centre[ 0 ])*(node->getX() - ele->centre[ 0 ])) +
                        ((node->getY() - ele->centre[ 1 ])*(node->getY() - ele->centre[ 1 ])) +
                        ((node->getZ() - ele->centre[ 2 ])*(node->getZ() - ele->centre[ 2 ]))
                );
                scalars[node_index] += ele->scalar *
                        (1 - dist / (sum_ele_dist[node_index] + sum_side_dist[node_index])) /
                        (sum_elements[node_index] + sum_sides[node_index] - 1);
            }
    }
    if (count_sides){
        FOR_SIDES(mesh_, side) {
            for (int li = 0; li < side->n_nodes; li++) {
                node = side->node[li];//!< get Node pointer from element */
                node_index = mesh_->node_vector.index(node); //!< get nod index from mesh */

                /**TODO - calculate it again or store it in prior pass*/
                dist = sqrt(
                        ((node->getX() - side->centre[ 0 ])*(node->getX() - side->centre[ 0 ])) +
                        ((node->getY() - side->centre[ 1 ])*(node->getY() - side->centre[ 1 ])) +
                        ((node->getZ() - side->centre[ 2 ])*(node->getZ() - side->centre[ 2 ]))
                );


                scalars[node_index] += side->scalar *
                        (1 - dist / (sum_ele_dist[node_index] + sum_side_dist[node_index])) /
                        (sum_sides[node_index] + sum_elements[node_index] - 1);
            }
        }
    }

//    xprintf(Msg, "**********************************************************************************************\n");
//    for (int i =0; i<n_nodes; i++){
//           xprintf(Msg, "make_node_scalar_param id: %i, %f\n", i, scalars[i]);
//       };
//    xprintf(Msg, "**********************************************************************************************\n");

    /** free memory */
    delete [] sum_elements;
    delete [] sum_sides;
    delete [] sum_ele_dist;
    delete [] sum_side_dist;

}



void DarcyFlowMHOutput::make_node_scalar() {
    int li;
    NodeIter nod;
    struct Side *sde;
    double dist;
    int max_side_id = 0;

    double **TED;
    double **TSD;


    TED = (double **) xmalloc((mesh_->element.size() + 1) * sizeof (double *));

    FOR_SIDES(mesh_, sde)
    if (max_side_id <= sde->id)
        max_side_id = sde->id;

    TSD = (double **) xmalloc((max_side_id + 1) * sizeof (double *));

    FOR_ELEMENTS(mesh_, ele)
        TED[ele.index()] = (double*) xmalloc(ele->n_nodes * sizeof (double));
    FOR_SIDES(mesh_, sde)
        TSD[sde->id] = (double*) xmalloc(sde->n_nodes * sizeof (double));

    FOR_NODES(mesh_, nod ) {
        nod->scalar = 0.0;
        nod->faux = 0.0;
        nod->aux = 0;
    }
    FOR_ELEMENTS(mesh_, ele)
        for (li = 0; li < ele->n_nodes; li++) {
            nod = ele->node[li];

            dist = sqrt(
                ((nod->getX() - ele->centre[ 0 ])*(nod->getX() - ele->centre[ 0 ])) +
                ((nod->getY() - ele->centre[ 1 ])*(nod->getY() - ele->centre[ 1 ])) +
                ((nod->getZ() - ele->centre[ 2 ])*(nod->getZ() - ele->centre[ 2 ]))
                );

            TED[ele.index()][li] = dist;
            nod->faux += dist; //       nod->faux += 1 / dist;
            nod->aux++;
        }

    FOR_SIDES(mesh_, sde) {
        for (li = 0; li < sde->n_nodes; li++) {
            nod = sde->node[li];

            dist = sqrt(
                    ((nod->getX() - sde->centre[ 0 ])*(nod->getX() - sde->centre[ 0 ])) +
                    ((nod->getY() - sde->centre[ 1 ])*(nod->getY() - sde->centre[ 1 ])) +
                    ((nod->getZ() - sde->centre[ 2 ])*(nod->getZ() - sde->centre[ 2 ]))
            );

            TSD[sde->id][li] = dist;
            nod->faux += dist; //      nod->faux += 1 / dist;
            nod->aux++;
        }
    }
    FOR_ELEMENTS(mesh_, ele)
        for (li = 0; li < ele->n_nodes; li++) {
            nod = ele->node[li];
            nod->scalar += ele->scalar * (1 - TED[ele.index()][li] / nod->faux)
                / (nod->aux - 1); // 1 / (dist * nod->faux);
        }
    FOR_SIDES(mesh_, sde)
        for (li = 0; li < sde->n_nodes; li++) {
            nod = sde->node[li];
            nod->scalar += sde->scalar * (1 - TSD[sde->id][li] / nod->faux)
                / (nod->aux - 1); // 1 / (dist * nod->faux);
        }
    xfree(TED);
    xfree(TSD);
}
/*
//=============================================================================
//
//=============================================================================
void make_node_vector(Mesh* mesh)
{
        int ni;
        int nei;
        TNode* nod;
        ElementIter ele;
        double cnt;

        for( ni = 0; ni < mesh->n_nodes; ni++ ) {
                nod = mesh->node + ni;
                nod->vector[ 0 ] = 0.0;
                nod->vector[ 1 ] = 0.0;
                cnt = 0.0;
                for( nei = 0; nei < nod->n_elements(); nei++ ) {
                        ele = nod->element[ nei ];
                        nod->vector[ 0 ] += ele->vector[ 0 ];
                        nod->vector[ 1 ] += ele->vector[ 1 ];
                        cnt += 1.0;
                }
                nod->vector[ 0 ] /= cnt;
                nod->vector[ 1 ] /= cnt;

        }
}
 */
//=============================================================================
// FILL TH "FLUX" FIELD FOR ALL VV NEIGHBOURS IN THE MESH
//=============================================================================

void DarcyFlowMHOutput::make_neighbour_flux() {
    struct Neighbour *ngh;

    FOR_NEIGHBOURS(mesh_, ngh) {
        if (ngh->type != VV_2E)
            continue;
        ngh->flux = ngh->sigma * ngh->geom_factor * (ngh->element[1]->scalar - ngh->element[0]->scalar);
        continue;
    }
}


//=============================================================================
//
//=============================================================================

void DarcyFlowMHOutput::water_balance() {
    F_ENTRY;

    double bal;
    int c_water;
    struct Boundary *bcd;

    xprintf(Msg, "Calculating balance of water by types...\n");

    std::vector<double> *bcd_balance = new std::vector<double>( mesh_->bcd_group_id.size(), 0.0 );

    FOR_BOUNDARIES(mesh_, bcd) (*bcd_balance)[bcd->group] += bcd->side->flux;
    for(int i=0; i < bcd_balance->size(); ++i)
        xprintf(Msg, "Boundary flux #%d\t%g\n", mesh_->bcd_group_id(i).id(), (*bcd_balance)[i]);

    delete bcd_balance;

    const FieldP0<double> *p_sources=darcy_flow->get_sources();
    if (p_sources != NULL) {
        xprintf(Msg, "Calculating sources of water by material types...\n");

        MaterialDatabase &mat_base = darcy_flow->get_mat_base();
        std::vector<double> *src_balance = new std::vector<double>( mat_base.size(), 0.0 ); // initialize by zero


        FOR_ELEMENTS(mesh_, elm) {
            (*src_balance)[mat_base.index(elm->material)] += elm->volume * p_sources->element_value(elm.index());
        }


        FOR_MATERIALS_IT(mat_base, mat) {
            xprintf(Msg, "Material flux #%d:\t% g\n", mat_base.get_id(mat), (*src_balance)[mat_base.index(mat)]);
        }


    }



}
//=============================================================================
//
//=============================================================================

double calc_water_balance(Mesh* mesh, int c_water) {
    double rc;
    struct Boundary *bcd;

    rc = 0.0;

    return rc;
}
