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

#include "flow/mh_fe_values.hh"
#include "flow/darcy_flow_mh.hh"
#include "flow/darcy_flow_mh_output.hh"
#include "field_p0.hh"
#include "system/system.hh"
#include "io_namehandler.hh"
#include <vector>

#include "io/output.h"


DarcyFlowMHOutput::DarcyFlowMHOutput(DarcyFlowMH *flow)
: darcy_flow(flow), mesh_(&darcy_flow->mesh()),
  ele_flux(mesh_->n_elements(),std::vector<double>(3,0.0)),
  balance_output_file(NULL),raw_output_file(NULL)
{
    // setup output
    string output_file = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Output", "Output_file", "\\"));
    DBGMSG("create output\n");
    output_writer = new OutputTime(mesh_, output_file);

    // allocate output containers
    ele_pressure.resize(mesh_->n_elements());
    node_pressure.resize(mesh_->node_vector.size());
    output_piezo_head=OptGetBool("Output","output_piezo_head","No");

    if (output_piezo_head) ele_piezo_head.resize(mesh_->n_elements());




    // set output time marks
    TimeMarks &marks = darcy_flow->time().marks();
    output_mark_type = darcy_flow->mark_type() | marks.type_fixed_time();
    marks.add_time_marks(0.0, OptGetDbl("Global", "Save_step", "1.0"), darcy_flow->time().end_time(), output_mark_type );
    DBGMSG("end create output\n");


    // temporary solution for balance output
    std::string balance_output_fname =
            IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Output", "balance_output", "water_balance"));
    balance_output_file = xfopen(balance_output_fname.c_str(), "wt");

    // optionally open raw output file
    std::string raw_output_fname=OptGetFileName("Output","raw_flow_output","//");
    if (raw_output_fname!= "//") {
        raw_output_fname=IONameHandler::get_instance()->get_output_file_name(raw_output_fname);
        raw_output_file = xfopen(raw_output_fname.c_str(), "wt");
    }

}

DarcyFlowMHOutput::~DarcyFlowMHOutput(){
    if (output_writer != NULL) delete output_writer;

    if (balance_output_file != NULL) xfclose(balance_output_file);
    if (raw_output_file != NULL) xfclose(raw_output_file);
};





//=============================================================================
// CONVERT SOLUTION, CALCULATE BALANCES, ETC...
//=============================================================================

void DarcyFlowMHOutput::postprocess() {

    //make_side_flux();

    /*  writes scalar values to mesh - cannot be moved to output!
     *  all other methods are moved to output
     */
    make_element_scalar();
//    make_element_scalar(ele_scalars);

//    make_element_vector();
//    make_sides_scalar();

    /* new version of make_node_scalar */
//    make_node_scalar_param(node_scalars);


//    make_neighbour_flux();
//    water_balance();
}

void DarcyFlowMHOutput::output()
{
    std::string eleVectorName = "velocity_elements";
    std::string eleVectorUnit = "L/T";

    unsigned int result = 0;

    if (darcy_flow->time().is_current(output_mark_type)) {

        make_element_vector();
        //make_sides_scalar();

        make_node_scalar_param(node_pressure);

        //make_neighbour_flux();

        water_balance();

        result = output_writer->register_node_data
                ("pressure_nodes","L", node_pressure);
        //xprintf(Msg, "Register_node_data - result: %i, node size: %i\n", result,  mesh_->node_vector.size());

        result = output_writer->register_elem_data
                ("pressure_elements","L",ele_pressure);
        //xprintf(Msg, "Register_elem_data scalars - result: %i\n", result);

        if (output_piezo_head) {
            result = output_writer->register_elem_data
                        ("piezo_head_elements","L",ele_piezo_head);
        }

        result = output_writer->register_elem_data("velocity_elements", "L/T", ele_flux);
        //xprintf(Msg, "Register_elem_data vectors - result: %i\n", result);

        //double time  = min(darcy_flow->solved_time(), 1.0E200);
        double time  = darcy_flow->solved_time();

        // Workaround for infinity time returned by steady solvers. Should be designed better. Maybe
        // consider begining of the interval of actual result as the output time. Or use
        // particular TimeMark. This can allow also interpolation and perform output even inside of time step interval.
        if (time == TimeGovernor::inf_time) time = 0.0;
        output_writer->write_data(time);

        output_internal_flow_data();
    }
}

//=============================================================================
// FILL TH "FLUX" FIELD FOR ALL SIDES IN THE MESH
//=============================================================================
/*
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
}*/
//=============================================================================
// FILL TH "SCALAR" FIELD FOR ALL ELEMENTS IN THE MESH
//=============================================================================

void DarcyFlowMHOutput::make_element_scalar() {
    unsigned int sol_size;
    double *sol;

    darcy_flow->get_solution_vector(sol, sol_size);
    unsigned int soi = mesh_->n_sides();
    unsigned int i = 0;
    FOR_ELEMENTS(mesh_,ele) {
        ele_pressure[i] = sol[ soi];
        if (output_piezo_head) ele_piezo_head[i] = sol[soi ] + ele->centre()[Mesh::z_coord];
        i++; soi++;
    }
}

/*
void DarcyFlowMHOutput::make_element_scalar(double* scalars) {
    int soi;
    unsigned int sol_size;
    double *sol;
    int ele_index = 0; //!< index of each element

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
*/

/****
 * compute Darcian velocity in centre of elements
 *
 */
void DarcyFlowMHOutput::make_element_vector() {
    const MH_DofHandler &dh = darcy_flow->get_mh_dofhandler();
    MHFEValues fe_values;

    int i_side=0;
    FOR_ELEMENTS(mesh_, ele) {
        arma::vec3 flux_in_centre;
        flux_in_centre.zeros();

        fe_values.update(ele);

        for (int li = 0; li < ele->n_sides(); li++) {
            flux_in_centre += dh.side_flux( *(ele->side( li ) ) )
                              * fe_values.RT0_value( ele, ele->centre(), li )
                              / ele->material->size;
        }

        for(int j=0;j<3;j++) ele_flux[i_side][j]=flux_in_centre[j];
        i_side++;
    }

}

//=============================================================================
//
//=============================================================================
/*
void DarcyFlowMHOutput::make_element_vector_line(ElementFullIter ele, arma::vec3 &vec) {
    const MH_DofHandler &dh = darcy_flow->get_mh_dofhandler();
    SideIter s1=ele->side(1);
    SideIter s0=ele->side(0);
    double darcy_vel =
            (dh.side_flux( *s1 ) - dh.side_flux( *s0 )) / 2.0 / ele->material->size;

    // normalize element vector [node 0, node 1]
    vec = ele->node[1]->point() - ele->node[0]->point();
    vec *= darcy_vel / arma::norm(vec, 2);
}*/
//=============================================================================
//
//=============================================================================
/*
void DarcyFlowMHOutput::make_element_vector_triangle(ElementFullIter ele, arma::vec3 &vec) {

    const MH_DofHandler & dh = darcy_flow->get_mh_dofhandler();
    double bas[ 3 ][ 3 ];
    double X[ 3 ];
    int i, li;

    // begin -- upravy Ji. -- prepocet vektoru do teziste elementu
    // 23.2.2007

    // make rotated coordinate system with triangle in plane XY, origin in A and axes X == AB
    arma::vec3 ex(ele->node[1]->point() - ele->node[0]->point());
    ex /= norm(ex,2);

    arma::vec3 ac(ele->node[2]->point() - ele->node[0]->point());
    arma::vec3 ez = cross(ex, ac);
    ez /= norm(ez,2);

    arma::vec3 ey = cross(ez,ex);
    ey /= norm(ey, 2);

    // compute barycenter in new coordinate system
    arma::vec3 u = ele->centre() - ele->node[0]->point();

    // compute flux form base functions
    ac.zeros();
    for (i = 0; i < 3; i++) {
        bas[ i ][ 0 ] = ele->bas_gama[ i ] * (dot(u,ex) - ele->bas_alfa[ i ]);
        bas[ i ][ 1 ] = ele->bas_gama[ i ] * (dot(u,ey) - ele->bas_beta[ i ]);
        bas[ i ][ 2 ] = 0;
        for (li = 0; li < 2; li++)
            ac[ li ] += dh.side_flux( *(ele->side( i ) ) ) * bas[ i ][ li ] / ele->material->size;
    }

    vec = ac[0] * ex + ac[1] * ey + ac[2] * ez;
}*/
//=============================================================================
//
//=============================================================================


//=============================================================================
// FILL TH "SCALAR" FIELD FOR ALL INTERNAL SIDES IN THE MESH
//=============================================================================
/*
void DarcyFlowMHOutput::make_sides_scalar() {
    double *sol;
    int soi, si;
    unsigned int sol_size;
    struct Side *sde;

    soi = mesh_->n_sides() + mesh_->n_elements();
    darcy_flow->get_solution_vector(sol, sol_size);

    FOR_EDGES(mesh_, edg) {
        for (si = 0; si < edg->n_sides; si++) {
            sde = edg->side[ si ];
            sde->scalar = sol[ soi ];
        }
        soi++;
    }
}
*/


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

void DarcyFlowMHOutput::make_node_scalar_param(std::vector<double> &scalars) {
    F_ENTRY;

    double dist; //!< tmp variable for storing particular distance node --> element, node --> side*/

    /** Iterators */
    const Node * node;
    ElementIter ele;
    //struct Side* side;

    int n_nodes = mesh_->node_vector.size(); //!< number of nodes in the mesh */
    int node_index = 0; //!< index of each node */

    int* sum_elements = new int [n_nodes]; //!< sum elements joined to node */
    int* sum_sides = new int [n_nodes]; //!< sum sides joined to node */
    double* sum_ele_dist = new double [n_nodes]; //!< sum distances to all joined elements */
    double* sum_side_dist = new double [n_nodes]; //!<  Sum distances to all joined sides */

    /** tmp variables, will be replaced by ini keys
     * TODO include them into ini file*/
    bool count_elements = true; //!< scalar is counted via elements*/
    bool count_sides = true; //!< scalar is counted via sides */


    const MH_DofHandler &dh = darcy_flow->get_mh_dofhandler();

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
            for (int li = 0; li < ele->n_nodes(); li++) {
                node = ele->node[li]; //!< get Node pointer from element */
                node_index = mesh_->node_vector.index(node); //!< get nod index from mesh */

                dist = sqrt(
                        ((node->getX() - ele->centre()[ 0 ])*(node->getX() - ele->centre()[ 0 ])) +
                        ((node->getY() - ele->centre()[ 1 ])*(node->getY() - ele->centre()[ 1 ])) +
                        ((node->getZ() - ele->centre()[ 2 ])*(node->getZ() - ele->centre()[ 2 ]))
                );
                sum_ele_dist[node_index] += dist;
                sum_elements[node_index]++;
            }
    }
    if (count_sides){
        FOR_SIDES(mesh_, side) {
            for (int li = 0; li < side->n_nodes(); li++) {
                node = side->node(li);//!< get Node pointer from element */
                node_index = mesh_->node_vector.index(node); //!< get nod index from mesh */
                dist = sqrt(
                        ((node->getX() - side->centre()[ 0 ])*(node->getX() - side->centre()[ 0 ])) +
                        ((node->getY() - side->centre()[ 1 ])*(node->getY() - side->centre()[ 1 ])) +
                        ((node->getZ() - side->centre()[ 2 ])*(node->getZ() - side->centre()[ 2 ]))
                );

                sum_side_dist[node_index] += dist;
                sum_sides[node_index]++;
            }
        }
    }

    /**second pass - calculate scalar  */
    if (count_elements){
        FOR_ELEMENTS(mesh_, ele)
            for (int li = 0; li < ele->n_nodes(); li++) {
                node = ele->node[li];//!< get Node pointer from element */
                node_index = mesh_->node_vector.index(node); //!< get nod index from mesh */

                /**TODO - calculate it again or store it in prior pass*/
                dist = sqrt(
                        ((node->getX() - ele->centre()[ 0 ])*(node->getX() - ele->centre()[ 0 ])) +
                        ((node->getY() - ele->centre()[ 1 ])*(node->getY() - ele->centre()[ 1 ])) +
                        ((node->getZ() - ele->centre()[ 2 ])*(node->getZ() - ele->centre()[ 2 ]))
                );
                scalars[node_index] += ele_pressure[ele.index()] *
                        (1 - dist / (sum_ele_dist[node_index] + sum_side_dist[node_index])) /
                        (sum_elements[node_index] + sum_sides[node_index] - 1);
            }
    }
    if (count_sides) {
        FOR_SIDES(mesh_, side) {
            for (int li = 0; li < side->n_nodes(); li++) {
                node = side->node(li);//!< get Node pointer from element */
                node_index = mesh_->node_vector.index(node); //!< get nod index from mesh */

                /**TODO - calculate it again or store it in prior pass*/
                dist = sqrt(
                        ((node->getX() - side->centre()[ 0 ])*(node->getX() - side->centre()[ 0 ])) +
                        ((node->getY() - side->centre()[ 1 ])*(node->getY() - side->centre()[ 1 ])) +
                        ((node->getZ() - side->centre()[ 2 ])*(node->getZ() - side->centre()[ 2 ]))
                );


                scalars[node_index] += dh.side_scalar( *side ) *
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


/*
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
                ((nod->getX() - ele->centre()[ 0 ])*(nod->getX() - ele->centre()[ 0 ])) +
                ((nod->getY() - ele->centre()[ 1 ])*(nod->getY() - ele->centre()[ 1 ])) +
                ((nod->getZ() - ele->centre()[ 2 ])*(nod->getZ() - ele->centre()[ 2 ]))
                );

            TED[ele.index()][li] = dist;
            nod->faux += dist; //       nod->faux += 1 / dist;
            nod->aux++;
        }

    FOR_SIDES(mesh_, sde) {
        for (li = 0; li < sde->n_nodes; li++) {
            nod = sde->node[li];

            dist = sqrt(
                    ((nod->getX() - sde->centre()[ 0 ])*(nod->getX() - sde->centre()[ 0 ])) +
                    ((nod->getY() - sde->centre()[ 1 ])*(nod->getY() - sde->centre()[ 1 ])) +
                    ((nod->getZ() - sde->centre()[ 2 ])*(nod->getZ() - sde->centre()[ 2 ]))
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
*/
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
/*
void DarcyFlowMHOutput::make_neighbour_flux() {
    struct Neighbour *ngh;

    FOR_NEIGHBOURS(mesh_, ngh) {
        if (ngh->type != VV_2E)
            continue;

        ngh->flux = ngh->sigma * ngh->geom_factor * (ngh->element[1]->scalar - ngh->element[0]->scalar);
        continue;
    }
}
*/

//=============================================================================
//
//=============================================================================

void DarcyFlowMHOutput::water_balance() {
    F_ENTRY;

    if (balance_output_file == NULL) return;
    const MH_DofHandler &dh = darcy_flow->get_mh_dofhandler();

    double bal;
    int c_water;
    struct Boundary *bcd;
    std::vector<double> *bcd_balance = new std::vector<double>( mesh_->bcd_group_id.size(), 0.0 );
    std::vector<double> *bcd_plus_balance = new std::vector<double>( mesh_->bcd_group_id.size(), 0.0 );
    std::vector<double> *bcd_minus_balance = new std::vector<double>( mesh_->bcd_group_id.size(), 0.0 );

    fprintf(balance_output_file,"********************************\n");
    fprintf(balance_output_file,"Boundary fluxes at time %f:\n",darcy_flow->time().t());
    fprintf(balance_output_file,"[total balance]    [total outflow]     [total inflow]\n");
    FOR_BOUNDARIES(mesh_, bcd) {
        double flux = dh.side_flux( *(bcd->side) );
        (*bcd_balance)[bcd->group] += flux;

        if (flux > 0) (*bcd_plus_balance)[bcd->group]+= flux;
        else (*bcd_minus_balance)[bcd->group]+= flux;
    }



    for(int i=0; i < bcd_balance->size(); ++i)
        fprintf(balance_output_file, "boundary #%d\t%g\t%g\t%g\n", mesh_->bcd_group_id(i).id(),
                (*bcd_balance)[i],(*bcd_plus_balance)[i],(*bcd_minus_balance)[i]);

    delete bcd_balance;
    delete bcd_plus_balance;
    delete bcd_minus_balance;

    const FieldP0<double> *p_sources=darcy_flow->get_sources();
    if (p_sources != NULL) {

        fprintf(balance_output_file,"\nSource fluxes over material subdomains:\n");
        MaterialDatabase &mat_base = darcy_flow->material_base();
        std::vector<double> *src_balance = new std::vector<double>( mat_base.size(), 0.0 ); // initialize by zero

        FOR_ELEMENTS(mesh_, elm) {
            (*src_balance)[mat_base.index(elm->material)] += elm->volume() * p_sources->element_value(elm.index());
        }


        FOR_MATERIALS_IT(mat_base, mat) {
            fprintf(balance_output_file, "material #%d:\t% g\n", mat_base.get_id(mat), (*src_balance)[mat_base.index(mat)]);
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

/*
 * Output of internal flow data.
 */

void DarcyFlowMHOutput::output_internal_flow_data()
{
    if (raw_output_file == NULL) return;

    char dbl_fmt[ 16 ]= "%.8g ";
    // header
    xfprintf( raw_output_file, "$FlowField\nT=");
    xfprintf( raw_output_file, dbl_fmt, darcy_flow->time().t());
    xfprintf( raw_output_file, "\n%d\n", mesh_->n_elements() );

    const MH_DofHandler &dh = darcy_flow->get_mh_dofhandler();

    int i;
    int cit = 0;
    FOR_ELEMENTS( mesh_,  ele ) {
        //xfprintf( raw_output_file, "%d ", cit);
        xfprintf( raw_output_file, "%d ", ele.id());
        xfprintf( raw_output_file, dbl_fmt, ele_pressure[cit]);
        for (i = 0; i < 3; i++)
            xfprintf( raw_output_file, dbl_fmt, ele_flux[cit][i]);

        xfprintf( raw_output_file, " %d ", ele->n_sides());
        for (i = 0; i < ele->n_sides(); i++)
            xfprintf( raw_output_file, dbl_fmt, dh.side_scalar( *(ele->side(i) ) ) );
        for (i = 0; i < ele->n_sides(); i++)
            xfprintf( raw_output_file, dbl_fmt, dh.side_flux( *(ele->side(i) ) ) );

        //xfprintf( raw_output_file, "%d ", ele->n_neighs_vv);
        //for (i = 0; i < ele->n_neighs_vv; i++)
        //    xfprintf( raw_output_file, "%d ", ele->neigh_vv[i]->id);

        xfprintf( raw_output_file, "\n" );
        cit ++;
    }
    xfprintf( raw_output_file, "$EndFlowField\n\n" );
}
