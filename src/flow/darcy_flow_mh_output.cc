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
#include "system/system.hh"
#include "system/sys_profiler.hh"

#include <vector>

#include "io/output.h"
#include <iostream>
#include <sstream>
#include <string>
#include "mesh/partitioning.hh"
namespace it = Input::Type;

it::Record DarcyFlowMHOutput::input_type
	= it::Record("DarcyMHOutput", "Parameters of MH output.")
	.declare_key("save_step", it::Double(0.0), it::Default("1.0"),
                    "Regular step between MH outputs.")
    .declare_key("output_stream", OutputTime::input_type, it::Default::obligatory(),
                    "Parameters of output stream.")
    .declare_key("velocity_p0", it::String(),
                    "Output stream for P0 approximation of the velocity field.")
    .declare_key("pressure_p0", it::String(),
                    "Output stream for P0 approximation of the pressure field.")
    .declare_key("pressure_p1", it::String(),
                    "Output stream for P1 approximation of the pressure field.")
    .declare_key("piezo_head_p0", it::String(),
                    "Output stream for P0 approximation of the piezometric head field.")
    .declare_key("balance_output", it::FileName::output(), it::Default("water_balance.txt"),
                    "Output file for water balance table.")
    .declare_key("subdomains", it::String(),
                    "Output stream for subdomain indices (partitioning of mesh elements) used by DarcyFlow module.")

    .declare_key("raw_flow_output", it::FileName::output(), it::Default::optional(),
                    "Output file with raw data form MH module.");



DarcyFlowMHOutput::DarcyFlowMHOutput(DarcyFlowMH *flow, Input::Record in_rec)
: darcy_flow(flow), mesh_(&darcy_flow->mesh()),
  in_rec_(in_rec),
  ele_flux(mesh_->n_elements(),std::vector<double>(3,0.0)),
  balance_output_file(NULL),raw_output_file(NULL)
{
    F_ENTRY;
    using namespace Input;
    unsigned int result = 0;
    
    // setup output
    // is created for every MPI process
    output_writer = OutputTime::output_stream(Record(in_rec).val<Record>("output_stream"));

    // allocate output containers
    ele_pressure.resize(mesh_->n_elements());
    node_pressure.resize(mesh_->node_vector.size());

    //local iterator it
    Iterator<string> it = in_rec.find<string>("piezo_head_p0");
    output_piezo_head=bool(it);
    DBGMSG("piezo set: %d \n", output_piezo_head);
      
    if (output_piezo_head) ele_piezo_head.resize(mesh_->n_elements());

    // set output time marks
    TimeMarks &marks = darcy_flow->time().marks();
    output_mark_type = darcy_flow->mark_type() | marks.type_fixed_time() | marks.type_output();
    marks.add_time_marks(0.0, in_rec.val<double>("save_step"),
          darcy_flow->time().end_time(), output_mark_type );
    DBGMSG("end create output\n");

      
    // testing MPI rank so that the files are opened only once
    /* When calling methods water_balance() and output_internal_flow_data()
     * it is tested if the file == NULL.
     * And it will be NULL for all processes except rank==0.
     */
    
    int ierr, rank;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    ASSERT(ierr == 0, "Error in MPI test of rank.");
    

	result = OutputTime::register_node_data
			(mesh_, "pressure_nodes", "L", in_rec.val<Input::Record>("output_stream"), node_pressure);

	result = OutputTime::register_elem_data
			(mesh_, "pressure_elements", "L", in_rec.val<Input::Record>("output_stream"), ele_pressure);

	if (output_piezo_head) {
		result = OutputTime::register_elem_data
				(mesh_, "piezo_head_elements", "L", in_rec.val<Input::Record>("output_stream"), ele_piezo_head);
	}

	result = OutputTime::register_elem_data
			(mesh_, "velocity_elements", "L/T", in_rec.val<Input::Record>("output_stream"), ele_flux);

	it = in_rec.find<string>("subdomains");
	if (bool(it)) {
		result = OutputTime::register_elem_data
				(mesh_, "subdomains", "", in_rec.val<Input::Record>("output_stream"), mesh_->get_part()->seq_output_partition() );
	}

	if (rank == 0) {
        // temporary solution for balance output
        balance_output_file = xfopen( in_rec.val<FilePath>("balance_output"), "wt");

        { // local iterator it
            // optionally open raw output file
            Iterator<FilePath> it = in_rec.find<FilePath>("raw_flow_output");

            if (it) {
                xprintf(Msg, "Opening raw output: %s\n", string(*it).c_str());
                raw_output_file = xfopen(*it, "wt");
            }
        }
    }
}



DarcyFlowMHOutput::~DarcyFlowMHOutput(){
    //if (output_writer != NULL) delete output_writer;

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
    START_TIMER("Darcy output");

    std::string eleVectorName = "velocity_elements";
    std::string eleVectorUnit = "L/T";

    //cout << "DMHO_output: rank: " << rank << "\t output_writer: " << output_writer << endl;
    
    // skip initial output for steady solver
    if (darcy_flow->time().is_steady() && darcy_flow->time().tlevel() ==0) return;

    if (darcy_flow->time().is_current(output_mark_type)) {

      make_element_vector();
      //make_sides_scalar();

      make_node_scalar_param(node_pressure);

      //make_neighbour_flux();

      DBGMSG("water_balance()\n");
      water_balance();

      //compute_l2_difference();

      double time  = darcy_flow->solved_time();

      // Workaround for infinity time returned by steady solvers. Should be designed better. Maybe
      // consider begining of the interval of actual result as the output time. Or use
      // particular TimeMark. This can allow also interpolation and perform output even inside of time step interval.
      if (time == TimeGovernor::inf_time) time = 0.0;

      if(output_writer) output_writer->write_data(time);
      
      output_internal_flow_data();
      
      //for synchronization when measuring time by Profiler
      MPI_Barrier(MPI_COMM_WORLD);
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

        fe_values.update(ele, darcy_flow->get_data().anisotropy, darcy_flow->get_data().cross_section, darcy_flow->get_data().conductivity );

        for (unsigned int li = 0; li < ele->n_sides(); li++) {
            flux_in_centre += dh.side_flux( *(ele->side( li ) ) )
                              * fe_values.RT0_value( ele, ele->centre(), li )
                              / darcy_flow->get_data().cross_section.value(ele->centre(), ele->element_accessor() );
        }

        for(unsigned int j=0; j<3; j++) 
            ele_flux[i_side][j]=flux_in_centre[j];
        
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
    //ElementIter ele;
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
            for (unsigned int li = 0; li < ele->n_nodes(); li++) {
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
            for (unsigned int li = 0; li < side->n_nodes(); li++) {
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
            for (unsigned int li = 0; li < ele->n_nodes(); li++) {
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
            for (unsigned int li = 0; li < side->n_nodes(); li++) {
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
}*/

//=============================================================================
//
//=============================================================================

void DarcyFlowMHOutput::water_balance() {
    F_ENTRY;
    const MH_DofHandler &dh = darcy_flow->get_mh_dofhandler();
    if (balance_output_file == NULL) return;

    //BOUNDARY
    //struct Boundary *bcd;
    std::vector<double> bcd_balance( mesh_->region_db().boundary_size(), 0.0 );
    std::vector<double> bcd_plus_balance( mesh_->region_db().boundary_size(), 0.0 );
    std::vector<double> bcd_minus_balance( mesh_->region_db().boundary_size(), 0.0 );

    using namespace std;
    //printing the head of water balance file
    unsigned int c = 5; //column number without label
    unsigned int w = 14;  //column width
    unsigned int wl = 2*(w-5)+7;  //label column width
    stringstream s; //helpful stringstream
    string bc_head_format = "# %-*s%-*s%-*s%-*s%-*s%-*s\n",
           bc_format = "%*s%-*d%-*s  %-*g%-*g%-*g%-*g\n",
           bc_total_format = "# %-*s%-*g%-*g%-*g\n\n\n";
    s << setw((w*c+wl-15)/2) << setfill('-') << "-"; //drawing half line
    fprintf(balance_output_file,"# %s WATER BALANCE %s\n",s.str().c_str(), s.str().c_str());
    fprintf(balance_output_file,"# Time of computed water balance: %f\n\n\n",darcy_flow->time().t());
    
    fprintf(balance_output_file,"# Boundary water balance:\n");
    fprintf(balance_output_file,bc_head_format.c_str(),w,"[boundary_id]",wl,"[label]",
                            w,"[total_balance]",w,"[total_outflow]",w,"[total_inflow]",w,"[time]");
    s.clear();
    s.str(std::string());
    s << setw(w*c+wl) << setfill('-') << "-"; 
    fprintf(balance_output_file,"# %s\n",s.str().c_str());  //drawing long line
    
    //computing water balance over boundaries
    FOR_BOUNDARIES(mesh_, bcd) {
        // !! there can be more sides per one boundary
        double flux = dh.side_flux( *(bcd->side()) );

        Region r = bcd->region();
        //DBGMSG("flux: %f side: %d %d reg: %s\n", flux, bcd->side()->element()->index(), bcd->side()->el_idx(), r.label().c_str() );
        if (! r.is_valid()) xprintf(Msg, "Invalid region, ele % d, edg: % d\n", bcd->bc_ele_idx_, bcd->edge_idx_);
        unsigned int bc_region_idx = r.boundary_idx();
        bcd_balance[bc_region_idx] += flux;

        if (flux > 0) bcd_plus_balance[bc_region_idx] += flux;
        else bcd_minus_balance[bc_region_idx] += flux;
    }
    //printing water balance over boundaries
    DBGMSG("DB[boundary] size: %u\n", mesh_->region_db().boundary_size());
    const RegionSet & b_set = mesh_->region_db().get_region_set("BOUNDARY");
    double total_balance = 0, // for computing total balance on boundary
           total_inflow = 0,
           total_outflow = 0; 
    for( RegionSet::const_iterator reg = b_set.begin(); reg != b_set.end(); ++reg) {
        //DBGMSG("writing reg->idx() and id() and boundary_idx(): %d\t%d\t%d\n", reg->idx(), reg->id(), reg->boundary_idx());
        total_balance += bcd_balance[reg->boundary_idx()];
        total_outflow += bcd_plus_balance[reg->boundary_idx()];
        total_inflow += bcd_minus_balance[reg->boundary_idx()];
        fprintf(balance_output_file, bc_format.c_str(),2,"",w,reg->id(),wl,reg->label().c_str(),
                w, bcd_balance[reg->boundary_idx()], w, bcd_plus_balance[reg->boundary_idx()],
                w, bcd_minus_balance[reg->boundary_idx()], w, darcy_flow->time().t());
    }
    //total boundary balance
    fprintf(balance_output_file,"# %s\n",s.str().c_str());  // drawing long line
    fprintf(balance_output_file, bc_total_format.c_str(),w+wl+2,"total boundary balance",
                w,total_balance, w, total_outflow, w, total_inflow);

    //SOURCES
    string src_head_format = "# %-*s%-*s%-*s%-*s%-*s\n",
           src_format = "%*s%-*d%-*s  %-*g%-*s%-*g\n",
           src_total_format = "# %-*s%-*g\n\n\n";
    //computing water balance of sources
    fprintf(balance_output_file,"# Source fluxes over material subdomains:\n");   //head
    fprintf(balance_output_file,src_head_format.c_str(),w,"[region_id]",wl,"[label]", 
                            w,"[total_balance]",2*w,"",w,"[time]");
    fprintf(balance_output_file,"# %s\n",s.str().c_str());  //long line
    std::vector<double> src_balance( mesh_->region_db().bulk_size(), 0.0 ); // initialize by zero
    FOR_ELEMENTS(mesh_, elm) {
      //DBGMSG("writing reg->idx() and id() and bulk_idx(): %d\t%d\t%d\n", 
      //       elm->element_accessor().region().idx(), 
      //       elm->element_accessor().region().id(),
      //       elm->element_accessor().region().bulk_idx());
      src_balance[elm->element_accessor().region().bulk_idx()] += elm->measure() * 
            darcy_flow->get_data().cross_section.value(elm->centre(), elm->element_accessor()) * 
            darcy_flow->get_data().water_source_density.value(elm->centre(), elm->element_accessor());
    }
  
    total_balance = 0;
    //printing water balance of sources
    DBGMSG("DB[bulk] size: %u\n", mesh_->region_db().bulk_size());
    const RegionSet & bulk_set = mesh_->region_db().get_region_set("BULK");
    for( RegionSet::const_iterator reg = bulk_set.begin(); reg != bulk_set.end(); ++reg)
      {
        total_balance += src_balance[reg->bulk_idx()];
        //"%*s%-*d%-*s  %-*g%-*s%-*g\n";
        fprintf(balance_output_file, src_format.c_str(), 2,"", w, reg->id(), wl,
                reg->label().c_str(), w, src_balance[reg->bulk_idx()],2*w,"", w,darcy_flow->time().t());
      }
    //total sources balance
    fprintf(balance_output_file,"# %s\n",s.str().c_str());  //drawing long line
    fprintf(balance_output_file, src_total_format.c_str(),w+wl+2,"total sources balance",
                w,total_balance);
}

//=============================================================================
// UNUSED FUNCTION
//=============================================================================
/*
double calc_water_balance(Mesh* mesh, int c_water) {
    double rc;
    struct Boundary *bcd;

    rc = 0.0;

    return rc;
}
*/

/*
 * Output of internal flow data.
 */

void DarcyFlowMHOutput::output_internal_flow_data()
{
    const MH_DofHandler &dh = darcy_flow->get_mh_dofhandler();

    if (raw_output_file == NULL) return;
    
    char dbl_fmt[ 16 ]= "%.8g ";
    // header
    xfprintf( raw_output_file, "// fields:\n//ele_id    ele_presure    flux_in_barycenter[3]    n_sides   side_pressures[n]    side_fluxes[n]\n");
    xfprintf( raw_output_file, "$FlowField\nT=");
    xfprintf( raw_output_file, dbl_fmt, darcy_flow->time().t());
    xfprintf( raw_output_file, "\n%d\n", mesh_->n_elements() );

    unsigned int i;
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

/*
#include "quadrature/quadrature_lib.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/mapping_p1.hh"
#include "functions/function_python.hh"
*/

/*
* Calculate approximation of L2 norm for:
 * 1) difference between regularized pressure and analytical solution (using FunctionPython)
 * 2) difference between RT velocities and analytical solution
 * 3) difference of divergence
 * */
/*
struct DiffData {
    double pressure_error[2], velocity_error[2], div_error[2];
    vector<double> pressure_diff;
    vector<double> velocity_diff;
    vector<double> div_diff;

    double * solution;
    const MH_DofHandler * dh;
    MHFEValues fe_values;

    std::vector< std::vector<double>  > *ele_flux;
};

template <int dim>
void l2_diff_local(ElementFullIter &ele, FEValues<dim,3> &fe_values, const FunctionPython<3> &anal_sol,  DiffData &result) {

    fe_values.reinit(ele);
    result.fe_values.update(ele);

    // get coefficients on the current element
    vector<double> fluxes(dim+1);
    vector<double> pressure_traces(dim+1);

    for (unsigned int li = 0; li < ele->n_sides(); li++) {
        fluxes[li] = result.dh->side_flux( *(ele->side( li ) ) );
        pressure_traces[li] = result.dh->side_scalar( *(ele->side( li ) ) );
    }
    double pressure_mean = result.dh->element_scalar(ele);

    vector<double> analytical(5); // values tuple (pressure, flux(dim), divergence)
    arma::vec3 flux_in_q_point;
    arma::vec3 anal_flux;

    double velocity_diff=0, divergence_diff=0, pressure_diff=0, diff;
    double resistance = ele->material->hydrodynamic_resistence[0];  // assumes isotropic medium

    // 1d:  mean_x_squared = 1/6 (v0^2 + v1^2 + v0*v1)
    // 2d:  mean_x_squared = 1/12 (v0^2 + v1^2 +v2^2 + v0*v1 + v0*v2 + v1*v2)
    double mean_x_squared=0;
    for(unsigned int i_node=0; i_node < ele->n_nodes(); i_node++ )
        for(unsigned int j_node=0; j_node < ele->n_nodes(); j_node++ )
        {
            mean_x_squared += (i_node == j_node ? 2.0 : 1.0) / ( 6 * dim )   // multiply by 2 on diagonal
                    * arma::dot( ele->node[i_node]->point(), ele->node[j_node]->point());
        }

    for(unsigned int i_point=0; i_point < fe_values.n_points(); i_point++) {
        arma::vec3 q_point = fe_values.point(i_point);

        anal_sol.vector_value(q_point, analytical);
        for(unsigned int i=0; i< 3; i++) anal_flux[i] = analytical[i+1];

        // compute postprocesed pressure
        diff = 0;
        for(unsigned int i_shape=0; i_shape < ele->n_sides(); i_shape++) {
            unsigned int oposite_node = (i_shape + dim) % (dim + 1);

            diff += fluxes[ i_shape ] *
                               (  arma::dot( q_point, q_point )/ 2
                                - mean_x_squared / 2
                                - arma::dot( q_point, ele->node[oposite_node]->point() )
                                + arma::dot( ele->centre(), ele->node[oposite_node]->point() )
                               );
        }

        diff = - resistance * diff / dim / ele->volume() + pressure_mean ;
        diff = ( diff - analytical[0]);
        pressure_diff += diff * diff * fe_values.JxW(i_point);


        // velocity difference
        flux_in_q_point.zeros();
        for(unsigned int i_shape=0; i_shape < ele->n_sides(); i_shape++) {
            flux_in_q_point += fluxes[ i_shape ]
                              * result.fe_values.RT0_value( ele, q_point, i_shape )
                              / ele->material->size;
        }

        flux_in_q_point -= anal_flux;
        velocity_diff += dot(flux_in_q_point, flux_in_q_point) * fe_values.JxW(i_point);

        // divergence diff
        diff = 0;
        for(unsigned int i_shape=0; i_shape < ele->n_sides(); i_shape++) diff += fluxes[ i_shape ];
        diff = ( diff / ele->volume() - analytical[4]);
        divergence_diff += diff * diff * fe_values.JxW(i_point);

    }


    result.velocity_diff[ele.index()] = velocity_diff;
    result.velocity_error[dim-1] += velocity_diff;

    result.pressure_diff[ele.index()] = pressure_diff;
    result.pressure_error[dim-1] += pressure_diff;

    result.div_diff[ele.index()] = divergence_diff;
    result.div_error[dim-1] += divergence_diff;

}





*/
void DarcyFlowMHOutput::compute_l2_difference() {
}
/*

    const unsigned int order = 4; // order of Gauss quadrature

    // we create trivial Dofhandler , for P0 elements, to get access to, FEValues on individual elements
    // this we use to integrate our own functions - difference of postprocessed pressure and analytical solution
    FE_P_disc<0,1,3> fe_1d;
    FE_P_disc<0,2,3> fe_2d;

    QGauss<1> quad_1d( order );
    QGauss<2> quad_2d( order );

    MappingP1<1,3> mapp_1d;
    MappingP1<2,3> mapp_2d;

    FEValues<1,3> fe_values_1d(mapp_1d, quad_1d,   fe_1d, update_JxW_values | update_quadrature_points);
    FEValues<2,3> fe_values_2d(mapp_2d, quad_2d,   fe_2d, update_JxW_values | update_quadrature_points);

    FunctionPython<3> anal_sol_1d(5);   // components: pressure, flux vector 3d, divergence
    anal_sol_1d.set_python_function_from_file("analytical_module.py", "all_values_1d");

    FunctionPython<3> anal_sol_2d(5);
    anal_sol_2d.set_python_function_from_file("analytical_module.py", "all_values_2d");


    static DiffData result;

    result.pressure_diff.resize( mesh_->n_elements() );
    result.velocity_diff.resize( mesh_->n_elements() );
    result.div_diff.resize( mesh_->n_elements() );

    result.pressure_error[0] = 0;
    result.velocity_error[0] = 0;
    result.div_error[0] = 0;
    result.pressure_error[1] = 0;
    result.velocity_error[1] = 0;
    result.div_error[1] = 0;

    result.ele_flux = &( ele_flux );

    output_writer->register_elem_data("pressure_diff","0",result.pressure_diff);
    output_writer->register_elem_data("velocity_diff","0",result.velocity_diff);
    output_writer->register_elem_data("div_diff","0",result.div_diff);

    unsigned int solution_size;
    darcy_flow->get_solution_vector(result.solution, solution_size);

    result.dh = &( darcy_flow->get_mh_dofhandler());

    FOR_ELEMENTS( mesh_, ele) {
        switch (ele->dim()) {
        case 1:
            l2_diff_local<1>( ele, fe_values_1d, anal_sol_1d, result);
            break;
        case 2:
            l2_diff_local<2>( ele, fe_values_2d, anal_sol_2d, result);
            break;
        }
    }

    xprintf(Msg,
            "\n"
            "pressure error 1d: %g\n"
            "pressure error 2d: %g\n"
            "velocity error 1d: %g\n"
            "velocity error 2d: %g\n"
            "div error: %g (1d) %g (2d)\n",
            sqrt(result.pressure_error[0]), sqrt(result.pressure_error[1]),
            sqrt(result.velocity_error[0]), sqrt(result.velocity_error[1]),
            sqrt(result.div_error[0]), sqrt(result.div_error[1])
            );
}
*/
