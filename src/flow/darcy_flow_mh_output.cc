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


#include <vector>
#include <iostream>
#include <sstream>
#include <string>

#include <system/global_defs.h>

#include "flow/mh_fe_values.hh"
#include "flow/darcy_flow_mh.hh"
#include "flow/darcy_flow_mh_output.hh"
#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/xio.h"

#include "fem/dofhandler.hh"
#include "fields/field_fe.hh"

#include "io/output.h"
#include "mesh/partitioning.hh"



namespace it = Input::Type;

it::Selection DarcyFlowMHOutput::OutputFields::output_selection
	= DarcyFlowMH::EqData().make_output_field_selection("DarcyMHOutput_Selection", "Selection of fields available for output.")
	.copy_values(OutputFields().make_output_field_selection("").close())
	.close();

it::Record DarcyFlowMHOutput::input_type
	= it::Record("DarcyMHOutput", "Parameters of MH output.")
    .declare_key("output_stream", OutputTime::input_type, it::Default::obligatory(),
                    "Parameters of output stream.")
    .declare_key("output_fields", it::Array(OutputFields::output_selection),
    		it::Default::obligatory(), "List of fields to write to output file.")
    .declare_key("balance_output", it::FileName::output(), it::Default("water_balance.txt"),
                    "Output file for water balance table.")
    .declare_key("compute_errors", it::Bool(), it::Default("false"),
    				"SPECIAL PURPOSE. Computing errors pro non-compatible coupling.")
    .declare_key("raw_flow_output", it::FileName::output(), it::Default::optional(),
                    "Output file with raw data form MH module.");


DarcyFlowMHOutput::OutputFields::OutputFields()
{
	*this += field_ele_pressure.name("pressure_p0").units(UnitSI().m());
	*this += field_node_pressure.name("pressure_p1").units(UnitSI().m());
	*this += field_ele_piezo_head.name("piezo_head_p0").units(UnitSI().m());
	*this += field_ele_flux.name("velocity_p0").units(UnitSI().m().s(-1));
	*this += subdomain.name("subdomain").units( UnitSI::dimensionless() );

	fields_for_output += *this;

	*this += pressure_diff.name("pressure_diff").units(UnitSI().m());
	*this += velocity_diff.name("velocity_diff").units(UnitSI().m().s(-1));
	*this += div_diff.name("div_diff").units(UnitSI().s(-1));

}


DarcyFlowMHOutput::DarcyFlowMHOutput(DarcyFlowMH_Steady *flow, Input::Record in_rec)
: darcy_flow(flow),
  mesh_(&darcy_flow->mesh()),
  in_rec_(in_rec),
  balance_output_file(NULL),
  raw_output_file(NULL)
{
    
    

	// we need to add data from the flow equation at this point, not in constructor of OutputFields
	output_fields.fields_for_output += darcy_flow->data();
	output_fields.set_mesh(*mesh_);

	// create shared pointer to a FieldElementwise and push this Field to output_field on all regions
	ele_pressure.resize(mesh_->n_elements());
	auto ele_pressure_ptr=make_shared< FieldElementwise<3, FieldValue<3>::Scalar> >(ele_pressure, 1);
	output_fields.field_ele_pressure.set_field(mesh_->region_db().get_region_set("ALL"), ele_pressure_ptr);

	dh = new DOFHandlerMultiDim(*mesh_);
	dh->distribute_dofs(fe1, fe2, fe3);
	corner_pressure.resize(dh->n_global_dofs());
	VecCreateSeqWithArray(PETSC_COMM_SELF, 1, dh->n_global_dofs(), &(corner_pressure[0]), &vec_corner_pressure);

	auto corner_ptr = make_shared< FieldFE<3, FieldValue<3>::Scalar> >();
	corner_ptr->set_fe_data(dh, &map1, &map2, &map3, &vec_corner_pressure);

	output_fields.field_node_pressure.set_field(mesh_->region_db().get_region_set("ALL"), corner_ptr);
	output_fields.field_node_pressure.output_type(OutputTime::NODE_DATA);

	ele_piezo_head.resize(mesh_->n_elements());
	output_fields.field_ele_piezo_head.set_field(mesh_->region_db().get_region_set("ALL"),
			make_shared< FieldElementwise<3, FieldValue<3>::Scalar> >(ele_piezo_head, 1));

	ele_flux.resize(3*mesh_->n_elements());
	output_fields.field_ele_flux.set_field(mesh_->region_db().get_region_set("ALL"),
			make_shared< FieldElementwise<3, FieldValue<3>::VectorFixed> >(ele_flux, 3));

	auto &vec_int_sub = mesh_->get_part()->seq_output_partition();
	subdomains.resize(vec_int_sub.size());
	for(unsigned int i=0; i<subdomains.size();i++)
		subdomains[i]=vec_int_sub[i];
	output_fields.subdomain.set_field(mesh_->region_db().get_region_set("ALL"),
			make_shared< FieldElementwise<3, FieldValue<3>::Integer> >(subdomains, 1));

	output_fields.set_limit_side(LimitSide::right);


	output_stream = OutputTime::create_output_stream(in_rec.val<Input::Record>("output_stream"));
	output_stream->add_admissible_field_names(in_rec.val<Input::Array>("output_fields"), OutputFields::output_selection);
	output_stream->mark_output_times(darcy_flow->time());

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
        // temporary solution for balance output
        balance_output_file = xfopen( in_rec.val<FilePath>("balance_output"), "wt");

        // optionally open raw output file
        FilePath raw_output_file_path;
        if (in_rec.opt_val("raw_flow_output", raw_output_file_path)) {
            xprintf(Msg, "Opening raw output: %s\n", string(raw_output_file_path).c_str() );
            raw_output_file = xfopen(raw_output_file_path, "wt");
        }

    }
}



DarcyFlowMHOutput::~DarcyFlowMHOutput(){

    if (balance_output_file != NULL) xfclose(balance_output_file);
    if (raw_output_file != NULL) xfclose(raw_output_file);
    VecDestroy(&vec_corner_pressure);
    delete output_stream;

    delete dh;
};





//=============================================================================
// CONVERT SOLUTION, CALCULATE BALANCES, ETC...
//=============================================================================


void DarcyFlowMHOutput::output()
{
    START_TIMER("Darcy output");

    if (darcy_flow->time().is_current( TimeGovernor::marks().type_output() )) {

      make_element_scalar();
      make_element_vector();

      make_node_scalar_param();

      water_balance();

      if (in_rec_.val<bool>("compute_errors")) compute_l2_difference();

	  output_fields.fields_for_output.set_time(darcy_flow->time());
	  output_fields.fields_for_output.output(output_stream);
	  output_stream->write_time_frame();

      output_internal_flow_data();

    }
    
}


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
        ele_piezo_head[i] = sol[soi ] + ele->centre()[Mesh::z_coord];
        i++; soi++;
    }
}



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

        fe_values.update(ele,
        		darcy_flow->data_.anisotropy,
        		darcy_flow->data_.cross_section,
        		darcy_flow->data_.conductivity );

        for (unsigned int li = 0; li < ele->n_sides(); li++) {
            flux_in_centre += dh.side_flux( *(ele->side( li ) ) )
                              * fe_values.RT0_value( ele, ele->centre(), li )
                              / darcy_flow->data_.cross_section.value(ele->centre(), ele->element_accessor() );
        }

        for(unsigned int j=0; j<3; j++) 
            ele_flux[3*i_side+j]=flux_in_centre[j];
        
        i_side++;
    }

}


void DarcyFlowMHOutput::make_corner_scalar(vector<double> &node_scalar)
{
	unsigned int ndofs = max(dh->fe<1>()->n_dofs(), max(dh->fe<2>()->n_dofs(), dh->fe<3>()->n_dofs()));
	unsigned int indices[ndofs];
	unsigned int i_node;
	FOR_ELEMENTS(mesh_, ele)
	{
		dh->get_dof_indices(ele, indices);
		FOR_ELEMENT_NODES(ele, i_node)
		{
			corner_pressure[indices[i_node]] = node_scalar[mesh_->node_vector.index(ele->node[i_node])];
		}
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

void DarcyFlowMHOutput::make_node_scalar_param() {

	vector<double> scalars(mesh_->n_nodes());

    double dist; //!< tmp variable for storing particular distance node --> element, node --> side*/

    /** Iterators */
    const Node * node;

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

    /** free memory */
    delete [] sum_elements;
    delete [] sum_sides;
    delete [] sum_ele_dist;
    delete [] sum_side_dist;

    make_corner_scalar(scalars);
}


//=============================================================================
//
//=============================================================================

void DarcyFlowMHOutput::water_balance() {
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
        if (! r.is_valid()) xprintf(Msg, "Invalid region, ele % d, edg: % d\n", bcd->bc_ele_idx_, bcd->edge_idx_);
        unsigned int bc_region_idx = r.boundary_idx();
        bcd_balance[bc_region_idx] += flux;

        if (flux > 0) bcd_plus_balance[bc_region_idx] += flux;
        else bcd_minus_balance[bc_region_idx] += flux;
    }
    //printing water balance over boundaries
    const RegionSet & b_set = mesh_->region_db().get_region_set("BOUNDARY");
    double total_balance = 0, // for computing total balance on boundary
           total_inflow = 0,
           total_outflow = 0; 
    for( RegionSet::const_iterator reg = b_set.begin(); reg != b_set.end(); ++reg) {
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
      src_balance[elm->element_accessor().region().bulk_idx()] += elm->measure() * 
            darcy_flow->data_.cross_section.value(elm->centre(), elm->element_accessor()) *
            darcy_flow->data_.water_source_density.value(elm->centre(), elm->element_accessor());
    }
  
    total_balance = 0;
    //printing water balance of sources
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
        xfprintf( raw_output_file, "%d ", ele.id());
        xfprintf( raw_output_file, dbl_fmt, ele_pressure[cit]);
        for (i = 0; i < 3; i++)
            xfprintf( raw_output_file, dbl_fmt, ele_flux[3*cit+i]);

        xfprintf( raw_output_file, " %d ", ele->n_sides());
        for (i = 0; i < ele->n_sides(); i++)
            xfprintf( raw_output_file, dbl_fmt, dh.side_scalar( *(ele->side(i) ) ) );
        for (i = 0; i < ele->n_sides(); i++)
            xfprintf( raw_output_file, dbl_fmt, dh.side_flux( *(ele->side(i) ) ) );

        xfprintf( raw_output_file, "\n" );
        cit ++;
    }
    xfprintf( raw_output_file, "$EndFlowField\n\n" );
}


#include "quadrature/quadrature_lib.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/mapping_p1.hh"
#include "fields/field_python.hh"
#include "fields/field_values.hh"

typedef FieldPython<3, FieldValue<3>::Vector > ExactSolution;

/*
* Calculate approximation of L2 norm for:
 * 1) difference between regularized pressure and analytical solution (using FunctionPython)
 * 2) difference between RT velocities and analytical solution
 * 3) difference of divergence
 * */

struct DiffData {
    double pressure_error[2], velocity_error[2], div_error[2];
    vector<double> pressure_diff;
    vector<double> velocity_diff;
    vector<double> div_diff;

    double * solution;
    const MH_DofHandler * dh;
    MHFEValues fe_values;

    DarcyFlowMH_Steady *darcy;
    DarcyFlowMH_Steady::EqData *data_;
};

template <int dim>
void l2_diff_local(ElementFullIter &ele, FEValues<dim,3> &fe_values, ExactSolution &anal_sol,  DiffData &result) {

    fe_values.reinit(ele);
    result.fe_values.update(ele,
    		result.data_->anisotropy,
    		result.data_->cross_section,
    		result.data_->conductivity);
    double conductivity = result.data_->conductivity.value(ele->centre(), ele->element_accessor() );
    double cross = result.data_->cross_section.value(ele->centre(), ele->element_accessor() );


    // get coefficients on the current element
    vector<double> fluxes(dim+1);
    vector<double> pressure_traces(dim+1);

    for (unsigned int li = 0; li < ele->n_sides(); li++) {
        fluxes[li] = result.dh->side_flux( *(ele->side( li ) ) );
        pressure_traces[li] = result.dh->side_scalar( *(ele->side( li ) ) );
    }
    double pressure_mean = result.dh->element_scalar(ele);

    arma::vec analytical(5);
    arma::vec3 flux_in_q_point;
    arma::vec3 anal_flux;

    double velocity_diff=0, divergence_diff=0, pressure_diff=0, diff;

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


        analytical = anal_sol.value(q_point, ele->element_accessor() );
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

        diff = - (1.0 / conductivity) * diff / dim / ele->measure() / cross + pressure_mean ;
        diff = ( diff - analytical[0]);
        pressure_diff += diff * diff * fe_values.JxW(i_point);


        // velocity difference
        flux_in_q_point.zeros();
        for(unsigned int i_shape=0; i_shape < ele->n_sides(); i_shape++) {
            flux_in_q_point += fluxes[ i_shape ]
                              * result.fe_values.RT0_value( ele, q_point, i_shape )
                              / cross;
        }

        flux_in_q_point -= anal_flux;
        velocity_diff += dot(flux_in_q_point, flux_in_q_point) * fe_values.JxW(i_point);

        // divergence diff
        diff = 0;
        for(unsigned int i_shape=0; i_shape < ele->n_sides(); i_shape++) diff += fluxes[ i_shape ];
        diff = ( diff / ele->measure() / cross - analytical[4]);
        divergence_diff += diff * diff * fe_values.JxW(i_point);

    }


    result.velocity_diff[ele.index()] = velocity_diff;
    result.velocity_error[dim-1] += velocity_diff;

    result.pressure_diff[ele.index()] = pressure_diff;
    result.pressure_error[dim-1] += pressure_diff;

    result.div_diff[ele.index()] = divergence_diff;
    result.div_error[dim-1] += divergence_diff;

}






void DarcyFlowMHOutput::compute_l2_difference() {
    DBGMSG("l2 norm output\n");
    ofstream os( FilePath("solution_error", FilePath::output_file) );

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

    FilePath source_file( "analytical_module.py", FilePath::input_file);
    ExactSolution  anal_sol_1d(5);   // components: pressure, flux vector 3d, divergence
    anal_sol_1d.set_python_field_from_file( source_file, "all_values_1d");

    ExactSolution anal_sol_2d(5);
    anal_sol_2d.set_python_field_from_file( source_file, "all_values_2d");


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


    auto vel_diff_ptr =	std::make_shared< FieldElementwise<3, FieldValue<3>::Scalar> >(&(result.velocity_diff[0]), 1, mesh_->n_elements());
    output_fields.velocity_diff.set_field(mesh_->region_db().get_region_set("ALL"), vel_diff_ptr, 0);
    auto pressure_diff_ptr =	std::make_shared< FieldElementwise<3, FieldValue<3>::Scalar> >(&(result.pressure_diff[0]), 1, mesh_->n_elements());
    output_fields.pressure_diff.set_field(mesh_->region_db().get_region_set("ALL"), pressure_diff_ptr, 0);
    auto div_diff_ptr =	std::make_shared< FieldElementwise<3, FieldValue<3>::Scalar> >(&(result.div_diff[0]), 1, mesh_->n_elements());
    output_fields.div_diff.set_field(mesh_->region_db().get_region_set("ALL"), div_diff_ptr, 0);

    output_fields.fields_for_output += output_fields;

    unsigned int solution_size;
    darcy_flow->get_solution_vector(result.solution, solution_size);

    result.dh = &( darcy_flow->get_mh_dofhandler());
    result.darcy = darcy_flow;
    result.data_ = &(darcy_flow->data_);

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

    os 	<< "l2 norm output\n\n"
    	<< "pressure error 1d: " << sqrt(result.pressure_error[0]) << endl
    	<< "pressure error 2d: " << sqrt(result.pressure_error[1]) << endl
    	<< "velocity error 1d: " << sqrt(result.velocity_error[0]) << endl
    	<< "velocity error 2d: " << sqrt(result.velocity_error[1]) << endl
    	<< "div error 1d: " << sqrt(result.div_error[0]) << endl
    	<< "div error 2d: " << sqrt(result.div_error[1]);
}

