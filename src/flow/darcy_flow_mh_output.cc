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
 * @file    darcy_flow_mh_output.cc
 * @ingroup flow
 * @brief   Output class for darcy_flow_mh model.
 * @author  Jan Brezina
 */

#include <vector>
#include <iostream>
#include <sstream>
#include <string>

#include <system/global_defs.h>

#include "flow/darcy_flow_mh.hh"
#include "flow/darcy_flow_assembly.hh"
#include "flow/darcy_flow_mh_output.hh"

#include "io/output_time.hh"
#include "io/observe.hh"
#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "fields/field_set.hh"
#include "fem/dofhandler.hh"
#include "fem/fe_values.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_values_views.hh"
#include "quadrature/quadrature_lib.hh"
#include "fields/field_fe.hh"
#include "fields/generic_field.hh"

#include "mesh/long_idx.hh"
#include "mesh/mesh.h"
#include "mesh/partitioning.hh"
#include "mesh/accessors.hh"
#include "mesh/node_accessor.hh"
#include "mesh/range_wrapper.hh"

// #include "coupling/balance.hh"
#include "intersection/mixed_mesh_intersections.hh"
#include "intersection/intersection_local.hh"

namespace it = Input::Type;


const it::Instance & DarcyFlowMHOutput::get_input_type() {
	OutputFields output_fields;
	DarcyMH::EqData eq_data;
	output_fields += eq_data;
	return output_fields.make_output_type("Flow_Darcy_MH", "");
}

const it::Instance & DarcyFlowMHOutput::get_input_type_specific_fields() {
	OutputSpecificFields output_fields;
	return output_fields.make_output_type("Flow_Darcy_MH_specific", "");
}

const it::Record & DarcyFlowMHOutput::get_input_type_specific() {
    return it::Record("Output_DarcyMHSpecific", "Specific Darcy flow MH output.")
        .declare_key("compute_errors", it::Bool(), it::Default("false"),
                        "SPECIAL PURPOSE. Computing errors pro non-compatible coupling.")
        .declare_key("python_solution", it::FileName::input(), it::Default::optional(),
                        "SPECIAL PURPOSE. Python script for computing analytical solution.")
        .declare_key("raw_flow_output", it::FileName::output(), it::Default::optional(),
                        "Output file with raw data from MH module.")
        .declare_key("output", DarcyFlowMHOutput::get_input_type_specific_fields(), it::Default::optional(),
                "Specific output fields.")
        .close();
}


DarcyFlowMHOutput::OutputFields::OutputFields()
: EquationOutput()
{

    *this += field_ele_pressure.name("pressure_p0").units(UnitSI().m());
    *this += field_node_pressure.name("pressure_p1").units(UnitSI().m());
	*this += field_ele_piezo_head.name("piezo_head_p0").units(UnitSI().m());
	*this += field_ele_flux.name("velocity_p0").units(UnitSI().m().s(-1));
	*this += subdomain.name("subdomain")
					  .units( UnitSI::dimensionless() )
					  .flags(FieldFlag::equation_external_output);
	*this += region_id.name("region_id")
	        .units( UnitSI::dimensionless())
	        .flags(FieldFlag::equation_external_output);
}


DarcyFlowMHOutput::OutputSpecificFields::OutputSpecificFields()
: EquationOutput()
{
    *this += field_ele_flux_enr.name("velocity_enr").units(UnitSI().m().s(-1));
    *this += field_ele_flux_reg.name("velocity_reg").units(UnitSI().m().s(-1));
    *this += pressure_diff.name("pressure_diff").units(UnitSI().m());
    *this += velocity_diff.name("velocity_diff").units(UnitSI().m().s(-1));
    *this += div_diff.name("div_diff").units(UnitSI().s(-1));
    *this += velocity_exact.name("velocity_exact").units(UnitSI().m().s(-1));
}

DarcyFlowMHOutput::DarcyFlowMHOutput(DarcyMH *flow)
: darcy_flow(flow),
  mesh_(&darcy_flow->mesh()),
  compute_errors_(false),
  is_output_specific_fields(false),
  fe1(1),
  fe2(1),
  fe3(1)
{}

void DarcyFlowMHOutput::initialize(Input::Record main_mh_in_rec)
{
    Input::Record in_rec_output = main_mh_in_rec.val<Input::Record>("output");
    
    output_stream = OutputTime::create_output_stream("flow",
                                                     main_mh_in_rec.val<Input::Record>("output_stream"),
                                                     darcy_flow->time().get_unit_string());
    prepare_output(in_rec_output);

    auto in_rec_specific = main_mh_in_rec.find<Input::Record>("output_specific");
    if (in_rec_specific) {
        in_rec_specific->opt_val("compute_errors", compute_errors_);
        if(compute_errors_){
            in_rec_specific->opt_val("python_solution", python_solution_filename_);
            ASSERT(python_solution_filename_.exists());
        }
        
        // raw output
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            // optionally open raw output file
            FilePath raw_output_file_path;
            if (in_rec_specific->opt_val("raw_flow_output", raw_output_file_path)) {
            	MessageOut() << "Opening raw flow output: " << raw_output_file_path << "\n";
            	try {
            		raw_output_file_path.open_stream(raw_output_file);
            	} INPUT_CATCH(FilePath::ExcFileOpen, FilePath::EI_Address_String, (*in_rec_specific))
            }
        }

        // possibly read the names of specific output fields
        auto in_rec_specific_output = in_rec_specific->find<Input::Record>("output");
        if(in_rec_specific_output){
            is_output_specific_fields = true;
            prepare_specific_output(*in_rec_specific_output);
        }
        else
            is_output_specific_fields = false;
    }
}

void DarcyFlowMHOutput::prepare_output(Input::Record in_rec)
{
  	// we need to add data from the flow equation at this point, not in constructor of OutputFields
	output_fields += darcy_flow->data();
	output_fields.set_mesh(*mesh_);
        
        all_element_idx_.resize(mesh_->n_elements());
	for(unsigned int i=0; i<all_element_idx_.size(); i++) all_element_idx_[i] = i;

	// create shared pointer to a FieldElementwise and push this Field to output_field on all regions
	ele_pressure.resize(mesh_->n_elements());
	auto ele_pressure_ptr=ele_pressure.create_field<3, FieldValue<3>::Scalar>(1);
	output_fields.field_ele_pressure.set_field(mesh_->region_db().get_region_set("ALL"), ele_pressure_ptr);

	dh_ = make_shared<DOFHandlerMultiDim>(*mesh_);
	dh_->distribute_dofs(fe1, fe2, fe3);
	corner_pressure.resize(dh_->n_global_dofs());

	auto corner_ptr = make_shared< FieldFE<3, FieldValue<3>::Scalar> >();
	corner_ptr->set_fe_data(dh_, &map1, &map2, &map3, &corner_pressure);

	output_fields.field_node_pressure.set_field(mesh_->region_db().get_region_set("ALL"), corner_ptr);
	output_fields.field_node_pressure.output_type(OutputTime::NODE_DATA);

	ele_piezo_head.resize(mesh_->n_elements());
	auto ele_piezo_head_ptr=ele_piezo_head.create_field<3, FieldValue<3>::Scalar>(1);
	output_fields.field_ele_piezo_head.set_field(mesh_->region_db().get_region_set("ALL"), ele_piezo_head_ptr);

	ele_flux.resize(3*mesh_->n_elements());
	auto ele_flux_ptr=ele_flux.create_field<3, FieldValue<3>::VectorFixed>(3);
	output_fields.field_ele_flux.set_field(mesh_->region_db().get_region_set("ALL"), ele_flux_ptr);

	output_fields.subdomain = GenericField<3>::subdomain(*mesh_);
	output_fields.region_id = GenericField<3>::region_id(*mesh_);

	//output_stream->add_admissible_field_names(in_rec_output.val<Input::Array>("fields"));
	//output_stream->mark_output_times(darcy_flow->time());
    output_fields.initialize(output_stream, mesh_, in_rec, darcy_flow->time() );
}

void DarcyFlowMHOutput::prepare_specific_output(Input::Record in_rec)
{
    diff_data.darcy = darcy_flow;
    diff_data.data_ = darcy_flow->data_.get();

    // mask 2d elements crossing 1d
    if (diff_data.data_->mortar_method_ != DarcyMH::NoMortar) {
        diff_data.velocity_mask.resize(mesh_->n_elements(),0);
        for(IntersectionLocal<1,2> & isec : mesh_->mixed_intersections().intersection_storage12_) {
            diff_data.velocity_mask[ isec.bulk_ele_idx() ]++;
        }
    }

    diff_data.pressure_diff.resize( mesh_->n_elements() );
    diff_data.velocity_diff.resize( mesh_->n_elements() );
    diff_data.div_diff.resize( mesh_->n_elements() );

    output_specific_fields.set_mesh(*mesh_);

    auto vel_diff_ptr =	diff_data.velocity_diff.create_field<3, FieldValue<3>::Scalar>(1);
    output_specific_fields.velocity_diff.set_field(mesh_->region_db().get_region_set("ALL"), vel_diff_ptr, 0);
    auto pressure_diff_ptr = diff_data.pressure_diff.create_field<3, FieldValue<3>::Scalar>(1);
    output_specific_fields.pressure_diff.set_field(mesh_->region_db().get_region_set("ALL"), pressure_diff_ptr, 0);
    auto div_diff_ptr =	diff_data.div_diff.create_field<3, FieldValue<3>::Scalar>(1);
    output_specific_fields.div_diff.set_field(mesh_->region_db().get_region_set("ALL"), div_diff_ptr, 0);
    
    /// Create empty fields
    std::shared_ptr<FieldElementwise<3, FieldValue<3>::VectorFixed> > field_ptr(new FieldElementwise<3, FieldValue<3>::VectorFixed>(3));
    output_specific_fields.field_ele_flux_enr.set_field(mesh_->region_db().get_region_set("ALL"), field_ptr);
    output_specific_fields.field_ele_flux_reg.set_field(mesh_->region_db().get_region_set("ALL"), field_ptr);
    output_specific_fields.velocity_exact.set_field(mesh_->region_db().get_region_set("ALL"), field_ptr);
    
    output_specific_fields.set_time(darcy_flow->time().step(), LimitSide::right);
    output_specific_fields.initialize(output_stream, mesh_, in_rec, darcy_flow->time() );
}

DarcyFlowMHOutput::~DarcyFlowMHOutput()
{};





//=============================================================================
// CONVERT SOLUTION, CALCULATE BALANCES, ETC...
//=============================================================================


void DarcyFlowMHOutput::output()
{
    START_TIMER("Darcy fields output");

    ElementSetRef observed_elements = output_stream->observe(mesh_)->observed_elements();
    {
        START_TIMER("post-process output fields");

        output_fields.set_time(darcy_flow->time().step(), LimitSide::right);

        if (output_fields.is_field_output_time(output_fields.field_ele_pressure,darcy_flow->time().step()) ||
            output_fields.is_field_output_time(output_fields.field_ele_piezo_head,darcy_flow->time().step()) )
                make_element_scalar(all_element_idx_);
        else
                make_element_scalar(observed_elements);

        if ( output_fields.is_field_output_time(output_fields.field_ele_flux,darcy_flow->time().step()) )
                make_element_vector(all_element_idx_);
        else
                make_element_vector(observed_elements);

        if ( output_fields.is_field_output_time(output_fields.field_node_pressure,darcy_flow->time().step()) )
                make_node_scalar_param(all_element_idx_);
        //else
        //        make_node_scalar_param(observed_elements);

        // Internal output only if both ele_pressure and ele_flux are output.
        if (output_fields.is_field_output_time(output_fields.field_ele_flux,darcy_flow->time().step()) &&
            output_fields.is_field_output_time(output_fields.field_ele_pressure,darcy_flow->time().step()) )
        {
                  output_internal_flow_data();
        }
    }

    {
        START_TIMER("evaluate output fields");
        output_fields.output(darcy_flow->time().step());
    }
    
    if (compute_errors_)
    {
        START_TIMER("compute specific output fields");
        compute_l2_difference();
    }
    
    if(is_output_specific_fields)
    {
        START_TIMER("evaluate output fields");
        output_specific_fields.set_time(darcy_flow->time().step(), LimitSide::right);
        output_specific_fields.output(darcy_flow->time().step());
    }

    {
        START_TIMER("write time frame");
        output_stream->write_time_frame();
    }

    
}


//=============================================================================
// FILL TH "SCALAR" FIELD FOR ALL ELEMENTS IN THE MESH
//=============================================================================

void DarcyFlowMHOutput::make_element_scalar(ElementSetRef element_indices)
{
    START_TIMER("DarcyFlowMHOutput::make_element_scalar");
    // need to call this to create mh solution vector
    darcy_flow->get_mh_dofhandler();
    
    for(unsigned int i_ele : element_indices) {
        ElementAccessor<3> ele = mesh_->element_accessor(i_ele);
        ele_pressure[i_ele] = darcy_flow->mh_dh.element_scalar(ele);
        ele_piezo_head[i_ele] = ele_pressure[i_ele] +
          - (darcy_flow->data_->gravity_[3] + arma::dot(darcy_flow->data_->gravity_vec_,ele.centre()));
    }
}



/****
 * compute Darcian velocity in centre of elements
 *
 */
void DarcyFlowMHOutput::make_element_vector(ElementSetRef element_indices) {
    START_TIMER("DarcyFlowMHOutput::make_element_vector");
    // need to call this to create mh solution vector
    darcy_flow->get_mh_dofhandler();

    auto multidim_assembler = AssemblyBase::create< AssemblyMH >(darcy_flow->data_);
    arma::vec3 flux_in_center;
    for(unsigned int i_ele : element_indices) {
    	ElementAccessor<3> ele = mesh_->element_accessor(i_ele);

        flux_in_center = multidim_assembler[ele->dim() -1]->make_element_vector(ele);

        // place it in the sequential vector
        for(unsigned int j=0; j<3; j++) ele_flux[3*i_ele + j]=flux_in_center[j];
    }
}


void DarcyFlowMHOutput::make_corner_scalar(vector<double> &node_scalar)
{
    START_TIMER("DarcyFlowMHOutput::make_corner_scalar");
	unsigned int ndofs = dh_->max_elem_dofs();
	std::vector<LongIdx> indices(ndofs);
	unsigned int i_node;
	for (auto ele : mesh_->bulk_elements_range()) {
		dh_->get_dof_indices(ele, indices);
		for (i_node=0; i_node<ele->n_nodes(); i_node++)
		{
			corner_pressure[indices[i_node]] = node_scalar[ ele.node_accessor(i_node).idx() ];
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

void DarcyFlowMHOutput::make_node_scalar_param(ElementSetRef element_indices) {
    START_TIMER("DarcyFlowMHOutput::make_node_scalar_param");

	vector<double> scalars(mesh_->n_nodes());

    double dist; //!< tmp variable for storing particular distance node --> element, node --> side*/

    /** Accessors */
    NodeAccessor<3> node;

    int n_nodes = mesh_->n_nodes(); //!< number of nodes in the mesh */
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
    	for (auto ele : mesh_->bulk_elements_range())
            for (unsigned int li = 0; li < ele->n_nodes(); li++) {
                node = ele.node_accessor(li); //!< get NodeAccessor from element */
                node_index = node.idx(); //!< get nod index from mesh */

                dist = sqrt(
                        ((node->getX() - ele.centre()[ 0 ])*(node->getX() - ele.centre()[ 0 ])) +
                        ((node->getY() - ele.centre()[ 1 ])*(node->getY() - ele.centre()[ 1 ])) +
                        ((node->getZ() - ele.centre()[ 2 ])*(node->getZ() - ele.centre()[ 2 ]))
                );
                sum_ele_dist[node_index] += dist;
                sum_elements[node_index]++;
            }
    }
    if (count_sides){
    	for (auto ele : mesh_->bulk_elements_range())
            for(SideIter side = ele.side(0); side->side_idx() < ele->n_sides(); ++side) {
                for (unsigned int li = 0; li < side->n_nodes(); li++) {
                    node = side->node(li);//!< get NodeAccessor from element */
                    node_index = node.idx(); //!< get nod index from mesh */
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
    	for (auto ele : mesh_->bulk_elements_range())
            for (unsigned int li = 0; li < ele->n_nodes(); li++) {
                node = ele.node_accessor(li);//!< get NodeAccessor from element */
                node_index = ele.node_accessor(li).idx(); //!< get nod index from mesh */

                /**TODO - calculate it again or store it in prior pass*/
                dist = sqrt(
                        ((node->getX() - ele.centre()[ 0 ])*(node->getX() - ele.centre()[ 0 ])) +
                        ((node->getY() - ele.centre()[ 1 ])*(node->getY() - ele.centre()[ 1 ])) +
                        ((node->getZ() - ele.centre()[ 2 ])*(node->getZ() - ele.centre()[ 2 ]))
                );
                scalars[node_index] += ele_pressure[ ele.idx() ] *
                        (1 - dist / (sum_ele_dist[node_index] + sum_side_dist[node_index])) /
                        (sum_elements[node_index] + sum_sides[node_index] - 1);
            }
    }
    if (count_sides) {
    	for (auto ele : mesh_->bulk_elements_range())
            for(SideIter side = ele.side(0); side->side_idx() < ele->n_sides(); ++side) {
                for (unsigned int li = 0; li < side->n_nodes(); li++) {
                    node = side->node(li);//!< get NodeAccessor from element */
                    node_index = node.idx(); //!< get nod index from mesh */

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


/*
 * Output of internal flow data.
 */
void DarcyFlowMHOutput::output_internal_flow_data()
{
    START_TIMER("DarcyFlowMHOutput::output_internal_flow_data");
    const MH_DofHandler &dh = darcy_flow->get_mh_dofhandler();

    if (! raw_output_file.is_open()) return;
    
    //char dbl_fmt[ 16 ]= "%.8g ";
    // header
    raw_output_file <<  "// fields:\n//ele_id    ele_presure    flux_in_barycenter[3]    n_sides   side_pressures[n]    side_fluxes[n]\n";
    raw_output_file <<  fmt::format("$FlowField\nT={}\n", darcy_flow->time().t());
    raw_output_file <<  fmt::format("{}\n" , mesh_->n_elements() );

    int cit = 0;
    for (auto ele : mesh_->bulk_elements_range()) {
        raw_output_file << fmt::format("{} {} ", ele.index(), ele_pressure[cit]);
        for (unsigned int i = 0; i < 3; i++)
            raw_output_file << ele_flux[3*cit+i] << " ";

        raw_output_file << ele->n_sides() << " ";

        for (unsigned int i = 0; i < ele->n_sides(); i++) {
            raw_output_file << dh.side_scalar( *(ele.side(i) ) ) << " ";
        }
        for (unsigned int i = 0; i < ele->n_sides(); i++) {
            raw_output_file << dh.side_flux( *(ele.side(i) ) ) << " ";
        }
        
        raw_output_file << endl;
        cit ++;
    }
    raw_output_file << "$EndFlowField\n" << endl;
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

template <int dim>
void DarcyFlowMHOutput::l2_diff_local(ElementAccessor<3> &ele,
                   FEValues<dim,3> &fe_values, FEValues<dim,3> &fv_rt,
                   ExactSolution &anal_sol,  DarcyFlowMHOutput::DiffData &result) {

    fv_rt.reinit(ele);
    fe_values.reinit(ele);
    
    double conductivity = result.data_->conductivity.value(ele.centre(), ele );
    double cross = result.data_->cross_section.value(ele.centre(), ele );


    // get coefficients on the current element
    vector<double> fluxes(dim+1);
//     vector<double> pressure_traces(dim+1);

    for (unsigned int li = 0; li < ele->n_sides(); li++) {
        fluxes[li] = result.dh->side_flux( *(ele.side( li ) ) );
//         pressure_traces[li] = result.dh->side_scalar( *(ele->side( li ) ) );
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
                    * arma::dot( ele.node(i_node)->point(), ele.node(j_node)->point());
        }

    for(unsigned int i_point=0; i_point < fe_values.n_points(); i_point++) {
        arma::vec3 q_point = fe_values.point(i_point);


        analytical = anal_sol.value(q_point, ele );
        for(unsigned int i=0; i< 3; i++) anal_flux[i] = analytical[i+1];

        // compute postprocesed pressure
        diff = 0;
        for(unsigned int i_shape=0; i_shape < ele->n_sides(); i_shape++) {
            unsigned int oposite_node = RefElement<dim>::oposite_node(i_shape);

            diff += fluxes[ i_shape ] *
                               (  arma::dot( q_point, q_point )/ 2
                                - mean_x_squared / 2
                                - arma::dot( q_point, ele.node(oposite_node)->point() )
                                + arma::dot( ele.centre(), ele.node(oposite_node)->point() )
                               );
        }

        diff = - (1.0 / conductivity) * diff / dim / ele.measure() / cross + pressure_mean ;
        diff = ( diff - analytical[0]);
        pressure_diff += diff * diff * fe_values.JxW(i_point);


        // velocity difference
        flux_in_q_point.zeros();
        for(unsigned int i_shape=0; i_shape < ele->n_sides(); i_shape++) {
            flux_in_q_point += fluxes[ i_shape ]
                              * fv_rt.vector_view(0).value(i_shape, i_point)
                              / cross;
        }

        flux_in_q_point -= anal_flux;
        velocity_diff += dot(flux_in_q_point, flux_in_q_point) * fe_values.JxW(i_point);

        // divergence diff
        diff = 0;
        for(unsigned int i_shape=0; i_shape < ele->n_sides(); i_shape++) diff += fluxes[ i_shape ];
        diff = ( diff / ele.measure() / cross - analytical[4]);
        divergence_diff += diff * diff * fe_values.JxW(i_point);

    }


    result.velocity_diff[ ele.idx() ] = velocity_diff;
    result.velocity_error[dim-1] += velocity_diff;
    if (dim == 2 && result.velocity_mask.size() != 0 ) {
    	result.mask_vel_error += (result.velocity_mask[ ele.idx() ])? 0 : velocity_diff;
    }

    result.pressure_diff[ ele.idx() ] = pressure_diff;
    result.pressure_error[dim-1] += pressure_diff;

    result.div_diff[ ele.idx() ] = divergence_diff;
    result.div_error[dim-1] += divergence_diff;

}

template
void DarcyFlowMHOutput::l2_diff_local<3>(ElementAccessor<3> &ele,
                   FEValues<3,3> &fe_values, FEValues<3,3> &fv_rt,
                   ExactSolution &anal_sol,  DarcyFlowMHOutput::DiffData &result);


void DarcyFlowMHOutput::compute_l2_difference() {
	DebugOut() << "l2 norm output\n";
    ofstream os( FilePath("solution_error", FilePath::output_file) );

    const unsigned int order = 4; // order of Gauss quadrature

    // we create trivial Dofhandler , for P0 elements, to get access to, FEValues on individual elements
    // this we use to integrate our own functions - difference of postprocessed pressure and analytical solution
    FE_P_disc<1> fe_1d(0);
    FE_P_disc<2> fe_2d(0);

    QGauss<1> quad_1d( order );
    QGauss<2> quad_2d( order );

    MappingP1<1,3> mapp_1d;
    MappingP1<2,3> mapp_2d;

    FEValues<1,3> fe_values_1d(mapp_1d, quad_1d,   fe_1d, update_JxW_values | update_quadrature_points);
    FEValues<2,3> fe_values_2d(mapp_2d, quad_2d,   fe_2d, update_JxW_values | update_quadrature_points);
    
    // FEValues for velocity.
    FE_RT0<1> fe_rt1d;
    FE_RT0<2> fe_rt2d;
    FEValues<1,3> fv_rt1d(mapp_1d,quad_1d, fe_rt1d, update_values | update_quadrature_points);
    FEValues<2,3> fv_rt2d(mapp_2d,quad_2d, fe_rt2d, update_values | update_quadrature_points);

    ASSERT(python_solution_filename_.exists());
    ExactSolution  anal_sol_1d(5);   // components: pressure, flux vector 3d, divergence
    anal_sol_1d.set_python_field_from_file( python_solution_filename_, "all_values_1d");

    ExactSolution anal_sol_2d(5);
    anal_sol_2d.set_python_field_from_file( python_solution_filename_, "all_values_2d");


    diff_data.dh = &( darcy_flow->get_mh_dofhandler());
    diff_data.pressure_error[0] = 0;
    diff_data.velocity_error[0] = 0;
    diff_data.div_error[0] = 0;
    diff_data.pressure_error[1] = 0;
    diff_data.velocity_error[1] = 0;
    diff_data.div_error[1] = 0;
    diff_data.mask_vel_error=0;

    //diff_data.ele_flux = &( ele_flux );
    
    unsigned int solution_size;
    darcy_flow->get_solution_vector(diff_data.solution, solution_size);


    for (auto ele : mesh_->bulk_elements_range()) {

    	switch (ele->dim()) {
        case 1:

            l2_diff_local<1>( ele, fe_values_1d, fv_rt1d, anal_sol_1d, diff_data);
            break;
        case 2:
            l2_diff_local<2>( ele, fe_values_2d, fv_rt2d, anal_sol_2d, diff_data);
            break;
        }
    }

    os 	<< "l2 norm output\n\n"
    	<< "pressure error 1d: " << sqrt(diff_data.pressure_error[0]) << endl
    	<< "pressure error 2d: " << sqrt(diff_data.pressure_error[1]) << endl
    	<< "velocity error 1d: " << sqrt(diff_data.velocity_error[0]) << endl
    	<< "velocity error 2d: " << sqrt(diff_data.velocity_error[1]) << endl
    	<< "masked velocity error 2d: " << sqrt(diff_data.mask_vel_error) <<endl
    	<< "div error 1d: " << sqrt(diff_data.div_error[0]) << endl
    	<< "div error 2d: " << sqrt(diff_data.div_error[1]);
}


