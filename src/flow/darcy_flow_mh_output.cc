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
#include "flow/darcy_flow_lmh.hh"
#include "flow/assembly_mh.hh"
#include "flow/assembly_lmh.hh"
#include "flow/darcy_flow_mh_output.hh"

#include "io/output_time.hh"
#include "io/observe.hh"
#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/index_types.hh"

#include "fields/field_set.hh"
#include "fem/dofhandler.hh"
#include "fem/fe_values.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_values_views.hh"
#include "quadrature/quadrature_lib.hh"
#include "fields/field_fe.hh"
#include "fields/fe_value_handler.hh"
#include "fields/generic_field.hh"

#include "mesh/mesh.h"
#include "mesh/partitioning.hh"
#include "mesh/accessors.hh"
#include "mesh/node_accessor.hh"
#include "mesh/range_wrapper.hh"

// #include "coupling/balance.hh"
#include "intersection/mixed_mesh_intersections.hh"
#include "intersection/intersection_local.hh"

namespace it = Input::Type;


const it::Instance & DarcyFlowMHOutput::get_input_type(FieldSet& eq_data, const std::string &equation_name) {
	OutputFields output_fields;
	output_fields += eq_data;
	return output_fields.make_output_type(equation_name, "");
}


const it::Instance & DarcyFlowMHOutput::get_input_type_specific() {
    
    static it::Record& rec = it::Record("Output_DarcyMHSpecific", "Specific Darcy flow MH output.")
        .copy_keys(OutputSpecificFields::get_input_type())
        .declare_key("compute_errors", it::Bool(), it::Default("false"),
                        "SPECIAL PURPOSE. Computes error norms of the solution, particulary suited for non-compatible coupling models.")
        .declare_key("raw_flow_output", it::FileName::output(), it::Default::optional(),
                        "Output file with raw data from MH module.")
        .close();
    
    OutputSpecificFields output_fields;
    return output_fields.make_output_type_from_record(rec,
                                                        "Flow_Darcy_MH_specific",
                                                        "");
}


DarcyFlowMHOutput::OutputFields::OutputFields()
: EquationOutput()
{

//     *this += field_ele_pressure.name("pressure_p0_old").units(UnitSI().m()) // TODO remove: obsolete field
//              .flags(FieldFlag::equation_result)
//              .description("Pressure solution - P0 interpolation.");
//     *this += field_node_pressure.name("pressure_p1").units(UnitSI().m())
//              .flags(FieldFlag::equation_result)
//              .description("Pressure solution - P1 interpolation.");
// 	*this += field_ele_piezo_head.name("piezo_head_p0_old").units(UnitSI().m()) // TODO remove: obsolete field
//              .flags(FieldFlag::equation_result)
//              .description("Piezo head solution - P0 interpolation.");
// 	*this += field_ele_flux.name("velocity_p0_old").units(UnitSI().m()) // TODO remove: obsolete field
//              .flags(FieldFlag::equation_result)
//              .description("Velocity solution - P0 interpolation.");
	*this += subdomain.name("subdomain")
					  .units( UnitSI::dimensionless() )
					  .flags(FieldFlag::equation_external_output)
                      .description("Subdomain ids of the domain decomposition.");
	*this += region_id.name("region_id")
	        .units( UnitSI::dimensionless())
	        .flags(FieldFlag::equation_external_output)
            .description("Region ids.");
}


DarcyFlowMHOutput::OutputSpecificFields::OutputSpecificFields()
: EquationOutput()
{
    *this += pressure_diff.name("pressure_diff").units(UnitSI().m())
             .flags(FieldFlag::equation_result) 
             .description("Error norm of the pressure solution. [Experimental]");
    *this += velocity_diff.name("velocity_diff").units(UnitSI().m().s(-1))
             .flags(FieldFlag::equation_result)
             .description("Error norm of the velocity solution. [Experimental]");
    *this += div_diff.name("div_diff").units(UnitSI().s(-1))
             .flags(FieldFlag::equation_result)
             .description("Error norm of the divergence of the velocity solution. [Experimental]");
}

DarcyFlowMHOutput::DarcyFlowMHOutput(DarcyFlowInterface *flow, Input::Record main_mh_in_rec)
: darcy_flow(flow),
  mesh_(&darcy_flow->mesh()),
  compute_errors_(false),
  is_output_specific_fields(false)
{
    Input::Record in_rec_output = main_mh_in_rec.val<Input::Record>("output");
    
    output_stream = OutputTime::create_output_stream("flow",
                                                     main_mh_in_rec.val<Input::Record>("output_stream"),
                                                     darcy_flow->time().get_unit_string());
    prepare_output(in_rec_output);

    auto in_rec_specific = main_mh_in_rec.find<Input::Record>("output_specific");
    if (in_rec_specific) {
        in_rec_specific->opt_val("compute_errors", compute_errors_);
        
        // raw output
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            
            // optionally open raw output file
            FilePath raw_output_file_path;
            if (in_rec_specific->opt_val("raw_flow_output", raw_output_file_path))
            {
                int mpi_size;
                MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
                if(mpi_size > 1)
                {
                    WarningOut() << "Raw output is not available in parallel computation. MPI size: " << mpi_size << "\n";
                }
                else
                {
                    MessageOut() << "Opening raw flow output: " << raw_output_file_path << "\n";
                    try {
                        raw_output_file_path.open_stream(raw_output_file);
                    } INPUT_CATCH(FilePath::ExcFileOpen, FilePath::EI_Address_String, (*in_rec_specific))
                }
            }
        }
        
        auto fields_array = in_rec_specific->val<Input::Array>("fields");
        if(fields_array.size() > 0){
            is_output_specific_fields = true;
            prepare_specific_output(*in_rec_specific);
        }
    }
}

void DarcyFlowMHOutput::prepare_output(Input::Record in_rec)
{
  	// we need to add data from the flow equation at this point, not in constructor of OutputFields
	output_fields += darcy_flow->data();
	output_fields.set_mesh(*mesh_);

	output_fields.subdomain = GenericField<3>::subdomain(*mesh_);
	output_fields.region_id = GenericField<3>::region_id(*mesh_);

	//output_stream->add_admissible_field_names(in_rec_output.val<Input::Array>("fields"));
	//output_stream->mark_output_times(darcy_flow->time());
    output_fields.initialize(output_stream, mesh_, in_rec, darcy_flow->time() );
}

void DarcyFlowMHOutput::prepare_specific_output(Input::Record in_rec)
{
    diff_data.data_ = nullptr;
    if(DarcyMH* d = dynamic_cast<DarcyMH*>(darcy_flow))
    {
        diff_data.data_ = d->data_.get();
    }
    else if(DarcyLMH* d = dynamic_cast<DarcyLMH*>(darcy_flow))
    {
        diff_data.data_ = d->data_.get();
    }
    ASSERT_PTR(diff_data.data_);

    { // init DOF handlers represents element DOFs
        uint p_elem_component = 1;
        diff_data.dh_ = std::make_shared<SubDOFHandlerMultiDim>(diff_data.data_->dh_, p_elem_component);
    }

    // mask 2d elements crossing 1d
    if (diff_data.data_->mortar_method_ != DarcyMH::NoMortar) {
        diff_data.velocity_mask.resize(mesh_->n_elements(),0);
        for(IntersectionLocal<1,2> & isec : mesh_->mixed_intersections().intersection_storage12_) {
            diff_data.velocity_mask[ isec.bulk_ele_idx() ]++;
        }
    }

    output_specific_fields.set_mesh(*mesh_);

    diff_data.vel_diff_ptr = std::make_shared< FieldFE<3, FieldValue<3>::Scalar> >();
    diff_data.vel_diff_ptr->set_fe_data(diff_data.dh_);
    output_specific_fields.velocity_diff.set_field(mesh_->region_db().get_region_set("ALL"), diff_data.vel_diff_ptr, 0);
    diff_data.pressure_diff_ptr = std::make_shared< FieldFE<3, FieldValue<3>::Scalar> >();
    diff_data.pressure_diff_ptr->set_fe_data(diff_data.dh_);
    output_specific_fields.pressure_diff.set_field(mesh_->region_db().get_region_set("ALL"), diff_data.pressure_diff_ptr, 0);
    diff_data.div_diff_ptr = std::make_shared< FieldFE<3, FieldValue<3>::Scalar> >();
    diff_data.div_diff_ptr->set_fe_data(diff_data.dh_);
    output_specific_fields.div_diff.set_field(mesh_->region_db().get_region_set("ALL"), diff_data.div_diff_ptr, 0);

    output_specific_fields.set_time(darcy_flow->time().step(), LimitSide::right);
    output_specific_fields.initialize(output_stream, mesh_, in_rec, darcy_flow->time() );
}

DarcyFlowMHOutput::~DarcyFlowMHOutput()
{}





//=============================================================================
// CONVERT SOLUTION, CALCULATE BALANCES, ETC...
//=============================================================================


void DarcyFlowMHOutput::output()
{
    START_TIMER("Darcy fields output");

    {
        START_TIMER("post-process output fields");

        {
                  output_internal_flow_data();
        }
    }

    {
        START_TIMER("evaluate output fields");
        output_fields.set_time(darcy_flow->time().step(), LimitSide::right);
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



/*
 * Output of internal flow data.
 */
void DarcyFlowMHOutput::output_internal_flow_data()
{
    START_TIMER("DarcyFlowMHOutput::output_internal_flow_data");

    if (! raw_output_file.is_open()) return;
    
    //char dbl_fmt[ 16 ]= "%.8g ";
    // header
    raw_output_file <<  "// fields:\n//ele_id    ele_presure    flux_in_barycenter[3]    n_sides   side_pressures[n]    side_fluxes[n]\n";
    raw_output_file <<  fmt::format("$FlowField\nT={}\n", darcy_flow->time().t());
    raw_output_file <<  fmt::format("{}\n" , mesh_->n_elements() );

    
    DarcyMH::EqData* data = nullptr;
    if(DarcyMH* d = dynamic_cast<DarcyMH*>(darcy_flow))
    {
        data = d->data_.get();
    }
    else if(DarcyLMH* d = dynamic_cast<DarcyLMH*>(darcy_flow))
    {
        data = d->data_.get();
    }
    ASSERT_PTR(data);
    
    arma::vec3 flux_in_center;
    
    int cit = 0;
    for ( DHCellAccessor dh_cell : data->dh_->own_range() ) {
        ElementAccessor<3> ele = dh_cell.elm();
        LocDofVec indices = dh_cell.get_loc_dof_indices();

        // pressure
        raw_output_file << fmt::format("{} {} ", dh_cell.elm().index(), data->full_solution[indices[ele->n_sides()]]);
        
        // velocity at element center
        flux_in_center = data->field_ele_velocity.value(ele.centre(), ele);
        for (unsigned int i = 0; i < 3; i++)
        	raw_output_file << flux_in_center[i] << " ";

        // number of sides
        raw_output_file << ele->n_sides() << " ";
        
        // pressure on edges
        unsigned int lid = ele->n_sides() + 1;
        for (unsigned int i = 0; i < ele->n_sides(); i++, lid++) {
            raw_output_file << data->full_solution[indices[lid]] << " ";
        }
        // fluxes on sides
        for (unsigned int i = 0; i < ele->n_sides(); i++) {
            raw_output_file << data->full_solution[indices[i]] << " ";
        }
        
        raw_output_file << endl;
        cit ++;
    }    
    
    raw_output_file << "$EndFlowField\n" << endl;
}


#include "quadrature/quadrature_lib.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fields/field_python.hh"
#include "fields/field_values.hh"

typedef FieldPython<3, FieldValue<3>::Vector > ExactSolution;

/*
* Calculate approximation of L2 norm for:
 * 1) difference between regularized pressure and analytical solution (using FunctionPython)
 * 2) difference between RT velocities and analytical solution
 * 3) difference of divergence
 * */

void DarcyFlowMHOutput::l2_diff_local(DHCellAccessor dh_cell,
                   FEValues<3> &fe_values, FEValues<3> &fv_rt,
                   ExactSolution &anal_sol,  DarcyFlowMHOutput::DiffData &result) {

    ASSERT_DBG( fe_values.dim() == fv_rt.dim());
    unsigned int dim = fe_values.dim();

    ElementAccessor<3> ele = dh_cell.elm();
    fv_rt.reinit(ele);
    fe_values.reinit(ele);
    
    double conductivity = result.data_->conductivity.value(ele.centre(), ele );
    double cross = result.data_->cross_section.value(ele.centre(), ele );


    // get coefficients on the current element
    vector<double> fluxes(dim+1);
//     vector<double> pressure_traces(dim+1);

    for (unsigned int li = 0; li < ele->n_sides(); li++) {
        fluxes[li] = diff_data.data_->full_solution[ dh_cell.get_loc_dof_indices()[li] ];
//         pressure_traces[li] = result.dh->side_scalar( *(ele->side( li ) ) );
    }
    const uint ndofs = dh_cell.n_dofs();
    // TODO: replace with DHCell getter when available for FESystem component
    double pressure_mean = diff_data.data_->full_solution[ dh_cell.get_loc_dof_indices()[ndofs/2] ];

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
                    * arma::dot( *ele.node(i_node), *ele.node(j_node));
        }

    for(unsigned int i_point=0; i_point < fe_values.n_points(); i_point++) {
        arma::vec3 q_point = fe_values.point(i_point);


        analytical = anal_sol.value(q_point, ele );
        for(unsigned int i=0; i< 3; i++) anal_flux[i] = analytical[i+1];

        // compute postprocesed pressure
        diff = 0;
        for(unsigned int i_shape=0; i_shape < ele->n_sides(); i_shape++) {
            unsigned int oposite_node = 0;
            switch (dim) {
                case 1: oposite_node =  RefElement<1>::oposite_node(i_shape); break;
                case 2: oposite_node =  RefElement<2>::oposite_node(i_shape); break;
                case 3: oposite_node =  RefElement<3>::oposite_node(i_shape); break;
                default: ASSERT(false)(dim).error("Unsupported FE dimension."); break;
            }

            diff += fluxes[ i_shape ] *
                               (  arma::dot( q_point, q_point )/ 2
                                - mean_x_squared / 2
                                - arma::dot( q_point, *ele.node(oposite_node) )
                                + arma::dot( ele.centre(), *ele.node(oposite_node) )
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

    // DHCell constructed with diff fields DH, get DOF indices of actual element
    DHCellAccessor sub_dh_cell = dh_cell.cell_with_other_dh(result.dh_.get());
    IntIdx idx = sub_dh_cell.get_loc_dof_indices()[0];

    auto velocity_data = result.vel_diff_ptr->get_data_vec();
    velocity_data[ idx ] = sqrt(velocity_diff);
    result.velocity_error[dim-1] += velocity_diff;
    if (dim == 2 && result.velocity_mask.size() != 0 ) {
    	result.mask_vel_error += (result.velocity_mask[ ele.idx() ])? 0 : velocity_diff;
    }

    auto pressure_data = result.pressure_diff_ptr->get_data_vec();
    pressure_data[ idx ] = sqrt(pressure_diff);
    result.pressure_error[dim-1] += pressure_diff;

    auto div_data = result.div_diff_ptr->get_data_vec();
    div_data[ idx ] = sqrt(divergence_diff);
    result.div_error[dim-1] += divergence_diff;

}

DarcyFlowMHOutput::FEData::FEData()
: order(4),
  quad(QGauss::make_array(order)),
  fe_p1(0), fe_p0(0),
  fe_rt( )
{
    UpdateFlags flags = update_values | update_JxW_values | update_quadrature_points;
    fe_values = mixed_fe_values(quad, fe_p0, flags);
    fv_rt = mixed_fe_values(quad, fe_rt, flags);
}


void DarcyFlowMHOutput::compute_l2_difference() {
	DebugOut() << "l2 norm output\n";
    ofstream os( FilePath("solution_error", FilePath::output_file) );

    FilePath source_file( "analytical_module.py", FilePath::input_file);
    ExactSolution  anal_sol_1d(5);   // components: pressure, flux vector 3d, divergence
    anal_sol_1d.set_python_field_from_file( source_file, "all_values_1d");

    ExactSolution anal_sol_2d(5);
    anal_sol_2d.set_python_field_from_file( source_file, "all_values_2d");

    ExactSolution anal_sol_3d(5);
    anal_sol_3d.set_python_field_from_file( source_file, "all_values_3d");

    diff_data.mask_vel_error=0;
    for(unsigned int j=0; j<3; j++){
        diff_data.pressure_error[j] = 0;
        diff_data.velocity_error[j] = 0;
        diff_data.div_error[j] = 0;
    }

    //diff_data.ele_flux = &( ele_flux );

    for (DHCellAccessor dh_cell : diff_data.data_->dh_->own_range()) {

    	switch (dh_cell.dim()) {
        case 1:
            l2_diff_local( dh_cell, fe_data.fe_values[1], fe_data.fv_rt[1], anal_sol_1d, diff_data);
            break;
        case 2:
            l2_diff_local( dh_cell, fe_data.fe_values[2], fe_data.fv_rt[2], anal_sol_2d, diff_data);
            break;
        case 3:
            l2_diff_local( dh_cell, fe_data.fe_values[3], fe_data.fv_rt[3], anal_sol_3d, diff_data);
            break;
        }
    }
    
    // square root for L2 norm
    for(unsigned int j=0; j<3; j++){
        diff_data.pressure_error[j] = sqrt(diff_data.pressure_error[j]);
        diff_data.velocity_error[j] = sqrt(diff_data.velocity_error[j]);
        diff_data.div_error[j] = sqrt(diff_data.div_error[j]);
    }
    diff_data.mask_vel_error = sqrt(diff_data.mask_vel_error);

    os 	<< "l2 norm output\n\n"
    	<< "pressure error 1d: " << diff_data.pressure_error[0] << endl
    	<< "pressure error 2d: " << diff_data.pressure_error[1] << endl
    	<< "pressure error 3d: " << diff_data.pressure_error[2] << endl
    	<< "velocity error 1d: " << diff_data.velocity_error[0] << endl
    	<< "velocity error 2d: " << diff_data.velocity_error[1] << endl
    	<< "velocity error 3d: " << diff_data.velocity_error[2] << endl
    	<< "masked velocity error 2d: " << diff_data.mask_vel_error <<endl
    	<< "div error 1d: " << diff_data.div_error[0] << endl
    	<< "div error 2d: " << diff_data.div_error[1] << endl
        << "div error 3d: " << diff_data.div_error[2];
}


