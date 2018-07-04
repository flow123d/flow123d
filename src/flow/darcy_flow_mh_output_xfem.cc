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
#include "darcy_flow_assembly_xfem.hh"
#include "flow/darcy_flow_mh_output_xfem.hh"
#include "flow/field_velocity.hh"

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
// #include "mesh/partitioning.hh"

// #include "coupling/balance.hh"
#include "intersection/mixed_mesh_intersections.hh"
#include "intersection/intersection_local.hh"

typedef FieldPython<3, FieldValue<3>::Vector > ExactSolution;
typedef FieldPython<3, FieldValue<3>::VectorFixed > ExactVelocity;


DarcyFlowMHOutputXFEM::DarcyFlowMHOutputXFEM(DarcyMH *flow)
: DarcyFlowMHOutput(flow),
fe_data_1d_xfem(true),
fe_data_2d_xfem(true),
fe_data_3d_xfem(true)
{}

void DarcyFlowMHOutputXFEM::prepare_output(Input::Record in_rec)
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
	dh_->distribute_dofs(fe_data_1d.fe_p1, fe_data_2d.fe_p1, fe_data_3d.fe_p1);
	corner_pressure.resize(dh_->n_global_dofs());

	auto corner_ptr = make_shared< FieldFE<3, FieldValue<3>::Scalar> >();
	corner_ptr->set_fe_data(dh_, &fe_data_1d.mapp, &fe_data_2d.mapp, &fe_data_3d.mapp, &corner_pressure);

	output_fields.field_node_pressure.set_field(mesh_->region_db().get_region_set("ALL"), corner_ptr);
	output_fields.field_node_pressure.output_type(OutputTime::NODE_DATA);

	ele_piezo_head.resize(mesh_->n_elements());
	auto ele_piezo_head_ptr=ele_piezo_head.create_field<3, FieldValue<3>::Scalar>(1);
	output_fields.field_ele_piezo_head.set_field(mesh_->region_db().get_region_set("ALL"), ele_piezo_head_ptr);

	ele_flux.resize(3*mesh_->n_elements());
        if(darcy_flow->mh_dh.enrich_velocity)
        {
            field_velocity = std::make_shared<FieldVelocity>(&darcy_flow->mh_dh, &darcy_flow->data_->cross_section, 
                                                            darcy_flow->mh_dh.enrich_velocity, true);
            output_fields.field_ele_flux.set_field(mesh_->region_db().get_region_set("ALL"), field_velocity);            
        }
        else{
            auto ele_flux_ptr=ele_flux.create_field<3, FieldValue<3>::VectorFixed>(3);
            output_fields.field_ele_flux.set_field(mesh_->region_db().get_region_set("ALL"), ele_flux_ptr);
        }

	output_fields.subdomain = GenericField<3>::subdomain(*mesh_);
	output_fields.region_id = GenericField<3>::region_id(*mesh_);

	//output_stream->add_admissible_field_names(in_rec_output.val<Input::Array>("fields"));
	//output_stream->mark_output_times(darcy_flow->time());
    output_fields.initialize(output_stream, mesh_, in_rec, darcy_flow->time() );
}

void DarcyFlowMHOutputXFEM::prepare_specific_output(Input::Record in_rec)
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

    field_velocity_enr_part = std::make_shared<FieldVelocity>(&darcy_flow->mh_dh, &darcy_flow->data_->cross_section,
                                                                darcy_flow->mh_dh.enrich_velocity, false);
    output_specific_fields.field_ele_flux_enr.set_field(mesh_->region_db().get_region_set("ALL"), field_velocity_enr_part);

    field_velocity_reg_part = std::make_shared<FieldVelocity>(&darcy_flow->mh_dh, &darcy_flow->data_->cross_section, false, true);
    output_specific_fields.field_ele_flux_reg.set_field(mesh_->region_db().get_region_set("ALL"), field_velocity_reg_part);

    std::shared_ptr<ExactVelocity> exact_vel_2d_ptr = std::make_shared<ExactVelocity>(3);
    if(python_solution_filename_.exists())
        exact_vel_2d_ptr->set_python_field_from_file( python_solution_filename_, "velocity_2d");
    output_specific_fields.velocity_exact.set_field(mesh_->region_db().get_region_set("ALL"), exact_vel_2d_ptr, 0);
    
    output_specific_fields.set_time(darcy_flow->time().step(), LimitSide::right);
    output_specific_fields.initialize(output_stream, mesh_, in_rec, darcy_flow->time() );
}

DarcyFlowMHOutputXFEM::~DarcyFlowMHOutputXFEM()
{};



/****
 * compute Darcian velocity in centre of elements
 *
 */
void DarcyFlowMHOutputXFEM::make_element_vector(ElementSetRef element_indices) {
    START_TIMER("DarcyFlowMHOutput::make_element_vector");
    // need to call this to create mh solution vector
    darcy_flow->get_mh_dofhandler();
    
    // create proper assembler
    AssemblyBase::MultidimAssembly multidim_assembler;
    if(darcy_flow->use_xfem && darcy_flow->mh_dh.enrich_velocity)
        multidim_assembler =  AssemblyBase::create< AssemblyMHXFEM >(darcy_flow->data_);
    else
        multidim_assembler =  AssemblyBase::create< AssemblyMH >(darcy_flow->data_);
    
    arma::vec3 flux_in_center;
    for(unsigned int i_ele : element_indices) {
        ElementAccessor<3> ele = mesh_->element_accessor(i_ele);

        unsigned int dim = ele.dim();
        flux_in_center = multidim_assembler[dim-1]->make_element_vector(ele);

        // place it in the sequential vector
        for(unsigned int j=0; j<3; j++) ele_flux[3*i_ele + j]=flux_in_center[j];
    }
}


#include "quadrature/quadrature_lib.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/mapping_p1.hh"
#include "fields/field_python.hh"
#include "fields/field_values.hh"


template <int dim>
void DarcyFlowMHOutputXFEM::l2_diff_local_xfem(LocalElementAccessorBase<3> &ele_ac,
                   XFEValues<dim,3> &fe_values, XFEValues<dim,3> &fv_rt,
                   ExactSolution &anal_sol,  DarcyFlowMHOutput::DiffData &result) {
//     DBGCOUT(<< "local diff\n");
    
    ElementAccessor<3> ele = ele_ac.element_accessor();
    
//     double conductivity = result.data_->conductivity.value(ele.centre(), ele );
    double cross = result.data_->cross_section.value(ele.centre(), ele );
    
    int dofs_vel[100];
    unsigned int ndofs_vel = ele_ac.get_dofs_vel(dofs_vel);
    
    // get values on the current element
    vector<double> dofs_vel_val(ndofs_vel);
    for (unsigned int i = 0; i < ndofs_vel; i++) {
        dofs_vel_val[i] = result.solution[dofs_vel[i]];
    }

    arma::vec analytical(5);
    arma::vec3 flux_in_q_point;
    arma::vec3 anal_flux;

    double velocity_diff=0, divergence_diff=0, pressure_diff=0, diff;

    double pressure_mean = result.dh->element_scalar(ele);
//     XFEMElementSingularData* xd = nullptr;
//     if(ele_ac.is_enriched())
//         xd = ele_ac.xfem_data_sing();
       
//     // 1d:  mean_x_squared = 1/6 (v0^2 + v1^2 + v0*v1)
//     // 2d:  mean_x_squared = 1/12 (v0^2 + v1^2 +v2^2 + v0*v1 + v0*v2 + v1*v2)
//     double mean_x_squared=0;
//     for(unsigned int i_node=0; i_node < ele->n_nodes(); i_node++ )
//         for(unsigned int j_node=0; j_node < ele->n_nodes(); j_node++ )
//         {
//             mean_x_squared += (i_node == j_node ? 2.0 : 1.0) / ( 6 * dim )   // multiply by 2 on diagonal
//                     * arma::dot( ele->node[i_node]->point(), ele->node[j_node]->point());
//         }

    auto velocity = fv_rt.vector_view(0);
    for(unsigned int i_point=0; i_point < fe_values.n_points(); i_point++) {
        arma::vec3 q_point = fe_values.point(i_point);

        analytical = anal_sol.value(q_point, ele );
        for(unsigned int i=0; i< 3; i++) anal_flux[i] = analytical[i+1];

        diff = pressure_mean;
//         if(xd){
//             arma::vec tmp_rho({3.33 + 0.03, 3.33, 0});
//             double aw = result.solution[ele_ac.sing_row(0)] / xd->enrichment_func(0)->value(tmp_rho);
//             double sr = xd->enrichment_func(0)->value(q_point);
// //             double srt = xd->enrichment_func(0)->value(ele->centre());
//             diff += aw*sr-;
//         }
        diff = diff - analytical[0];
        pressure_diff += diff * diff * fe_values.JxW(i_point);
//         // compute postprocesed pressure
//         diff = 0;
//         for(unsigned int i_shape=0; i_shape < ele->n_sides(); i_shape++) {
//             unsigned int oposite_node = RefElement<dim>::oposite_node(i_shape);
// 
//             diff += fluxes[ i_shape ] *
//                                (  arma::dot( q_point, q_point )/ 2
//                                 - mean_x_squared / 2
//                                 - arma::dot( q_point, ele->node[oposite_node]->point() )
//                                 + arma::dot( ele->centre(), ele->node[oposite_node]->point() )
//                                );
//         }
// 
//         diff = - (1.0 / conductivity) * diff / dim / ele->measure() / cross + pressure_mean ;
//         diff = ( diff - analytical[0]);
//         pressure_diff += diff * diff * fe_values.JxW(i_point);


        // velocity difference
        flux_in_q_point.zeros();
        for(unsigned int i_shape=0; i_shape < ndofs_vel; i_shape++) {
            flux_in_q_point += dofs_vel_val[ i_shape ] * velocity.value(i_shape,i_point);
        }
        flux_in_q_point = flux_in_q_point / cross;

        flux_in_q_point -= anal_flux;
        velocity_diff += dot(flux_in_q_point, flux_in_q_point) * fe_values.JxW(i_point);
//         velocity_diff += arma::norm(anal_flux,1) * fe_values.JxW(i_point);
        
        // divergence diff - flow over sides (RT dofs)
//         diff = 0;
//         for(unsigned int i_shape=0; i_shape < ele->n_sides(); i_shape++) diff += dofs_vel_val[ i_shape ];
//         diff = ( diff / ele->measure() / cross - analytical[4]);
//         divergence_diff += diff * diff * fe_values.JxW(i_point);

    }


//     DBGVAR(velocity_diff);
    result.velocity_diff[ele.idx()] = std::sqrt(velocity_diff);
    result.velocity_error[dim-1] += velocity_diff;
//     if (dim == 2 && result.velocity_mask.size() != 0 ) {
//         result.mask_vel_error += (result.velocity_mask[ ele.idx() ])? 0 : velocity_diff;
//     }

    result.pressure_diff[ele.idx()] = std::sqrt(pressure_diff);
    result.pressure_error[dim-1] += pressure_diff;

//     result.div_diff[ele.idx()] = std::sqrt(divergence_diff);
//     result.div_error[dim-1] += divergence_diff;

}

template<int dim>
DarcyFlowMHOutputXFEM::FEDataXFEM<dim>::FEDataXFEM(bool single_enr)
: DarcyFlowMHOutput::FEData<dim>(),
    single_enr(single_enr),
    qfactory(13-2*dim),
    fv_rt_xfem(this->mapp, this->fe_rt, this->fe_p0, update_values |
                  update_JxW_values | update_jacobians |
                  update_inverse_jacobians | update_quadrature_points
                  | update_divergence),
    fv_p0_xfem(this->mapp, this->fe_p0, this->fe_p0,
                  update_JxW_values | update_jacobians |
                  update_inverse_jacobians | update_quadrature_points)
{}
    
template<int dim>
void DarcyFlowMHOutputXFEM::FEDataXFEM<dim>::prepare_xfem(LocalElementAccessorBase<3> ele_ac)
{
    ElementAccessor<3> ele = ele_ac.element_accessor();
    XFEMElementSingularData<dim> * xdata = ele_ac.xfem_data_sing<dim>();
    
    qxfem = qfactory.create_singular(xdata->sing_vec(), ele);
    
    fv_rt_xfem.reinit(ele, *xdata, *qxfem);
    fv_p0_xfem.reinit(ele, *xdata, *qxfem);
}

    
void DarcyFlowMHOutputXFEM::compute_l2_difference() {
	DebugOut() << "l2 norm output\n";
    ofstream os( FilePath("solution_error", FilePath::output_file) );
    
    ASSERT(python_solution_filename_.exists());
    ExactSolution  anal_sol_1d(5);   // components: pressure, flux vector 3d, divergence
    anal_sol_1d.set_python_field_from_file( python_solution_filename_, "all_values_1d");

    ExactSolution anal_sol_2d(5);
    anal_sol_2d.set_python_field_from_file( python_solution_filename_, "all_values_2d");

    ExactSolution anal_sol_3d(5);
    anal_sol_3d.set_python_field_from_file( python_solution_filename_, "all_values_3d");

    diff_data.dh = &( darcy_flow->get_mh_dofhandler());
    fe_data_1d_xfem.single_enr = diff_data.dh->single_enr;
    fe_data_2d_xfem.single_enr = diff_data.dh->single_enr;
    fe_data_3d_xfem.single_enr = diff_data.dh->single_enr;
    diff_data.mask_vel_error=0;
    for(unsigned int j=0; j<3; j++){
        diff_data.pressure_error[j] = 0;
        diff_data.velocity_error[j] = 0;
        diff_data.div_error[j] = 0;
    }
    //diff_data.ele_flux = &( ele_flux );
    
    unsigned int solution_size;
    darcy_flow->get_solution_vector(diff_data.solution, solution_size);


    for (unsigned int i_loc = 0; i_loc < diff_data.dh->el_ds->lsize(); i_loc++) {
        DBGVAR(i_loc);
        auto ele_ac = const_cast<MH_DofHandler*>(diff_data.dh)->accessor(i_loc);
        ElementAccessor<3> ele = ele_ac.element_accessor();
        unsigned int dim = ele_ac.dim();
        
        switch (dim) {
//             case 1:
//                 l2_diff_local<1>( ele, fe_data_1d.fe_values, fe_data_1d.fv_rt, anal_sol_1d, diff_data);
//                 break;
        case 2:
            if(ele_ac.is_enriched()){
                fe_data_2d_xfem.prepare_xfem(ele_ac);
                l2_diff_local_xfem<2>(ele_ac, fe_data_2d_xfem.fv_p0_xfem, fe_data_2d_xfem.fv_rt_xfem,
                                      anal_sol_2d, diff_data);
            }
            else
                l2_diff_local<2>(ele, fe_data_2d.fe_values, fe_data_2d.fv_rt,
                                  anal_sol_2d, diff_data);
            break;
        case 3:
            if(ele_ac.is_enriched()){
                fe_data_3d_xfem.prepare_xfem(ele_ac);
                l2_diff_local_xfem<3>(ele_ac, fe_data_3d_xfem.fv_p0_xfem, fe_data_3d_xfem.fv_rt_xfem,
                                      anal_sol_3d, diff_data);
            }
            else
                l2_diff_local<3>(ele, fe_data_3d.fe_values, fe_data_3d.fv_rt,
                                  anal_sol_3d, diff_data);
            break;
        }
    }

    DebugOut() << "l2 norm output end\n";
    
    os 	<< "l2 norm output\n\n"
    	<< "pressure error 1d: " << sqrt(diff_data.pressure_error[0]) << endl
    	<< "pressure error 2d: " << sqrt(diff_data.pressure_error[1]) << endl
    	<< "pressure error 3d: " << sqrt(diff_data.pressure_error[2]) << endl
    	<< "velocity error 1d: " << sqrt(diff_data.velocity_error[0]) << endl
    	<< "velocity error 2d: " << sqrt(diff_data.velocity_error[1]) << endl
    	<< "velocity error 3d: " << sqrt(diff_data.velocity_error[2]) << endl
    	<< "masked vel error 2d: " << sqrt(diff_data.mask_vel_error) <<endl
    	<< "div error 1d: " << sqrt(diff_data.div_error[0]) << endl
    	<< "div error 2d: " << sqrt(diff_data.div_error[1]) << endl
    	<< "div error 3d: " << sqrt(diff_data.div_error[2]);
    
    if(darcy_flow->use_xfem){
        os << endl
        << "enr vel dof: " << diff_data.dh->mh_solution[mesh_->n_sides_] << endl
        << "sing LP: " << diff_data.dh->mh_solution[diff_data.dh->row_4_sing[0]];
    }
}


