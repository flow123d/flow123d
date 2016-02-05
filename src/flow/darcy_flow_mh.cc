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
 * @file    darcy_flow_mh.cc
 * @ingroup flow
 * @brief   Setup and solve linear system of mixed-hybrid discretization of the linear
 *          porous media flow with possible preferential flow in fractures and chanels.
 */

#include "petscmat.h"
#include "petscviewer.h"
#include "petscerror.h"
#include <armadillo>
#include <boost/foreach.hpp>


#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "input/factory.hh"

#include "mesh/mesh.h"
#include "mesh/intersection.hh"
#include "mesh/partitioning.hh"
#include "la/distribution.hh"
#include "la/linsys.hh"
#include "la/linsys_PETSC.hh"
#include "la/linsys_BDDC.hh"
#include "la/schur.hh"
#include "la/sparse_graph.hh"
#include "la/local_to_global_map.hh"

#include "system/file_path.hh"
#include "flow/darcy_flow_mh.hh"

#include "flow/darcy_flow_mh_output.hh"

#include "fem/mapping_p1.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include <fem/fe_rt.hh>
#include "quadrature/quadrature_lib.hh"

#include <limits>
#include <set>
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>

#include "tools/time_governor.hh"
#include "fields/field_algo_base.hh"
#include "fields/field.hh"
#include "fields/field_values.hh"
#include "fields/field_add_potential.hh"


#include "coupling/balance.hh"

#include "fields/vec_seq_double.hh"


FLOW123D_FORCE_LINK_IN_CHILD(darcy_flow_mh);




namespace it = Input::Type;

const it::Selection & DarcyFlowMH::get_mh_mortar_selection() {
	return it::Selection("MH_MortarMethod")
		.add_value(NoMortar, "None", "Mortar space: P0 on elements of lower dimension.")
		.add_value(MortarP0, "P0", "Mortar space: P0 on elements of lower dimension.")
		.add_value(MortarP1, "P1", "Mortar space: P1 on intersections, using non-conforming pressures.")
		.close();
}


const it::Selection & DarcyFlowMH::EqData::get_bc_type_selection() {
	return it::Selection("DarcyFlow_BC_Type")
               .add_value(none, "none", "Homogeneous Neumann boundary condition. Zero flux")
               .add_value(dirichlet, "dirichlet",
                       "Dirichlet boundary condition. "
                       "Specify the pressure head through the 'bc_pressure' field "
                       "or the piezometric head through the 'bc_piezo_head' field.")
//                .add_value(neumann, "neumann", "Neumann boundary condition. Prescribe water inflow by the 'bc_flux' field.")
//                .add_value(robin, "robin", "Robin boundary condition. Water inflow equal to (($\\sigma (h^R - h)$)). "
//                        "Specify the transition coefficient by 'bc_sigma' and the reference pressure head or pieaozmetric head "
//                        "through 'bc_pressure' and 'bc_piezo_head' respectively.")
               .add_value(total_flux, "total_flux", "Flux boundary condition (combines Neumann and Robin type). "
		       "Water inflow equal to (($q^N + \\sigma (h^R - h)$)). "
                        "Specify the water inflow by the 'bc_flux' field, the transition coefficient by 'bc_robin_sigma' "
			"and the reference pressure head or pieozmetric head through 'bc_pressure' or 'bc_piezo_head' respectively.")
			   .close();
}


const it::Record & DarcyFlowMH_Steady::get_input_type() {
    return it::Record("SteadyDarcy_MH", "Mixed-Hybrid  solver for STEADY saturated Darcy flow.")
		.derive_from(DarcyFlowInterface::get_input_type())
		.declare_key("n_schurs", it::Integer(0,2), it::Default("2"),
				"Number of Schur complements to perform when solving MH system.")
		.declare_key("solver", LinSys::get_input_type(), it::Default::obligatory(),
				"Linear solver for MH problem.")
		.declare_key("output",
		        DarcyFlowMHOutput::get_input_type(),
		        it::Default::obligatory(),
				"Parameters of output form MH module.")
		.declare_key("mortar_method", DarcyFlowMH::get_mh_mortar_selection(), it::Default("\"None\""),
				"Method for coupling Darcy flow between dimensions." )
		.declare_key("balance", Balance::get_input_type(), it::Default::obligatory(),
				"Settings for computing mass balance.")
		.declare_key("bc_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(),
				"Boundary condition for pressure as piezometric head." )
		.declare_key("init_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(),
				"Initial condition for pressure as piezometric head." )
		.declare_key("input_fields", it::Array(
					it::Record("DarcyFlowMH_Data", FieldCommon::field_descriptor_record_description("DarcyFlowMH_Data") )
					.copy_keys( DarcyFlowMH_Steady::EqData()
                        .make_field_descriptor_type("DarcyFlowMH_Data_aux") )
					.declare_key("bc_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "Boundary condition for pressure as piezometric head." )
					.declare_key("init_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "Initial condition for pressure as piezometric head." )
					.close()
					), it::Default::obligatory(), ""  )
		.close();
}


const int DarcyFlowMH_Steady::registrar =
		Input::register_class< DarcyFlowMH_Steady, Mesh &, const Input::Record >("SteadyDarcy_MH") +
		DarcyFlowMH_Steady::get_input_type().size();


const it::Record & DarcyFlowMH_Unsteady::get_input_type() {
	return it::Record("UnsteadyDarcy_MH", "Mixed-Hybrid solver for unsteady saturated Darcy flow.")
		.derive_from(DarcyFlowInterface::get_input_type())
		.copy_keys(DarcyFlowMH_Steady::get_input_type())
		.declare_key("time", TimeGovernor::get_input_type(), it::Default::obligatory(),
					 "Time governor setting for the unsteady Darcy flow model.")
		.close();
}


const int DarcyFlowMH_Unsteady::registrar =
		Input::register_class< DarcyFlowMH_Unsteady, Mesh &, const Input::Record >("UnsteadyDarcy_MH") +
		DarcyFlowMH_Unsteady::get_input_type().size();



DarcyFlowMH::EqData::EqData()
{
    ADD_FIELD(anisotropy, "Anisotropy of the conductivity tensor.", "1.0" );
    	anisotropy.units( UnitSI::dimensionless() );

    ADD_FIELD(cross_section, "Complement dimension parameter (cross section for 1D, thickness for 2D).", "1.0" );
    	cross_section.units( UnitSI().m(3).md() );

    ADD_FIELD(conductivity, "Isotropic conductivity scalar.", "1.0" );
    	conductivity.units( UnitSI().m().s(-1) );

    ADD_FIELD(sigma, "Transition coefficient between dimensions.", "1.0");
    	sigma.units( UnitSI::dimensionless() );

    ADD_FIELD(water_source_density, "Water source density.", "0.0");
    	water_source_density.units( UnitSI().s(-1) );
    
    ADD_FIELD(bc_type,"Boundary condition type, possible values:", "\"none\"" );
    	// TODO: temporary solution, we should try to get rid of pointer to the selection after having generic types
        bc_type.input_selection( &get_bc_type_selection() );
        bc_type.units( UnitSI::dimensionless() );

    ADD_FIELD(bc_pressure,"Prescribed pressure value for bc_type=\"dirichlet\" or reference pressure for bc_type=\"total_flux\".", "0.0");
    	bc_pressure.disable_where(bc_type, {none/*, neumann*/} );
        bc_pressure.units( UnitSI().m() );

    ADD_FIELD(bc_flux,"Incoming flux in total flux boundary condition.", "0.0");
    	bc_flux.disable_where(bc_type, {none, dirichlet/*, robin*/} );
        bc_flux.units( UnitSI().m(4).s(-1).md() );

    ADD_FIELD(bc_robin_sigma,"Conductivity coefficient in total flux boundary condition.", "0.0");
    	bc_robin_sigma.disable_where(bc_type, {none, dirichlet/*, neumann*/} );
        bc_robin_sigma.units( UnitSI().m(3).s(-1).md() );

    //these are for unsteady
    ADD_FIELD(init_pressure, "Initial condition as pressure", "0.0" );
    	init_pressure.units( UnitSI().m() );

    ADD_FIELD(storativity,"Storativity.", "1.0" );
    	storativity.units( UnitSI().m(-1) );

    //time_term_fields = this->subset({"storativity"});
    //main_matrix_fields = this->subset({"anisotropy", "conductivity", "cross_section", "sigma", "bc_type", "bc_robin_sigma"});
    //rhs_fields = this->subset({"water_source_density", "bc_pressure", "bc_flux"});

}


template<unsigned int dim>
DarcyFlowMH_Steady::Assembly<dim>::Assembly(AssemblyData ad)
: quad_(3),
  fe_values_(map_, quad_, fe_rt_, 
            update_values | update_gradients | update_JxW_values | update_quadrature_points),

  side_quad_(1),
  fe_side_values_(map_, side_quad_, fe_p_disc_, update_normal_vectors),

  velocity_interpolation_quad_(0), // veloctiy values in barycenter
  velocity_interpolation_fv_(map_,velocity_interpolation_quad_, fe_rt_, update_values | update_quadrature_points),

  d(ad)
{
}

template<unsigned int dim>
DarcyFlowMH_Steady::Assembly<dim>::~Assembly()
{
}


//=============================================================================
// CREATE AND FILL GLOBAL MH MATRIX OF THE WATER MODEL
// - do it in parallel:
//   - initial distribution of elements, edges
//
/*! @brief CREATE AND FILL GLOBAL MH MATRIX OF THE WATER MODEL
 *
 * Parameters {Solver,NSchurs} number of performed Schur
 * complements (0,1,2) for water flow MH-system
 *
 */
//=============================================================================
DarcyFlowMH_Steady::DarcyFlowMH_Steady(Mesh &mesh_in, const Input::Record in_rec, bool make_tg )
: DarcyFlowMH(mesh_in, in_rec)

{
    using namespace Input;
    START_TIMER("Darcy constructor");

    this->eq_data_ = &data_;

    //connecting data fields with mesh
    START_TIMER("data init");
    data_.set_mesh(mesh_in);
    //data_.gravity_ = arma::vec4( in_rec.val<std::string>("gravity") );
    data_.gravity_ =  arma::vec4(" 0 0 -1 0");
    data_.bc_pressure.add_factory(
    		std::make_shared<FieldAddPotential<3, FieldValue<3>::Scalar>::FieldFactory>
    		(data_.gravity_, "bc_piezo_head") );

    data_.set_input_list( in_rec.val<Input::Array>("input_fields") );


    
    END_TIMER("data init");

    
    size = mesh_->n_elements() + mesh_->n_sides() + mesh_->n_edges();
    n_schur_compls = in_rec.val<int>("n_schurs");
    
    solution = NULL;
    schur0   = NULL;

    
    mortar_method_= in_rec.val<MortarMethod>("mortar_method");
    if (mortar_method_ != NoMortar) {
        mesh_->make_intersec_elements();
    }

    mh_dh.reinit(mesh_);
    
    AssemblyData assembly_data;
    assembly_data.data = &data_;
    assembly_data.mesh = mesh_;
    assembly_data.mh_dh = &mh_dh;
    assembly_.push_back(new Assembly<1>(assembly_data));
    assembly_.push_back(new Assembly<2>(assembly_data));
    assembly_.push_back(new Assembly<3>(assembly_data));

    
    prepare_parallel(in_rec.val<AbstractRecord>("solver"));

    //side_ds->view( std::cout );
    //el_ds->view( std::cout );
    //edge_ds->view( std::cout );
    //rows_ds->view( std::cout );
    

    // initialization of balance object
    Input::Iterator<Input::Record> it = in_rec.find<Input::Record>("balance");
    if (it->val<bool>("balance_on"))
    {
        balance_ = boost::make_shared<Balance>("water", mesh_, *it);
        water_balance_idx_ = balance_->add_quantity("water_volume");
        balance_->allocate(rows_ds->lsize(), 1);
        balance_->units(UnitSI().m(3));
    }


    if (make_tg) {
    	// steady time governor
    	time_ = new TimeGovernor();
    	data_.mark_input_times(this->mark_type());
    	data_.set_time(time_->step(), LimitSide::right);

    	create_linear_system();
    	output_object = new DarcyFlowMHOutput(this, in_rec.val<Input::Record>("output"));
    }
}



//=============================================================================
// COMPOSE and SOLVE WATER MH System possibly through Schur complements
//=============================================================================
void DarcyFlowMH_Steady::update_solution() {
    START_TIMER("Solving MH system");


    if (time_->is_end()) return;

    if (! time_->is_steady()) time_->next_time();
    


    assembly_linear_system();

    int convergedReason = schur0->solve();

    xprintf(MsgLog, "Linear solver ended with reason: %d \n", convergedReason );
    ASSERT( convergedReason >= 0, "Linear solver failed to converge. Convergence reason %d \n", convergedReason );

    this -> postprocess();

    solution_changed_for_scatter=true;

    output_data();

    if (time_->is_steady()) time_->next_time();
}

void DarcyFlowMH_Steady::postprocess() 
{
    START_TIMER("postprocess");
    //ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);

    // modify side fluxes in parallel
    // for every local edge take time term on digonal and add it to the corresponding flux
    /*
    for (unsigned int i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {
        ele = mesh_->element(el_4_loc[i_loc]);
        FOR_ELEMENT_SIDES(ele,i) {
            side_rows[i] = side_row_4_id[ mh_dh.side_dof( ele->side(i) ) ];
            values[i] = -1.0 * ele->measure() *
              data.cross_section.value(ele->centre(), ele->element_accessor()) *
              data.water_source_density.value(ele->centre(), ele->element_accessor()) /
              ele->n_sides();
        }
        VecSetValues(schur0->get_solution(), ele->n_sides(), side_rows, values, ADD_VALUES);
    }
    VecAssemblyBegin(schur0->get_solution());
    VecAssemblyEnd(schur0->get_solution());
    */
}


void DarcyFlowMH_Steady::output_data() {
    START_TIMER("Darcy output data");
    time_->view("DARCY"); //time governor information output
	this->output_object->output();


	if (balance_ != nullptr)
	{
	    START_TIMER("Darcy balance output");
		if (balance_->cumulative() && time_->tlevel() > 0)
		{
			balance_->calculate_cumulative_sources(water_balance_idx_, schur0->get_solution(), time_->dt());
			balance_->calculate_cumulative_fluxes(water_balance_idx_, schur0->get_solution(), time_->dt());
		}

		if (time_->is_current( TimeGovernor::marks().type_output() ))
		{
			balance_->calculate_mass(water_balance_idx_, schur0->get_solution());
			balance_->calculate_source(water_balance_idx_, schur0->get_solution());
			balance_->calculate_flux(water_balance_idx_, schur0->get_solution());
			balance_->output(time_->t());
		}
	}
}

double DarcyFlowMH_Steady::solution_precision() const
{
	return schur0->get_solution_precision();
}


void  DarcyFlowMH_Steady::get_solution_vector(double * &vec, unsigned int &vec_size)
{
    // TODO: make class for vectors (wrapper for PETSC or other) derived from LazyDependency
    // and use its mechanism to manage dependency between vectors
    if (solution_changed_for_scatter) {

        // scatter solution to all procs
        VecScatterBegin(par_to_all, schur0->get_solution(), sol_vec, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(  par_to_all, schur0->get_solution(), sol_vec, INSERT_VALUES, SCATTER_FORWARD);
        solution_changed_for_scatter=false;
    }

    vec = solution;
    vec_size = this->size;
    ASSERT(vec != NULL, "Requested solution is not allocated!\n");
}

void  DarcyFlowMH_Steady::get_parallel_solution_vector(Vec &vec)
{
    vec=schur0->get_solution();
    ASSERT(vec != NULL, "Requested solution is not allocated!\n");
}


template<unsigned int dim>
void DarcyFlowMH_Steady::Assembly<dim>::assembly_local_matrix(arma::mat& local_matrix, ElementFullIter ele)
{
    //START_TIMER("Assembly<dim>::assembly_local_matrix");
    fe_values_.reinit(ele);
    const unsigned int ndofs = fe_values_.get_fe()->n_dofs(), qsize = fe_values_.get_quadrature()->size();
    local_matrix.zeros(ndofs, ndofs);

    double scale = 1
                   / d.data->conductivity.value( ele->centre(), ele->element_accessor() )
                   / d.data->cross_section.value( ele->centre(), ele->element_accessor() );
                           
    for (unsigned int k=0; k<qsize; k++)
    {
        for (unsigned int i=0; i<ndofs; i++)
        {
             for (unsigned int j=0; j<ndofs; j++)
                local_matrix[i*ndofs+j] += 
                        scale
                        * arma::dot(fe_values_.shape_vector(i,k),
                                    (d.data->anisotropy.value(ele->centre(), ele->element_accessor() )).i()
                                     * fe_values_.shape_vector(j,k)
                                   ) 
                        * fe_values_.JxW(k);
        }
    }
}

template<unsigned int dim>
void DarcyFlowMH_Steady::Assembly<dim>::assembly_local_vb(double* local_vb, ElementFullIter ele, Neighbour *ngh)
{
    //START_TIMER("Assembly<dim>::assembly_local_vb");
    // compute normal vector to side
    arma::vec3 nv;
    ElementFullIter ele_higher = d.mesh->element.full_iter(ngh->side()->element());
    fe_side_values_.reinit(ele_higher, ngh->side()->el_idx());
    nv = fe_side_values_.normal_vector(0);

    double value = d.data->sigma.value( ele->centre(), ele->element_accessor()) *
                    2*d.data->conductivity.value( ele->centre(), ele->element_accessor()) *
                    arma::dot(d.data->anisotropy.value( ele->centre(), ele->element_accessor())*nv, nv) *
                    d.data->cross_section.value( ngh->side()->centre(), ele_higher->element_accessor() ) * // cross-section of higher dim. (2d)
                    d.data->cross_section.value( ngh->side()->centre(), ele_higher->element_accessor() ) /
                    d.data->cross_section.value( ele->centre(), ele->element_accessor() ) *      // crossection of lower dim.
                    ngh->side()->measure();

    local_vb[0] = -value;   local_vb[1] = value;
    local_vb[2] = value;    local_vb[3] = -value;
}



template<unsigned int dim>
arma::vec3 DarcyFlowMH_Steady::Assembly<dim>::make_element_vector(ElementFullIter ele)
{
    //START_TIMER("Assembly<dim>::make_element_vector");
    arma::vec3 flux_in_center;
    flux_in_center.zeros();
    
    velocity_interpolation_fv_.reinit(ele);
    for (unsigned int li = 0; li < ele->n_sides(); li++) {
        flux_in_center += d.mh_dh->side_flux( *(ele->side( li ) ) )
                  * velocity_interpolation_fv_.shape_vector(li,0);
    }

    flux_in_center /= d.data->cross_section.value(ele->centre(), ele->element_accessor() );
    return flux_in_center;
}


// ===========================================================================================
//
//   MATRIX ASSEMBLY - we use abstract assembly routine, where  LS Mat/Vec SetValues
//   are in fact pointers to allocating or filling functions - this is governed by Linsystem roitunes
//
// =======================================================================================

void DarcyFlowMH_Steady::assembly_steady_mh_matrix()
{
    START_TIMER("DarcyFlowMH_Steady::assembly_steady_mh_matrix");
    
    LinSys *ls = schur0;
    ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);

    class Boundary *bcd;
    class Neighbour *ngh;

    bool fill_matrix = ls->is_preallocated();
    int el_row, side_row, edge_row, loc_b = 0;
    int tmp_rows[100];
    int side_rows[4], edge_rows[4]; // rows for sides and edges of one element
    double local_vb[4]; // 2x2 matrix

    // to make space for second schur complement, max. 10 neighbour edges of one el.
    double zeros[1000];
    for(int i=0; i<1000; i++) zeros[i]=0.0;

    double minus_ones[4] = { -1.0, -1.0, -1.0, -1.0 };
    double loc_side_rhs[4];

    if (balance_ != nullptr)
        balance_->start_flux_assembly(water_balance_idx_);

    for (unsigned int i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {
        arma::mat local_matrix;
        ele = mesh_->element(el_4_loc[i_loc]);
        el_row = row_4_el[el_4_loc[i_loc]];
        unsigned int nsides = ele->n_sides();
        
        if (fill_matrix) 
            assembly_[ele->dim()-1]->assembly_local_matrix(local_matrix, ele);
        
        double cross_section = data_.cross_section.value(ele->centre(), ele->element_accessor());

        for (unsigned int i = 0; i < nsides; i++) {
            side_row = side_rows[i] = side_row_4_id[ mh_dh.side_dof( ele->side(i) ) ];
            edge_row = edge_rows[i] = row_4_edge[ele->side(i)->edge_idx()];
            bcd=ele->side(i)->cond();

            // gravity term on RHS
            loc_side_rhs[i] = (ele->centre()[ 2 ] - ele->side(i)->centre()[ 2 ]);

            // set block C and C': side-edge, edge-side
            double c_val = 1.0;

            if (bcd) {
                ElementAccessor<3> b_ele = bcd->element_accessor();
                EqData::BC_Type type = (EqData::BC_Type)data_.bc_type.value(b_ele.centre(), b_ele);
                if ( type == EqData::none) {
                    // homogeneous neumann
                } else if ( type == EqData::dirichlet ) {
                    c_val = 0.0;
                    double bc_pressure = data_.bc_pressure.value(b_ele.centre(), b_ele);
                    loc_side_rhs[i] -= bc_pressure;
                    ls->rhs_set_value(edge_row, -bc_pressure);
                    ls->mat_set_value(edge_row, edge_row, -1.0);

//                } else if ( type == EqData::neumann) {
//                    double bc_flux = -data_.bc_flux.value(b_ele.centre(), b_ele);
//                    ls->rhs_set_value(edge_row, bc_flux * bcd->element()->measure() * cross_section);
//
//                } else if ( type == EqData::robin) {
//                    double bc_pressure = data_.bc_pressure.value(b_ele.centre(), b_ele);
//                    double bc_sigma = data_.bc_robin_sigma.value(b_ele.centre(), b_ele);
//                    ls->rhs_set_value(edge_row, -bcd->element()->measure() * bc_sigma * bc_pressure * cross_section );
//                    ls->mat_set_value(edge_row, edge_row, -bcd->element()->measure() * bc_sigma * cross_section );
//
		} else if ( type == EqData::total_flux) {
		    // internally we work with outward flux
		    double bc_flux = -data_.bc_flux.value(b_ele.centre(), b_ele);
		    double bc_pressure = data_.bc_pressure.value(b_ele.centre(), b_ele);
                    double bc_sigma = data_.bc_robin_sigma.value(b_ele.centre(), b_ele);
                    ls->mat_set_value(edge_row, edge_row, -bcd->element()->measure() * bc_sigma * cross_section );
                    ls->rhs_set_value(edge_row, (bc_flux - bc_sigma * bc_pressure) * bcd->element()->measure() * cross_section);

                } else {
                    xprintf(UsrErr, "BC type not supported.\n");
                }

                if (balance_ != nullptr)
                {
                    balance_->add_flux_matrix_values(water_balance_idx_, loc_b, {side_row}, {1});
                }
                ++loc_b;
            }
            ls->mat_set_value(side_row, edge_row, c_val);
            ls->mat_set_value(edge_row, side_row, c_val);

            // assemble matrix for weights in BDDCML
            // approximation to diagonal of 
            // S = -C - B*inv(A)*B'
            // as 
            // diag(S) ~ - diag(C) - 1./diag(A)
            // the weights form a partition of unity to average a discontinuous solution from neighbouring subdomains
            // to a continuous one
            // it is important to scale the effect - if conductivity is low for one subdomain and high for the other,
            // trust more the one with low conductivity - it will be closer to the truth than an arithmetic average
            if ( typeid(*ls) == typeid(LinSys_BDDC) ) {
               double val_side =  local_matrix(i,i);
               double val_edge =  -1./local_matrix(i,i);

               static_cast<LinSys_BDDC*>(ls)->diagonal_weights_set_value( side_row, val_side );
               static_cast<LinSys_BDDC*>(ls)->diagonal_weights_set_value( edge_row, val_edge );
            }
        }

        ls->rhs_set_values(nsides, side_rows, loc_side_rhs);

        
        // set block A: side-side on one element - block diagonal matrix
        ls->mat_set_values(nsides, side_rows, nsides, side_rows, local_matrix.memptr());
        // set block B, B': element-side, side-element
        ls->mat_set_values(1, &el_row, nsides, side_rows, minus_ones);
        ls->mat_set_values(nsides, side_rows, 1, &el_row, minus_ones);


        // D block: non-compatible conections and diagonal: element-element

        ls->mat_set_value(el_row, el_row, 0.0);         // maybe this should be in virtual block for schur preallocation

        if ( typeid(*ls) == typeid(LinSys_BDDC) ) {
           double val_ele =  1.;
           static_cast<LinSys_BDDC*>(ls)->diagonal_weights_set_value( el_row, val_ele );
        }

        // D, E',E block: compatible connections: element-edge
        
        for (unsigned int i = 0; i < ele->n_neighs_vb; i++) {
            // every compatible connection adds a 2x2 matrix involving
            // current element pressure  and a connected edge pressure
            ngh= ele->neigh_vb[i];
            tmp_rows[0]=el_row;
            tmp_rows[1]=row_4_edge[ ngh->edge_idx() ];
            
            assembly_[ngh->side()->dim()]->assembly_local_vb(local_vb, ele, ngh);
            
            ls->mat_set_values(2, tmp_rows, 2, tmp_rows, local_vb);

            // update matrix for weights in BDDCML
            if ( typeid(*ls) == typeid(LinSys_BDDC) ) {
               int ind = tmp_rows[1];
               // there is -value on diagonal in block C!
               double new_val = local_vb[0];
               static_cast<LinSys_BDDC*>(ls)->diagonal_weights_set_value( ind, new_val );
            }

            if (n_schur_compls == 2) {
                // for 2. Schur: N dim edge is conected with N dim element =>
                // there are nz between N dim edge and N-1 dim edges of the element

                ls->mat_set_values(nsides, edge_rows, 1, tmp_rows+1, zeros);
                ls->mat_set_values(1, tmp_rows+1, nsides, edge_rows, zeros);

                // save all global edge indices to higher positions
                tmp_rows[2+i] = tmp_rows[1];
            }
        }


        // add virtual values for schur complement allocation
        switch (n_schur_compls) {
        case 2:
            // Connections between edges of N+1 dim. elements neighboring with actual N dim element 'ele'
            ASSERT(ele->n_neighs_vb*ele->n_neighs_vb<1000, "Too many values in E block.");
            ls->mat_set_values(ele->n_neighs_vb, tmp_rows+2,
                               ele->n_neighs_vb, tmp_rows+2, zeros);

        case 1: // included also for case 2
            // -(C')*(A-)*B block and its transpose conect edge with its elements
            ls->mat_set_values(1, &el_row, nsides, edge_rows, zeros);
            ls->mat_set_values(nsides, edge_rows, 1, &el_row, zeros);
            // -(C')*(A-)*C block conect all edges of every element
            ls->mat_set_values(nsides, edge_rows, nsides, edge_rows, zeros);
        }
    }    
    
    if (balance_ != nullptr)
        balance_->finish_flux_assembly(water_balance_idx_);

    assembly_source_term();


    if (mortar_method_ == MortarP0) {
        P0_CouplingAssembler(*this).assembly(*ls);
    } else if (mortar_method_ == MortarP1) {
        P1_CouplingAssembler(*this).assembly(*ls);
    }  
}


void DarcyFlowMH_Steady::assembly_source_term()
{
    START_TIMER("assembly source term");
    if (balance_ != nullptr)
    	balance_->start_source_assembly(water_balance_idx_);

    for (unsigned int i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {

        ElementFullIter ele = mesh_->element(el_4_loc[i_loc]);
        int el_row = row_4_el[el_4_loc[i_loc]];

        // set sources
        double source = ele->measure() *
                data_.cross_section.value(ele->centre(), ele->element_accessor()) *
                data_.water_source_density.value(ele->centre(), ele->element_accessor());
        schur0->rhs_set_value(el_row, -1.0 * source );

        if (balance_ != nullptr)
        	balance_->add_source_rhs_values(water_balance_idx_, ele->region().bulk_idx(), {el_row}, {source});
    }

    if (balance_ != nullptr)
    	balance_->finish_source_assembly(water_balance_idx_);
}


void P0_CouplingAssembler::pressure_diff(int i_ele,
		vector<int> &dofs, unsigned int &ele_type, double &delta, arma::vec &dirichlet) {

	const Element *ele;

	if (i_ele == (int)(ml_it_->size()) ) { // master element .. 1D
		ele_type = 0;
		delta = -delta_0;
		ele=master_;
	} else {
		ele_type = 1;
		const Intersection &isect=intersections_[ (*ml_it_)[i_ele] ];
		delta = isect.intersection_true_size();
		ele = isect.slave_iter();
	}

	dofs.resize(ele->n_sides());
        dirichlet.resize(ele->n_sides());
        dirichlet.zeros();

	for(unsigned int i_side=0; i_side < ele->n_sides(); i_side++ ) {
		dofs[i_side]=darcy_.row_4_edge[ele->side(i_side)->edge_idx()];
		Boundary * bcd = ele->side(i_side)->cond();
		if (bcd) {			
			ElementAccessor<3> b_ele = bcd->element_accessor();
			DarcyFlowMH::EqData::BC_Type type = (DarcyFlowMH::EqData::BC_Type)darcy_.data_.bc_type.value(b_ele.centre(), b_ele);
			//DBGMSG("bcd id: %d sidx: %d type: %d\n", ele->id(), i_side, type);
			if (type == DarcyFlowMH::EqData::dirichlet) {
				//DBGMSG("Dirichlet: %d\n", ele->index());
				dofs[i_side] = -dofs[i_side];
				double bc_pressure = darcy_.data_.bc_pressure.value(b_ele.centre(), b_ele);
				dirichlet[i_side] = bc_pressure;
			}
		} 
	}

}

/**
 * Works well but there is large error next to the boundary.
 */
 void P0_CouplingAssembler::assembly(LinSys &ls) {
	double delta_i, delta_j;
	arma::mat product;
	arma::vec dirichlet_i, dirichlet_j;
	unsigned int ele_type_i, ele_type_j;

	unsigned int i,j;
	vector<int> dofs_i,dofs_j;

	for(ml_it_ = master_list_.begin(); ml_it_ != master_list_.end(); ++ml_it_) {

    	if (ml_it_->size() == 0) continue; // skip empty masters


		// on the intersection element we consider
		// * intersection dofs for master and slave
		//   those are dofs of the space into which we interpolate
		//   base functions from individual master and slave elements
		//   For the master dofs both are usualy eqivalent.
		// * original dofs - for master same as intersection dofs, for slave
		//   all dofs of slave elements

		// form list of intersection dofs, in this case pressures in barycenters
		// but we do not use those form MH system in order to allow second schur somplement. We rather map these
		// dofs to pressure traces, i.e. we use averages of traces as barycentric values.


		master_ = intersections_[ml_it_->front()].master_iter();
		delta_0 = master_->measure();

		double master_sigma=darcy_.data_.sigma.value( master_->centre(), master_->element_accessor());

		// rows
		double check_delta_sum=0;
		for(i = 0; i <= ml_it_->size(); ++i) {
			pressure_diff(i, dofs_i, ele_type_i, delta_i, dirichlet_i);
			check_delta_sum+=delta_i;
			//columns
			for (j = 0; j <= ml_it_->size(); ++j) {
				pressure_diff(j, dofs_j, ele_type_j, delta_j, dirichlet_j);

				double scale =  -master_sigma * delta_i * delta_j / delta_0;
				product = scale * tensor_average[ele_type_i][ele_type_j];

				arma::vec rhs(dofs_i.size());
				rhs.zeros();
				ls.set_values( dofs_i, dofs_j, product, rhs, dirichlet_i, dirichlet_j);
				auto dofs_i_cp=dofs_i;
				auto dofs_j_cp=dofs_j;
				ls.set_values( dofs_i_cp, dofs_j_cp, product, rhs, dirichlet_i, dirichlet_j);
			}
		}
                ASSERT(check_delta_sum < 1E-5*delta_0, "sum err %f > 0\n", check_delta_sum/delta_0);
    } // loop over master elements
}



 void P1_CouplingAssembler::add_sides(const Element * ele, unsigned int shift, vector<int> &dofs, vector<double> &dirichlet)
 {

		for(unsigned int i_side=0; i_side < ele->n_sides(); i_side++ ) {
			dofs[shift+i_side] =  darcy_.row_4_edge[ele->side(i_side)->edge_idx()];
			Boundary * bcd = ele->side(i_side)->cond();

			if (bcd) {
				ElementAccessor<3> b_ele = bcd->element_accessor();
				DarcyFlowMH::EqData::BC_Type type = (DarcyFlowMH::EqData::BC_Type)darcy_.data_.bc_type.value(b_ele.centre(), b_ele);
				//DBGMSG("bcd id: %d sidx: %d type: %d\n", ele->id(), i_side, type);
				if (type == DarcyFlowMH::EqData::dirichlet) {
					//DBGMSG("Dirichlet: %d\n", ele->index());
					dofs[shift + i_side] = -dofs[shift + i_side];
					double bc_pressure = darcy_.data_.bc_pressure.value(b_ele.centre(), b_ele);
					dirichlet[shift + i_side] = bc_pressure;
				}
			}
		}
 }


/**
 * P1 connection of different dimensions
 *
 * - 20.11. 2014 - very poor convergence, big error in pressure even at internal part of the fracture
 */

void P1_CouplingAssembler::assembly(LinSys &ls) {

	for (const Intersection &intersec : intersections_) {
    	const Element * master = intersec.master_iter();
       	const Element * slave = intersec.slave_iter();

	add_sides(slave, 0, dofs, dirichlet);
       	add_sides(master, 3, dofs, dirichlet);
       	
		double master_sigma=darcy_.data_.sigma.value( master->centre(), master->element_accessor());

/*
 * Local coordinates on 1D
 *         t0
 * node 0: 0.0
 * node 1: 1.0
 *
 * base fce points
 * t0 = 0.0    on side 0 node 0
 * t0 = 1.0    on side 1 node 1
 *
 * Local coordinates on 2D
 *         t0  t1
 * node 0: 0.0 0.0
 * node 1: 1.0 0.0
 * node 2: 0.0 1.0
 *
 * base fce points
 * t0=0.5, t1=0.0        on side 0 nodes (0,1)
 * t0=0.5, t1=0.5        on side 1 nodes (1,2)
 * t0=0.0, t1=0.5        on side 2 nodes (2,0)
 */



        arma::vec point_Y(1);
        point_Y.fill(1.0);
        arma::vec point_2D_Y(intersec.map_to_slave(point_Y)); // local coordinates of  Y on slave (1, t0, t1)
        arma::vec point_1D_Y(intersec.map_to_master(point_Y)); //  local coordinates of  Y on master (1, t0)

        arma::vec point_X(1);
        point_X.fill(0.0);
        arma::vec point_2D_X(intersec.map_to_slave(point_X)); // local coordinates of  X on slave (1, t0, t1)
        arma::vec point_1D_X(intersec.map_to_master(point_X)); // local coordinates of  X on master (1, t0)

        arma::mat base_2D(3, 3);
        // basis functions are numbered as sides
        // TODO:
        // Use RT finite element to evaluate following matrices.

        // Ravirat - Thomas base functions evaluated in points (0,0), (1,0), (0,1)
        // 2D RT_i(t0, t1) = a0 + a1*t0 + a2*t1
        //         a0     a1      a2
        base_2D << 1.0 << 0.0 << -2.0 << arma::endr // RT for side 0
                << 1.0 << -2.0 << 0.0 << arma::endr // RT for side 1
                << -1.0 << 2.0 << 2.0 << arma::endr;// RT for side 2
                

        arma::mat base_1D(2, 2);
        // Ravirat - Thomas base functions evaluated in points (0,0), (1,0), (0,1)
        // 1D RT_i(t0) =   a0 + a1 * t0
        //          a0     a1
        base_1D << 1.0 << -1.0 << arma::endr // RT for side 0,
                << 0.0 << 1.0 << arma::endr; // RT for side 1,



        // Consider both 2D and 1D value are defined for the test function
        // related to the each of 5 DOFs involved in the intersection.
        // One of these values is always zero.
        // Compute difference of the 2D and 1D value for every DOF.
        // Compute value of this difference in both endpoints X,Y of the intersection.

        arma::vec difference_in_Y(5);
        arma::vec difference_in_X(5);
        // slave sides 0,1,2
        difference_in_Y.subvec(0, 2) = -base_2D * point_2D_Y;
        difference_in_X.subvec(0, 2) = -base_2D * point_2D_X;
        // master sides 3,4
        difference_in_Y.subvec(3, 4) = base_1D * point_1D_Y;
        difference_in_X.subvec(3, 4) = base_1D * point_1D_X;

        // applying the Simpson's rule
        // to the product of two linear functions f, g we get
        // (b-a)/6 * ( 3*f(a)*g(a) + 3*f(b)*g(b) + 2*f(a)*g(b) + 2*f(b)*g(a) )
        arma::mat A(5, 5);
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                A(i, j) = -master_sigma * intersec.intersection_true_size() *
                        ( difference_in_Y[i] * difference_in_Y[j]
                          + difference_in_Y[i] * difference_in_X[j]/2
                          + difference_in_X[i] * difference_in_Y[j]/2
                          + difference_in_X[i] * difference_in_X[j]
                        ) * (1.0 / 3.0);

            }
        }
        auto dofs_cp=dofs;
        ls.set_values( dofs_cp, dofs_cp, A, rhs, dirichlet, dirichlet);

    }
}



/*******************************************************************************
 * COMPOSE WATER MH MATRIX WITHOUT SCHUR COMPLEMENT
 ******************************************************************************/

void DarcyFlowMH_Steady::create_linear_system() {
  
    START_TIMER("preallocation");

    auto in_rec = this->input_record_.val<Input::AbstractRecord>("solver");

    if (schur0 == NULL) { // create Linear System for MH matrix
       
    	if (in_rec.type() == LinSys_BDDC::get_input_type()) {
#ifdef FLOW123D_HAVE_BDDCML
            xprintf(Warn, "For BDDC is using no Schur complements.");
            n_schur_compls = 0;
            LinSys_BDDC *ls = new LinSys_BDDC(global_row_4_sub_row->size(), &(*rows_ds),
                    3,  // 3 == la::BddcmlWrapper::SPD_VIA_SYMMETRICGENERAL
                    1,  // 1 == number of subdomains per process
                    true); // swap signs of matrix and rhs to make the matrix SPD
            ls->set_from_input(in_rec);
            ls->set_solution( NULL );
            // possible initialization particular to BDDC
            START_TIMER("BDDC set mesh data");
            set_mesh_data_for_bddc(ls);
            schur0=ls;
            END_TIMER("BDDC set mesh data");
#else
            xprintf(Err, "Flow123d was not build with BDDCML support.\n");
#endif // FLOW123D_HAVE_BDDCML
        } 
        else if (in_rec.type() == LinSys_PETSC::get_input_type()) {
        // use PETSC for serial case even when user wants BDDC
            if (n_schur_compls > 2) {
                xprintf(Warn, "Invalid number of Schur Complements. Using 2.");
                n_schur_compls = 2;
            }

            LinSys_PETSC *schur1, *schur2;

            if (n_schur_compls == 0) {
                LinSys_PETSC *ls = new LinSys_PETSC( &(*rows_ds) );

                // temporary solution; we have to set precision also for sequantial case of BDDC
                // final solution should be probably call of direct solver for oneproc case
                if (in_rec.type() != LinSys_BDDC::get_input_type()) ls->set_from_input(in_rec);
                else {
                    ls->LinSys::set_from_input(in_rec); // get only common options
                }

                ls->set_solution( NULL );
                schur0=ls;
            } else {
                IS is;
                ISCreateStride(PETSC_COMM_WORLD, side_ds->lsize(), rows_ds->begin(), 1, &is);
                //ASSERT(err == 0,"Error in ISCreateStride.");

                SchurComplement *ls = new SchurComplement(is, &(*rows_ds));
                ls->set_from_input(in_rec);
                ls->set_solution( NULL );

                // make schur1
                Distribution *ds = ls->make_complement_distribution();
                if (n_schur_compls==1) {
                    schur1 = new LinSys_PETSC(ds);
                    schur1->set_positive_definite();
                } else {
                    IS is;
                    ISCreateStride(PETSC_COMM_WORLD, el_ds->lsize(), ls->get_distribution()->begin(), 1, &is);
                    //ASSERT(err == 0,"Error in ISCreateStride.");
                    SchurComplement *ls1 = new SchurComplement(is, ds); // is is deallocated by SchurComplement
                    ls1->set_negative_definite();

                    // make schur2
                    schur2 = new LinSys_PETSC( ls1->make_complement_distribution() );
                    schur2->set_positive_definite();
                    ls1->set_complement( schur2 );
                    schur1 = ls1;
                }
                ls->set_complement( schur1 );
                schur0=ls;
            }

            START_TIMER("PETSC PREALLOCATION");
            schur0->set_symmetric();
            schur0->start_allocation();
            assembly_steady_mh_matrix(); // preallocation
    	    VecZeroEntries(schur0->get_solution());
            END_TIMER("PETSC PREALLOCATION");
        }
        else {
            xprintf(Err, "Unknown solver type. Internal error.\n");
        }
    }


    END_TIMER("preallocation");

    make_serial_scatter();

}




void DarcyFlowMH_Steady::assembly_linear_system() {
    START_TIMER("DarcyFlowMH_Steady::assembly_linear_system");

    {
        data_.set_time(time_->step(), LimitSide::right);
    }
	//DBGMSG("Assembly linear system\n");
	if (data_.changed()) {
		//DBGMSG("  Data changed\n");
		// currently we have no optimization for cases when just time term data or RHS data are changed
	    START_TIMER("full assembly");
        if (typeid(*schur0) != typeid(LinSys_BDDC)) {
            schur0->start_add_assembly(); // finish allocation and create matrix
        }
	    schur0->mat_zero_entries();
	    schur0->rhs_zero_entries();
	    assembly_steady_mh_matrix(); // fill matrix
	    schur0->finish_assembly();
	    schur0->set_matrix_changed();
            //MatView( *const_cast<Mat*>(schur0->get_matrix()), PETSC_VIEWER_STDOUT_WORLD  );
            //VecView( *const_cast<Vec*>(schur0->get_rhs()),   PETSC_VIEWER_STDOUT_WORLD);

	    if (!time_->is_steady()) {
	        START_TIMER("fix time term");
	    	//DBGMSG("    setup time term\n");
	    	// assembly time term and rhs
	    	setup_time_term();
	    	modify_system();
	    }
	    else if (balance_ != nullptr)
	    {
	    	balance_->start_mass_assembly(water_balance_idx_);
	    	balance_->finish_mass_assembly(water_balance_idx_);
	    }
	    END_TIMER("full assembly");
	} else {
		START_TIMER("modify system");
		if (!time_->is_steady()) {
			modify_system();
		} else {
			xprintf(PrgErr, "Planned computation time for steady solver, but data are not changed.\n");
		}
		END_TIMER("modiffy system");
	}

}



void DarcyFlowMH_Steady::set_mesh_data_for_bddc(LinSys_BDDC * bddc_ls) {
    START_TIMER("DarcyFlowMH_Steady::set_mesh_data_for_bddc");
    // prepare mesh for BDDCML
    // initialize arrays
    // auxiliary map for creating coordinates of local dofs and global-to-local numbering
    std::map<int,arma::vec3> localDofMap;
    // connectivity for the subdomain, i.e. global dof numbers on element, stored element-by-element
    // Indices of Nodes on Elements
    std::vector<int> inet;
    // number of degrees of freedom on elements - determines elementwise chunks of INET array
    // Number of Nodes on Elements
    std::vector<int> nnet;
    // Indices of Subdomain Elements in Global Numbering - for local elements, their global indices
    std::vector<int> isegn;

    // This array is currently not used in BDDCML, it was used as an interface scaling alternative to scaling
    // by diagonal. It corresponds to the rho-scaling.
    std::vector<double> element_permeability;

    // maximal and minimal dimension of elements
    int elDimMax = 1;
    int elDimMin = 3;
    for ( unsigned int i_loc = 0; i_loc < el_ds->lsize(); i_loc++ ) {
        // for each element, create local numbering of dofs as fluxes (sides), pressure (element centre), Lagrange multipliers (edges), compatible connections
        ElementFullIter el = mesh_->element(el_4_loc[i_loc]);
        int e_idx = el.index();

        int elDim = el->dim();
        elDimMax = std::max( elDimMax, elDim );
        elDimMin = std::min( elDimMin, elDim );

        isegn.push_back( e_idx );
        int nne = 0;

        FOR_ELEMENT_SIDES(el,si) {
            // insert local side dof
            int side_row = side_row_4_id[ mh_dh.side_dof( el->side(si) ) ];
            arma::vec3 coord = el->side(si)->centre();

            localDofMap.insert( std::make_pair( side_row, coord ) );
            inet.push_back( side_row );
            nne++;
        }

        // insert local pressure dof
        int el_row  = row_4_el[ el_4_loc[i_loc] ];
        arma::vec3 coord = el->centre();
        localDofMap.insert( std::make_pair( el_row, coord ) );
        inet.push_back( el_row );
        nne++;

        FOR_ELEMENT_SIDES(el,si) {
            // insert local edge dof
            int edge_row = row_4_edge[ el->side(si)->edge_idx() ];
            arma::vec3 coord = el->side(si)->centre();

            localDofMap.insert( std::make_pair( edge_row, coord ) );
            inet.push_back( edge_row );
            nne++;
        }

        // insert dofs related to compatible connections
        for ( unsigned int i_neigh = 0; i_neigh < el->n_neighs_vb; i_neigh++) {
            int edge_row = row_4_edge[ el->neigh_vb[i_neigh]->edge_idx()  ];
            arma::vec3 coord = el->neigh_vb[i_neigh]->edge()->side(0)->centre();

            localDofMap.insert( std::make_pair( edge_row, coord ) );
            inet.push_back( edge_row );
            nne++;
        }

        nnet.push_back( nne );

        // version for rho scaling
        // trace computation
        arma::vec3 centre = el->centre();
        double conduct = data_.conductivity.value( centre , el->element_accessor() );
        arma::mat33 aniso = data_.anisotropy.value( centre, el->element_accessor() );

        // compute mean on the diagonal
        double coef = 0.;
        for ( int i = 0; i < 3; i++) {
            coef = coef + aniso.at(i,i);
        }
        // Maybe divide by cs
        coef = conduct*coef / 3;

        ASSERT( coef > 0.,
                "Zero coefficient of hydrodynamic resistance %f . \n ", coef );
        element_permeability.push_back( 1. / coef );
    }
    //convert set of dofs to vectors
    // number of nodes (= dofs) on the subdomain
    int numNodeSub = localDofMap.size();
    ASSERT_EQUAL( (unsigned int)numNodeSub, global_row_4_sub_row->size() );
    // Indices of Subdomain Nodes in Global Numbering - for local nodes, their global indices
    std::vector<int> isngn( numNodeSub );
    // pseudo-coordinates of local nodes (i.e. dofs)
    // they need not be exact, they are used just for some geometrical considerations in BDDCML, 
    // such as selection of corners maximizing area of a triangle, bounding boxes fro subdomains to 
    // find candidate neighbours etc.
    std::vector<double> xyz( numNodeSub * 3 ) ;
    int ind = 0;
    std::map<int,arma::vec3>::iterator itB = localDofMap.begin();
    for ( ; itB != localDofMap.end(); ++itB ) {
        isngn[ind] = itB -> first;

        arma::vec3 coord = itB -> second;
        for ( int j = 0; j < 3; j++ ) {
            xyz[ j*numNodeSub + ind ] = coord[j];
        }

        ind++;
    }
    localDofMap.clear();

    // Number of Nodal Degrees of Freedom
    // nndf is trivially one - dofs coincide with nodes
    std::vector<int> nndf( numNodeSub, 1 );

    // prepare auxiliary map for renumbering nodes
    typedef std::map<int,int> Global2LocalMap_; //! type for storage of global to local map
    Global2LocalMap_ global2LocalNodeMap;
    for ( unsigned ind = 0; ind < isngn.size(); ++ind ) {
        global2LocalNodeMap.insert( std::make_pair( static_cast<unsigned>( isngn[ind] ), ind ) );
    }

    // renumber nodes in the inet array to locals
    int indInet = 0;
    for ( unsigned int iEle = 0; iEle < isegn.size(); iEle++ ) {
        int nne = nnet[ iEle ];
        for ( int ien = 0; ien < nne; ien++ ) {

            int indGlob = inet[indInet];
            // map it to local node
            Global2LocalMap_::iterator pos = global2LocalNodeMap.find( indGlob );
            ASSERT( pos != global2LocalNodeMap.end(),
                    "Cannot remap node index %d to local indices. \n ", indGlob );
            int indLoc = static_cast<int> ( pos -> second );

            // store the node
            inet[ indInet++ ] = indLoc;
        }
    }

    int numNodes    = size;
    int numDofsInt  = size;
    int spaceDim    = 3;    // TODO: what is the proper value here?
    int meshDim     = elDimMax;

    bddc_ls -> load_mesh( spaceDim, numNodes, numDofsInt, inet, nnet, nndf, isegn, isngn, isngn, xyz, element_permeability, meshDim );
}




//=============================================================================
// DESTROY WATER MH SYSTEM STRUCTURE
//=============================================================================
DarcyFlowMH_Steady::~DarcyFlowMH_Steady() {
    if (schur0 != NULL) delete schur0;

	delete edge_ds;
//	delete el_ds;
	delete side_ds;

//	xfree(el_4_loc);
	xfree(row_4_el);
	xfree(side_id_4_loc);
	xfree(side_row_4_id);
	xfree(edge_4_loc);
	xfree(row_4_edge);

	if (solution != NULL) {
		VecDestroy(&sol_vec);
		xfree(solution);
	}

	delete output_object;

	VecScatterDestroy(&par_to_all);
    
}


// ================================================
// PARALLLEL PART
//

// ========================================================================
// to finish row_4_id arrays we have to convert individual numberings of
// sides/els/edges to whole numbering of rows. To this end we count shifts
// for sides/els/edges on each proc and then we apply them on row_4_id
// arrays.
// we employ macros to avoid code redundancy
// =======================================================================
void DarcyFlowMH_Steady::make_row_numberings() {
    int i, shift;
    int np = edge_ds->np();
    int edge_shift[np], el_shift[np], side_shift[np];
    unsigned int rows_starts[np];

    int edge_n_id = mesh_->n_edges(),
            el_n_id = mesh_->element.size(),
            side_n_id = mesh_->n_sides();

    // compute shifts on every proc
    shift = 0; // in this var. we count new starts of arrays chunks
    for (i = 0; i < np; i++) {
        side_shift[i] = shift - (side_ds->begin(i)); // subtract actual start of the chunk
        shift += side_ds->lsize(i);
        el_shift[i] = shift - (el_ds->begin(i));
        shift += el_ds->lsize(i);
        edge_shift[i] = shift - (edge_ds->begin(i));
        shift += edge_ds->lsize(i);
        rows_starts[i] = shift;
    }
    // apply shifts
    for (i = 0; i < side_n_id; i++) {
        int &what = side_row_4_id[i];
        if (what >= 0)
            what += side_shift[side_ds->get_proc(what)];
    }
    for (i = 0; i < el_n_id; i++) {
        int &what = row_4_el[i];
        if (what >= 0)
            what += el_shift[el_ds->get_proc(what)];

    }
    for (i = 0; i < edge_n_id; i++) {
        int &what = row_4_edge[i];
        if (what >= 0)
            what += edge_shift[edge_ds->get_proc(what)];
    }
    // make distribution of rows
    for (i = np - 1; i > 0; i--)
        rows_starts[i] -= rows_starts[i - 1];

    rows_ds = boost::make_shared<Distribution>(&(rows_starts[0]), PETSC_COMM_WORLD);
}

void DarcyFlowMH_Steady::make_serial_scatter() {
    START_TIMER("prepare scatter");
    // prepare Scatter form parallel to sequantial in original numbering
    {
            IS is_loc;
            int i, *loc_idx; //, si;

            // create local solution vector
            solution = (double *) xmalloc(size * sizeof(double));
            VecCreateSeqWithArray(PETSC_COMM_SELF,1, size, solution,
                    &(sol_vec));

            // create seq. IS to scatter par solutin to seq. vec. in original order
            // use essentialy row_4_id arrays
            loc_idx = (int *) xmalloc(size * sizeof(int));
            i = 0;
            FOR_ELEMENTS(mesh_, ele) {
                FOR_ELEMENT_SIDES(ele,si) {
                    loc_idx[i++] = side_row_4_id[ mh_dh.side_dof( ele->side(si) ) ];
                }
            }
            FOR_ELEMENTS(mesh_, ele) {
                loc_idx[i++] = row_4_el[ele.index()];
            }
            for(unsigned int i_edg=0; i_edg < mesh_->n_edges(); i_edg++) {
                loc_idx[i++] = row_4_edge[i_edg];
            }
            ASSERT( i==size,"Size of array does not match number of fills.\n");
            //DBGPRINT_INT("loc_idx",size,loc_idx);
            ISCreateGeneral(PETSC_COMM_SELF, size, loc_idx, PETSC_COPY_VALUES, &(is_loc));
            xfree(loc_idx);
            VecScatterCreate(schur0->get_solution(), is_loc, sol_vec,
                    PETSC_NULL, &par_to_all);
            ISDestroy(&(is_loc));
    }
    solution_changed_for_scatter=true;

    END_TIMER("prepare scatter");

}

// ====================================================================================
// - compute optimal edge partitioning
// - compute appropriate partitioning of elements and sides
// - make arrays: *_id_4_loc and *_row_4_id to allow parallel assembly of the MH matrix
// ====================================================================================
void DarcyFlowMH_Steady::prepare_parallel( const Input::AbstractRecord in_rec) {
    
    START_TIMER("prepare parallel");
    
    int *loc_part; // optimal (edge,el) partitioning (local chunk)
    int *id_4_old; // map from old idx to ids (edge,el)
    int loc_i;

    int e_idx;

    
    //ierr = MPI_Barrier(PETSC_COMM_WORLD);
    //ASSERT(ierr == 0, "Error in MPI_Barrier.");

    // row_4_el will be modified so we make a copy of the array from mesh
    row_4_el = new int[mesh_->n_elements()];
    std::copy(mesh_->get_row_4_el(), mesh_->get_row_4_el()+mesh_->n_elements(), row_4_el);
    el_4_loc = mesh_->get_el_4_loc();
    el_ds = mesh_->get_el_ds();

        //optimal element part; loc. els. id-> new el. numbering
        Distribution init_edge_ds(DistributionLocalized(), mesh_->n_edges(), PETSC_COMM_WORLD);
        // partitioning of edges, edge belongs to the proc of his first element
        // this is not optimal but simple
        loc_part = new int[init_edge_ds.lsize()];
        id_4_old = new int[mesh_->n_edges()];
        {
            loc_i = 0;
            FOR_EDGES(mesh_, edg ) {
                unsigned int i_edg = edg - mesh_->edges.begin();
                // partition
                e_idx = mesh_->element.index(edg->side(0)->element());
                if (init_edge_ds.is_local(i_edg)) {
                    // find (new) proc of the first element of the edge
                    loc_part[loc_i++] = el_ds->get_proc(row_4_el[e_idx]);
                }
                // id array
                id_4_old[i_edg] = i_edg;
            }
        }

        Partitioning::id_maps(mesh_->n_edges(), id_4_old, init_edge_ds, loc_part, edge_ds, edge_4_loc, row_4_edge);
        delete[] loc_part;
        delete[] id_4_old;

    //optimal side part; loc. sides; id-> new side numbering
    Distribution init_side_ds(DistributionBlock(), mesh_->n_sides(), PETSC_COMM_WORLD);
    // partitioning of sides follows elements
    loc_part = new int[init_side_ds.lsize()];
    id_4_old = new int[mesh_->n_sides()];
    {
        int is = 0;
        loc_i = 0;
        FOR_SIDES(mesh_, side ) {
            // partition
            if (init_side_ds.is_local(is)) {
                // find (new) proc of the element of the side
                loc_part[loc_i++] = el_ds->get_proc(
                        row_4_el[mesh_->element.index(side->element())]);
            }
            // id array
            id_4_old[is++] = mh_dh.side_dof( side );
        }
    }

    Partitioning::id_maps(mesh_->n_sides(), id_4_old, init_side_ds, loc_part, side_ds,
            side_id_4_loc, side_row_4_id);
    delete [] loc_part;
    delete [] id_4_old;

    // convert row_4_id arrays from separate numberings to global numbering of rows
    make_row_numberings();

    // prepare global_row_4_sub_row

#ifdef FLOW123D_HAVE_BDDCML
    if (in_rec.type() ==  LinSys_BDDC::get_input_type()) {
        // auxiliary
        Element *el;
        int side_row, edge_row;

        global_row_4_sub_row = boost::make_shared<LocalToGlobalMap>(rows_ds);

        //
        // ordering of dofs
        // for each subdomain:
        // | velocities (at sides) | pressures (at elements) | L. mult. (at edges) |
        for (unsigned int i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {
            el = mesh_->element(el_4_loc[i_loc]);
            int el_row = row_4_el[el_4_loc[i_loc]];

            global_row_4_sub_row->insert( el_row );

            unsigned int nsides = el->n_sides();
            for (unsigned int i = 0; i < nsides; i++) {
                side_row = side_row_4_id[ mh_dh.side_dof( el->side(i) ) ];
		        edge_row = row_4_edge[el->side(i)->edge_idx()];

		        global_row_4_sub_row->insert( side_row );
		        global_row_4_sub_row->insert( edge_row );
            }

            for (unsigned int i_neigh = 0; i_neigh < el->n_neighs_vb; i_neigh++) {
                // mark this edge
                edge_row = row_4_edge[el->neigh_vb[i_neigh]->edge_idx() ];
                global_row_4_sub_row->insert( edge_row );
            }
        }
        global_row_4_sub_row->finalize();
    }
#endif // FLOW123D_HAVE_BDDCML

}



void mat_count_off_proc_values(Mat m, Vec v) {
    int n, first, last;
    const PetscInt *cols;
    Distribution distr(v);

    int n_off = 0;
    int n_on = 0;
    int n_off_rows = 0;
    MatGetOwnershipRange(m, &first, &last);
    for (int row = first; row < last; row++) {
        MatGetRow(m, row, &n, &cols, PETSC_NULL);
        bool exists_off = false;
        for (int i = 0; i < n; i++)
            if (distr.get_proc(cols[i]) != distr.myp())
                n_off++, exists_off = true;
            else
                n_on++;
        if (exists_off)
            n_off_rows++;
        MatRestoreRow(m, row, &n, &cols, PETSC_NULL);
    }
}












// ========================
// unsteady

DarcyFlowMH_Unsteady::DarcyFlowMH_Unsteady(Mesh &mesh_in, const Input::Record in_rec)
    : DarcyFlowMH_Steady(mesh_in, in_rec, false)
{
    time_ = new TimeGovernor(in_rec.val<Input::Record>("time"));
	data_.mark_input_times(this->mark_type());
	data_.set_time(time_->step(), LimitSide::right);

	output_object = new DarcyFlowMHOutput(this, in_rec.val<Input::Record>("output"));
	//balance_->units(output_object->get_output_fields().field_ele_pressure.units()*data_.cross_section.units()*data_.storativity.units());

	//time_->fix_dt_until_mark();
	create_linear_system();

	VecDuplicate(schur0->get_solution(), &previous_solution);
    VecCreateMPI(PETSC_COMM_WORLD,rows_ds->lsize(),PETSC_DETERMINE,&(steady_diagonal));
    VecDuplicate(steady_diagonal,& new_diagonal);
    VecZeroEntries(new_diagonal);
    VecDuplicate(*( schur0->get_rhs()), &steady_rhs);

    assembly_linear_system();
	read_init_condition();

    output_data();
}

void DarcyFlowMH_Unsteady::read_init_condition()
{

	// read inital condition
	VecZeroEntries(schur0->get_solution());

	double *local_sol = schur0->get_solution_array();

	// cycle over local element rows
	ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);

	DBGMSG("Setup with dt: %f\n",time_->dt());
	for (unsigned int i_loc_el = 0; i_loc_el < el_ds->lsize(); i_loc_el++) {
		ele = mesh_->element(el_4_loc[i_loc_el]);
		int i_loc_row = i_loc_el + side_ds->lsize();

		// set initial condition
		local_sol[i_loc_row] = data_.init_pressure.value(ele->centre(),ele->element_accessor());
	}

	solution_changed_for_scatter=true;

}

void DarcyFlowMH_Unsteady::setup_time_term() {
    // save diagonal of steady matrix
    MatGetDiagonal(*( schur0->get_matrix() ), steady_diagonal);
    // save RHS
    VecCopy(*( schur0->get_rhs()), steady_rhs);


    PetscScalar *local_diagonal;
    VecGetArray(new_diagonal,& local_diagonal);

    ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);
    DBGMSG("Setup with dt: %f\n",time_->dt());

    if (balance_ != nullptr)
    	balance_->start_mass_assembly(water_balance_idx_);

    for (unsigned int i_loc_el = 0; i_loc_el < el_ds->lsize(); i_loc_el++) {
        ele = mesh_->element(el_4_loc[i_loc_el]);
        int i_loc_row = i_loc_el + side_ds->lsize();

        // set new diagonal
        double diagonal_coeff = data_.cross_section.value(ele->centre(), ele->element_accessor())
        		* data_.storativity.value(ele->centre(), ele->element_accessor())
				* ele->measure();
        local_diagonal[i_loc_row]= - diagonal_coeff / time_->dt();

        if (balance_ != nullptr)
        	balance_->add_mass_matrix_values(water_balance_idx_, ele->region().bulk_idx(), {i_loc_row}, {diagonal_coeff});
    }
    VecRestoreArray(new_diagonal,& local_diagonal);
    MatDiagonalSet(*( schur0->get_matrix() ), new_diagonal, ADD_VALUES);

    solution_changed_for_scatter=true;
    schur0->set_matrix_changed();

    if (balance_ != nullptr)
    	balance_->finish_mass_assembly(water_balance_idx_);
}

void DarcyFlowMH_Unsteady::modify_system() {
	START_TIMER("modify system");
	if (time_->is_changed_dt() && time_->step().index()>0) {
        double scale_factor=time_->step(-2).length()/time_->step().length();
        if (scale_factor != 1.0) {
            // if time step has changed and setup_time_term not called
            MatDiagonalSet(*( schur0->get_matrix() ),steady_diagonal, INSERT_VALUES);

            VecScale(new_diagonal, time_->last_dt()/time_->dt());
            MatDiagonalSet(*( schur0->get_matrix() ),new_diagonal, ADD_VALUES);
            schur0->set_matrix_changed();
        }
	}

    // modify RHS - add previous solution
    VecPointwiseMult(*( schur0->get_rhs()), new_diagonal, schur0->get_solution());
    VecAXPY(*( schur0->get_rhs()), 1.0, steady_rhs);
    schur0->set_rhs_changed();

    // swap solutions
    VecSwap(previous_solution, schur0->get_solution());
}


//-----------------------------------------------------------------------------
// vim: set cindent:
