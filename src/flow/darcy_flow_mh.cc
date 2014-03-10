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
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 *
 * @file
 * @ingroup flow
 * @brief  Setup and solve linear system of mixed-hybrid discretization of the linear
 * porous media flow with possible preferential flow in fractures and chanels.
 *
 */

#include "petscmat.h"
#include "petscviewer.h"
#include "petscerror.h"
#include <armadillo>
#include <boost/foreach.hpp>


#include "system/system.hh"

#include "mesh/mesh.h"
#include "mesh/intersection.hh"
#include "mesh/partitioning.hh"
#include "la/distribution.hh"
#include "la/linsys.hh"
#include "la/linsys_PETSC.hh"
#include "la/linsys_BDDC.hh"
//#include "la/solve.h"
#include "la/schur.hh"
#include "la/sparse_graph.hh"
#include "la/local_to_global_map.hh"

#include "system/file_path.hh"
#include "flow/mh_fe_values.hh"
#include "flow/darcy_flow_mh.hh"

#include "flow/darcy_flow_mh_output.hh"

#include "fem/mapping_p1.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "quadrature/quadrature_lib.hh"

#include <limits>
#include <set>
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>

#include "coupling/time_governor.hh"

#include "fields/field_base.hh"
#include "fields/field_values.hh"
#include "system/sys_profiler.hh"




namespace it = Input::Type;

it::Selection DarcyFlowMH::mh_mortar_selection
	= it::Selection("MH_MortarMethod")
	.add_value(NoMortar, "None", "Mortar space: P0 on elements of lower dimension.")
	.add_value(MortarP0, "P0", "Mortar space: P0 on elements of lower dimension.")
	.add_value(MortarP1, "P1", "Mortar space: P1 on intersections, using non-conforming pressures.");


it::Selection DarcyFlowMH::EqData::bc_type_selection =
              it::Selection("EqData_bc_Type")
               .add_value(none, "none", "Homogeneous Neoumann BC.")
               .add_value(dirichlet, "dirichlet")
               .add_value(neumann, "neumann")
               .add_value(robin, "robin")
               .add_value(total_flux, "total_flux");

//new input type with FIELDS
it::AbstractRecord DarcyFlowMH::input_type=
        it::AbstractRecord("DarcyFlowMH", "Mixed-Hybrid  solver for saturated Darcy flow.")
        .declare_key("n_schurs", it::Integer(0,2), it::Default("2"),
                "Number of Schur complements to perform when solving MH sytem.")
        .declare_key("solver", LinSys::input_type, it::Default::obligatory(),
                "Linear solver for MH problem.")
        .declare_key("output", DarcyFlowMHOutput::input_type, it::Default::obligatory(),
                "Parameters of output form MH module.")
        .declare_key("mortar_method", mh_mortar_selection, it::Default("None"),
                "Method for coupling Darcy flow between dimensions." );


it::Record DarcyFlowMH_Steady::input_type
    = it::Record("Steady_MH", "Mixed-Hybrid  solver for STEADY saturated Darcy flow.")
    .derive_from(DarcyFlowMH::input_type)
    .declare_key("bc_data", it::Array(
                DarcyFlowMH_Steady::EqData().boundary_input_type()
                .declare_key("bc_piezo_head", FieldBase< 3, FieldValue<3>::Scalar >::get_input_type(), "Boundary condition for pressure as piezometric head." )
                .declare_key("flow_old_bcd_file", it::FileName::input(), "")
                ), it::Default::obligatory(), ""  )
    .declare_key("bulk_data", it::Array(
                DarcyFlowMH_Steady::EqData().bulk_input_type() 
                .declare_key("init_piezo_head", FieldBase< 3, FieldValue<3>::Scalar >::get_input_type(), "Initial condition for pressure as piezometric head." )
                ), it::Default::obligatory(), "");


it::Record DarcyFlowMH_Unsteady::input_type
	= it::Record("Unsteady_MH", "Mixed-Hybrid solver for unsteady saturated Darcy flow.")
	.derive_from(DarcyFlowMH::input_type)
	.declare_key("time", TimeGovernor::input_type, it::Default::obligatory(),
                 "Time governor setting for the unsteady Darcy flow model.")
  .declare_key("bc_data", it::Array(
                DarcyFlowMH_Unsteady::EqData().boundary_input_type()
                .declare_key("bc_piezo_head", FieldBase< 3, FieldValue<3>::Scalar >::get_input_type(), "Boundary condition for piezometric head." )
                .declare_key("flow_old_bcd_file", it::FileName::input(), "")
                ), it::Default::obligatory(), ""  )
  .declare_key("bulk_data", it::Array(
                DarcyFlowMH_Unsteady::EqData().bulk_input_type()
                .declare_key("init_piezo_head", FieldBase< 3, FieldValue<3>::Scalar >::get_input_type(), "Initial piezometric head." )
                ), it::Default::obligatory(), "");


it::Record DarcyFlowLMH_Unsteady::input_type
    = it::Record("Unsteady_LMH", "Lumped Mixed-Hybrid solver for unsteady saturated Darcy flow.")
    .derive_from(DarcyFlowMH::input_type)
    .declare_key("time",         TimeGovernor::input_type, it::Default::obligatory(),
                                "Time governor setting for the unsteady Darcy flow model.")
    .declare_key("bc_data", it::Array(
                DarcyFlowLMH_Unsteady::EqData().boundary_input_type()
                .declare_key("bc_piezo_head", FieldBase< 3, FieldValue<3>::Scalar >::get_input_type(), "Boundary condition for piezometric head." )
                .declare_key("flow_old_bcd_file", it::FileName::input(), "")
                ), it::Default::obligatory(), ""  )
    .declare_key("bulk_data", it::Array(
                DarcyFlowLMH_Unsteady::EqData().bulk_input_type()
                .declare_key("init_piezo_head", FieldBase< 3, FieldValue<3>::Scalar >::get_input_type(), "Initial piezometric head." )
                ), it::Default::obligatory(), "");
    







DarcyFlowMH::EqData::EqData(const std::string &name)
: EqDataBase(name)
{
    gravity_ = arma::vec4("0 0 -1 0"); // gravity vector + constant shift of values

    ADD_FIELD(anisotropy, "Anisotropy of the conductivity tensor.", Input::Type::Default("1.0"));
    ADD_FIELD(cross_section, "Complement dimension parameter (cross section for 1D, thickness for 2D).", Input::Type::Default("1.0") );
    ADD_FIELD(conductivity, "Isotropic conductivity scalar.", Input::Type::Default("1.0") );
    ADD_FIELD(sigma, "Transition coefficient between dimensions.", Input::Type::Default("1.0"));
    ADD_FIELD(water_source_density, "Water source density.", Input::Type::Default("0.0"));
    
    ADD_FIELD(bc_type,"Boundary condition type, possible values:", it::Default("none") );
              bc_type.set_selection(&bc_type_selection);

    ADD_FIELD(bc_pressure,"Dirichlet BC condition value for pressure.");
    std::vector<FieldEnum> list; list.push_back(none); list.push_back(neumann);
    bc_pressure.disable_where(& bc_type, list );

    ADD_FIELD(bc_flux,"Flux in Neumman or Robin boundary condition.");
    list.clear(); list.push_back(none); list.push_back(dirichlet); list.push_back(robin);
    bc_flux.disable_where(& bc_type, list );

    ADD_FIELD(bc_robin_sigma,"Conductivity coefficient in Robin boundary condition.");
    list.clear(); list.push_back(none); list.push_back(dirichlet); list.push_back(neumann);
    bc_robin_sigma.disable_where(& bc_type, list );
    
    //these are for unsteady
    ADD_FIELD(init_pressure, "Initial condition as pressure", it::Default("0.0") );
    ADD_FIELD(storativity,"Storativity.", it::Default("1.0") );
}



RegionSet DarcyFlowMH::EqData::read_boundary_list_item(Input::Record rec) {
    RegionSet domain=EqDataBase::read_boundary_list_item(rec);
    Input::AbstractRecord field_a_rec;
    if (rec.opt_val("bc_piezo_head", field_a_rec))
                bc_pressure.set_field(domain, boost::make_shared< FieldAddPotential<3, FieldValue<3>::Scalar > >( this->gravity_, field_a_rec) );
    FilePath flow_bcd_file;
    if (rec.opt_val("flow_old_bcd_file", flow_bcd_file) ) {
        OldBcdInput::instance()->read_flow(flow_bcd_file, bc_type, bc_pressure, bc_flux, bc_robin_sigma);
    }
    return domain;
}

RegionSet DarcyFlowMH::EqData::read_bulk_list_item(Input::Record rec) {
    RegionSet domain=EqDataBase::read_bulk_list_item(rec);
    Input::AbstractRecord field_a_rec;
    if (rec.opt_val("init_piezo_head", field_a_rec)) {
                init_pressure.set_field(domain, boost::make_shared< FieldAddPotential<3, FieldValue<3>::Scalar > >( this->gravity_, field_a_rec) );
    }
    return domain;
}

/* TODO: this can be applied when Unstedy is no longer descendant from Steady
DarcyFlowMH_Unsteady::EqData::EqData()
  : DarcyFlowMH::EqData("DarcyFlowMH_Unsteady") 
{
    ADD_FIELD(init_pressure, "Initial condition as pressure");
    ADD_FIELD(storativity,"Storativity.");
}

DarcyFlowLMH_Unsteady::EqData::EqData()
  : DarcyFlowMH::EqData("DarcyFlowLMH_Unsteady") 
 {
    ADD_FIELD(init_pressure, "Initial condition as pressure");
    ADD_FIELD(storativity,"Storativity.");
}
//*/













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
DarcyFlowMH_Steady::DarcyFlowMH_Steady(Mesh &mesh_in, const Input::Record in_rec)
: DarcyFlowMH(mesh_in, in_rec)

{
    using namespace Input;
    F_ENTRY;
    START_TIMER("Darcy constructor");

    //connecting data fields with mesh
    START_TIMER("data init");
    data.set_mesh(&mesh_in);
    data.init_from_input( in_rec.val<Input::Array>("bulk_data"), in_rec.val<Input::Array>("bc_data") );
        
    // steady time governor
    time_ = new TimeGovernor();
    
    //initializing data fields at the beginning (time = 0)
    data.set_time(*time_);
    END_TIMER("data init");
    
    int ierr;

    size = mesh_->n_elements() + mesh_->n_sides() + mesh_->n_edges();
    n_schur_compls = in_rec.val<int>("n_schurs");
    

    solution = NULL;
    schur0   = NULL;
    schur1   = NULL;
    schur2   = NULL;

    
    
    // init paralel structures
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &(myp));
    ierr += MPI_Comm_size(PETSC_COMM_WORLD, &(np));
    if (ierr)
        xprintf(Err, "Some error in MPI.\n");


    
    mortar_method_= in_rec.val<MortarMethod>("mortar_method");
    if (mortar_method_ != NoMortar) {
        mesh_->make_intersec_elements();
    }


    mh_dh.reinit(mesh_);

    prepare_parallel(in_rec.val<AbstractRecord>("solver"));

    //side_ds->view( std::cout );
    //el_ds->view( std::cout );
    //edge_ds->view( std::cout );
    //rows_ds->view( std::cout );
    
    make_schurs(in_rec.val<AbstractRecord>("solver"));

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




//=============================================================================
// COMPOSE and SOLVE WATER MH System possibly through Schur complements
//=============================================================================
void DarcyFlowMH_Steady::update_solution() {
    START_TIMER("Solving MH system");
    F_ENTRY;

    if (time_->is_end()) return;

    time_->next_time();
    
    START_TIMER("data reinit");
    //reinitializing data fields after time step
    //TODO: workaround for the steady problem
    //if (time_->t() != TimeGovernor::inf_time) //this test cannot be here due to (mainly implicit) transport - the fields are not neccesary (or cannot) to be read again but the time must be set to infinity
    //the problem of time==infinity shows up in field_elementwise and field_interpolatedP0 where a gmsh file is read and there is no such data at infinity
    //temporarily solved directly in field_elementwise and field_interpolatedP0
    data.set_time(*time_);
    END_TIMER("data reinit");

    //xprintf(Msg, "DARCY:  t: %f  dt: %f\n",time_->t(), time_->dt());
    time_->view("DARCY"); //time governor information output
    
    modify_system(); // hack for unsteady model
    int convergedReason = schur0->solve();

    xprintf(MsgLog, "Linear solver ended with reason: %d \n", convergedReason );
    ASSERT( convergedReason >= 0, "Linear solver failed to converge. Convergence reason %d \n", convergedReason );

    this -> postprocess();

    //int rank;
    //MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
    //if ( rank == 0 ) {
    //    PetscViewer solViewer;
    //    PetscViewerASCIIOpen( PETSC_COMM_SELF, "sol.m", &solViewer );
    //    PetscViewerSetFormat(solViewer,PETSC_VIEWER_ASCII_MATLAB);
    //    VecView( sol_vec, solViewer );
    //    PetscViewerDestroy(solViewer);
    //}
    solution_changed_for_scatter=true;
    DBGMSG("solution updated\n");
}

void DarcyFlowMH_Steady::postprocess() 
{
    START_TIMER("postprocess");
    int side_rows[4];
    double values[4];
    ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);

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


double DarcyFlowMH_Steady::solution_precision() const
{
	return schur0->get_solution_precision();
}


void  DarcyFlowMH_Steady::get_solution_vector(double * &vec, unsigned int &vec_size)
{
    // TODO: make class for vectors (wrapper for PETSC or other) derived from LazyDependency
    // and use its mechanism to manage dependency between vectors
    if (solution_changed_for_scatter) {
    //    // scatter solution to all procs
    //    VecScatterBegin(par_to_all, schur0->get_solution(), sol_vec,
    //            INSERT_VALUES, SCATTER_FORWARD);
    //    VecScatterEnd(par_to_all, schur0->get_solution(), sol_vec,
    //            INSERT_VALUES, SCATTER_FORWARD);
    //    solution_changed_for_scatter=false;

        std::vector<double> sol_disordered(this->size);
        schur0 -> get_whole_solution( sol_disordered );

        // reorder solution to application ordering
        if ( solution_.empty() ) solution_.resize( this->size, 0. );
        for ( int i = 0; i < this->size; i++ ) {
            solution_[i] = sol_disordered[solver_indices_[i]];
        }
    }

    vec=&(solution_[0]);
    vec_size = solution_.size();
    ASSERT(vec != NULL, "Requested solution is not allocated!\n");
}

void  DarcyFlowMH_Steady::get_partitioning_vector(int * &elem_part, unsigned &lelem_part)
{
  
//    elem_part=&(element_part[0]);
//    lelem_part = element_part.size();
//    ASSERT(elem_part != NULL, "Requested vector is not allocated!\n");
}

void  DarcyFlowMH_Steady::get_parallel_solution_vector(Vec &vec)
{
    vec=schur0->get_solution();
    ASSERT(vec != NULL, "Requested solution is not allocated!\n");
}



// ===========================================================================================
//
//   MATRIX ASSEMBLY - we use abstract assembly routine, where  LS Mat/Vec SetValues
//   are in fact pointers to allocating or filling functions - this is governed by Linsystem roitunes
//
// =======================================================================================


// ******************************************
// ABSTRACT ASSEMBLY OF MH matrix
// TODO: matice by se mela sestavovat zvlast pro kazdou dimenzi (objem, pukliny, pruseciky puklin)
//       konekce by se mely sestavovat cyklem pres konekce, konekce by mely byt paralelizovany podle
//       distribuce elementu nizssi dimenze
//       k tomuto je treba nejdriv spojit s JK verzi, aby se vedelo co se deje v transportu a
//       predelat mesh a neigbouring
// *****************************************
void DarcyFlowMH_Steady::assembly_steady_mh_matrix() {
    LinSys *ls = schur0;
    ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);
    MHFEValues fe_values;

    // We use FESideValues for calculating normal vectors.
    // For initialization of FESideValues some auxiliary objects are needed.
    MappingP1<1,3> map1;
    MappingP1<2,3> map2;
    MappingP1<3,3> map3;
    QGauss<0> q0(1);
    QGauss<1> q1(1);
    QGauss<2> q2(1);
    FE_P_disc<1,1,3> fe1;
    FE_P_disc<0,2,3> fe2;
    FE_P_disc<0,3,3> fe3;
    FESideValues<1,3> fe_side_values1(map1, q0, fe1, update_normal_vectors);
    FESideValues<2,3> fe_side_values2(map2, q1, fe2, update_normal_vectors);
    FESideValues<3,3> fe_side_values3(map3, q2, fe3, update_normal_vectors);

    class Boundary *bcd;
    class Neighbour *ngh;

    bool fill_matrix = schur0->is_preallocated();
    DBGMSG("fill_matrix: %d\n", fill_matrix);
    int el_row, side_row, edge_row;
    int tmp_rows[100];
    //int  nsides;
    int side_rows[4], edge_rows[4]; // rows for sides and edges of one element
    double local_vb[4]; // 2x2 matrix
    double zeros[1000]; // to make space for second schur complement, max. 10 neighbour edges of one el.
    double minus_ones[4] = { -1.0, -1.0, -1.0, -1.0 };
    double loc_side_rhs[4];
    std::map<int,double> subdomain_diagonal_map;
    F_ENTRY;

    //DBGPRINT_INT("side_row_4_id",mesh->max_side_id+1,side_row_4_id);
    //DBGPRINT_INT("el_row_4_id",mesh->max_elm_id+1,el_row_4_id);
    //DBGPRINT_INT("edge_row_4_id",mesh->max_edg_id+1,edge_row_4_id);
    //DBGPRINT_INT("el_id_4_loc",el_ds->lsize(),el_id_4_loc);

    SET_ARRAY_ZERO(zeros,1000);
    for (unsigned int i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {

        ele = mesh_->element(el_4_loc[i_loc]);
        el_row = row_4_el[el_4_loc[i_loc]];
        unsigned int nsides = ele->n_sides();
        if (fill_matrix) fe_values.update(ele, data.anisotropy, data.cross_section, data.conductivity);
        double cross_section = data.cross_section.value(ele->centre(), ele->element_accessor());

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
                EqData::BC_Type type = (EqData::BC_Type)data.bc_type.value(b_ele.centre(), b_ele);
                //DBGMSG("BCD id: %d sidx: %d type: %d\n", ele->id(), i, type);
                if ( type == EqData::none) {
                    // homogeneous neumann
                } else if ( type == EqData::dirichlet ) {
                    c_val = 0.0;
                    double bc_pressure = data.bc_pressure.value(b_ele.centre(), b_ele);
                    loc_side_rhs[i] -= bc_pressure;
                    ls->rhs_set_value(edge_row, -bc_pressure);
                    ls->mat_set_value(edge_row, edge_row, -1.0);

                } else if ( type == EqData::neumann) {
                    double bc_flux = data.bc_flux.value(b_ele.centre(), b_ele);
                    ls->rhs_set_value(edge_row, bc_flux * bcd->element()->measure() * cross_section);
                    //DBGMSG("neumann edge_row, ele_index,el_idx: %d \t %d \t %d\n", edge_row, ele->index(), ele->side(i)->el_idx());

                } else if ( type == EqData::robin) {
                    double bc_pressure = data.bc_pressure.value(b_ele.centre(), b_ele);
                    double bc_sigma = data.bc_robin_sigma.value(b_ele.centre(), b_ele);
                    ls->rhs_set_value(edge_row, -bcd->element()->measure() * bc_sigma * bc_pressure * cross_section );
                    ls->mat_set_value(edge_row, edge_row, -bcd->element()->measure() * bc_sigma * cross_section );

                } else {
                    xprintf(UsrErr, "BC type not supported.\n");
                }
            }
            ls->mat_set_value(side_row, edge_row, c_val);
            ls->mat_set_value(edge_row, side_row, c_val);

            // update matrix for weights in BDDCML
            if ( typeid(*ls) == typeid(LinSys_BDDC) ) {
               double val_side =  (fe_values.local_matrix())[i*nsides+i];
               double val_edge =  -1./ (fe_values.local_matrix())[i*nsides+i];
               subdomain_diagonal_map.insert( std::make_pair( side_row, val_side ) );
               subdomain_diagonal_map.insert( std::make_pair( edge_row, val_edge ) );
            }
        }

        ls->rhs_set_values(nsides, side_rows, loc_side_rhs);


        // set block A: side-side on one element - block diagonal matrix

        //std::cout << "subMat in flow" << std::endl;
        //for ( unsigned i = 0; i < nsides; i++) {
        //    for ( unsigned j = 0; j < nsides; j++) {
        //        std::cout << fe_values.local_matrix()[i*nsides+j]  << " ";
        //    }
        //    std::cout << std::endl;
        //}

        ls->mat_set_values(nsides, side_rows, nsides, side_rows, fe_values.local_matrix() );
        // set block B, B': element-side, side-element
        ls->mat_set_values(1, &el_row, nsides, side_rows, minus_ones);
        ls->mat_set_values(nsides, side_rows, 1, &el_row, minus_ones);


        // set sources
        ls->rhs_set_value(el_row, -1.0 * ele->measure() *
                          data.cross_section.value(ele->centre(), ele->element_accessor()) * 
                          data.water_source_density.value(ele->centre(), ele->element_accessor()) );
        
        
        // D block: non-compatible conections and diagonal: element-element

        ls->mat_set_value(el_row, el_row, 0.0);         // maybe this should be in virtual block for schur preallocation

        // D, E',E block: compatible connections: element-edge
        
        for (unsigned int i = 0; i < ele->n_neighs_vb; i++) {
            // every compatible connection adds a 2x2 matrix involving
            // current element pressure  and a connected edge pressure
            ngh= ele->neigh_vb[i];
            tmp_rows[0]=el_row;
            tmp_rows[1]=row_4_edge[ ngh->edge_idx() ];

            // compute normal vector to side
            arma::vec3 nv;
            ElementFullIter ele_higher = mesh_->element.full_iter(ngh->side()->element());
            switch (ele_higher->dim()) {
            case 1:
            	fe_side_values1.reinit(ele_higher, ngh->side()->el_idx());
            	nv = fe_side_values1.normal_vector(0);
            	break;
            case 2:
            	fe_side_values2.reinit(ele_higher, ngh->side()->el_idx());
            	nv = fe_side_values2.normal_vector(0);
            	break;
            case 3:
            	fe_side_values3.reinit(ele_higher, ngh->side()->el_idx());
            	nv = fe_side_values3.normal_vector(0);
            	break;
            }

            double value = data.sigma.value( ele->centre(), ele->element_accessor()) *
            		2*data.conductivity.value( ele->centre(), ele->element_accessor()) *
            		arma::dot(data.anisotropy.value( ele->centre(), ele->element_accessor())*nv, nv) *
                    data.cross_section.value( ngh->side()->centre(), ele_higher->element_accessor() ) * // cross-section of higher dim. (2d)
                    data.cross_section.value( ngh->side()->centre(), ele_higher->element_accessor() ) /
                    data.cross_section.value( ele->centre(), ele->element_accessor() ) *      // crossection of lower dim.
                    ngh->side()->measure();


            local_vb[0] = -value;   local_vb[1] = value;
            local_vb[2] = value;    local_vb[3] = -value;

            ls->mat_set_values(2, tmp_rows, 2, tmp_rows, local_vb);

            // update matrix for weights in BDDCML
            if ( typeid(*ls) == typeid(LinSys_BDDC) ) {
               int ind = tmp_rows[1];
               double new_val = value;
               std::map<int,double>::iterator it = subdomain_diagonal_map.find( ind );
               if ( it != subdomain_diagonal_map.end() ) {
                  new_val = new_val + it->second;
               }
               subdomain_diagonal_map.insert( std::make_pair( ind, new_val ) );
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

            // ls->mat_set_values(nsides, edge_rows, ele->e_row_count, tmp_rows, zeros);
            // ls->mat_set_values(ele->e_row_count, tmp_rows, nsides, edge_rows, zeros);
        case 1: // included also for case 2
            // -(C')*(A-)*B block and its transpose conect edge with its elements
            ls->mat_set_values(1, &el_row, nsides, edge_rows, zeros);
            ls->mat_set_values(nsides, edge_rows, 1, &el_row, zeros);
            // -(C')*(A-)*C block conect all edges of every element
            ls->mat_set_values(nsides, edge_rows, nsides, edge_rows, zeros);
        }
    }

    if ( typeid(*ls) == typeid(LinSys_BDDC) ) {
       ls->load_diagonal( subdomain_diagonal_map );
    }

    //if (! mtx->ins_mod == ALLOCATE ) {
    //    MatAssemblyBegin(mtx->A,MAT_FINAL_ASSEMBLY);
    //    MatAssemblyEnd(mtx->A,MAT_FINAL_ASSEMBLY);
    // }
    // set block F - diagonal: edge-edge from Newton BC
    // also Dirichlet BC
    /*
    for (i_loc = 0; i_loc < edge_ds->lsize(); i_loc++) {
        edge_row = row_4_edge[edge_4_loc[i_loc]];
        EdgeFullIter edg = mesh_->edge(edge_4_loc[i_loc]);

        //xprintf(Msg,"F: %d %f\n",old_4_new[edge_row],edg->f_val);
        //ls->mat_set_value(edge_row, edge_row, edg->f_val);
        //ls->rhs_set_value(edge_row, edg->f_rhs);
    }*/


    if (mortar_method_ == MortarP0) {
    	P0_CouplingAssembler(*this).assembly(*ls);
    } else if (mortar_method_ == MortarP1) {
        P1_CouplingAssembler(*this).assembly(*ls);
    }  
}

void P0_CouplingAssembler::pressure_diff(int i_ele,
		vector<int> &dofs, unsigned int &ele_type, double &delta, arma::vec &dirichlet) {

	const Element *ele;

	if (i_ele == ml_it_->size() ) { // master element .. 1D
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

	for(unsigned int i_side=0; i_side < ele->n_sides(); i_side++ ) {
		dofs[i_side]=darcy_.row_4_edge[ele->side(i_side)->edge_idx()];
		Boundary * bcd = ele->side(i_side)->cond();
		if (bcd) {
			dirichlet.resize(ele->n_sides());
			ElementAccessor<3> b_ele = bcd->element_accessor();
			DarcyFlowMH::EqData::BC_Type type = (DarcyFlowMH::EqData::BC_Type)darcy_.data.bc_type.value(b_ele.centre(), b_ele);
			//DBGMSG("bcd id: %d sidx: %d type: %d\n", ele->id(), i_side, type);
			if (type == DarcyFlowMH::EqData::dirichlet) {
				//DBGMSG("Dirichlet: %d\n", ele->index());
				dofs[i_side] = -dofs[i_side];
				double bc_pressure = darcy_.data.bc_pressure.value(b_ele.centre(), b_ele);
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

		double master_sigma=darcy_.data.sigma.value( master_->centre(), master_->element_accessor());

		// rows
		for(i = 0; i <= ml_it_->size(); ++i) {
			//cout << "Intersection:" << master_->index() << ", "
			//		<<  intersections_[ (*ml_it_)[i] ].slave_iter()->index();
			pressure_diff(i, dofs_i, ele_type_i, delta_i, dirichlet_i);
			//columns
			for (j = 0; j <= ml_it_->size(); ++j) {
				pressure_diff(j, dofs_j, ele_type_j, delta_j, dirichlet_j);

				double scale =  -master_sigma * delta_i * delta_j / delta_0;
				product = scale * tensor_average[ele_type_i][ele_type_j];


				//cout << product << endl;


				arma::vec rhs(dofs_i.size());
				rhs.zeros();
				ls.set_values( dofs_i, dofs_j, product, rhs, dirichlet_i, dirichlet_j);
			}
		}
    } // loop over master elements
}



 void P1_CouplingAssembler::add_sides(const Element * ele, unsigned int shift, vector<int> &dofs, vector<double> &dirichlet)
 {

		for(unsigned int i_side=0; i_side < ele->n_sides(); i_side++ ) {
			dofs[shift+i_side] =  darcy_.row_4_edge[ele->side(i_side)->edge_idx()];
			Boundary * bcd = ele->side(i_side)->cond();

			if (bcd) {
				ElementAccessor<3> b_ele = bcd->element_accessor();
				DarcyFlowMH::EqData::BC_Type type = (DarcyFlowMH::EqData::BC_Type)darcy_.data.bc_type.value(b_ele.centre(), b_ele);
				//DBGMSG("bcd id: %d sidx: %d type: %d\n", ele->id(), i_side, type);
				if (type == DarcyFlowMH::EqData::dirichlet) {
					//DBGMSG("Dirichlet: %d\n", ele->index());
					dofs[shift + i_side] = -dofs[shift + i_side];
					double bc_pressure = darcy_.data.bc_pressure.value(b_ele.centre(), b_ele);
					dirichlet[shift + i_side] = bc_pressure;
				}
			}
		}
 }


/**
 * P1 coonection of different dimensions
 * - demonstrated convergence, but still major open questions:
 * ? in all test cases the error on the fracture is less on the left and greater on the right
 *   with incresing trend
 *   tried:
 *   - changed order of 1d elements (no change)
 *   - higher precision of ngh output and linear solver (no change)
 * ? seems that there should be some factor 6.0 in the communication term
 * ? in the case of infinite k2 -> simplest 1d-constant communication, the biggest difference on borders,
 *   numerical solution greater then analytical
 *
 * TODO:
 * * full implementation of Dirichlet BC ( also for 2d sides)
 */

void P1_CouplingAssembler::assembly(LinSys &ls) {

	for (const Intersection &intersec : intersections_) {
    	const Element * master = intersec.master_iter();
       	const Element * slave = intersec.slave_iter();

       	add_sides(master, 0, dofs, dirichlet);
       	add_sides(slave, 2, dofs, dirichlet);
		double master_sigma=darcy_.data.sigma.value( master->centre(), master->element_accessor());

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
        arma::vec point_2D_Y(intersec.map_to_slave(point_Y)); // local coordinates of  Y on slave
        arma::vec point_1D_Y(intersec.map_to_master(point_Y)); //  local coordinates of  Y on master

        arma::vec point_X(1);
        point_X.fill(0.0);
        arma::vec point_2D_X(intersec.map_to_slave(point_X)); // local coordinates of  X on slave
        arma::vec point_1D_X(intersec.map_to_master(point_X)); // local coordinates of  X on master

        arma::mat base_2D(3, 3);
        // base fce = a0 + a1*t0 + a2*t1
        //         a0     a1      a2
        base_2D << 1.0 << 0.0 << -2.0 << arma::endr //point on side 0
                << -1.0 << 2.0 << 2.0 << arma::endr // point on side 1
                << 1.0 << -2.0 << 0.0 << arma::endr;// point on side 2

        arma::mat base_1D(2, 2);
        //    base fce =   a0 + a1 * t0
        //          a0     a1
        base_1D << 1.0 << -1.0 << arma::endr // point on side 0,
                << 0.0 << 1.0 << arma::endr; // point on side 1,


        //base_1D.print("base_1D: ");
        //base_2D.print("base_2D: ");
        //point_2D_X.print("point_2D_X: ");
        //point_2D_Y.print("point_2D_Y: ");

        arma::vec difference_in_Y(5);
        arma::vec difference_in_X(5);

        // slave sides 0,1,2
        difference_in_Y.subvec(0, 2) = -base_2D * point_2D_Y;
        difference_in_X.subvec(0, 2) = -base_2D * point_2D_X;
        // master sides 3,4
        difference_in_Y.subvec(3, 4) = base_1D * point_1D_Y;
        difference_in_X.subvec(3, 4) = base_1D * point_1D_X;

        //difference_in_X.print("difference_in_X: ");
        //difference_in_Y.print("difference_in_Y: ");

        //prvky matice A[i,j]
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
        //for(int i=0;i<5;i++) DBGMSG("idx %d: %d\n",i, dofs[i]);
        //A.print("A:");


		ls.set_values( dofs, dofs, A, rhs, dirichlet, dirichlet);

    }
}



/*******************************************************************************
 * COMPOSE WATER MH MATRIX WITHOUT SCHUR COMPLEMENT
 ******************************************************************************/

void DarcyFlowMH_Steady::make_schurs( const Input::AbstractRecord in_rec) {
  
    START_TIMER("preallocation");
    int i_loc, el_row;
    Element *ele;
    Vec aux;
    PetscErrorCode err;

    F_ENTRY;

    //xprintf(Msg,"****************** problem statistics \n");
    //xprintf(Msg,"edges: %d \n",mesh_->n_edges());
    //xprintf(Msg,"sides: %d \n",mesh_->n_sides());
    //xprintf(Msg,"elements: %d \n",mesh_->n_elements());
    //xprintf(Msg,"************************************* \n");
    //xprintf(Msg,"problem size: %d \n",this->size);
    //xprintf(Msg,"****************** problem statistics \n");

    if (schur0 == NULL) { // create Linear System for MH matrix
       
    	if (in_rec.type() == LinSys_BDDC::input_type) {
    		xprintf(Warn, "For BDDC is using no Schur complements.");
            n_schur_compls = 0;
    	} else if (n_schur_compls > 2) {
            xprintf(Warn, "Invalid number of Schur Complements. Using 2.");
            n_schur_compls = 2;
        }
        if (in_rec.type() == LinSys_BDDC::input_type && rows_ds->np() > 1) {
#ifdef HAVE_BDDCML
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
#endif
        }


        // use PETSC for serial case even when user want BDDC
        if (in_rec.type() == LinSys_PETSC::input_type || schur0==NULL) {
        	LinSys_PETSC *schur1, *schur2;

        	if (n_schur_compls == 0) {
                LinSys_PETSC *ls = new LinSys_PETSC( &(*rows_ds) );

                // temporary solution; we have to set precision also for sequantial case of BDDC
                // final solution should be probably call of direct solver for oneproc case
                if (in_rec.type() != LinSys_BDDC::input_type) ls->set_from_input(in_rec);
                else {
                    ls->LinSys::set_from_input(in_rec); // get only common options
                }

                ls->set_solution( NULL );
                schur0=ls;
        	} else {
        		IS is;
        		err = ISCreateStride(PETSC_COMM_WORLD, side_ds->lsize(), rows_ds->begin(), 1, &is);
        		ASSERT(err == 0,"Error in ISCreateStride.");

        		SchurComplement *ls = new SchurComplement(is, &(*rows_ds));
        		ls->set_from_input(in_rec);
        		ls->set_solution( NULL );
        		ls->set_positive_definite();

        		// make schur1
            	Distribution *ds = ls->make_complement_distribution();
            	if (n_schur_compls==1) {
            		schur1 = new LinSys_PETSC(ds);
            	} else {
        			IS is;
        			err = ISCreateStride(PETSC_COMM_WORLD, el_ds->lsize(), ls->get_distribution()->begin(), 1, &is);
        			ASSERT(err == 0,"Error in ISCreateStride.");
        			SchurComplement *ls1 = new SchurComplement(is, ds); // is is deallocated by SchurComplement
        			ls1->set_negative_definite();

        			// make schur2
        			schur2 = new LinSys_PETSC( ls1->make_complement_distribution() );
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
            schur0->start_add_assembly(); // finish allocation and create matrix

            END_TIMER("PETSC PREALLOCATION");

            VecZeroEntries(schur0->get_solution());
        }

        if (schur0==NULL) {
            xprintf(Err, "Unknown solver type. Internal error.\n");
        }
    }


    END_TIMER("preallocation");
    
    START_TIMER("assembly");

    assembly_steady_mh_matrix(); // fill matrix
    schur0->finish_assembly();

    END_TIMER("assembly");
    //schur0->view_local_matrix();
    //PetscViewer myViewer;
    //PetscViewerASCIIOpen(PETSC_COMM_WORLD,"matis.m",&myViewer);
    //PetscViewerSetFormat(myViewer,PETSC_VIEWER_ASCII_MATLAB);

    //MatView( schur0->get_matrix(),PETSC_VIEWER_STDOUT_WORLD  );
    //DBGMSG("RHS\n");
    //VecView(schur0->get_rhs(),   PETSC_VIEWER_STDOUT_WORLD);
    //VecView(schur0->get_solution(),   PETSC_VIEWER_STDOUT_WORLD);

    //PetscViewerDestroy(myViewer);


    // add time term
}



void DarcyFlowMH_Steady::set_mesh_data_for_bddc(LinSys_BDDC * bddc_ls) {
    // prepare mesh for BDDCML
    // initialize arrays
    std::map<int,arma::vec3> localDofMap;
    std::vector<int> inet;
    std::vector<int> nnet;
    std::vector<int> isegn;

    std::vector<double> element_permeability;

    // maximal and minimal dimension of elements
    int elDimMax = 1;
    int elDimMin = 3;
    for ( int i_loc = 0; i_loc < el_ds->lsize(); i_loc++ ) {
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
            Edge *edg=el->side(si)->edge();

            // insert local edge dof
            int edge_row = row_4_edge[ el->side(si)->edge_idx() ];
            arma::vec3 coord = el->side(si)->centre();

            localDofMap.insert( std::make_pair( edge_row, coord ) );
            inet.push_back( edge_row );
            nne++;
        }

        // insert dofs related to compatible connections
        for ( int i_neigh = 0; i_neigh < el->n_neighs_vb; i_neigh++) {
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
        double conduct = data.conductivity.value( centre , el->element_accessor() );
        double cs = data.cross_section.value( centre, el->element_accessor() );
        arma::mat33 aniso = data.anisotropy.value( centre, el->element_accessor() );

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
    int numNodeSub = localDofMap.size();
    ASSERT_EQUAL( numNodeSub, global_row_4_sub_row->size() );
    std::vector<int> isngn( numNodeSub );
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

    // nndf is trivially one
    std::vector<int> nndf( numNodeSub, 1 );

    // prepare auxiliary map for renumbering nodes
    typedef std::map<int,int> Global2LocalMap_; //! type for storage of global to local map
    Global2LocalMap_ global2LocalNodeMap;
    for ( unsigned ind = 0; ind < isngn.size(); ++ind ) {
        global2LocalNodeMap.insert( std::make_pair( static_cast<unsigned>( isngn[ind] ), ind ) );
    }

    //std::cout << "INET: \n";
    //std::copy( inet.begin(), inet.end(), std::ostream_iterator<int>( std::cout, " " ) );
    //std::cout << std::endl;
    //std::cout << "ISNGN: \n";
    //std::copy( isngn.begin(), isngn.end(), std::ostream_iterator<int>( std::cout, " " ) );
    //std::cout << std::endl << std::flush;
    //std::cout << "ISEGN: \n";
    //std::copy( isegn.begin(), isegn.end(), std::ostream_iterator<int>( std::cout, " " ) );
    //std::cout << std::endl << std::flush;
    //MPI_Barrier( PETSC_COMM_WORLD );

    // renumber nodes in the inet array to locals
    int indInet = 0;
    for ( int iEle = 0; iEle < isegn.size(); iEle++ ) {
        int nne = nnet[ iEle ];
        for ( unsigned ien = 0; ien < nne; ien++ ) {

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
    //std::cout << "I have identified following dimensions: max " << elDimMax << ", min " << elDimMin << std::endl;

    bddc_ls -> load_mesh( spaceDim, numNodes, numDofsInt, inet, nnet, nndf, isegn, isngn, isngn, xyz, element_permeability, meshDim );
}




//=============================================================================
// DESTROY WATER MH SYSTEM STRUCTURE
//=============================================================================
DarcyFlowMH_Steady::~DarcyFlowMH_Steady() {
    //if (schur2 != NULL) {
    //	delete schur2;
    	//ISDestroy(&IS2);
    //}
    //if (schur1 != NULL) {
    //	delete schur1;
    	//ISDestroy(&IS1);
    //}
    if (schur0 != NULL) delete schur0;

	delete edge_ds;
	delete el_ds;
	delete side_ds;

	xfree(el_4_loc);
	xfree(row_4_el);
	xfree(side_id_4_loc);
	xfree(side_row_4_id);
	xfree(edge_4_loc);
	xfree(row_4_edge);

	if (solution != NULL) {
		VecDestroy(&sol_vec);
		xfree(solution);
	}

	VecScatterDestroy(&par_to_all);

}


// ================================================
// PARALLLEL PART
//

/**
 * Make connectivity graph of the second Schur complement and compute optimal partitioning.
 * This verison assign 1D and 2D edges to one processor, represent them as
 * a weighted vertex, and 2D-3D neghbourings as weighted edges. 3D edges are
 * then distributed over all procesors.
 */
/*
void make_edge_conection_graph(Mesh *mesh, SparseGraph * &graph) {

    Distribution edistr = graph->get_distr();
    //Edge *edg;
    Element *ele;
    int li, eid, , i_edg;
    unsigned int si, i_neigh;
    int e_weight;

    int edge_dim_weights[3] = { 100, 10, 1 };
    F_ENTRY;

    i_edg=0;
    FOR_EDGES(mesh, edg) {

        // skip non-local edges
        if (!edistr.is_local(i_edg))   continue;

        e_weight = edge_dim_weights[edg->side(0)->element()->dim() - 1];
        // for all connected elements
        FOR_EDGE_SIDES( edg, li ) {
            ASSERT( edg->side(li)->valid(),"NULL side of edge.");
            ele = edg->side(li)->element();
            ASSERT(NONULL(ele),"NULL element of side.");

            // for sides of connected element, excluding edge itself
            for(si=0; si<ele->n_sides(); si++) {
                eid = ele->side(si)->edge_idx();
                if (eid != i_edg)
                    graph->set_edge(i_edg, eid, e_weight);
            }

            // include connections from lower dim. edge
            // to the higher dimension
            for (i_neigh = 0; i_neigh < ele->n_neighs_vb; i_neigh++) {
                eid = ele->neigh_vb[i_neigh]->edge_idx();
                graph->set_edge(i_edg, eid, e_weight);
                graph->set_edge(eid, i_edg, e_weight);
            }
        }
        i_edg++;
    }

    graph->finalize();
}
*/
/**
 * Make connectivity graph of elements of mesh - dual graph: elements vertices of graph.
 */
/*
void make_element_connection_graph(Mesh *mesh, SparseGraph * &graph, bool neigh_on) {

    Distribution edistr = graph->get_distr();

    Edge *edg;
    int li, e_idx, i_neigh;
    int i_s, n_s;
    F_ENTRY;

    //int elDimMax = 1;
    //int elDimMin = 3;
    //FOR_ELEMENTS(mesh, ele) {
    //    //xprintf(Msg,"Element id %d , its index %d.\n",ele.id(), i_ele);
    //    int elDim = ele->dim();
    //    elDimMax = std::max( elDimMax, elDim );
    //    elDimMin = std::min( elDimMin, elDim );
    //}
    //std::cout << "max and min element dimensions: " << elDimMax << " " << elDimMin << std::endl;

    FOR_ELEMENTS(mesh, ele) {
        //xprintf(Msg,"Element id %d , its index %d.\n",ele.id(), i_ele);

        // skip non-local elements
        if (!edistr.is_local(ele.index()))
            continue;

        // for all connected elements
        FOR_ELEMENT_SIDES( ele, si ) {
            edg = ele->side(si)->edge();

            FOR_EDGE_SIDES( edg, li ) {
                ASSERT(edg->side(li)->valid(),"NULL side of edge.");
                e_idx = ELEMENT_FULL_ITER(mesh, edg->side(li)->element()).index();

                // for elements of connected elements, excluding element itself
                if (e_idx != ele.index()) {
                    graph->set_edge(ele.index(), e_idx);
                }
            }
        }

        // include connections from lower dim. edge
        // to the higher dimension
        if (neigh_on) {
            for (i_neigh = 0; i_neigh < ele->n_neighs_vb; i_neigh++) {
               n_s = ele->neigh_vb[i_neigh]->edge()->n_sides;
                for (i_s = 0; i_s < n_s; i_s++) {
                   e_idx=ELEMENT_FULL_ITER(mesh, ele->neigh_vb[i_neigh]->edge()->side(i_s)->element()).index();
                    graph->set_edge(ele.index(), e_idx);
                    graph->set_edge(e_idx, ele.index());
                }
            }
        }
    }
    graph->finalize();
}*/


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
    //DBGPRINT_INT("side_shift",np,side_shift);
    //DBGPRINT_INT("el_shift",np,el_shift);
    //DBGPRINT_INT("edge_shift",np,edge_shift);
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

// ====================================================================================
// - compute optimal edge partitioning
// - compute appropriate partitioning of elements and sides
// - make arrays: *_id_4_loc and *_row_4_id to allow parallel assembly of the MH matrix
// ====================================================================================
void DarcyFlowMH_Steady::prepare_parallel( const Input::AbstractRecord in_rec) {
    
    START_TIMER("prepare parallel");
    
    int *loc_part; // optimal (edge,el) partitioning (local chunk)
    int *id_4_old; // map from old idx to ids (edge,el)
    int i, loc_i;

    int e_idx;

    
    F_ENTRY;

    //ierr = MPI_Barrier(PETSC_COMM_WORLD);
    //ASSERT(ierr == 0, "Error in MPI_Barrier.");

        id_4_old = new int[mesh_->n_elements()];
        i = 0;
        FOR_ELEMENTS(mesh_, el) id_4_old[i++] = el.index();

        mesh_->get_part()->id_maps(mesh_->element.size(), id_4_old, el_ds, el_4_loc, row_4_el);
        //DBGPRINT_INT("el_4_loc",el_ds->lsize(),el_4_loc);
        //xprintf(Msg,"Number of elements in subdomain %d \n",el_ds->lsize());
        delete[] id_4_old;

        //el_ds->view( std::cout );
        //
        //DBGMSG("Compute appropriate edge partitioning ...\n");
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
                //xprintf(Msg,"Index of edge: %d first element: %d \n",edgid,e_idx);
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

    //edge_ds->view( std::cout );
    //side_ds->view( std::cout );

    /*
     DBGPRINT_INT("edge_id_4_loc",edge_ds->lsize,edge_id_4_loc);
     DBGPRINT_INT("el_4_loc",el_ds->lsize,el_4_loc);
     DBGPRINT_INT("side_id_4_loc",side_ds->lsize,side_id_4_loc);
     DBGPRINT_INT("edge_row_4_id",mesh_->n_edges,edge_row_4_id);
     DBGPRINT_INT("el_row_4_id",mesh_->max_elm_id+1,el_row_4_id);
     DBGPRINT_INT("side_row_4_id",mesh_->max_side_id+1,side_row_4_id);
     */
    // convert row_4_id arrays from separate numberings to global numbering of rows
    make_row_numberings();
    //DBGPRINT_INT("edge_row_4_id",mesh_->n_edges,edge_row_4_id);
    //DBGPRINT_INT("el_row_4_id",mesh_->max_elm_id+1,el_row_4_id);
    //DBGPRINT_INT("side_row_4_id",mesh_->max_side_id+1,side_row_4_id);

    lsize = side_ds->lsize() + el_ds->lsize() + edge_ds->lsize();


    // prepare global_row_4_sub_row

#ifdef HAVE_BDDCML
    if (in_rec.type() ==  LinSys_BDDC::input_type) {
        // auxiliary
        Element *el;
        int side_row, edge_row;

        //xprintf(Msg,"Compute mapping of local subdomain rows to global rows.\n");

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
                int side_row = side_row_4_id[ mh_dh.side_dof( el->side(i) ) ];
		        int edge_row = row_4_edge[el->side(i)->edge_idx()];

		        global_row_4_sub_row->insert( side_row );
		        global_row_4_sub_row->insert( edge_row );

                // edge neighbouring overlap
                //if (edg->neigh_vb != NULL) {
		//	int neigh_el_row=row_4_el[mesh_->element.index(edg->neigh_vb->element[0])];
                //}
            }

            for (unsigned int i_neigh = 0; i_neigh < el->n_neighs_vb; i_neigh++) {
                // mark this edge
                int edge_row = row_4_edge[el->neigh_vb[i_neigh]->edge_idx() ];
                global_row_4_sub_row->insert( edge_row );
            }
        }
        global_row_4_sub_row->finalize();
    }
#endif // HAVE_BDDCML

    // common to both solvers - create renumbering of unknowns
    solver_indices_.reserve(size);
    FOR_ELEMENTS(mesh_, ele) {
        FOR_ELEMENT_SIDES(ele,si) {
            solver_indices_.push_back( side_row_4_id[ mh_dh.side_dof( ele->side(si) ) ] );
        }
    }
    FOR_ELEMENTS(mesh_, ele) {
        solver_indices_.push_back( row_4_el[ele.index()] );
    }

    unsigned int i_edg=0;
    FOR_EDGES(mesh_, edg) {
        solver_indices_.push_back( row_4_edge[i_edg++] );
    }
    ASSERT( solver_indices_.size() == size, "Size of array does not match number of fills.\n" );
    //std::cout << "Solve rindices:" << std::endl;
    //std::copy( solver_indices_.begin(), solver_indices_.end(), std::ostream_iterator<int>( std::cout, " " ) );
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
    //printf("[%d] rows: %d off_rows: %d on: %d off: %d\n",distr.myp(),last-first,n_off_rows,n_on,n_off);
}












// ========================
// unsteady

DarcyFlowMH_Unsteady::DarcyFlowMH_Unsteady(Mesh &mesh_in, const Input::Record in_rec)
    : DarcyFlowMH_Steady(mesh_in, in_rec)
{
    delete time_; // delete steady TG
 
    time_ = new TimeGovernor(in_rec.val<Input::Record>("time"), equation_mark_type_);
    
    time_->fix_dt_until_mark();
    setup_time_term();
}



void DarcyFlowMH_Unsteady::setup_time_term() {

    // have created full steady linear system
    // save diagonal of steady matrix
    VecCreateMPI(PETSC_COMM_WORLD, rows_ds->lsize(), PETSC_DETERMINE, &(steady_diagonal));
    MatGetDiagonal(schur0->get_matrix(), steady_diagonal);

    // read inital condition
    VecZeroEntries(schur0->get_solution());

    double *local_sol = schur0->get_solution_array();

    PetscScalar *local_diagonal;
    VecDuplicate(steady_diagonal, &new_diagonal);
    VecZeroEntries(new_diagonal);
    VecGetArray(new_diagonal,& local_diagonal);

    // apply initial condition and modify matrix diagonal
    // cycle over local element rows
    ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);

    DBGMSG("Setup with dt: %f\n",time_->dt());
    for (unsigned int i_loc_el = 0; i_loc_el < el_ds->lsize(); i_loc_el++) {
        ele = mesh_->element(el_4_loc[i_loc_el]);
        int i_loc_row = i_loc_el + side_ds->lsize();

        // set initial condition
        local_sol[i_loc_row] = data.init_pressure.value(ele->centre(),ele->element_accessor());
        // set new diagonal
        local_diagonal[i_loc_row]= - data.storativity.value(ele->centre(), ele->element_accessor()) * 
                                  ele->measure() / time_->dt();
    }
    VecRestoreArray(new_diagonal,& local_diagonal);
    MatDiagonalSet(schur0->get_matrix(), new_diagonal, ADD_VALUES);

    // set previous solution as copy of initial condition
    VecDuplicate(schur0->get_solution(), &previous_solution);
    VecCopy(schur0->get_solution(), previous_solution);

    // save RHS
    VecDuplicate(schur0->get_rhs(), &steady_rhs);
    VecCopy(schur0->get_rhs(), steady_rhs);

    solution_changed_for_scatter=true;
    schur0->set_matrix_changed();
}

void DarcyFlowMH_Unsteady::modify_system() {
  START_TIMER("modify system");
  if (time_->is_changed_dt()) {
      MatDiagonalSet(schur0->get_matrix(),steady_diagonal, INSERT_VALUES);

      VecScale(new_diagonal, time_->last_dt()/time_->dt());
      MatDiagonalSet(schur0->get_matrix(),new_diagonal, ADD_VALUES);
      schur0->set_matrix_changed();
  }

    // modify RHS - add previous solution
    VecPointwiseMult(schur0->get_rhs(), new_diagonal, schur0->get_solution());
    VecAXPY(schur0->get_rhs(), 1.0, steady_rhs);
    schur0->set_rhs_changed();

    // swap solutions
    VecSwap(previous_solution, schur0->get_solution());
}

// ========================
// unsteady

DarcyFlowLMH_Unsteady::DarcyFlowLMH_Unsteady(Mesh &mesh_in, const  Input::Record in_rec)
    : DarcyFlowMH_Steady(mesh_in, in_rec)
{
    delete time_; // delete steady TG
 
    time_ = new TimeGovernor(in_rec.val<Input::Record>("time"), equation_mark_type_);

    time_->fix_dt_until_mark();
    setup_time_term();
}





void DarcyFlowLMH_Unsteady::setup_time_term()
{
    // have created full steady linear system
     // save diagonal of steady matrix
     VecCreateMPI(PETSC_COMM_WORLD,rows_ds->lsize(),PETSC_DETERMINE,&(steady_diagonal));
     MatGetDiagonal(schur0->get_matrix(), steady_diagonal);

     // read inital condition
     VecZeroEntries(schur0->get_solution());

     VecDuplicate(steady_diagonal,& new_diagonal);

     // apply initial condition and modify matrix diagonal
     // cycle over local element rows

     ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);
     double init_value;

     for (unsigned int i_loc_el = 0; i_loc_el < el_ds->lsize(); i_loc_el++) {
         ele = mesh_->element(el_4_loc[i_loc_el]);

         init_value = data.init_pressure.value(ele->centre(), ele->element_accessor());

         FOR_ELEMENT_SIDES(ele,i) {
             int edge_row = row_4_edge[ele->side(i)->edge_idx()];
             // set new diagonal
             VecSetValue(new_diagonal,edge_row, - ele->measure() *
                          data.storativity.value(ele->centre(), ele->element_accessor()) *
                          data.cross_section.value(ele->centre(), ele->element_accessor()) /
                          time_->dt() / ele->n_sides(),ADD_VALUES);
             // set initial condition
             VecSetValue(schur0->get_solution(),edge_row,init_value/ele->n_sides(),ADD_VALUES);
         }
     }
     VecAssemblyBegin(new_diagonal);
     VecAssemblyBegin(schur0->get_solution());
     VecAssemblyEnd(new_diagonal);
     VecAssemblyEnd(schur0->get_solution());

     MatDiagonalSet(schur0->get_matrix(),new_diagonal, ADD_VALUES);

     // set previous solution as copy of initial condition
     VecDuplicate(schur0->get_solution(), &previous_solution);
     VecCopy(schur0->get_solution(), previous_solution);

     // save RHS
     VecDuplicate(schur0->get_rhs(), &steady_rhs);
     VecCopy(schur0->get_rhs(),steady_rhs);

     // auxiliary vector for time term
     VecDuplicate(schur0->get_rhs(), &time_term);

     solution_changed_for_scatter=true;
     schur0->set_matrix_changed();

}

void DarcyFlowLMH_Unsteady::modify_system() {
    START_TIMER("modify system");
    if (time_->is_changed_dt()) {
        MatDiagonalSet(schur0->get_matrix(),steady_diagonal, INSERT_VALUES);
        VecScale(new_diagonal, time_->last_dt()/time_->dt());
        MatDiagonalSet(schur0->get_matrix(),new_diagonal, ADD_VALUES);
        schur0->set_matrix_changed();
    }

    // modify RHS - add previous solution
    VecPointwiseMult(schur0->get_rhs(), new_diagonal, schur0->get_solution());
    VecAXPY(schur0->get_rhs(), 1.0, steady_rhs);
    schur0->set_rhs_changed();

    // swap solutions
    VecSwap(previous_solution, schur0->get_solution());
}

// TODO: make this operating on parallel solution
// i.e. access from elements to edge values (possibly by constructing specific matrix)

// is it really necessary what is natural value of element pressures ?
// Since
void DarcyFlowLMH_Unsteady::postprocess() {
    int i_loc, side_row, loc_edge_row, i;
    Edge* edg;
    ElementIter ele;
    double new_pressure, old_pressure, time_coef;

    PetscScalar *loc_prev_sol;
    VecGetArray(previous_solution, &loc_prev_sol);

    // modify side fluxes in parallel
    // for every local edge take time term on digonal and add it to the corresponding flux
    for (i_loc = 0; i_loc < edge_ds->lsize(); i_loc++) {

        edg = &( mesh_->edges[ edge_4_loc[i_loc] ] );
        loc_edge_row = side_ds->lsize() + el_ds->lsize() + i_loc;

        new_pressure = (schur0->get_solution_array())[loc_edge_row];
        old_pressure = loc_prev_sol[loc_edge_row];
        FOR_EDGE_SIDES(edg,i) {
          ele = edg->side(i)->element();
          side_row = side_row_4_id[ mh_dh.side_dof( edg->side(i) ) ];
          time_coef = - ele->measure() *
              data.cross_section.value(ele->centre(), ele->element_accessor()) *
              data.storativity.value(ele->centre(), ele->element_accessor()) /
              time_->dt() / ele->n_sides();
            VecSetValue(schur0->get_solution(), side_row, time_coef * (new_pressure - old_pressure), ADD_VALUES);
        }
    }
  VecRestoreArray(previous_solution, &loc_prev_sol);

    VecAssemblyBegin(schur0->get_solution());
    VecAssemblyEnd(schur0->get_solution());

  //DarcyFlowMH_Steady::postprocess();

    int side_rows[4];
    double values[4];
    //ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);

  // modify side fluxes in parallel
  // for every local edge take time term on digonal and add it to the corresponding flux

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
}

//-----------------------------------------------------------------------------
// vim: set cindent:
