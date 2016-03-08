/*
 * richards_lmh.cc
 *
 *  Created on: Sep 16, 2015
 *      Author: jb
 */

#include "input/input_type.hh"
#include "input/factory.hh"
#include "flow/richards_lmh.hh"
#include "flow/darcy_flow_mh_output.hh"
#include "tools/time_governor.hh"

#include "petscmat.h"
#include "petscviewer.h"
#include "petscerror.h"
#include <armadillo>

#include "system/global_defs.h"
#include "system/sys_profiler.hh"
#include "la/schur.hh"

#include "coupling/balance.hh"

#include "fields/vec_seq_double.hh"

FLOW123D_FORCE_LINK_IN_CHILD(richards_lmh);


namespace it=Input::Type;
const it::Record & DarcyFlowLMH_Unsteady::get_input_type() {
    return it::Record("UnsteadyDarcy_LMH", "Lumped Mixed-Hybrid solver for unsteady saturated Darcy flow.")
        .derive_from(DarcyFlowInterface::get_input_type())
        .copy_keys(DarcyFlowMH_Steady::get_input_type())
        .declare_key("time",         TimeGovernor::get_input_type(), it::Default::obligatory(),
                                    "Time governor setting for the unsteady Darcy flow model.")
        .close();
}


const int DarcyFlowLMH_Unsteady::registrar =
        Input::register_class< DarcyFlowLMH_Unsteady, Mesh &, const Input::Record >("UnsteadyDarcy_LMH") +
        DarcyFlowLMH_Unsteady::get_input_type().size();



DarcyFlowLMH_Unsteady::DarcyFlowLMH_Unsteady(Mesh &mesh_in, const  Input::Record in_rec)
    : DarcyFlowMH_Steady(mesh_in, in_rec,false)
{
    time_ = new TimeGovernor(in_rec.val<Input::Record>("time"));
    data_.mark_input_times(this->time());

    data_.set_time(time_->step(), LimitSide::right);

    output_object = new DarcyFlowMHOutput(this, in_rec.val<Input::Record>("output"));
    //balance_->units(output_object->get_output_fields().field_ele_pressure.units()*data_.cross_section.units()*data_.storativity.units());

    //time_->fix_dt_until_mark();
    create_linear_system();
    if ( typeid(SchurComplement) != typeid(*(this->schur0)) ) {
        DBGMSG( "%s != %s\n", typeid(LinSys_PETSC).name(),typeid(*(this->schur0)).name() );
        THROW( ExcMessage() << EI_Message("Only SchurComplement linear solver currently allowed for DarcyFlowLMH.") );
    }

    VecDuplicate(schur0->get_solution(), &previous_solution);
    VecCreateMPI(PETSC_COMM_WORLD,rows_ds->lsize(),PETSC_DETERMINE,&(steady_diagonal));
    VecDuplicate(steady_diagonal,& new_diagonal);
    VecDuplicate(*( schur0->get_rhs()), &steady_rhs);

    assembly_linear_system();
    read_init_condition();
    output_data();
}




void DarcyFlowLMH_Unsteady::read_init_condition()
{
    VecZeroEntries(schur0->get_solution());

    // apply initial condition
    // cycle over local element rows

    ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);
    double init_value;

    for (unsigned int i_loc_el = 0; i_loc_el < el_ds->lsize(); i_loc_el++) {
     ele = mesh_->element(el_4_loc[i_loc_el]);

     init_value = data_.init_pressure.value(ele->centre(), ele->element_accessor());

     FOR_ELEMENT_SIDES(ele,i) {
         int edge_row = row_4_edge[ele->side(i)->edge_idx()];
         VecSetValue(schur0->get_solution(),edge_row,init_value/ele->n_sides(),ADD_VALUES);
     }
    }
    VecAssemblyBegin(schur0->get_solution());
    VecAssemblyEnd(schur0->get_solution());

    solution_changed_for_scatter=true;
}



void DarcyFlowLMH_Unsteady::setup_time_term()
{
    // save diagonal of steady matrix
    MatGetDiagonal(*( schur0->get_matrix() ), steady_diagonal);
    // save RHS
    VecCopy(*( schur0->get_rhs()),steady_rhs);

    VecZeroEntries(new_diagonal);

    // modify matrix diagonal
    // cycle over local element rows
    ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);
    DBGMSG("setup time term with dt: %f\n", time_->dt());

    if (balance_ != nullptr)
        balance_->start_mass_assembly(water_balance_idx_);

    for (unsigned int i_loc_el = 0; i_loc_el < el_ds->lsize(); i_loc_el++) {
        ele = mesh_->element(el_4_loc[i_loc_el]);

        data_.init_pressure.value(ele->centre(), ele->element_accessor());

        FOR_ELEMENT_SIDES(ele,i) {
            int edge_row = row_4_edge[ele->side(i)->edge_idx()];
            // set new diagonal
            double diagonal_coef = ele->measure() *
                      data_.storativity.value(ele->centre(), ele->element_accessor()) *
                      data_.cross_section.value(ele->centre(), ele->element_accessor())
                      / ele->n_sides();
            VecSetValue(new_diagonal, edge_row, -diagonal_coef / time_->dt(), ADD_VALUES);

            if (balance_ != nullptr)
                balance_->add_mass_matrix_values(water_balance_idx_, ele->region().bulk_idx(), {edge_row}, {diagonal_coef});


        }
    }
    VecAssemblyBegin(new_diagonal);
    VecAssemblyEnd(new_diagonal);

    MatDiagonalSet(*( schur0->get_matrix() ),new_diagonal, ADD_VALUES);

    solution_changed_for_scatter=true;
    schur0->set_matrix_changed();

    if (balance_ != nullptr)
        balance_->finish_mass_assembly(water_balance_idx_);
}

void DarcyFlowLMH_Unsteady::modify_system() {
    START_TIMER("modify system");
    //if (time_->step().index()>0)
    //    DBGMSG("dt: %f dt-1: %f indexchanged: %d matrix: %d\n", time_->step().length(), time_->step(-1).length(), time_->is_changed_dt(), schur0->is_matrix_changed() );

    if (time_->is_changed_dt() && time_->step().index()>0) {
        // if time step has changed and setup_time_term not called

        double scale_factor=time_->step(-2).length()/time_->step().length();
        if (scale_factor != 1.0) {
            MatDiagonalSet(*( schur0->get_matrix() ),steady_diagonal, INSERT_VALUES);

            //DBGMSG("Scale factor: %f\n",scale_factor);
            VecScale(new_diagonal, scale_factor);
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


void DarcyFlowLMH_Unsteady::assembly_source_term()
{
    if (balance_ != nullptr)
        balance_->start_source_assembly(water_balance_idx_);

    for (unsigned int i_loc = 0; i_loc < el_ds->lsize(); i_loc++)
    {
        ElementFullIter ele = mesh_->element(el_4_loc[i_loc]);

        // set lumped source
        double diagonal_coef = ele->measure()
                  * data_.cross_section.value(ele->centre(), ele->element_accessor())
                  * data_.water_source_density.value(ele->centre(), ele->element_accessor())
                  / ele->n_sides();

        FOR_ELEMENT_SIDES(ele,i)
        {
            int edge_row = row_4_edge[ele->side(i)->edge_idx()];

            schur0->rhs_set_value(edge_row, -diagonal_coef);

            if (balance_ != nullptr)
                balance_->add_source_rhs_values(water_balance_idx_, ele->region().bulk_idx(), {edge_row}, {diagonal_coef});
        }
    }

    if (balance_ != nullptr)
        balance_->finish_source_assembly(water_balance_idx_);
}


void DarcyFlowLMH_Unsteady::postprocess() {
    int side_row, loc_edge_row, i;
    Edge* edg;
    ElementIter ele;
    double new_pressure, old_pressure, time_coef;

    PetscScalar *loc_prev_sol;
    VecGetArray(previous_solution, &loc_prev_sol);

    // modify side fluxes in parallel
    // for every local edge take time term on diagonal and add it to the corresponding flux
    for (unsigned int i_loc = 0; i_loc < edge_ds->lsize(); i_loc++) {

        edg = &( mesh_->edges[ edge_4_loc[i_loc] ] );
        loc_edge_row = side_ds->lsize() + el_ds->lsize() + i_loc;

        new_pressure = (schur0->get_solution_array())[loc_edge_row];
        old_pressure = loc_prev_sol[loc_edge_row];
        FOR_EDGE_SIDES(edg,i) {
          ele = edg->side(i)->element();
          side_row = side_row_4_id[ mh_dh.side_dof( edg->side(i) ) ];
          time_coef = - ele->measure() *
              data_.cross_section.value(ele->centre(), ele->element_accessor()) *
              data_.storativity.value(ele->centre(), ele->element_accessor()) /
              time_->dt() / ele->n_sides();
            VecSetValue(schur0->get_solution(), side_row, time_coef * (new_pressure - old_pressure), ADD_VALUES);
        }
    }
  VecRestoreArray(previous_solution, &loc_prev_sol);

    VecAssemblyBegin(schur0->get_solution());
    VecAssemblyEnd(schur0->get_solution());

    int side_rows[4];
    double values[4];

  // modify side fluxes in parallel
  // for every local edge take time term on digonal and add it to the corresponding flux

  for (unsigned int i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {
      ele = mesh_->element(el_4_loc[i_loc]);
      FOR_ELEMENT_SIDES(ele,i) {
          side_rows[i] = side_row_4_id[ mh_dh.side_dof( ele->side(i) ) ];
          values[i] = 1.0 * ele->measure() *
            data_.cross_section.value(ele->centre(), ele->element_accessor()) *
            data_.water_source_density.value(ele->centre(), ele->element_accessor()) /
            ele->n_sides();
      }
      VecSetValues(schur0->get_solution(), ele->n_sides(), side_rows, values, ADD_VALUES);
  }
  VecAssemblyBegin(schur0->get_solution());
  VecAssemblyEnd(schur0->get_solution());
}

