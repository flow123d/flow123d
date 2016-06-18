/*
 * richards_lmh.cc
 *
 *  Created on: Sep 16, 2015
 *      Author: jb
 */


#include "system/global_defs.h"
#include "system/sys_profiler.hh"


#include "input/input_type.hh"
#include "input/factory.hh"
#include "flow/richards_lmh.hh"
#include "flow/darcy_flow_mh_output.hh"
#include "tools/time_governor.hh"

#include "petscmat.h"
#include "petscviewer.h"
#include "petscerror.h"
#include <armadillo>

#include "la/schur.hh"

#include "coupling/balance.hh"

#include "fields/vec_seq_double.hh"

// in the third_party/FADBAD++ dir, namespace "fadbad"
#include "fadbad.h"
#include "badiff.h"
#include "fadiff.h"

#include "flow/assembly_lmh.hh"


FLOW123D_FORCE_LINK_IN_CHILD(richards_lmh);


namespace it=Input::Type;


RichardsLMH::EqData::EqData()
{

    ADD_FIELD(water_content_saturated,
            "Saturated water content (($ \theta_s $)).\n"
            "Relative volume of the water in a reference volume of a saturated porous media.", "0.0");
        water_content_saturated.units( UnitSI::dimensionless() );

    ADD_FIELD(water_content_residual,
            "Residual water content (($ \theta_r $)).\n"
            "Relative volume of the water in a reference volume of an ideally dry porous media.", "0.0");
        water_content_residual.units( UnitSI::dimensionless() );

    ADD_FIELD(genuchten_p_head_scale,
            "The van Genuchten pressure head scaling parameter (($ \alpha $)).\n"
            "The parameter of the van Genuchten's model to scale the pressure head."
            "Related to the inverse of the air entry pressure, i.e. the pressure where the relative water content starts to decrease below 1.", "0.0");
        genuchten_p_head_scale.units( UnitSI().m(-1) );

    ADD_FIELD(genuchten_n_exponent,
            "The van Genuchten exponent parameter (($ n $)).\n", "2.0");
        genuchten_n_exponent.units( UnitSI::dimensionless() );

}


const it::Record & RichardsLMH::get_input_type() {
    it::Record field_descriptor = it::Record("RichardsLMH_Data",FieldCommon::field_descriptor_record_description("RichardsLMH_Data"))
    .copy_keys( DarcyMH::type_field_descriptor() )
    .copy_keys( RichardsLMH::EqData().make_field_descriptor_type("RichardsLMH_Data_aux") )
    .close();

    return it::Record("Flow_Richards_LMH", "Lumped Mixed-Hybrid solver for unsteady saturated Darcy flow.")
        .derive_from(DarcyFlowInterface::get_input_type())
        .copy_keys(DarcyMH::get_input_type())
        .declare_key("input_fields", it::Array( field_descriptor ), it::Default::obligatory(),
                "Input data for Darcy flow model.")

        .close();
}


const int RichardsLMH::registrar =
        Input::register_class< RichardsLMH, Mesh &, const Input::Record >("Flow_Richards_LMH") +
        RichardsLMH::get_input_type().size();



RichardsLMH::RichardsLMH(Mesh &mesh_in, const  Input::Record in_rec)
    : DarcyMH(mesh_in, in_rec)
{
    data_ = make_shared<EqData>();
    DarcyMH::data_ = data_;
    EquationBase::eq_data_ = data_.get();
    //data_->edge_new_local_4_mesh_idx_ = &(this->edge_new_local_4_mesh_idx_);
}

void RichardsLMH::initialize_specific() {

    // create edge vectors
    unsigned int n_local_edges = mh_dh.edge_new_local_4_mesh_idx_.size();
    unsigned int n_local_sides = mh_dh.side_ds->lsize();
    data_->phead_edge_.resize( n_local_edges);
    data_->water_content_previous_it.resize(n_local_sides);
    data_->water_content_previous_time.resize(n_local_sides);
    data_->capacity.resize(n_local_sides);

    Distribution ds_split_edges(n_local_edges, PETSC_COMM_WORLD);
    vector<int> local_edge_rows(n_local_edges);

    IS is_loc;
    for(auto  item : mh_dh.edge_new_local_4_mesh_idx_) {
        local_edge_rows[item.second]=mh_dh.row_4_edge[item.first];
    }
    ISCreateGeneral(PETSC_COMM_SELF, local_edge_rows.size(),
            &(local_edge_rows[0]), PETSC_COPY_VALUES, &(is_loc));

    VecScatterCreate(schur0->get_solution(), is_loc,
            data_->phead_edge_.petsc_vec(), PETSC_NULL, &solution_2_edge_scatter_);
    ISDestroy(&is_loc);


    // test the scatter
    /*
    vector<unsigned int> loc_to_glob(n_local_edges);
    for(auto item : edge_new_local_4_mesh_idx_)
        loc_to_glob[item.second] = row_4_edge[item.first];

    VectorMPI tmp_solution(rows_ds->lsize());
    for(unsigned int i=0; i< rows_ds->lsize(); i++) tmp_solution[i] = i + rows_ds->begin();
    VecScatterBegin(solution_2_edge_scatter_, tmp_solution.petsc_vec(), phead_edge_.petsc_vec() , INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(solution_2_edge_scatter_, tmp_solution.petsc_vec(), phead_edge_.petsc_vec() , INSERT_VALUES, SCATTER_FORWARD);
    for(unsigned int i=0; i< phead_edge_.data().size(); i++)
        cout << "p: " << el_ds->myp() << "i: " << i
             << "phead: " << phead_edge_[i] << "check: " << loc_to_glob[i] << endl;
    */
}

/*
void RichardsLMH::local_assembly_specific(LocalAssemblyData &local_data)
{

}
*/

void RichardsLMH::read_initial_condition()
{
    // apply initial condition
    // cycle over local element rows
    double init_value;

    for (unsigned int i_loc_el = 0; i_loc_el < mh_dh.el_ds->lsize(); i_loc_el++) {
         auto ele_ac = mh_dh.accessor(i_loc_el);

         init_value = data_->init_pressure.value(ele_ac.centre(), ele_ac.element_accessor());

         FOR_ELEMENT_SIDES(ele_ac.full_iter(),i) {
             int edge_row = ele_ac.edge_row(i);
             uint n_sides_of_edge =  ele_ac.full_iter()->side(i)->edge()->n_sides;
             VecSetValue(schur0->get_solution(),edge_row, init_value/n_sides_of_edge, ADD_VALUES);
         }
         VecSetValue(schur0->get_solution(),ele_ac.ele_row(), init_value,ADD_VALUES);
    }
    VecAssemblyBegin(schur0->get_solution());
    VecAssemblyEnd(schur0->get_solution());

    // set water_content
    // pretty ugly since postprocess change fluxes, which cause bad balance, so we must set them back
    VecCopy(schur0->get_solution(), previous_solution); // store solution vector
    postprocess();
    VecSwap(schur0->get_solution(), previous_solution); // restore solution vector

    //DBGMSG("init sol:\n");
    //VecView( schur0->get_solution(),   PETSC_VIEWER_STDOUT_WORLD);
    //DBGMSG("init water content:\n");
    //VecView( data_->water_content_previous_it.petsc_vec(),   PETSC_VIEWER_STDOUT_WORLD);

    solution_changed_for_scatter=true;
}


void RichardsLMH::prepare_new_time_step()
{
    VecCopy(schur0->get_solution(), previous_solution);
    data_->water_content_previous_time.copy(data_->water_content_previous_it);
    //VecCopy(schur0->get_solution(), previous_solution);
}

bool RichardsLMH::zero_time_term() {
    return data_->storativity.field_result(mesh_->region_db().get_region_set("BULK")) == result_zeros &&
           data_->genuchten_p_head_scale.field_result(mesh_->region_db().get_region_set("BULK")) == result_zeros;
}


void RichardsLMH::assembly_linear_system()
{

    START_TIMER("RicharsLMH::assembly_linear_system");

    VecScatterBegin(solution_2_edge_scatter_, schur0->get_solution(), data_->phead_edge_.petsc_vec() , INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(solution_2_edge_scatter_, schur0->get_solution(), data_->phead_edge_.petsc_vec() , INSERT_VALUES, SCATTER_FORWARD);

    is_linear_ = data_->genuchten_p_head_scale.field_result(mesh_->region_db().get_region_set("BULK")) == result_zeros;

    bool is_steady = zero_time_term();
    //DBGMSG("Assembly linear system\n");
        START_TIMER("full assembly");
        if (typeid(*schur0) != typeid(LinSys_BDDC)) {
            schur0->start_add_assembly(); // finish allocation and create matrix
        }
        data_->time_step_ = time_->dt();
        auto multidim_assembler = AssemblyBase::create< AssemblyLMH >(data_);


        schur0->mat_zero_entries();
        schur0->rhs_zero_entries();

        if (balance_ != nullptr) {
            balance_->start_source_assembly(water_balance_idx_);
            balance_->start_mass_assembly(water_balance_idx_);
        }

        assembly_mh_matrix( multidim_assembler ); // fill matrix

        if (balance_ != nullptr) {
            balance_->finish_source_assembly(water_balance_idx_);
            balance_->finish_mass_assembly(water_balance_idx_);
        }
            //MatView( *const_cast<Mat*>(schur0->get_matrix()), PETSC_VIEWER_STDOUT_WORLD  );
            //VecView( *const_cast<Vec*>(schur0->get_rhs()),   PETSC_VIEWER_STDOUT_WORLD);

        schur0->finish_assembly();
        schur0->set_matrix_changed();


        if (! is_steady) {
            START_TIMER("fix time term");
            //DBGMSG("    setup time term\n");
            // assembly time term and rhs
            solution_changed_for_scatter=true;


            //VecPointwiseMult(*( schur0->get_rhs()), new_diagonal, schur0->get_solution());
            //VecAXPY(*( schur0->get_rhs()), 1.0, steady_rhs);
            //schur0->set_rhs_changed();

            // swap solutions
            //VecSwap(previous_solution, schur0->get_solution());
        }
}

/*
void RichardsLMH::compute_per_element_nonlinearities() {

}
*/


void RichardsLMH::setup_time_term()
{
    FEAL_ASSERT(false).error("Shold not be called.");
    // save diagonal of steady matrix
    //MatGetDiagonal(*( schur0->get_matrix() ), steady_diagonal);
    // save RHS
    //VecCopy(*( schur0->get_rhs()),steady_rhs);
/*
    VecZeroEntries(new_diagonal);

    // modify matrix diagonal
    // cycle over local element rows
    ElementFullIter ele = ELEMENT_FULL_ITER(mesh_, NULL);
    DBGMSG("setup time term with dt: %f\n", time_->dt());


    for (unsigned int i_loc_el = 0; i_loc_el < el_ds->lsize(); i_loc_el++) {
        ele = mesh_->element(el_4_loc[i_loc_el]);

        //data_->init_pressure.value(ele_ac.centre(), ele_ac.element_accessor());

        FOR_ELEMENT_SIDES(ele,i) {
            int edge_row = row_4_edge[ele_ac.side(i)->edge_idx()];
            // set new diagonal
            double diagonal_coef = ele_ac.measure() *
                      data_->storativity.value(ele_ac.centre(), ele_ac.element_accessor()) *
                      data_->cross_section.value(ele_ac.centre(), ele_ac.element_accessor())
                      / ele_ac.n_sides();
            VecSetValue(new_diagonal, edge_row, -diagonal_coef / time_->dt(), ADD_VALUES);

            if (balance_ != nullptr)
                balance_->add_mass_matrix_values(water_balance_idx_, ele_ac.region().bulk_idx(), {edge_row}, {diagonal_coef});


        }
    }
    VecAssemblyBegin(new_diagonal);
    VecAssemblyEnd(new_diagonal);

    MatDiagonalSet(*( schur0->get_matrix() ),new_diagonal, ADD_VALUES);
*/
    /*
    solution_changed_for_scatter=true;
    schur0->set_matrix_changed();


    VecPointwiseMult(*( schur0->get_rhs()), new_diagonal, schur0->get_solution());
    VecAXPY(*( schur0->get_rhs()), 1.0, steady_rhs);
    schur0->set_rhs_changed();

    // swap solutions
    VecSwap(previous_solution, schur0->get_solution());
*/
}



void RichardsLMH::assembly_source_term()
{

}

/*
void RichardsLMH::prepare_new_time_step()
{
    DBGMSG("sol swap\n");
    //VecSwap(previous_solution, schur0->get_solution());
    data_->water_content_previous_time.copy(data_->water_content_previous_it);
}*/

void RichardsLMH::postprocess() {
    //int side_row, loc_edge_row, i;
    //Edge* edg;
    //ElementIter ele;

/*
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
          time_coef =
            VecSetValue(schur0->get_solution(), side_row, time_coef * (new_pressure - old_pressure), ADD_VALUES);
        }
    }

    VecAssemblyBegin(schur0->get_solution());
    VecAssemblyEnd(schur0->get_solution());
*/
    int side_rows[4];
    double values[4];

    VecScatterBegin(solution_2_edge_scatter_, schur0->get_solution(), data_->phead_edge_.petsc_vec() , INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(solution_2_edge_scatter_, schur0->get_solution(), data_->phead_edge_.petsc_vec() , INSERT_VALUES, SCATTER_FORWARD);


  // modify side fluxes in parallel
  // for every local edge take time term on digonal and add it to the corresponding flux
    //PetscScalar *loc_prev_sol;
    auto multidim_assembler = AssemblyBase::create< AssemblyLMH >(data_);

    //VecGetArray(previous_solution, &loc_prev_sol);

    for (unsigned int i_loc = 0; i_loc < mh_dh.el_ds->lsize(); i_loc++) {
      auto ele_ac = mh_dh.accessor(i_loc);
      multidim_assembler[ele_ac.dim()]->update_water_content(ele_ac);

      double ele_scale = ele_ac.measure() *
              data_->cross_section.value(ele_ac.centre(), ele_ac.element_accessor()) / ele_ac.n_sides();
      double ele_source = data_->water_source_density.value(ele_ac.centre(), ele_ac.element_accessor());
      //double storativity = data_->storativity.value(ele_ac.centre(), ele_ac.element_accessor());

      FOR_ELEMENT_SIDES(ele_ac.full_iter(),i) {
          //unsigned int loc_edge_row = ele_ac.edge_local_row(i);
          side_rows[i] = ele_ac.side_row(i);
          double water_content = data_->water_content_previous_it[ele_ac.side_local_idx(i)];
          double water_content_previous_time = data_->water_content_previous_time[ele_ac.side_local_idx(i)];

          values[i] = ele_scale * ele_source - ele_scale * (water_content - water_content_previous_time) / time_->dt();
      }
      VecSetValues(schur0->get_solution(), ele_ac.n_sides(), side_rows, values, ADD_VALUES);
    }
    VecAssemblyBegin(schur0->get_solution());
    //VecRestoreArray(previous_solution, &loc_prev_sol);
    VecAssemblyEnd(schur0->get_solution());
}

