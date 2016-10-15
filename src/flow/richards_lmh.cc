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
    string desc;
    
    desc = R"(
Saturated water content (($ \theta_s $)).
relative volume of the water in a reference volume of a saturated porous media.
)" ;
    ADD_FIELD(water_content_saturated, desc, "0.0");
    water_content_saturated.units( UnitSI::dimensionless() );

    desc = R"(
Residual water content (($ \theta_r $)).
Relative volume of the water in a reference volume of an ideally dry porous media.
)";
    ADD_FIELD(water_content_residual, desc, "0.0");
    water_content_residual.units( UnitSI::dimensionless() );

    desc = R"(
The van Genuchten pressure head scaling parameter (($ \alpha $)).
The parameter of the van Genuchten's model to scale the pressure head.
Related to the inverse of the air entry pressure, i.e. the pressure where the relative water content starts to decrease below 1.
)";
    ADD_FIELD(genuchten_p_head_scale, desc, "0.0");
    genuchten_p_head_scale.units( UnitSI().m(-1) );

    ADD_FIELD(genuchten_n_exponent,
            "The van Genuchten exponent parameter (($ n $)).", "2.0");
    genuchten_n_exponent.units( UnitSI::dimensionless() );

}


const it::Record & RichardsLMH::get_input_type() {
    it::Record field_descriptor = it::Record("RichardsLMH_Data",FieldCommon::field_descriptor_record_description("RichardsLMH_Data"))
    .copy_keys( DarcyMH::type_field_descriptor() )
    .copy_keys( RichardsLMH::EqData().make_field_descriptor_type("RichardsLMH_Data_aux") )
    .close();

    auto model_selection = it::Selection("Flow_Darcy_BC_Type")
            .add_value(SoilModelBase::van_genuchten, "van_genuchten", "Van Genuchten soil model with cutting near zero.")
            .add_value(SoilModelBase::irmay, "irmay", "Irmay model for conductivity, Van Genuchten model for the water content. Suitable for bentonite.")
            .close();


    return it::Record("Flow_Richards_LMH", "Lumped Mixed-Hybrid solver for unsteady saturated Darcy flow.")
        .derive_from(DarcyFlowInterface::get_input_type())
        .copy_keys(DarcyMH::get_input_type())
        .declare_key("input_fields", it::Array( field_descriptor ), it::Default::obligatory(),
                "Input data for Darcy flow model.")
        .declare_key("soil_model", model_selection, it::Default("\"van_genuchten\""),
                "Selection of the globally applied soil model. In future we replace this key by a field for selection of the model."
                "That will allow usage of different soil model in a single simulation.")
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

    auto model_type = input_record_.val<SoilModelBase::SoilModelType>("soil_model");
    if (model_type == SoilModelBase::van_genuchten)
        data_->soil_model_ = std::make_shared<SoilModel_VanGenuchten>();
    else if (model_type == SoilModelBase::irmay)
        data_->soil_model_ = std::make_shared<SoilModel_Irmay>();
    else
        ASSERT(false);

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

}


void RichardsLMH::assembly_source_term()
{

}


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

    //DebugOut() << "init sol:\n";
    //VecView( schur0->get_solution(),   PETSC_VIEWER_STDOUT_WORLD);
    //DebugOut() << "init water content:\n";
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
    //DebugOut() << "Assembly linear system\n";
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
            //DebugOut() << "setup time term\n";
            // assembly time term and rhs
            solution_changed_for_scatter=true;
        }
}



void RichardsLMH::setup_time_term()
{
    FEAL_ASSERT(false).error("Shold not be called.");
}





void RichardsLMH::postprocess() {

    // update structures for balance of water volume
    if (balance_ != nullptr)
      assembly_linear_system();



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
      multidim_assembler[ele_ac.dim()-1]->update_water_content(ele_ac);

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
