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
 * @file    dual_porosity.cc
 * @brief   
 */

#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "reaction/dual_porosity.hh"
#include "reaction/reaction_term.hh"
#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "la/distribution.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "mesh/region.hh"
#include "fields/field_elementwise.hh" 

#include "reaction/sorption.hh"
#include "reaction/first_order_reaction.hh"
#include "reaction/radioactive_decay.hh"
#include "semchem/semchem_interface.hh"
#include "input/factory.hh"

FLOW123D_FORCE_LINK_IN_CHILD(dualPorosity)


using namespace Input::Type;



const Record & DualPorosity::get_input_type() {
    return Record("DualPorosity",
            "Dual porosity model in transport problems.\n"
            "Provides computing the concentration of substances in mobile and immobile zone.\n"
            )
		.derive_from(ReactionTerm::get_input_type())
		.declare_key("input_fields", Array(DualPorosity::EqData().make_field_descriptor_type("DualPorosity")), Default::obligatory(),
						"Containes region specific data necessary to construct dual porosity model.")
		.declare_key("scheme_tolerance", Double(0.0), Default("1e-3"),
					 "Tolerance according to which the explicit Euler scheme is used or not."
					 "Set 0.0 to use analytic formula only (can be slower).")

		.declare_key("reaction_mobile", ReactionTerm::get_input_type(), Default::optional(), "Reaction model in mobile zone.")
		.declare_key("reaction_immobile", ReactionTerm::get_input_type(), Default::optional(), "Reaction model in immobile zone.")
		.declare_key("output",
		                    EqData().output_fields.make_output_type("DualPorosity", ""),
		                    IT::Default("{ \"fields\": [ \"conc_immobile\" ] }"),
		                    "Setting of the fields output.")
		.close();
}
    
const int DualPorosity::registrar =
		Input::register_class< DualPorosity, Mesh &, Input::Record >("DualPorosity") +
		DualPorosity::get_input_type().size();

DualPorosity::EqData::EqData()
{
  *this += diffusion_rate_immobile
           .name("diffusion_rate_immobile")
           .description("Diffusion coefficient of non-equilibrium linear exchange between mobile and immobile zone.")
           .input_default("0")
           .units( UnitSI().s(-1) );
  
  *this += porosity_immobile
          .name("porosity_immobile")
          .description("Porosity of the immobile zone.")
          .input_default("0")
          .units( UnitSI::dimensionless() );

  *this += init_conc_immobile
          .name("init_conc_immobile")
          .description("Initial concentration of substances in the immobile zone.")
          .units( UnitSI().kg().m(-3) );

  //creating field for porosity that is set later from the governing equation (transport)
  *this +=porosity
        .name("porosity")
        .units( UnitSI::dimensionless() )
        .flags( FieldFlag::input_copy );

  *this += conc_immobile
          .name("conc_immobile")
          .units( UnitSI().kg().m(-3) )
          .flags(FieldFlag::equation_result);

  output_fields += *this;

}

DualPorosity::DualPorosity(Mesh &init_mesh, Input::Record in_rec)
	: ReactionTerm(init_mesh, in_rec)
{
  //set pointer to equation data fieldset
  this->eq_data_ = &data_;
  
  //reads input and creates possibly other reaction terms
  make_reactions();
  //read value from input
  scheme_tolerance_ = input_record_.val<double>("scheme_tolerance");
}

DualPorosity::~DualPorosity(void)
{
  VecScatterDestroy(&(vconc_out_scatter));


  for (unsigned int sbi = 0; sbi < substances_.size(); sbi++)
  {

      //no mpi vectors
      VecDestroy(&(vconc_immobile[sbi]));
      delete [] conc_immobile[sbi];
  }

  delete [] vconc_immobile;
  delete [] conc_immobile;
}


void DualPorosity::make_reactions() {
    Input::Iterator<Input::AbstractRecord> reactions_it = input_record_.find<Input::AbstractRecord>("reaction_mobile");
    if ( reactions_it )
    {
      // TODO: allowed instances in this case are only
      // FirstOrderReaction, RadioactiveDecay and SorptionMob
      reaction_mobile = (*reactions_it).factory< ReactionTerm, Mesh &, Input::Record >(*mesh_, *reactions_it);
    } else
    {
      reaction_mobile = nullptr;
    }

    reactions_it = input_record_.find<Input::AbstractRecord>("reaction_immobile");
    if ( reactions_it )
    {
      // TODO: allowed instances in this case are only
      // FirstOrderReaction, RadioactiveDecay and SorptionImmob
      reaction_immobile = (*reactions_it).factory< ReactionTerm, Mesh &, Input::Record >(*mesh_, *reactions_it);
    } else
    {
      reaction_immobile = nullptr;
    }

}

void DualPorosity::initialize()
{
  OLD_ASSERT(distribution_ != nullptr, "Distribution has not been set yet.\n");
  OLD_ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  OLD_ASSERT(output_stream_,"Null output stream.");
  OLD_ASSERT_LESS(0, substances_.size());
  
  //allocating memory for immobile concentration matrix
  conc_immobile = new double* [substances_.size()];
  conc_immobile_out.clear();
  conc_immobile_out.resize( substances_.size() );
  for (unsigned int sbi = 0; sbi < substances_.size(); sbi++)
  {
    conc_immobile[sbi] = new double [distribution_->lsize()];
    conc_immobile_out[sbi].resize( distribution_->size() );
  }
  allocate_output_mpi();
  
  initialize_fields();

  if(reaction_mobile)
  {
    reaction_mobile->substances(substances_)
                .output_stream(output_stream_)
                .concentration_matrix(concentration_matrix_, distribution_, el_4_loc_, row_4_el_)
                .set_time_governor(*time_);
    reaction_mobile->initialize();
  }

  if(reaction_immobile)
  {
    reaction_immobile->substances(substances_)
                .output_stream(output_stream_)
                .concentration_matrix(conc_immobile, distribution_, el_4_loc_, row_4_el_)
                .set_time_governor(*time_);
    reaction_immobile->initialize();
  }

}

void DualPorosity::initialize_fields()
{
  data_.set_components(substances_.names());
  //setting fields that are set from input file
  input_data_set_+=data_;
  input_data_set_.set_input_list(input_record_.val<Input::Array>("input_fields"));

  //setting fields in data
  data_.set_mesh(*mesh_);
  
  //initialization of output
  //output_array = input_record_.val<Input::Array>("output_fields");
  data_.output_fields.set_components(substances_.names());
  data_.output_fields.set_mesh(*mesh_);
  data_.output_fields.output_type(OutputTime::ELEM_DATA);
  data_.conc_immobile.setup_components();
  for (unsigned int sbi=0; sbi<substances_.size(); sbi++)
  {
    // create shared pointer to a FieldElementwise and push this Field to output_field on all regions
	auto output_field_ptr = conc_immobile_out[sbi].create_field<3, FieldValue<3>::Scalar>(substances_.size());
    data_.conc_immobile[sbi].set_field(mesh_->region_db().get_region_set("ALL"), output_field_ptr, 0);
  }
  //output_stream_->add_admissible_field_names(output_array);
  data_.output_fields.initialize(output_stream_, input_record_.val<Input::Record>("output"),time());
}


void DualPorosity::zero_time_step()
{
  OLD_ASSERT(distribution_ != nullptr, "Distribution has not been set yet.\n");
  OLD_ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  OLD_ASSERT(output_stream_,"Null output stream.");
  OLD_ASSERT_LESS(0, substances_.size());
 
  //coupling - passing fields
  if(reaction_mobile)
  if (typeid(*reaction_mobile) == typeid(SorptionMob))
  {
	  reaction_mobile->data().set_field("porosity", data_["porosity"]);
	  reaction_mobile->data().set_field("porosity_immobile", data_["porosity_immobile"]);
  }
  if(reaction_immobile)
  if (typeid(*reaction_immobile) == typeid(SorptionImmob))
  {
	  reaction_immobile->data().set_field("porosity", data_["porosity"]);
	  reaction_immobile->data().set_field("porosity_immobile", data_["porosity_immobile"]);
  }
  
  data_.set_time(time_->step(0), LimitSide::right);
  set_initial_condition();
  
  // write initial condition
  //output_vector_gather();
  //data_.output_fields.set_time(time_->step(0), LimitSide::right);
  //data_.output_fields.output(output_stream_);

  output_data();
  
  if(reaction_mobile)
    reaction_mobile->zero_time_step();

  if(reaction_immobile)
    reaction_immobile->zero_time_step();

}

void DualPorosity::set_initial_condition()
{
  //setting initial condition for immobile concentration matrix
  for (unsigned int loc_el = 0; loc_el < distribution_->lsize(); loc_el++)
  {
    unsigned int index = el_4_loc_[loc_el];
    ElementAccessor<3> ele_acc = mesh_->element_accessor(index);
    arma::vec value = data_.init_conc_immobile.value(ele_acc.centre(), ele_acc);
        
    for (unsigned int sbi=0; sbi < substances_.size(); sbi++)
    {
      conc_immobile[sbi][loc_el] = value(sbi);
    }
  }
}

void DualPorosity::update_solution(void) 
{
  data_.set_time(time_->step(-2), LimitSide::right);
 
  START_TIMER("dual_por_exchange_step");
  for (unsigned int loc_el = 0; loc_el < distribution_->lsize(); loc_el++) 
  {
    compute_reaction(conc_immobile, loc_el);
  }
  END_TIMER("dual_por_exchange_step");
  
  if(reaction_mobile) reaction_mobile->update_solution();
  if(reaction_immobile) reaction_immobile->update_solution();
}


double **DualPorosity::compute_reaction(double **concentrations, int loc_el) 
{
  unsigned int sbi;
  double conc_average, // weighted (by porosity) average of concentration
         conc_mob, conc_immob,  // new mobile and immobile concentration
         previous_conc_mob, previous_conc_immob, // mobile and immobile concentration in previous time step
         conc_max, //difference between concentration and average concentration
         por_mob, por_immob; // mobile and immobile porosity
   
  // get data from fields
  ElementFullIter ele = mesh_->element(el_4_loc_[loc_el]);
  por_mob = data_.porosity.value(ele->centre(),ele->element_accessor());
  por_immob = data_.porosity_immobile.value(ele->centre(),ele->element_accessor());
  arma::Col<double> diff_vec = data_.diffusion_rate_immobile.value(ele->centre(), ele->element_accessor());
 
    // if porosity_immobile == 0 then mobile concentration stays the same 
    // and immobile concentration cannot change
    if (por_immob == 0.0) return conc_immobile;
    
    double exponent,
           temp_exponent = (por_mob + por_immob) / (por_mob * por_immob) * time_->dt();
  
    for (sbi = 0; sbi < substances_.size(); sbi++) //over all substances
    {
        exponent = diff_vec[sbi] * temp_exponent;
        //previous values
        previous_conc_mob = concentration_matrix_[sbi][loc_el];
        previous_conc_immob = conc_immobile[sbi][loc_el];
        
        // ---compute average concentration------------------------------------------
        conc_average = ((por_mob * previous_conc_mob) + (por_immob * previous_conc_immob)) 
                       / (por_mob + por_immob);
        
        conc_max = std::max(previous_conc_mob-conc_average, previous_conc_immob-conc_average);
        
        // the following 2 conditions guarantee:
        // 1) stability of forward Euler's method
        // 2) the error of forward Euler's method will not be large
        if(time_->dt() <= por_mob*por_immob/(max(diff_vec)*(por_mob+por_immob)) &&
           conc_max <= (2*scheme_tolerance_/(exponent*exponent)*conc_average))               // forward euler
        {
            double temp = diff_vec[sbi]*(previous_conc_immob - previous_conc_mob) * time_->dt();
            // ---compute concentration in mobile area
            conc_mob = temp / por_mob + previous_conc_mob;

            // ---compute concentration in immobile area
            conc_immob = -temp / por_immob + previous_conc_immob;
        }
        else                                                        //analytic solution
        {
            double temp = exp(-exponent);
            // ---compute concentration in mobile area
            conc_mob = (previous_conc_mob - conc_average) * temp + conc_average;

            // ---compute concentration in immobile area
            conc_immob = (previous_conc_immob - conc_average) * temp + conc_average;
        }
        
        concentration_matrix_[sbi][loc_el] = conc_mob;
        conc_immobile[sbi][loc_el] = conc_immob;
    }
    
  return conc_immobile;
}


void DualPorosity::allocate_output_mpi(void )
{
    int sbi, n_subst;
    n_subst = substances_.size();

    vconc_immobile = new Vec [n_subst];


    for (sbi = 0; sbi < n_subst; sbi++) {
        VecCreateMPIWithArray(PETSC_COMM_WORLD,1, distribution_->lsize(), mesh_->n_elements(), conc_immobile[sbi],
                &vconc_immobile[sbi]);
        VecZeroEntries(vconc_immobile[sbi]);

        //  if(rank == 0)
        VecZeroEntries(conc_immobile_out[sbi].get_data_petsc());
    }
    
    // create output vector scatter
    IS is;
    ISCreateGeneral(PETSC_COMM_SELF, mesh_->n_elements(), row_4_el_, PETSC_COPY_VALUES, &is); //WithArray
    VecScatterCreate(vconc_immobile[0], is, conc_immobile_out[0].get_data_petsc(), PETSC_NULL, &vconc_out_scatter);
    ISDestroy(&(is));
}


void DualPorosity::output_vector_gather() 
{
    unsigned int sbi;

    for (sbi = 0; sbi < substances_.size(); sbi++) {
        VecScatterBegin(vconc_out_scatter, vconc_immobile[sbi], conc_immobile_out[sbi].get_data_petsc(), INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(vconc_out_scatter, vconc_immobile[sbi], conc_immobile_out[sbi].get_data_petsc(), INSERT_VALUES, SCATTER_FORWARD);
    }
}


void DualPorosity::output_data(void )
{
    data_.output_fields.set_time(time_->step(), LimitSide::right);
    if ( data_.output_fields.is_field_output_time(data_.conc_immobile, time().step()) ) {
        output_vector_gather();
    }

    // Register fresh output data
    data_.output_fields.output(time_->step());

    if (time_->tlevel() !=0) {
        // zero_time_step call zero_time_Step of subreactions which performs its own output
        if (reaction_mobile) reaction_mobile->output_data();
        if (reaction_immobile) reaction_immobile->output_data();
    }
}


bool DualPorosity::evaluate_time_constraint(double &time_constraint)
{
    bool cfl_changed = false;
    if (reaction_mobile)
    {
        if (reaction_mobile->evaluate_time_constraint(time_constraint))
            cfl_changed = true;
    }
    if (reaction_immobile)
    {
        double cfl_immobile;
        if (reaction_immobile->evaluate_time_constraint(cfl_immobile))
        {
            time_constraint = std::min(time_constraint, cfl_immobile);
            cfl_changed = true;
        }
    }
    
    return cfl_changed;
}
