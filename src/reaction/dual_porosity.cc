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
#include "mesh/accessors.hh"
#include "fields/field_fe.hh"
#include "fields/fe_value_handler.hh"

#include "reaction/sorption.hh"
#include "reaction/first_order_reaction.hh"
#include "reaction/radioactive_decay.hh"
#include "input/factory.hh"

FLOW123D_FORCE_LINK_IN_CHILD(dualPorosity)


using namespace Input::Type;



const Record & DualPorosity::get_input_type() {
    return Record("DualPorosity",
            "Dual porosity model in transport problems.\n"
            "Provides computing the concentration of substances in mobile and immobile zone.\n"
            )
		.derive_from(ReactionTerm::it_abstract_term())
		.declare_key("input_fields", Array(DualPorosity::EqFields().make_field_descriptor_type("DualPorosity")), Default::obligatory(),
						"Containes region specific data necessary to construct dual porosity model.")
		.declare_key("scheme_tolerance", Double(0.0), Default("1e-3"),
					 "Tolerance according to which the explicit Euler scheme is used or not."
					 "Set 0.0 to use analytic formula only (can be slower).")

		.declare_key("reaction_mobile", ReactionTerm::it_abstract_mobile_term(), Default::optional(), "Reaction model in mobile zone.")
		.declare_key("reaction_immobile", ReactionTerm::it_abstract_immobile_term(), Default::optional(), "Reaction model in immobile zone.")
		.declare_key("output",
		                    EqFields().output_fields.make_output_type("DualPorosity", ""),
		                    IT::Default("{ \"fields\": [ \"conc_immobile\" ] }"),
		                    "Setting of the fields output.")
		.close();
}
    
const int DualPorosity::registrar =
		Input::register_class< DualPorosity, Mesh &, Input::Record >("DualPorosity") +
		DualPorosity::get_input_type().size();

DualPorosity::EqFields::EqFields()
: ReactionTerm::EqFields()
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
          .units( UnitSI::dimensionless() )
		  .set_limits(0.0);

  *this += init_conc_immobile
          .name("init_conc_immobile")
          .description("Initial concentration of substances in the immobile zone.")
          .units( UnitSI().kg().m(-3) );

  //creating field for porosity that is set later from the governing equation (transport)
  *this +=porosity
        .name("porosity")
        .description("Concentration solution in the mobile phase.")
        .units( UnitSI::dimensionless() )
        .flags( FieldFlag::input_copy )
		.set_limits(0.0);

  *this += conc_immobile
          .name("conc_immobile")
          .units( UnitSI().kg().m(-3) )
          .flags(FieldFlag::equation_result);

  output_fields += *this;

}

DualPorosity::EqData::EqData()
{}

DualPorosity::DualPorosity(Mesh &init_mesh, Input::Record in_rec)
	: ReactionTerm(init_mesh, in_rec)
{
    eq_fields_ = std::make_shared<EqFields>();
    eq_fields_->add_coords_field();
    //set pointer to equation data fieldset
    this->eq_fieldset_ = eq_fields_.get();
    this->eq_fields_base_ = std::static_pointer_cast<ReactionTerm::EqFields>(eq_fields_);

    //reads input and creates possibly other reaction terms
    make_reactions();
    //read value from input
    scheme_tolerance_ = input_record_.val<double>("scheme_tolerance");
}

DualPorosity::~DualPorosity(void)
{
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
  ASSERT(time_ != nullptr).error("Time governor has not been set yet.\n");
  ASSERT_LT(0, eq_data_->substances_.size()).error("No substances for rection term.\n");
  ASSERT(output_stream_ != nullptr).error("Null output stream.\n");
  
  initialize_fields();

  if(reaction_mobile)
  {
    reaction_mobile->substances(eq_data_->substances_)
                .output_stream(output_stream_)
                .concentration_fields(eq_fields_->conc_mobile_fe)
                .set_time_governor(*time_);
    reaction_mobile->initialize();
  }

  if(reaction_immobile)
  {
    reaction_immobile->substances(eq_data_->substances_)
                .output_stream(output_stream_)
                .concentration_fields(eq_fields_->conc_immobile_fe)
                .set_time_governor(*time_);
    reaction_immobile->initialize();
  }

}

void DualPorosity::initialize_fields()
{
  eq_fields_->set_components(eq_data_->substances_.names());
  //setting fields that are set from input file
  input_field_set_ += *eq_fields_;
  input_field_set_.set_input_list(input_record_.val<Input::Array>("input_fields"), *time_);

  //setting fields in data
  eq_fields_->set_mesh(*mesh_);
  
  //initialization of output
  eq_fields_->output_fields.set_components(eq_data_->substances_.names());
  eq_fields_->output_fields.set_mesh(*mesh_);
  eq_fields_->output_fields.output_type(OutputTime::ELEM_DATA);
  eq_fields_->conc_immobile.setup_components();


  //creating field fe and output multifield for sorbed concentrations
  eq_fields_->conc_immobile_fe.resize(eq_data_->substances_.size());
  for (unsigned int sbi = 0; sbi < eq_data_->substances_.size(); sbi++)
  {
      // create shared pointer to a FieldFE and push this Field to output_field on all regions
	  eq_fields_->conc_immobile_fe[sbi] = create_field_fe< 3, FieldValue<3>::Scalar >(eq_data_->dof_handler_);
	  eq_fields_->conc_immobile[sbi].set(eq_fields_->conc_immobile_fe[sbi], 0);
  }

  eq_fields_->output_fields.initialize(output_stream_, mesh_, input_record_.val<Input::Record>("output"),time());
}


void DualPorosity::zero_time_step()
{
  ASSERT(time_ != nullptr).error("Time governor has not been set yet.\n");
  ASSERT_LT(0, eq_data_->substances_.size()).error("No substances for rection term.\n");
  ASSERT(output_stream_ != nullptr).error("Null output stream.\n");
  
  //coupling - passing fields
  if(reaction_mobile)
  if (typeid(*reaction_mobile) == typeid(SorptionMob))
  {
	  reaction_mobile->eq_fieldset().set_field("porosity", (*eq_fields_)["porosity"]);
	  reaction_mobile->eq_fieldset().set_field("porosity_immobile", (*eq_fields_)["porosity_immobile"]);
  }
  if(reaction_immobile)
  if (typeid(*reaction_immobile) == typeid(SorptionImmob))
  {
	  reaction_immobile->eq_fieldset().set_field("porosity", (*eq_fields_)["porosity"]);
	  reaction_immobile->eq_fieldset().set_field("porosity_immobile", (*eq_fields_)["porosity_immobile"]);
  }
  
  eq_fields_->set_time(time_->step(0), LimitSide::right);
  std::stringstream ss; // print warning message with table of uninitialized fields
  if ( FieldCommon::print_message_table(ss, "dual porosity") ) {
      WarningOut() << ss.str();
  }
  set_initial_condition();

  output_data();
  
  if(reaction_mobile)
    reaction_mobile->zero_time_step();

  if(reaction_immobile)
    reaction_immobile->zero_time_step();

}

void DualPorosity::set_initial_condition()
{
    for ( DHCellAccessor dh_cell : eq_data_->dof_handler_->own_range() ) {
        IntIdx dof_p0 = dh_cell.get_loc_dof_indices()[0];
        const ElementAccessor<3> ele = dh_cell.elm();

        //setting initial solid concentration for substances involved in adsorption
        for (unsigned int sbi = 0; sbi < eq_data_->substances_.size(); sbi++)
        {
            eq_fields_->conc_immobile_fe[sbi]->vec().set( dof_p0, eq_fields_->init_conc_immobile[sbi].value(ele.centre(), ele) );
        }
    }
}

void DualPorosity::update_solution(void) 
{
  eq_fields_->set_time(time_->step(-2), LimitSide::right);
 
  START_TIMER("dual_por_exchange_step");
  for ( DHCellAccessor dh_cell : eq_data_->dof_handler_->own_range() )
  {
      compute_reaction(dh_cell);
  }
  END_TIMER("dual_por_exchange_step");
  
  if(reaction_mobile) reaction_mobile->update_solution();
  if(reaction_immobile) reaction_immobile->update_solution();
}


void DualPorosity::compute_reaction(const DHCellAccessor& dh_cell)
{
    unsigned int sbi;
    double conc_average, // weighted (by porosity) average of concentration
          conc_mob, conc_immob,  // new mobile and immobile concentration
          previous_conc_mob, previous_conc_immob, // mobile and immobile concentration in previous time step
          conc_max, //difference between concentration and average concentration
          por_mob, por_immob; // mobile and immobile porosity
    
    // get data from fields
    ElementAccessor<3> ele = dh_cell.elm();
    IntIdx dof_p0 = dh_cell.get_loc_dof_indices()[0];
    por_mob = eq_fields_->porosity.value(ele.centre(),ele);
    por_immob = eq_fields_->porosity_immobile.value(ele.centre(),ele);
    arma::Col<double> diff_vec(eq_data_->substances_.size());
    for (sbi=0; sbi<eq_data_->substances_.size(); sbi++) // Optimize: SWAP LOOPS
        diff_vec[sbi] = eq_fields_->diffusion_rate_immobile[sbi].value(ele.centre(), ele);
 
    // if porosity_immobile == 0 then mobile concentration stays the same 
    // and immobile concentration cannot change
    if (por_immob == 0.0) return;
    
    double exponent,
           temp_exponent = (por_mob + por_immob) / (por_mob * por_immob) * time_->dt();
  
    for (sbi = 0; sbi < eq_data_->substances_.size(); sbi++) //over all substances
    {
        exponent = diff_vec[sbi] * temp_exponent;
        //previous values
        previous_conc_mob = eq_fields_->conc_mobile_fe[sbi]->vec().get(dof_p0);
        previous_conc_immob = eq_fields_->conc_immobile_fe[sbi]->vec().get(dof_p0);
        
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
        
        eq_fields_->conc_mobile_fe[sbi]->vec().set(dof_p0, conc_mob);
        eq_fields_->conc_immobile_fe[sbi]->vec().set(dof_p0, conc_immob);
    }
}


void DualPorosity::output_data(void )
{
    eq_fields_->output_fields.set_time(time_->step(), LimitSide::right);

    // Register fresh output data
    eq_fields_->output_fields.output(time_->step());

    if (time_->tlevel() !=0) {
        // zero_time_step call zero_time_Step of subreactions which performs its own output
        if (reaction_mobile) reaction_mobile->output_data();
        if (reaction_immobile) reaction_immobile->output_data();
    }
}
