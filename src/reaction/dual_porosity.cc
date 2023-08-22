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

#include "reaction/sorption.hh"
#include "reaction/first_order_reaction.hh"
#include "reaction/radioactive_decay.hh"
#include "reaction/assembly_reaction.hh"
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
          .input_default("0")
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
: ReactionTerm::EqData() {}

DualPorosity::DualPorosity(Mesh &init_mesh, Input::Record in_rec)
	: ReactionTerm(init_mesh, in_rec),
	  init_condition_assembly_(nullptr),
	  reaction_assembly_(nullptr)
{
    eq_fields_ = std::make_shared<EqFields>();
    eq_fields_->add_coords_field();
    //set pointer to equation data fieldset
    this->eq_fieldset_ = eq_fields_;
    this->eq_fields_base_ = std::static_pointer_cast<ReactionTerm::EqFields>(eq_fields_);

    eq_data_ = std::make_shared<EqData>();
    this->eq_data_base_ = std::static_pointer_cast<ReactionTerm::EqData>(eq_data_);

    //reads input and creates possibly other reaction terms
    make_reactions();
    //read value from input
    eq_data_->scheme_tolerance_ = input_record_.val<double>("scheme_tolerance");
}

DualPorosity::~DualPorosity(void)
{
    if (init_condition_assembly_!=nullptr) delete init_condition_assembly_;
    if (reaction_assembly_!=nullptr) delete reaction_assembly_;
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
  ASSERT_PERMANENT(time_ != nullptr).error("Time governor has not been set yet.\n");
  ASSERT_PERMANENT_LT(0, eq_data_->substances_.size()).error("No substances for rection term.\n");
  ASSERT_PERMANENT(output_stream_ != nullptr).error("Null output stream.\n");
  
  eq_data_->time_ = this->time_;

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

  init_condition_assembly_ = new GenericAssembly< InitConditionAssemblyDp >(eq_fields_.get(), eq_data_.get());
  reaction_assembly_ = new GenericAssembly< ReactionAssemblyDp >(eq_fields_.get(), eq_data_.get());

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
  init_condition_assembly_->assemble(eq_data_->dof_handler_);

  output_data();
  
  if(reaction_mobile)
    reaction_mobile->zero_time_step();

  if(reaction_immobile)
    reaction_immobile->zero_time_step();

}

void DualPorosity::update_solution(void) 
{
  eq_fields_->set_time(time_->step(-2), LimitSide::right);
 
  START_TIMER("dual_por_exchange_step");
  reaction_assembly_->assemble(eq_data_->dof_handler_);
  END_TIMER("dual_por_exchange_step");
  
  if(reaction_mobile) reaction_mobile->update_solution();
  if(reaction_immobile) reaction_immobile->update_solution();
}


void DualPorosity::compute_reaction(FMT_UNUSED const DHCellAccessor& dh_cell)
{
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
