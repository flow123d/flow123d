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
 * @file    transport_operator_splitting.cc
 * @brief   
 */

#include <iostream>
#include <iomanip>

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/index_types.hh"

#include "transport/transport_operator_splitting.hh"
#include <petscmat.h>

#include "io/output_time.hh"
#include "tools/time_governor.hh"
#include "coupling/equation.hh"
#include "coupling/balance.hh"
#include "transport/transport.h"
#include "mesh/mesh.h"

#include "reaction/reaction_term.hh"
#include "reaction/first_order_reaction.hh"
#include "reaction/radioactive_decay.hh"
#include "reaction/sorption.hh"
#include "reaction/dual_porosity.hh"

#include "la/distribution.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/factory.hh"
#include "input/flow_attribute_lib.hh"
#include "fem/fe_p.hh"

FLOW123D_FORCE_LINK_IN_CHILD(transportOperatorSplitting)

FLOW123D_FORCE_LINK_IN_PARENT(firstOrderReaction)
FLOW123D_FORCE_LINK_IN_PARENT(radioactiveDecay)
FLOW123D_FORCE_LINK_IN_PARENT(dualPorosity)
FLOW123D_FORCE_LINK_IN_PARENT(sorptionMobile)
FLOW123D_FORCE_LINK_IN_PARENT(sorptionImmobile)
FLOW123D_FORCE_LINK_IN_PARENT(sorption)


using namespace Input::Type;



Abstract & ConcentrationTransportBase::get_input_type() {
	return Abstract("Solute",
			"Transport of soluted  substances.")
			.close();
}


const Record & TransportOperatorSplitting::get_input_type() {
	return Record("Coupling_OperatorSplitting",
            "Transport by convection and/or diffusion\n"
            "coupled with reaction and adsorption model (ODE per element)\n"
            " via operator splitting.")
		.derive_from(AdvectionProcessBase::get_input_type())
		.add_attribute( FlowAttribute::subfields_address(), "\"/problem/solute_equation/substances/*/name\"")
		.declare_key("time", TimeGovernor::get_input_type(), Default::obligatory(),
				"Time governor settings for the transport equation.")
		.declare_key("balance", Balance::get_input_type(), Default("{}"),
				"Settings for computing mass balance.")
		.declare_key("output_stream", OutputTime::get_input_type(), Default("{}"),
				"Output stream settings.\n Specify file format, precision etc.")
		.declare_key("substances", Array( Substance::get_input_type(), 1), Default::obligatory(),
				"Specification of transported substances.")
				// input data
		.declare_key("transport", ConcentrationTransportBase::get_input_type(), Default::obligatory(),
				"Type of the numerical method for the transport equation.")
		.declare_key("reaction_term", ReactionTerm::it_abstract_term(), Default::optional(),
					"Reaction model involved in transport.")
/*
		.declare_key("output_fields", Array(ConvectionTransport::get_output_selection()),
				Default("\"conc\""),
				"List of fields to write to output file.")*/
		.close();
}


const int TransportOperatorSplitting::registrar =
		Input::register_class< TransportOperatorSplitting, Mesh &, const Input::Record>("Coupling_OperatorSplitting") +
		TransportOperatorSplitting::get_input_type().size();







TransportEqData::TransportEqData()
{
    *this += porosity.name("porosity")
            .description("Porosity of the mobile phase.")
            .input_default("1.0")
            .units( UnitSI::dimensionless() )
            .flags_add(in_main_matrix & in_rhs);

    *this += water_content.name("water_content")
            .description("INTERNAL. Water content passed from unsaturated Darcy flow model.")
            .units( UnitSI::dimensionless() )
            .flags_add(input_copy & in_time_term & in_main_matrix & in_rhs);

    *this += cross_section.name("cross_section")
               .flags( FieldFlag::input_copy )
               .flags_add(in_time_term & in_main_matrix & in_rhs);

    *this += sources_density.name("sources_density")
            .description("Density of concentration sources.")
            .input_default("0.0")
            .units( UnitSI().kg().m(-3).s(-1) )
            .flags_add(in_rhs);

    *this += sources_sigma.name("sources_sigma")
            .description("Concentration flux.")
            .input_default("0.0")
            .units( UnitSI().s(-1) )
            .flags_add(in_main_matrix & in_rhs);

    *this += sources_conc.name("sources_conc")
            .description("Concentration sources threshold.")
            .input_default("0.0")
            .units( UnitSI().kg().m(-3) )
            .flags_add(in_rhs);
}







/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



TransportOperatorSplitting::TransportOperatorSplitting(Mesh &init_mesh, const Input::Record in_rec)
: AdvectionProcessBase(init_mesh, in_rec),
//  Semchem_reactions(NULL),
  cfl_convection(numeric_limits<double>::max()),
  cfl_reaction(numeric_limits<double>::max())
{
	START_TIMER("TransportOperatorSpliting");

	Distribution *el_distribution;
	LongIdx *el_4_loc;

	Input::AbstractRecord trans = in_rec.val<Input::AbstractRecord>("transport");
	convection = trans.factory< ConcentrationTransportBase, Mesh &, const Input::Record >(init_mesh, trans);

	time_ = new TimeGovernor(in_rec.val<Input::Record>("time"), TimeMark::none_type, false);
	convection->set_time_governor(time());

	// Initialize list of substances.
	convection->substances().initialize(in_rec.val<Input::Array>("substances"));

	// Initialize output stream.
    convection->set_output_stream(OutputTime::create_output_stream("solute", in_rec.val<Input::Record>("output_stream"), time().get_unit_string()));


    // initialization of balance object

    balance_ = make_shared<Balance>("mass", mesh_);
    balance_->init_from_input(in_rec.val<Input::Record>("balance"), this->time());

    balance_->units(UnitSI().kg(1));
    convection->set_balance_object(balance_);

	convection->initialize(); //




	time_ = new TimeGovernor(in_rec.val<Input::Record>("time"), convection->mark_type());

    this->eq_data_ = &(convection->data());

    convection->get_par_info(el_4_loc, el_distribution);
    Input::Iterator<Input::AbstractRecord> reactions_it = in_rec.find<Input::AbstractRecord>("reaction_term");
	if ( reactions_it ) {
		// TODO: allowed instances in this case are only
		// FirstOrderReaction, RadioactiveDecay, SorptionSimple and DualPorosity
		reaction = (*reactions_it).factory< ReactionTerm, Mesh &, Input::Record >(init_mesh, *reactions_it);

        //initialization of DOF handler
        static MixedPtr<FE_P_disc> fe(0);
        shared_ptr<DOFHandlerMultiDim> dof_handler = make_shared<DOFHandlerMultiDim>(*mesh_);
        shared_ptr<DiscreteSpace> ds = make_shared<EqualOrderDiscreteSpace>( mesh_, fe);
        dof_handler->distribute_dofs(ds);

        reaction->substances(convection->substances())
                    .concentration_matrix(convection->get_concentration_matrix(),
						el_distribution, el_4_loc, convection->get_row_4_el())
				.output_stream(convection->output_stream())
				.set_dh(dof_handler)
				.set_time_governor((TimeGovernor &)convection->time());

		reaction->initialize();

	} else {
		reaction = nullptr;
		//Semchem_reactions = nullptr;
	}
}

TransportOperatorSplitting::~TransportOperatorSplitting()
{
    //delete field_output;
    //if (Semchem_reactions) delete Semchem_reactions;
    delete time_;
}


void TransportOperatorSplitting::initialize()
{
    //coupling - passing fields
  if(reaction)
  if( typeid(*reaction) == typeid(SorptionSimple) ||
          typeid(*reaction) == typeid(DualPorosity)
        )
  {
    reaction->data().set_field("porosity", convection->data()["porosity"]);
  }
}




void TransportOperatorSplitting::output_data(){

        
        START_TIMER("TOS-output data");


        convection->output_data();
        if(reaction) reaction->output_data(); // do not perform write_time_frame
        convection->output_stream()->write_time_frame();

        END_TIMER("TOS-output data");
}


void TransportOperatorSplitting::zero_time_step()
{
    //DebugOut() << "tos ZERO TIME STEP.\n";
    convection->zero_time_step();
    convection->calculate_concentration_matrix();   // due to reading of init_conc in reactions
    if(reaction)
    {
      reaction->zero_time_step();
      reaction->output_data(); // do not perform write_time_frame
    }
    convection->output_stream()->write_time_frame();

}



void TransportOperatorSplitting::update_solution() {

	vector<double> source(convection->n_substances()), region_mass(mesh_->region_db().bulk_size());

    time_->next_time();
    time_->view("TOS");    //show time governor

    convection->set_target_time(time_->t());
    
    START_TIMER("TOS-one step");
    int steps=0;
    while ( convection->time().step().lt(time_->t()) )
    {
        steps++;
	    // one internal step
        // we call evaluate_time_constraint() of convection and reaction separately to
        // make sure that both routines are executed.
        bool cfl_convection_changed =  convection->evaluate_time_constraint(cfl_convection);
        bool cfl_reaction_changed = (reaction?reaction->evaluate_time_constraint(cfl_reaction):0);
        bool cfl_changed = cfl_convection_changed || cfl_reaction_changed;
        
        if (steps == 1 || cfl_changed)
        {
            convection->time().set_upper_constraint(cfl_convection, "Time step constrained by transport CFL condition (including both flow and sources).");
            convection->time().set_upper_constraint(cfl_reaction, "Time step constrained by reaction CFL condition.");
            
            // fix step with new constraint
            convection->time().fix_dt_until_mark();
        
            convection->time().view("Convection");   // write TG only once on change
        }
        
	    convection->update_solution();
        

	    if (balance_->cumulative())
	    {
	    	START_TIMER("TOS-balance");

			// save mass after transport step
	    	for (unsigned int sbi=0; sbi<convection->n_substances(); sbi++)
	    	{
	    		balance_->calculate_mass(convection->get_subst_idx()[sbi], convection->get_solution(sbi), region_mass);
	    		source[sbi] = 0;
	    		for (unsigned int ri=0; ri<mesh_->region_db().bulk_size(); ri++)
	    			source[sbi] -= region_mass[ri];
	    	}

	    	END_TIMER("TOS-balance");
	    }

        if(reaction) {
        	convection->calculate_concentration_matrix();
        	reaction->update_solution();
        	convection->update_after_reactions(true);
        }
        else
        	convection->update_after_reactions(false);

	    //if(Semchem_reactions) Semchem_reactions->update_solution();



	    if (balance_->cumulative())
	    {
	    	START_TIMER("TOS-balance");

	    	for (unsigned int sbi=0; sbi<convection->n_substances(); sbi++)
	    	{
	    		// compute mass difference due to reactions
	    		balance_->calculate_mass(convection->get_subst_idx()[sbi], convection->get_solution(sbi), region_mass);
	    		for (unsigned int ri=0; ri<mesh_->region_db().bulk_size(); ri++)
	    			source[sbi] += region_mass[ri];

	    		// update balance of sources due to reactions
	    		balance_->add_cumulative_source(sbi, source[sbi]);
	    	}

	    	END_TIMER("TOS-balance");
	    }
	}

    LogOut().fmt("CONVECTION: steps: {}\n", steps);
}





void TransportOperatorSplitting::set_velocity_field(std::shared_ptr<FieldFE<3, FieldValue<3>::VectorFixed>> flux_field)
{
	convection->set_velocity_field( flux_field );
}




