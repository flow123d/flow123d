/*
 * transport_operator_splitting.cc
 *
 *  Created on: May 21, 2011
 *      Author: jiri
 */

#include <iostream>
#include <iomanip>

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/xio.h"

#include "transport/transport_operator_splitting.hh"
#include <petscmat.h>

#include "tools/time_governor.hh"
#include "system/sys_vector.hh"
#include "coupling/equation.hh"
#include "transport/transport.h"
#include "mesh/mesh.h"
#include "flow/old_bcd.hh"

#include "reaction/reaction_term.hh"
#include "reaction/first_order_reaction.hh"
#include "reaction/radioactive_decay.hh"
#include "reaction/sorption.hh"
#include "reaction/dual_porosity.hh"

#include "semchem/semchem_interface.hh"

#include "la/distribution.hh"
#include "io/output_data.hh"

#include "input/input_type.hh"
#include "input/accessors.hh"

using namespace Input::Type;

AbstractRecord AdvectionProcessBase::input_type
	= AbstractRecord("Transport", "Secondary equation for transport of substances.")
	.declare_key("time", TimeGovernor::input_type, Default::obligatory(),
			"Time governor setting for the secondary equation.")
	.declare_key("balance", Balance::input_type, Default::obligatory(),
			"Settings for computing balance.")
	.declare_key("output_stream", OutputTime::input_type, Default::obligatory(),
			"Parameters of output stream.");


Record TransportBase::input_type_output_record
	= Record("TransportOutput", "Output setting for transport equations.")
	.declare_key("output_stream", OutputTime::input_type, Default::obligatory(),
			"Parameters of output stream.");


Record TransportOperatorSplitting::input_type
	= Record("TransportOperatorSplitting",
            "Explicit FVM transport (no diffusion)\n"
            "coupled with reaction and adsorption model (ODE per element)\n"
            " via operator splitting.")
    .derive_from(AdvectionProcessBase::input_type)
    .declare_key("substances", Array(Substance::input_type), Default::obligatory(),
    		"Specification of transported substances.")
    	    // input data
    .declare_key("reaction_term", ReactionTerm::input_type, Default::optional(),
                "Reaction model involved in transport.")

    .declare_key("input_fields", Array(
    		ConvectionTransport::EqData().make_field_descriptor_type("TransportOperatorSplitting")
    		.declare_key(OldBcdInput::transport_old_bcd_file_key(), IT::FileName::input(), "File with mesh dependent boundary conditions (obsolete).")
    		), IT::Default::obligatory(), "")
    .declare_key("output_fields", Array(ConvectionTransport::EqData::output_selection),
    		Default("conc"),
       		"List of fields to write to output file.");







TransportBase::TransportEqData::TransportEqData()
{

	ADD_FIELD(porosity, "Mobile porosity", "1");
	porosity.units( UnitSI::dimensionless() );

	ADD_FIELD(cross_section, "");
	cross_section.flags( FieldFlag::input_copy );

	ADD_FIELD(sources_density, "Density of concentration sources.", "0");
	sources_density.units( UnitSI().kg().m(-3).s(-1) );

	ADD_FIELD(sources_sigma, "Concentration flux.", "0");
	sources_sigma.units( UnitSI().s(-1) );

	ADD_FIELD(sources_conc, "Concentration sources threshold.", "0");
	sources_conc.units( UnitSI().kg().m(-3) );
}


TransportBase::TransportBase(Mesh &mesh, const Input::Record in_rec)
: AdvectionProcessBase(mesh, in_rec ),
  mh_dh(nullptr)
{
}

TransportBase::~TransportBase()
{
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



TransportOperatorSplitting::TransportOperatorSplitting(Mesh &init_mesh, const Input::Record &in_rec)
: TransportBase(init_mesh, in_rec),
  convection(NULL),
  Semchem_reactions(NULL)
{
	START_TIMER("TransportOperatorSpliting");

	Distribution *el_distribution;
	int *el_4_loc;

	// Initialize list of substances.
	substances_.initialize(in_rec.val<Input::Array>("substances"));
    n_subst_ = substances_.size();

	convection = new ConvectionTransport(*mesh_, in_rec);
	this->eq_data_ = &(convection->data());

    time_ = new TimeGovernor(in_rec.val<Input::Record>("time"), convection->mark_type() );


    convection->get_par_info(el_4_loc, el_distribution);
    Input::Iterator<Input::AbstractRecord> reactions_it = in_rec.find<Input::AbstractRecord>("reaction_term");
	if ( reactions_it ) {
		// TODO: allowed instances in this case are only
		// FirstOrderReaction, RadioactiveDecay, SorptionSimple and DualPorosity
		reaction = (*reactions_it).factory< ReactionTerm, Mesh &, Input::Record >(init_mesh, *reactions_it);
		if (reactions_it->type() == FirstOrderReaction::input_type ) {}
		else if (reactions_it->type() == RadioactiveDecay::input_type) {}
		else if (reactions_it->type() == DualPorosity::input_type ) {}
		else if (reactions_it->type() == SorptionMob::input_type ) {}
		else if (reactions_it->type() == SorptionImmob::input_type ) {}
		else if (reactions_it->type() == SorptionSimple::input_type ) {}
		//temporary, until new mass balance considering reaction term is created
		xprintf(Warn, "The mass balance is not computed correctly when reaction term is present. "
					  "Only the mass flux over boundaries is correct.\n");

		reaction->substances(substances_)
                    .concentration_matrix(convection->get_concentration_matrix(),
						el_distribution, el_4_loc, convection->get_row_4_el())
				.output_stream(*(convection->output_stream()))
				.set_time_governor(*(convection->time_));

		reaction->initialize();

	} else {
		reaction = nullptr;
		Semchem_reactions = nullptr;
	}
        
  //coupling - passing fields
  if(reaction)
  if( typeid(*reaction) == typeid(SorptionSimple) ||
		  typeid(*reaction) == typeid(DualPorosity)
		)
  {
	reaction->data().set_field("porosity", convection->data()["porosity"]);
  }

  // initialization of balance object
  Input::Iterator<Input::Record> it = in_rec.find<Input::Record>("balance");
  if (it->val<bool>("balance_on"))
  {
//	  convection->get_par_info(el_4_loc, el_distribution);

	  balance_ = boost::make_shared<Balance>("mass", mesh_, el_distribution, el_4_loc, *it);

	  convection->set_balance_object(balance_);

	  balance_->allocate(el_distribution->lsize(), 1);
  }
}

TransportOperatorSplitting::~TransportOperatorSplitting()
{
    //delete field_output;
    delete convection;
    if (Semchem_reactions) delete Semchem_reactions;
    delete time_;
}




void TransportOperatorSplitting::output_data(){

        
        START_TIMER("TOS-output data");

        time_->view("TOS");    //show time governor
        convection->output_data();
        if(reaction) reaction->output_data(); // do not perform write_time_frame
        convection->output_stream_->write_time_frame();

        if (balance_ != nullptr && time_->is_current( time_->marks().type_output() ))
        {
        	START_TIMER("TOS-balance");
        	convection->calculate_instant_balance();
        	balance_->output(time_->t());
        	END_TIMER("TOS-balance");
        }

        END_TIMER("TOS-output data");
}


void TransportOperatorSplitting::zero_time_step()
{
  
    convection->zero_time_step();
    if(reaction) reaction->zero_time_step();
    convection->output_stream_->write_time_frame();
    if (balance_ != nullptr)
    {
    	balance_->units(
    	        convection->data_.cross_section.units()*UnitSI().md(1)
    	        *convection->data_.porosity.units()
    	        *convection->data_.conc_mobile.units());
    	balance_->output(time_->t());
    }

}



void TransportOperatorSplitting::update_solution() {

	vector<double> source(n_substances()), region_mass(mesh_->region_db().bulk_size());

    time_->next_time();

    convection->set_target_time(time_->t());
    convection->time_->estimate_dt();

    START_TIMER("TOS-one step");
    int steps=0;
    while ( convection->time().lt(time_->t()) )
    {
        steps++;
	    // one internal step
	    convection->update_solution();

	    if (balance_ != nullptr && balance_->cumulative())
	    {
			// save mass after transport step
	    	for (unsigned int sbi=0; sbi<n_substances(); sbi++)
	    	{
	    		balance_->calculate_mass(convection->get_subst_idx()[sbi], convection->get_concentration_vector()[sbi], region_mass);
	    		source[sbi] = 0;
	    		for (unsigned int ri=0; ri<mesh_->region_db().bulk_size(); ri++)
	    			source[sbi] -= region_mass[ri];
	    	}
	    }

        if(reaction) reaction->update_solution();
	    if(Semchem_reactions) Semchem_reactions->update_solution();

//	    if (convection->mass_balance() != NULL)
//	    	convection->mass_balance()->calculate(convection->time().t());

	    if (balance_ != nullptr && balance_->cumulative())
	    {
	    	START_TIMER("TOS-balance");
	    	convection->calculate_cumulative_balance();

	    	for (unsigned int sbi=0; sbi<n_substances(); sbi++)
	    	{
	    		// compute mass difference due to reactions
	    		balance_->calculate_mass(convection->get_subst_idx()[sbi], convection->get_concentration_vector()[sbi], region_mass);
	    		for (unsigned int ri=0; ri<mesh_->region_db().bulk_size(); ri++)
	    			source[sbi] += region_mass[ri];

	    		// update balance of sources due to reactions
	    		balance_->add_cumulative_source(sbi, source[sbi]);
	    	}

	    	END_TIMER("TOS-balance");
	    }
	}

    xprintf( MsgLog, "    CONVECTION: steps: %d\n",steps);
}





void TransportOperatorSplitting::set_velocity_field(const MH_DofHandler &dh)
{
    mh_dh = &dh;
	convection->set_velocity_field( dh );
};







