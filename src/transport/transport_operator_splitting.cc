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

#include "io/output_time.hh"
#include "tools/time_governor.hh"
#include "system/sys_vector.hh"
#include "coupling/equation.hh"
#include "coupling/balance.hh"
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
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/factory.hh"

FLOW123D_FORCE_LINK_IN_CHILD(transportOperatorSplitting);

FLOW123D_FORCE_LINK_IN_PARENT(firstOrderReaction);
FLOW123D_FORCE_LINK_IN_PARENT(radioactiveDecay);
FLOW123D_FORCE_LINK_IN_PARENT(dualPorosity);
FLOW123D_FORCE_LINK_IN_PARENT(sorptionMobile);
FLOW123D_FORCE_LINK_IN_PARENT(sorptionImmobile);
FLOW123D_FORCE_LINK_IN_PARENT(sorption);


using namespace Input::Type;

AbstractRecord & AdvectionProcessBase::get_input_type() {
	return AbstractRecord("Transport",
			"Secondary equation for transport of substances.")
			.close();
}


const Record & TransportOperatorSplitting::get_input_type() {
	return Record("TransportOperatorSplitting",
            "Explicit FVM transport (no diffusion)\n"
            "coupled with reaction and adsorption model (ODE per element)\n"
            " via operator splitting.")
		.derive_from(AdvectionProcessBase::get_input_type())
		.declare_key("time", TimeGovernor::get_input_type(), Default::obligatory(),
				"Time governor setting for the secondary equation.")
		.declare_key("balance", Balance::get_input_type(), Default::obligatory(),
				"Settings for computing balance.")
		.declare_key("output_stream", OutputTime::get_input_type(), Default::obligatory(),
				"Parameters of output stream.")
		.declare_key("substances", Array( Substance::get_input_type() ), Default::obligatory(),
				"Specification of transported substances.")
				// input data
		.declare_key("reaction_term", ReactionTerm::get_input_type(), Default::optional(),
					"Reaction model involved in transport.")

		.declare_key("input_fields", Array(
				Record("TransportOperatorSplitting_Data", FieldCommon::field_descriptor_record_decsription("TransportOperatorSplitting_Data") )
				.copy_keys( ConvectionTransport::EqData().make_field_descriptor_type("TransportOperatorSplitting") )
				.declare_key(OldBcdInput::transport_old_bcd_file_key(), IT::FileName::input(), "File with mesh dependent boundary conditions (obsolete).")
				.close()
				), IT::Default::obligatory(), "")
		.declare_key("output_fields", Array(ConvectionTransport::EqData::get_output_selection()),
				Default("conc"),
				"List of fields to write to output file.")
		.close();
}


const int TransportOperatorSplitting::registrar =
		Input::register_class< TransportOperatorSplitting, Mesh &, const Input::Record & >("TransportOperatorSplitting") +
		TransportOperatorSplitting::get_input_type().size();







TransportBase::TransportEqData::TransportEqData()
{

	ADD_FIELD(porosity, "Mobile porosity", "1");
	porosity.units( UnitSI::dimensionless() ).flags_add(in_time_term & in_main_matrix & in_rhs);

	ADD_FIELD(cross_section, "");
	cross_section.flags( FieldFlag::input_copy ).flags_add(in_time_term & in_main_matrix & in_rhs);

	ADD_FIELD(sources_density, "Density of concentration sources.", "0");
	sources_density.units( UnitSI().kg().m(-3).s(-1) )
			.flags_add(in_rhs);

	ADD_FIELD(sources_sigma, "Concentration flux.", "0");
	sources_sigma.units( UnitSI().s(-1) )
			.flags_add(in_main_matrix & in_rhs);

	ADD_FIELD(sources_conc, "Concentration sources threshold.", "0");
	sources_conc.units( UnitSI().kg().m(-3) )
			.flags_add(in_rhs);
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

		reaction->substances(substances_)
                    .concentration_matrix(convection->get_concentration_matrix(),
						el_distribution, el_4_loc, convection->get_row_4_el())
				.output_stream(convection->output_stream())
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
    time_->view("TOS");    //show time governor

    convection->set_target_time(time_->t());
    convection->time_->estimate_dt();
    convection->time_->view("Convection");    //show time governor

    START_TIMER("TOS-one step");
    int steps=0;
    while ( convection->time().step().lt(time_->t()) )
    {
        steps++;
	    // one internal step
	    convection->update_solution();

	    if (balance_ != nullptr && balance_->cumulative())
	    {
	    	START_TIMER("TOS-balance");

			// save mass after transport step
	    	for (unsigned int sbi=0; sbi<n_substances(); sbi++)
	    	{
	    		balance_->calculate_mass(convection->get_subst_idx()[sbi], convection->get_concentration_vector()[sbi], region_mass);
	    		source[sbi] = 0;
	    		for (unsigned int ri=0; ri<mesh_->region_db().bulk_size(); ri++)
	    			source[sbi] -= region_mass[ri];
	    	}

	    	convection->calculate_cumulative_balance();

	    	END_TIMER("TOS-balance");
	    }

        if(reaction) reaction->update_solution();
	    if(Semchem_reactions) Semchem_reactions->update_solution();

	    if (balance_ != nullptr && balance_->cumulative())
	    {
	    	START_TIMER("TOS-balance");

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




