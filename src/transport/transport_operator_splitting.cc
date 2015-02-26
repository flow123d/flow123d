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
#include "system/sys_vector.hh"
#include "coupling/time_governor.hh"
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
	.declare_key("output_stream", OutputTime::input_type, Default::obligatory(),
			"Parameters of output stream.")
	.declare_key("mass_balance", MassBalance::input_type, Default::optional(), "Settings for computing mass balance.");


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
  mh_dh(nullptr),
  mass_balance_(nullptr)
{
}

TransportBase::~TransportBase()
{
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



TransportOperatorSplitting::TransportOperatorSplitting(Mesh &init_mesh, const Input::Record &in_rec)
: TransportBase(init_mesh, in_rec),
  convection(NULL),
  reaction(nullptr),
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
		if (reactions_it->type() == FirstOrderReaction::input_type ) {
			reaction =  new FirstOrderReaction(init_mesh, *reactions_it);
		} else
		if (reactions_it->type() == RadioactiveDecay::input_type) {
			reaction = new RadioactiveDecay(init_mesh, *reactions_it);
		} else
		if (reactions_it->type() == SorptionSimple::input_type ) {
			reaction =  new SorptionSimple(init_mesh, *reactions_it);
		} else
		if (reactions_it->type() == DualPorosity::input_type ) {
			reaction =  new DualPorosity(init_mesh, *reactions_it);
		} else
		if (reactions_it->type() == Semchem_interface::input_type ) {
// 			Semchem_reactions = new Semchem_interface(0.0, mesh_, n_subst_, false); //false instead of convection->get_dual_porosity
// 			Semchem_reactions->set_el_4_loc(el_4_loc);
//                 //Semchem works with phases 0-3; this is not supported no more!
//                 semchem_conc_ptr = new double**[1];
//                 semchem_conc_ptr[0] = convection->get_concentration_matrix();
//                 Semchem_reactions->set_concentration_matrix(semchem_conc_ptr, el_distribution, el_4_loc);
            THROW( ReactionTerm::ExcWrongDescendantModel() 
                << ReactionTerm::EI_Model((*reactions_it).type().type_name())
                << EI_Message("This model is not currently supported!") 
                << (*reactions_it).ei_address());

		} else {
			//This point cannot be reached. The TYPE_selection will throw an error first. 
            THROW( ExcMessage() 
                << EI_Message("Descending model type selection failed (SHOULD NEVER HAPPEN).") 
                << (*reactions_it).ei_address());
		}
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
}

TransportOperatorSplitting::~TransportOperatorSplitting()
{
    //delete field_output;
    delete convection;
    if (reaction) delete reaction;
    if (Semchem_reactions) delete Semchem_reactions;
    delete time_;
}




void TransportOperatorSplitting::output_data(){

        
        START_TIMER("TOS-output data");

        time_->view("TOS");    //show time governor
        convection->output_data();
        if(reaction) reaction->output_data(); // do not perform write_time_frame
        convection->output_stream_->write_time_frame();

}


void TransportOperatorSplitting::zero_time_step()
{
  
    convection->zero_time_step();
    if(reaction) reaction->zero_time_step();
    convection->output_stream_->write_time_frame();

}



void TransportOperatorSplitting::update_solution() {



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
            if(reaction) reaction->update_solution();
	    if(Semchem_reactions) Semchem_reactions->update_solution();
	    if (convection->mass_balance() != NULL)
	    	convection->mass_balance()->calculate(convection->time().t());

	}
    END_TIMER("TOS-one step");


    
    xprintf( MsgLog, "    CONVECTION: steps: %d\n",steps);
}





void TransportOperatorSplitting::set_velocity_field(const MH_DofHandler &dh)
{
    mh_dh = &dh;
	convection->set_velocity_field( dh );
};


void TransportOperatorSplitting::calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance)
{}

void TransportOperatorSplitting::calc_elem_sources(vector<vector<double> > &mass, vector<vector<double> > &src_balance)
{}





