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

#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "reaction/sorption_base.hh"
#include "reaction/sorption.hh"
#include "reaction/dual_por_exchange.hh"

#include "semchem/semchem_interface.hh"

#include "la/distribution.hh"
#include "io/output.h"

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
    .declare_key("substances", Array(String()), Default::obligatory(),
    		"Names of transported substances.")
    	    // input data
    .declare_key("reaction_term", ReactionTerm::input_type, Default::optional(),
                "Reaction model involved in transport.")

    .declare_key("data", Array(
    		ConvectionTransport::EqData().make_field_descriptor_type("TransportOperatorSplitting")
    		.declare_key(OldBcdInput::transport_old_bcd_file_key(), IT::FileName::input(), "File with mesh dependent boundary conditions (obsolete).")
    		), IT::Default::obligatory(), "")
    .declare_key("output_fields", Array(ConvectionTransport::EqData::output_selection),
    		Default("conc"),
       		"List of fields to write to output file.");


TransportBase::TransportEqData::TransportEqData()
{

	ADD_FIELD(porosity, "Mobile porosity", "1");

	ADD_FIELD(sources_density, "Density of concentration sources.", "0");
	ADD_FIELD(sources_sigma, "Concentration flux.", "0");
	ADD_FIELD(sources_conc, "Concentration sources threshold.", "0");

}


TransportBase::TransportBase(Mesh &mesh, const Input::Record in_rec)
: AdvectionProcessBase(mesh, in_rec ),
  mh_dh(NULL),
  mass_balance_(NULL)
{
}

TransportBase::~TransportBase()
{
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



TransportOperatorSplitting::TransportOperatorSplitting(Mesh &init_mesh, const Input::Record &in_rec)
: TransportBase(init_mesh, in_rec),
  reaction(nullptr),
  convection(NULL),
  Semchem_reactions(NULL)
{
	Distribution *el_distribution;
	int *el_4_loc;

    // double problem_save_step = OptGetDbl("Global", "Save_step", "1.0");

    in_rec.val<Input::Array>("substances").copy_to(subst_names_);
    n_subst_ = subst_names_.size();
	convection = new ConvectionTransport(*mesh_, in_rec);

	output_mark_type = convection->mark_type() | TimeGovernor::marks().type_fixed_time() | TimeGovernor::marks().type_output();
    time_ = new TimeGovernor(in_rec.val<Input::Record>("time"), output_mark_type );


    convection->get_par_info(el_4_loc, el_distribution);
    Input::Iterator<Input::AbstractRecord> reactions_it = in_rec.find<Input::AbstractRecord>("reaction_term");
        if ( reactions_it ) {
            if (reactions_it->type() == Linear_reaction::input_type ) {
                reaction =  new Linear_reaction(init_mesh, *reactions_it, subst_names_);
                
            } else
            if (reactions_it->type() == Pade_approximant::input_type) {
                reaction = new Pade_approximant(init_mesh, *reactions_it, subst_names_ );
              
            } else
            if (reactions_it->type() == SorptionSimple::input_type ) {
                reaction =  new SorptionSimple(init_mesh, *reactions_it, subst_names_);
                
                static_cast<SorptionSimple *> (reaction) -> init_from_input(*reactions_it);
                static_cast<SorptionSimple *> (reaction) -> set_porosity(convection->get_data()->porosity);
                
            } else
            if (reactions_it->type() == DualPorosity::input_type ) {
                reaction =  new DualPorosity(init_mesh, *reactions_it, subst_names_);
                
                static_cast<DualPorosity *> (reaction) -> set_porosity(convection->get_data()->porosity);
                
            } else
            if (reactions_it->type() == Semchem_interface::input_type ) {
                Semchem_reactions = new Semchem_interface(0.0, mesh_, n_subst_, false); //false instead of convection->get_dual_porosity
                Semchem_reactions->set_el_4_loc(el_4_loc);
                Semchem_reactions->set_concentration_matrix(convection->get_concentration_matrix(), el_distribution, el_4_loc);

            } else {
                xprintf(UsrErr, "Wrong reaction type.\n");
            }
            //temporary, until new mass balance considering reaction term is created
            xprintf(Warn, "The mass balance is not computed correctly when reaction term is present. "
                          "Only the mass flux over boundaries is correct.\n");
            reaction->set_time_governor(*(convection->time_));
            reaction->set_concentration_matrix(convection->get_concentration_matrix()[MOBILE], el_distribution, el_4_loc, convection->get_row_4_el());
            reaction->initialize(convection->output_stream());
            
        } else {
            reaction = nullptr;
            Semchem_reactions = nullptr;
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

    if (time_->is_current(output_mark_type)) {
        
        START_TIMER("TOS-output data");
        DBGMSG("\nTOS: output time: %f\n", time_->t());

        convection->output_data();
        if(reaction) reaction->output_data();
    }
}



void TransportOperatorSplitting::update_solution() {

	static bool first_time_call = true;

	if (first_time_call)
	{
		if (convection->mass_balance() != NULL)
			convection->mass_balance()->output(convection->time().t());
		first_time_call = false;
	}

    time_->next_time();
#ifdef DEBUG_MESSAGES
    time_->view("TOS");    //show time governor
#endif
    
    convection->set_target_time(time_->t());
    convection->time_->estimate_dt();
        
    xprintf( Msg, "TOS: time: %f        CONVECTION: time: %f      dt_estimate: %f\n", 
             time_->t(), convection->time().t(), convection->time().estimate_dt() );
    
    START_TIMER("TOS-one step");
    int steps=0;
    while ( convection->time().lt(time_->t()) )
    {
        steps++;
	    // one internal step
	    convection->compute_one_step();
            if(reaction) reaction->update_solution();
	    if(Semchem_reactions) Semchem_reactions->update_solution();
	    if (convection->mass_balance() != NULL)
	    	convection->mass_balance()->calculate(convection->time().t());
	}
    END_TIMER("TOS-one step");
    
    xprintf( Msg, "CONVECTION: steps: %d\n",steps);
}





void TransportOperatorSplitting::set_velocity_field(const MH_DofHandler &dh)
{
    mh_dh = &dh;
	convection->set_velocity_field( dh );
};



void TransportOperatorSplitting::get_parallel_solution_vector(Vec &vec){
	convection->compute_one_step();
};



void TransportOperatorSplitting::get_solution_vector(double * &x, unsigned int &a){
	convection->compute_one_step();
};



void TransportOperatorSplitting::set_cross_section_field(Field< 3, FieldValue<3>::Scalar >* cross_section)
{
    convection->set_cross_section_field(cross_section);

    /*if (Semchem_reactions != NULL) {
        Semchem_reactions->set_cross_section(cross_section);
        Semchem_reactions->set_sorption_fields(&convection->get_data()->porosity, &convection->get_data()->por_imm, &convection->get_data()->phi);
    }
  if (sorptions != NULL)
  {
	  sorptions->set_porosity(&(convection->get_data()->porosity),&(convection->get_data()->por_imm));
  }*/
}



void TransportOperatorSplitting::calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance)
{}

void TransportOperatorSplitting::calc_elem_sources(vector<vector<double> > &mass, vector<vector<double> > &src_balance)
{}





