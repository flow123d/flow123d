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
#include "reaction/sorption.hh"
//#include "reaction/dual_por_exchange.hh"

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
	.declare_key("mass_balance", MassBalance::input_type, Default::optional(), "Settings for computing mass balance.")
	.declare_key("output", TransportBase::input_type_output_record, Default::obligatory(),
    		"Parameters of output stream.");


Record TransportBase::input_type_output_record
	= Record("TransportOutput", "Output setting for transport equations.")
//	.declare_key("output_stream", OutputTime::input_type, Default::obligatory(),
//			"Parameters of output stream.")
	.declare_key("save_step", Double(0.0), Default::obligatory(),
			"Interval between outputs.")
	.declare_key("output_times", Array(Double(0.0)),
			"Explicit array of output times (can be combined with 'save_step'.")
	.declare_key("conc_mobile_p0", String(),
			"Name of output stream for P0 approximation of the concentration in mobile phase.")
	.declare_key("conc_immobile_p0", String(),
			"Name of output stream for P0 approximation of the concentration in immobile phase.")
	.declare_key("conc_mobile_sorbed_p0", String(),
			"Name of output stream for P0 approximation of the surface concentration of sorbed mobile phase.")
	.declare_key("conc_immobile_sorbed_p0", String(),
			"Name of output stream for P0 approximation of the surface concentration of sorbed immobile phase.");


Record TransportOperatorSplitting::input_type
	= Record("TransportOperatorSplitting",
            "Explicit FVM transport (no diffusion)\n"
            "coupled with reaction and adsorption model (ODE per element)\n"
            " via operator splitting.")
    .derive_from(AdvectionProcessBase::input_type)
    .declare_key("substances", Array(String()), Default::obligatory(),
    		"Names of transported substances.")
    	    // input data
    .declare_key("sorption_enable", Bool(), Default("false"),
    		"Model of sorption.")
    .declare_key("dual_porosity", Bool(), Default("false"),
    		"Dual porosity model.")
	.declare_key("reactions", Reaction::input_type, Default::optional(),
                "Initialization of per element reactions.")
    .declare_key("adsorptions", Sorption::input_type, Default::optional(),
    			"Initialization of per element sorptions.")
    .declare_key("bc_data", Array(ConvectionTransport::EqData().boundary_input_type()), IT::Default::obligatory(), "")
    .declare_key("bulk_data", Array(ConvectionTransport::EqData().bulk_input_type()),
    		IT::Default::obligatory(), "");


TransportBase::TransportEqData::TransportEqData(const std::string& eq_name)
: EqDataBase(eq_name)
{

	ADD_FIELD(por_m, "Mobile porosity", Default("1"));

	ADD_FIELD(sources_density, "Density of concentration sources.", Default("0"));
	ADD_FIELD(sources_sigma, "Concentration flux.", Default("0"));
	ADD_FIELD(sources_conc, "Concentration sources threshold.", Default("0"));

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
  convection(NULL),
  decayRad(NULL),
  sorptions(NULL),
  sorptions_immob(NULL),
  Semchem_reactions(NULL)
{
	Distribution *el_distribution;
	int *el_4_loc;

    // double problem_save_step = OptGetDbl("Global", "Save_step", "1.0");

    in_rec.val<Input::Array>("substances").copy_to(subst_names_);
    n_subst_ = subst_names_.size();
	convection = new ConvectionTransport(*mesh_, in_rec);

	Input::Iterator<Input::AbstractRecord> reactions_it = in_rec.find<Input::AbstractRecord>("reactions");
	if ( reactions_it ) {
		if (reactions_it->type() == Linear_reaction::input_type ) {
	        decayRad =  new Linear_reaction(init_mesh, *reactions_it, subst_names_);
	        convection->get_par_info(el_4_loc, el_distribution);
	        decayRad->set_dual_porosity(convection->get_dual_porosity());
	        static_cast<Linear_reaction *> (decayRad) -> modify_reaction_matrix();
	        decayRad->set_concentration_matrix(convection->get_concentration_matrix(), el_distribution, el_4_loc);

	        //Supresses possibility to combine reactions
	        /*Semchem_reactions = NULL;
	        sorptions = NULL;*/
		} else
	    if (reactions_it->type() == Pade_approximant::input_type) {
            decayRad = new Pade_approximant(init_mesh, *reactions_it, subst_names_ );
	        convection->get_par_info(el_4_loc, el_distribution);
	        decayRad->set_dual_porosity(convection->get_dual_porosity());
	        static_cast<Pade_approximant *> (decayRad) -> modify_reaction_matrix();
	        decayRad->set_concentration_matrix(convection->get_concentration_matrix(), el_distribution, el_4_loc);

	        //Supresses possibility to combine reactions
	        /*Semchem_reactions = NULL;
	        sorptions = NULL;*/
	    } else
	    if (reactions_it->type() == Semchem_interface::input_type ) {
	        Semchem_reactions = new Semchem_interface(0.0, mesh_, n_subst_, convection->get_dual_porosity()); //(mesh->n_elements(),convection->get_concentration_matrix(), mesh);
	        Semchem_reactions->set_el_4_loc(el_4_loc);
	        Semchem_reactions->set_concentration_matrix(convection->get_concentration_matrix(), el_distribution, el_4_loc);

	        /*decayRad = NULL;
	        sorptions = NULL;*/
	    } else {
	        xprintf(UsrErr, "Wrong reaction type.\n");
	    }
	} else {
	    decayRad = NULL;
	    Semchem_reactions = NULL;
	}

	Input::Iterator<Input::Record> sorptions_it = in_rec.find<Input::Record>("adsorptions");
	if (sorptions_it){
        // Part for mobile zone description follows.
	    sorptions = new Sorption(init_mesh, *sorptions_it, subst_names_);
        convection->get_par_info(el_4_loc, el_distribution);
	    sorptions->set_dual_porosity(convection->get_dual_porosity());
	    //xprintf(Msg,"sorption->set_dual_porosity() finished successfuly.\n");
	    sorptions->set_porosity(&(convection->get_data()->por_m), &(convection->get_data()->por_imm)); //, &(convection->get_data()->por_imm));
	    sorptions->set_phi(&(convection->get_data()->phi));
	    //xprintf(Msg,"sorption->set_phi() finished successfuly.\n");
	    sorptions->prepare_inputs(*sorptions_it, MOBILE);
	    //xprintf(Msg,"sorption->prepare_inputs() finished successfuly.\n");
	    double ***conc_matrix = convection->get_concentration_matrix();
	    sorptions->set_concentration_matrix(conc_matrix[MOBILE], el_distribution, el_4_loc);
	    sorptions->set_sorb_conc_array(el_distribution->lsize());

	    if(convection->get_dual_porosity()){
	    	sorptions_immob = new Sorption(init_mesh, *sorptions_it, subst_names_);
	    	//dual_por_exchange = new Dual_por_exchange(init_mesh, *sorptions_it, subst_names_);
		    sorptions_immob->set_dual_porosity(convection->get_dual_porosity());
	    	sorptions_immob->set_porosity(&(convection->get_data()->por_m), &(convection->get_data()->por_imm));
	    	sorptions_immob->set_phi(&(convection->get_data()->phi));
		    sorptions_immob->prepare_inputs(*sorptions_it, IMMOBILE);
		    sorptions_immob->set_concentration_matrix(conc_matrix[MOBILE], el_distribution, el_4_loc);
		    sorptions_immob->set_immob_concentration_matrix(conc_matrix[IMMOBILE], el_distribution, el_4_loc);
		    sorptions_immob->set_sorb_conc_array(el_distribution->lsize());
	    }else{
		    sorptions_immob = NULL;
	    }
	  } else{
	    sorptions = NULL;
	    sorptions_immob = NULL;
	}
	
	output_mark_type = convection->mark_type() | TimeGovernor::marks().type_fixed_time() | TimeGovernor::marks().type_output();
    time_ = new TimeGovernor(in_rec.val<Input::Record>("time"), output_mark_type );

}

TransportOperatorSplitting::~TransportOperatorSplitting()
{
    //delete field_output;
    delete convection;
    if (decayRad) delete decayRad;
    if (sorptions) delete sorptions;
    if (sorptions_immob) delete sorptions_immob;
    if (Semchem_reactions) delete Semchem_reactions;
    delete time_;
}




void TransportOperatorSplitting::output_data(){

    if (time_->is_current(output_mark_type)) {
        
        START_TIMER("TOS-output data");
        DBGMSG("\nTOS: output time: %f\n", time_->t());

        convection->output_vector_gather();
        
        convection->output_data();
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
    time_->view("TOS");    //show time governor
    
    convection->set_target_time(time_->t());
	if (decayRad) decayRad->set_time_step(convection->time().estimate_dt());
	if (sorptions) sorptions->set_time_step(convection->time().estimate_dt());
	if (sorptions_immob) sorptions_immob->set_time_step(convection->time().estimate_dt());
	//if (dual_por_exchange) dual_por_exchange->set_time_step(convection->time().estimate_dt());
	// TODO: update Semchem time step here!!
	if (Semchem_reactions) Semchem_reactions->set_timestep(convection->time().estimate_dt());

        
    xprintf( Msg, "TOS: time: %f        CONVECTION: time: %f      dt_estimate: %f\n", 
             time_->t(), convection->time().t(), convection->time().estimate_dt() );
    
    START_TIMER("TOS-one step");
    int steps=0;
    while ( convection->time().lt(time_->t()) )
    {
        steps++;
	    // one internal step
	    convection->compute_one_step();
		//Just temporarly commented.
	    //if (dual_por_exchange) dual_por_exchange->compute_one_step();
	    if(decayRad) decayRad->compute_one_step();
	    if(Semchem_reactions) Semchem_reactions->compute_one_step();
	    if(sorptions) sorptions->compute_one_step();//equilibrial sorption at the end of simulated time-step
	    if(sorptions_immob) sorptions_immob->compute_one_step();
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
        Semchem_reactions->set_sorption_fields(&convection->get_data()->por_m, &convection->get_data()->por_imm, &convection->get_data()->phi);
    }
  if (sorptions != NULL)
  {
	  sorptions->set_porosity(&(convection->get_data()->por_m),&(convection->get_data()->por_imm));
  }*/
}



void TransportOperatorSplitting::calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance)
{}

void TransportOperatorSplitting::calc_elem_sources(vector<vector<double> > &mass, vector<vector<double> > &src_balance)
{}





