/*
 * transport_operator_splitting.cc
 *
 *  Created on: May 21, 2011
 *      Author: jiri
 */

#include "transport/transport_operator_splitting.hh"
#include <petscmat.h>
#include "system/sys_vector.hh"
#include <time_governor.hh>
#include <materials.hh>
#include "coupling/equation.hh"
#include "transport/transport.h"
#include "transport/transport_dg.hh"
#include "mesh/mesh.h"

#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"

#include "semchem/semchem_interface.hh"

#include "la/distribution.hh"
#include "io/output.h"

#include "input/input_type.hh"
#include "input/accessors.hh"


Input::Type::AbstractRecord & TransportBase::get_input_type()
{
	using namespace Input::Type;
	static AbstractRecord rec("Transport", "Secondary equation for transport of substances.");

	if (!rec.is_finished()) {
	    rec.declare_key("time", TimeGovernor::get_input_type(), Default::obligatory(),
	                    "Time governor setting for the transport model.");
        rec.declare_key("substances", Array(String()), Default::obligatory(),
                        "Names of transported substances.");

	    // input data
	    rec.declare_key("sorption_enable", Bool(), Default("false"),
						                "Model of sorption.");
		rec.declare_key("dual_porosity", Bool(), Default("false"),
						                "Dual porosity model.");

		rec.declare_key("initial_file", FileName::input(), Default::obligatory(),
						                "Input file with initial concentrations.");
		rec.declare_key("boundary_file", FileName::input(), Default::obligatory(),
						                "Input file with boundary conditions.");
		rec.declare_key("bc_times", Array(Double()), Default::optional(),
				 	 	                "Times for changing the boundary conditions.");
		rec.declare_key("sources_file", FileName::input(), Default::optional(),
		                                "File with data for the source term in the transport equation.");
        rec.declare_key("output", TransportBase::get_input_type_output_record(), Default::obligatory(),
                                        "Parameters of output stream.");

		rec.finish();

		TransportOperatorSplitting::get_input_type();
		TransportDG::get_input_type();

		rec.no_more_descendants();
	}
	return rec;
}



Input::Type::Record & TransportBase::get_input_type_output_record()
{
    using namespace Input::Type;
    static Record out_rec("TransportOutput", "Output setting for transport equations.");

    if (!out_rec.is_finished()) {
        out_rec.declare_key("output_stream", OutputTime::get_input_type(), Default::obligatory(),
                "Parameters of output stream.");

        out_rec.declare_key("save_step", Double(0.0), Default::obligatory(),
                "Interval between outputs.");

        out_rec.declare_key("output_times", Array(Double(0.0)),
                        "Explicit array of output times (can be combined with 'save_step'.");

        out_rec.declare_key("conc_mobile_p0", String(),
                        "Name of output stream for P0 approximation of the concentration in mobile phase.");
        out_rec.declare_key("conc_immobile_p0", String(),
                                "Name of output stream for P0 approximation of the concentration in immobile phase.");
        out_rec.declare_key("conc_mobile_sorbed_p0", String(),
                                "Name of output stream for P0 approximation of the surface concentration of sorbed mobile phase.");
        out_rec.declare_key("conc_immobile_sorbed_p0", String(),
                                "Name of output stream for P0 approximation of the surface concentration of sorbed immobile phase.");

        out_rec.finish();

    }
    return out_rec;
}


TransportOperatorSplitting::TransportOperatorSplitting(TimeMarks &marks, Mesh &init_mesh, MaterialDatabase &material_database, const Input::Record &in_rec)
: TransportBase(marks, init_mesh, material_database, in_rec)
{
	Distribution *el_distribution;
	int *el_4_loc;

    // double problem_save_step = OptGetDbl("Global", "Save_step", "1.0");

	convection = new ConvectionTransport(marks, *mesh_, *mat_base, in_rec);

	Input::Iterator<Input::AbstractRecord> reactions_it = in_rec.find<Input::AbstractRecord>("reactions");
	if ( reactions_it ) {
		if (reactions_it->type() == Linear_reaction::get_input_type() ) {
	        decayRad =  new Linear_reaction(marks, init_mesh, material_database, *reactions_it,
	                                        convection->get_substance_names());
	        convection->get_par_info(el_4_loc, el_distribution);
	        decayRad->set_dual_porosity(convection->get_dual_porosity());
	        static_cast<Linear_reaction *> (decayRad) -> modify_reaction_matrix();
	        decayRad->set_concentration_matrix(convection->get_prev_concentration_matrix(), el_distribution, el_4_loc);

	        Semchem_reactions = NULL;
		} else
	    if (reactions_it->type() == Pade_approximant::get_input_type() ) {
	        decayRad = new Pade_approximant(marks, init_mesh, material_database, in_rec,
	                                        convection->get_substance_names());
	        Semchem_reactions = NULL;
	    } else
	    if (reactions_it->type() == Semchem_interface::get_input_type() ) {
	        Semchem_reactions = new Semchem_interface(0.0, mesh_, convection->get_n_substances(), convection->get_dual_porosity()); //(mesh->n_elements(),convection->get_concentration_matrix(), mesh);
	        Semchem_reactions->set_el_4_loc(el_4_loc);
	        Semchem_reactions->set_concentration_matrix(convection->get_prev_concentration_matrix(), el_distribution, el_4_loc);
	    } else {
	        xprintf(UsrErr, "Wrong reaction type.\n");
	    }
	} else {
	    decayRad = NULL;
	    Semchem_reactions = NULL;
	}

	time_ = new TimeGovernor(in_rec.val<Input::Record>("time"), *time_marks, this->mark_type());
	output_mark_type = this->mark_type() | time_marks->type_fixed_time() | time_marks->type_output();

    time_marks->add_time_marks(0.0,
            in_rec.val<Input::Record>("output").val<double>("save_step"),
            time_->end_time(), output_mark_type );
	// TODO: this has to be set after construction of transport matrix !!


	// register output vectors from convection
	double ***out_conc = convection->get_out_conc();
	vector<string> substance_name = convection->get_substance_names();

	// TODO: Add corresponding record to the in_rec
	Input::Record output_rec = in_rec.val<Input::Record>("output");

	field_output = new OutputTime(mesh_, output_rec.val<Input::Record>("output_stream"));


    for(int subst_id=0; subst_id < convection->get_n_substances(); subst_id++) {
         // TODO: What about output also other "phases", IMMOBILE and so on.
         std::string subst_name = substance_name[subst_id] + "_mobile";
         double *data = out_conc[MOBILE][subst_id];
         field_output->register_elem_data<double>(subst_name, "M/L^3", data , mesh_->n_elements());
    }
    // write initial condition
    convection->output_vector_gather();
    field_output->write_data(time_->t());

}

TransportOperatorSplitting::~TransportOperatorSplitting()
{
    delete field_output;
    delete convection;
    if (decayRad) delete decayRad;
    if (Semchem_reactions) delete Semchem_reactions;
    delete time_;
}


Input::Type::Record &TransportOperatorSplitting::get_input_type()
{
    using namespace Input::Type;
    static Record rec("TransportOperatorSplitting",
            "Explicit FVM transport (no diffusion)\n"
            "coupled with reaction and sorption model (ODE per element)\n"
            " via. operator splitting.");

    if (!rec.is_finished()) {
        rec.derive_from(TransportBase::get_input_type());
        rec.declare_key("reactions", Reaction::get_input_type(), Default::optional(),
                "Initialization of per element reactions.");

        rec.finish();
    }
    return rec;
}



void TransportOperatorSplitting::output_data(){

    if (time_->is_current(output_mark_type)) {
        convection->output_vector_gather();
        field_output->write_data(time_->t());
    }
}

void TransportOperatorSplitting::read_simulation_step(double sim_step) {
	time_->set_constrain(sim_step);
}


void TransportOperatorSplitting::update_solution() {


    time_->next_time();
	convection->set_target_time(time_->t());

	if (decayRad) decayRad->set_time_step(convection->time().estimate_dt());
	//cout << "recent time step value is " << decayRad->get_time_step() << endl;
	// TODO: update Semchem time step here!!
	if (Semchem_reactions) Semchem_reactions->set_timestep(convection->time().estimate_dt());

    xprintf( Msg, "t: %f (TOS)                  cfl_dt: %f ", convection->time().t(), convection->time().estimate_dt() );
    START_TIMER("transport_steps");
    int steps=0;
    while ( convection->time().lt(time_->t()) )
    {
        steps++;
	    // one internal step
	    //xprintf( Msg, "Time : %f\n", convection->time().t() );
	    convection->compute_one_step();
	    // Calling linear reactions and Semchem, temporarily commented
	    if(decayRad) decayRad->compute_one_step();
	    if (Semchem_reactions) Semchem_reactions->compute_one_step();
	}
    END_TIMER("transport_steps");
    //DBGMSG("conv time: %f TOS time: %f\n", convection->time().t(), time_->t());
    xprintf( Msg, " steps: %d\n",steps);
}


void TransportOperatorSplitting::set_velocity_field(const MH_DofHandler &dh)
{
    mh_dh = &dh;
	convection->set_flow_field_vector( dh );
};


void TransportOperatorSplitting::get_parallel_solution_vector(Vec &vec){
	convection->compute_one_step();
};

void TransportOperatorSplitting::get_solution_vector(double * &x, unsigned int &a){
	convection->compute_one_step();
};

