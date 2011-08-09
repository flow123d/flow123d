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
#include "equation.hh"
#include "transport/transport.h"
#include "mesh/mesh.h"
#include "reaction/linear_reaction.hh"
#include "semchem/semchem_interface.hh"
#include "system/par_distribution.hh"


TransportOperatorSplitting::TransportOperatorSplitting(TimeMarks &marks, Mesh &init_mesh, MaterialDatabase &material_database )
: TransportBase(marks, init_mesh, material_database)
{
	Distribution *el_distribution;
	int *el_4_loc;

    double problem_save_step = OptGetDbl("Global", "Save_step", "1.0");
    double problem_stop_time = OptGetDbl("Global", "Stop_time", "1.0");

	convection = new ConvectionTransport(marks, *mesh_, *mat_base);

	// Chemistry initialization
	decayRad = new Linear_reaction(0.0, mesh_, convection->get_n_substances(), convection->get_dual_porosity());
	convection->get_par_info(el_4_loc, el_distribution);
	//decayRad->release_reaction_matrix();
	decayRad->set_concentration_matrix(convection->get_concentration_matrix(), el_distribution, el_4_loc);
	Semchem_reactions = new Semchem_interface(0.0, mesh_, convection->get_n_substances(), convection->get_dual_porosity()); //(mesh->n_elements(),convection->get_concentration_matrix(), mesh);
	Semchem_reactions->set_el_4_loc(el_4_loc);
	Semchem_reactions->set_concentration_matrix(convection->get_concentration_matrix(), el_distribution, el_4_loc);


	time_ = new TimeGovernor(0.0, problem_stop_time, *time_marks);
    output_mark_type = time_marks->new_strict_mark_type();
    time_marks->add_time_marks(0.0, OptGetDbl("Global", "Save_step", "1.0"), time_->end_time(), output_mark_type );
	// TOdO: this has to be set after construction of transport matrix !!

	solved = true;

	// register output vectors from convection
	double ***out_conc = convection->get_out_conc();
	char    **substance_name = convection->get_substance_names();

	string output_file = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Transport", "Transport_out", "\\"));
	DBGMSG("create output\n");
	field_output = new OutputTime(mesh_, output_file);

	/*
    transport_out_fname = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Transport", "Transport_out", "\\"));
    transport_out_im_fname = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Transport", "Transport_out_im", "\\"));
    transport_out_sorp_fname = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Transport", "Transport_out_sorp", "\\"));
    transport_out_im_sorp_fname = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Transport", "Transport_out_im_sorp", "\\"));
    */

    for(int subst_id=0; subst_id < convection->get_n_substances(); subst_id++) {
         // TODO: What about output also other "phases", IMMOBILE and so on.
         std::string subst_name = std::string(substance_name[subst_id]) + "_mobile";
         double *data = out_conc[MOBILE][subst_id];
         field_output->register_elem_data<double>(subst_name, "M/L^3", data , mesh_->n_elements());
    }
    // write initial condition
    field_output->write_data(time_->t());

}

TransportOperatorSplitting::~TransportOperatorSplitting()
{
    delete field_output;
    delete convection;
    delete decayRad;
    delete Semchem_reactions;
    delete time_;
}

void TransportOperatorSplitting::output_data(){

    if (time_->is_current(output_mark_type)) {
        field_output->write_data(time_->t());
    }
}

void TransportOperatorSplitting::read_simulation_step(double sim_step) {
	time_->set_constrain(sim_step);
}


void TransportOperatorSplitting::update_solution() {

    // setup convection matrix for actual CFL time step
	//double cfl_dt =  convection->get_cfl_time_constrain();
	//DBGMSG("cfl: %f dt: %f\n",cfl_dt, );
	//int steps = (int) ceil(time_->dt() / cfl_dt);
	//cfl_dt = time_->dt() / steps;

	convection->set_target_time(time_->t());
	// TODO: update linear reaciton marix here !!
	decayRad->set_time_step(convection->time().dt());

    START_TIMER("transport_steps");
	while ( convection->time().lt(time_->t()) ) {
	    // one internal step
	    xprintf( Msg, "Time : %f\n", convection->time().t() );
	    convection->compute_one_step();
	    // Calling linear reactions and Semchem
	    decayRad->compute_one_step();
	    Semchem_reactions->compute_one_step();
	}
    END_TIMER("transport_steps");
    xprintf( Msg, "O.K.\n");
	solved=true;
}


void TransportOperatorSplitting::set_velocity_field(Vec &vec)
{
	convection->set_flow_field_vector(&vec);
};


void TransportOperatorSplitting::get_parallel_solution_vector(Vec &vec){
	convection->compute_one_step();
};

void TransportOperatorSplitting::get_solution_vector(double * &x, unsigned int &a){
	convection->compute_one_step();
};

