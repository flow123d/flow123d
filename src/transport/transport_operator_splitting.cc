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
	Distribution *distribution;
	int *el_4_loc;

    double problem_save_step = OptGetDbl("Global", "Save_step", "1.0");
    double problem_stop_time = OptGetDbl("Global", "Stop_time", "1.0");

	convection = new ConvectionTransport(mat_base, mesh_);

	// Chemistry initialization
	decayRad = new Linear_reaction(convection->get_cfl_time_constrain(), mesh_, convection->get_n_substances(), convection->get_dual_porosity());
	convection->get_par_info(el_4_loc, distribution);
	decayRad->set_concentration_matrix(convection->get_concentration_matrix(), distribution, el_4_loc);
	Semchem_reactions = new Semchem_interface(convection->get_cfl_time_constrain(), mesh_, convection->get_n_substances(), convection->get_dual_porosity()); //(mesh->n_elements(),convection->get_concentration_matrix(), mesh);
	Semchem_reactions->set_el_4_loc(el_4_loc);
	Semchem_reactions->set_concentration_matrix(convection->get_concentration_matrix(), distribution, el_4_loc);


	time_ = new TimeGovernor(0.0, problem_stop_time, *time_marks);
    TimeMark::Type output_mark_type = time_marks->new_strict_mark_type();
    time_marks->add_time_marks(0.0, OptGetDbl("Global", "Save_step", "1.0"), time_->end_time(), output_mark_type );
	// TOdO: this has to be set after construction of transport matrix !!

	solved = true;

	// register output vectors from convection
	double ***out_conc = convection->get_out_conc();
	char    **substance_name = convection->get_substance_names();

	string output_file = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Output", "Output_file", "\\"));
	DBGMSG("create output\n");
	output_time = new OutputTime(mesh_, output_file);

    for(int subst_id=0; subst_id < convection->get_n_substances(); subst_id++) {
         // TODO: What about output also other "phases", IMMOBILE and so on.
         std::string subst_name = std::string(substance_name[subst_id]);
         double *data = out_conc[MOBILE][subst_id];
         output_time->register_elem_data<double>(subst_name, "", data , mesh_->n_elements());
    }
    // write initial condition
    output_time->write_data(time_->t());

}
void TransportOperatorSplitting::output_data(){


	output_time->write_data(time_->t());
}

void TransportOperatorSplitting::read_simulation_step(double sim_step) {
	time_->set_constrain(sim_step);
}


void TransportOperatorSplitting::update_solution() {


	convection->convection();
	steps = (int) ceil(time_->dt() / convection->get_cfl_time_constrain());

    START_TIMER("transport_steps");
	for(int i=0;i < steps;i++)
		compute_internal_step();
    END_TIMER("transport_steps");
    xprintf( Msg, "O.K.\n");
	solved=true;
}

void TransportOperatorSplitting::compute_internal_step(){

	convection->compute_one_step();
    // Calling linear reactions and Semchem
	decayRad->compute_one_step();
	Semchem_reactions->compute_one_step();
}

void TransportOperatorSplitting::set_velocity_field(Vec &vec)
{
	convection->read_flow_field_vector(&vec);
};


void TransportOperatorSplitting::get_parallel_solution_vector(Vec &vec){
	convection->compute_one_step();
};

void TransportOperatorSplitting::get_solution_vector(double * &x, unsigned int &a){
	convection->compute_one_step();
};

