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
	// TOdO: this has to be set after construction of transport matrix !!
	//time_->set_constrain(convection->get_cfl_time_constrain());

	solved = true;
}

void TransportOperatorSplitting::update_solution() {
    convection->convection();
	//convection->compute_one_step();
    // Calling linear reactions and Semchem
	decayRad->compute_one_step();
	Semchem_reactions->compute_one_step();

	solved=true;
}

void TransportOperatorSplitting::compute_until_save_time(){

	//while(time->output_time())
		compute_one_step();

	//call output

}

void TransportOperatorSplitting::set_velocity_field(Vec &vec)
{
	//convection->read_flow_field_vector(&vec);
};


void TransportOperatorSplitting::get_parallel_solution_vector(Vec &vec){
	convection->compute_one_step();
};

void TransportOperatorSplitting::get_solution_vector(double * &x, unsigned int &a){
	convection->compute_one_step();
};

