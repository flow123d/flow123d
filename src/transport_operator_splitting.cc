/*
 * transport_operator_splitting.cc
 *
 *  Created on: May 21, 2011
 *      Author: jiri
 */

#include "transport_operator_splitting.hh"
#include <petscmat.h>
#include "system/sys_vector.hh"
#include <time_governor.hh>
#include <materials.hh>
//#include "equation.hh"
#include "transport.h"
#include "mesh/mesh.h"


TransportOperatorSplitting::TransportOperatorSplitting(MaterialDatabase *material_database, Mesh *init_mesh)
{
	mat_base = material_database;
	mesh = init_mesh;
    ConvectionTransport *convection;


    double problem_save_step = OptGetDbl("Global", "Save_step", "1.0");
    double problem_stop_time = OptGetDbl("Global", "Stop_time", "1.0");


	convection = new ConvectionTransport(mat_base, mesh);

	/*
	if(OptGetBool("Semchem_module", "Compute_reactions", "no") == yes)
		chemistry = new Chemistry;
	*/



	time = new TimeGovernor(0.0,0.0,problem_stop_time,problem_stop_time);
	time->constrain_dt(convection->cfl_time_constrain());

}

void TransportOperatorSplitting::compute_one_step(){

	convection->compute_one_step();
	//chemistry->compute_one_step();
	time->next_time();
}

void TransportOperatorSplitting::compute_until_save_time(){

	//while(time->output_time())
		compute_one_step();

	//call output

}




void TransportOperatorSplitting::get_parallel_solution_vector(Vec &vec){
	convection->compute_one_step();
};

void TransportOperatorSplitting::get_solution_vector(double * &x, unsigned int &a){
	convection->compute_one_step();
};

