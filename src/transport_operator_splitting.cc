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

	convection = new ConvectionTransport(mat_base, mesh);

	/*
	if(OptGetBool("Semchem_module", "Compute_reactions", "no") == yes)
		chemistry = new Chemistry;
	*/
//	time = new TimeGovernor();
	//time->

}

void TransportOperatorSplitting::compute_one_step(){

	convection->compute_one_step();
	/*

	chemistry->compute_one_step();

	 */

	// TimeGovernor
}


void TransportOperatorSplitting::get_parallel_solution_vector(Vec &vec){
	convection->compute_one_step();
};

void TransportOperatorSplitting::get_solution_vector(double * &x, unsigned int &a){
	convection->compute_one_step();
};

