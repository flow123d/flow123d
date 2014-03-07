#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>

#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "semchem/semchem_interface.hh"

#include "system/system.hh"
//#include "system/par_distribution.hh"
#include "mesh/mesh.h"

#include "input/accessors.hh"


using namespace Input::Type;

AbstractRecord Reaction::input_type
	= AbstractRecord("Reactions", "Equation for reading information about simple chemical reactions.");
//		rec.declare_key("substances", Array(String()), Default::obligatory(),
//								"Names of transported chemical species.");


using namespace std;

Reaction::Reaction(Mesh &init_mesh, Input::Record in_rec, const  vector<string> &names) //(double timeStep, Mesh * mesh, int nrOfSpecies, bool dualPorosity) //(double timestep, int nrOfElements, double ***ConvectionMatrix)
    : EquationBase(init_mesh, in_rec),
//      dual_porosity_on(false), 
//      prev_conc(NULL),
      names_(names)
{
//	prev_conc = new double[ n_substances() ];
	//if(timeStep < 1e-12) this->set_time_step(timeStep); else this->set_time_step(0.5); // temporary solution
}

Reaction::~Reaction()
{
/*
	if(prev_conc != NULL){
		delete[](prev_conc);
		prev_conc = NULL;
	}
*/
}


double **Reaction::compute_reaction(double **concentrations, int loc_el) //multiplication of concentrations array by reaction matrix
{
    cout << "double **Reaction::compute_reaction(double **concentrations, int loc_el) needs to be re-implemented in ancestors." << endl;
	return concentrations;
}

void Reaction::set_concentration_matrix(double ***ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc_)
{
	concentration_matrix = ConcentrationMatrix;
	distribution = conc_distr;
	el_4_loc = el_4_loc_;
	return;
}


void Reaction::get_parallel_solution_vector(Vec &vec){
	cout << "Reaction.get_parallel_solution_vector(Vec &vec) is not implemented." << endl; //convection->compute_one_step();
}

void Reaction::get_solution_vector(double * &x, unsigned int &a){
	cout << "Reaction.get_solution_vector(double * &x, unsigned int &a) is not implemented." << endl; //convection->compute_one_step();
}

void Reaction::update_solution(void)
{
	cout << "Reaction::update_solution() is not implemented." << endl;
}

void Reaction::choose_next_time(void)
{
	cout << "Reaction::choose_next_time() is not implemented." << endl;
}

void Reaction::set_time_step_constrain(double dt)
{
	cout << "Reaction::choose_time_step_constrain(double dt) is not implemented." << endl;
}


unsigned int Reaction::find_subst_name(const string &name)
{

    unsigned int k=0;
	for(; k < names_.size(); k++)
		if (name == names_[k]) return k;

	return k;
}

void Reaction::set_mesh(Mesh &mesh)
{
	mesh_ = &mesh;
}

void Reaction::set_names(const std::vector<string> &names)
{
	names_ = names;
}

Element * Reaction::get_element_for_dof_index(unsigned int idx)
{
	return &( mesh_->element[idx] );
}
