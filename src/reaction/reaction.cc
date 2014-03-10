#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>

#include "reaction/reaction.hh"
#include "system/system.hh"
#include "mesh/mesh.h"
#include "input/accessors.hh"

using namespace Input::Type;
using namespace std;

AbstractRecord Reaction::input_type
    = AbstractRecord("Reactions", "Equation for reading information about simple chemical reactions.")
        .declare_key("species", Array(String()), Default::obligatory(),
                     "Names of the species that take part in the reaction model.");


Reaction::Reaction(Mesh &init_mesh, Input::Record in_rec, const  vector<string> &names) //(double timeStep, Mesh * mesh, int nrOfSpecies, bool dualPorosity) //(double timestep, int nrOfElements, double ***ConvectionMatrix)
    : EquationBase(init_mesh, in_rec),
      names_(names),
      n_all_substances_ (names.size())
{
  initialize_substance_ids(names, in_rec);
}

Reaction::~Reaction()
{
}

void Reaction::initialize_substance_ids(const vector< string >& names, Input::Record in_rec)
{
  Input::Array species_array = in_rec.val<Input::Array>("species");
  unsigned int k, idx, i_spec = 0;
  
  for(Input::Iterator<string> spec_iter = species_array.begin<string>(); spec_iter != species_array.end(); ++spec_iter, i_spec++)
  {
    //finding name in the global array of names
    for(k = 0; k < names.size(); k++)
    {
      if (*spec_iter == names[k]) 
      {
        idx = k;
        break;
      }
    }
    
    if ((idx < names.size()) && (idx >= 0)) 
    {
      substance_id[i_spec] = idx;       //mapping - if not found, it creates new map
    }
      else    xprintf(UsrErr,"Wrong name of %d-th adsorbing specie.\n", i_spec);
    }
    n_substances_ = substance_id.size();
}

double **Reaction::compute_reaction(double **concentrations, int loc_el) //multiplication of concentrations array by reaction matrix
{
    cout << "double **Reaction::compute_reaction(double **concentrations, int loc_el) needs to be re-implemented in ancestors." << endl;
	return concentrations;
}

void Reaction::set_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc_)
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
