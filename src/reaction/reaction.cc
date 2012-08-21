#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>

#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"

#include "system/system.hh"
#include "materials.hh"
#include "transport/transport.h"
//#include "system/par_distribution.hh"
#include "mesh/mesh.h"

#include "input/accessors.hh"

Input::Type::AbstractRecord & Reaction::get_input_type()
{
	using namespace Input::Type;
	static AbstractRecord rec("Reactions", "Equation for reading information about simple chemical reactions.");

	if (!rec.is_finished()) {
		rec.declare_key("substances", Array(String()), Default::obligatory(),
								"Names of transported chemical species.");
		rec.declare_key("decays", Array( Linear_reaction::get_one_decay_substep() ), Default::optional(),
				"Description of particular decay chain substeps.");
		rec.declare_key("general_reaction", Array( Linear_reaction::get_one_decay_substep() ), Default::optional(),
				"Description of general chemical reactions.");
		rec.finish();

		Linear_reaction::get_input_type();

		rec.no_more_descendants();
	}
	return rec;
}

using namespace std;

Reaction::Reaction(TimeMarks &marks, Mesh &init_mesh, MaterialDatabase &material_database, Input::Record in_rec)//(double timeStep, Mesh * mesh, int nrOfSpecies, bool dualPorosity) //(double timestep, int nrOfElements, double ***ConvectionMatrix)
	: EquationBase(marks, init_mesh, material_database, Input::Record() ), dual_porosity_on(false), prev_conc(NULL)
{
	Input::Array decay_array = in_rec.val<Input::Array>("decays"); //Input::Array decay_array = in_rec.find<Input::Array>("decays");
	nr_of_decays = decay_array.size();
	//Input::Iterator<Input::Array> dec_it = in_rec.find<Input::Array>("decays");
	//Asi nebude treba rozlisovat kineticke reakce a rozpady
	/*nr_of_decays = 0;
	nr_of_FoR = 0;
	int i = 0;
	for(Input::Iterator<Input::Array> dec_it = decay_array.find<Input::Array>("decays"); dec_it != decay_array.end(); ++dec_it)
	{
		i++;
		Input::Iterator<double> it_hl = it->find<double>("half_life");
		if (it_hl) {
		    nr_of_decays++;
		} else {
		    it_hl = it->find<double>("kinetic");
		    if (it_hl) {
		   	   nr_of_FoR++;
		    } else {
		       xprintf(Msg, "You did not specify either the half life nor kinetic konstant for %d substep of decay chain.\n", i);
		       exit(1);
		    }
		}
	}*/

	Input::Array subst_array = in_rec.val<Input::Array>("substance_names"); //Number of all the substances contained in groundwater.
	nr_of_species = subst_array.size();

	//proverit nize uvedene predani dat
	this->dual_porosity_on = in_rec.val<bool>("dual_porosity");
	this->time_step = in_rec.val<double>("time_step");
	set_nr_of_elements(mesh_->n_elements());
	if(prev_conc != NULL){
		free(prev_conc);
		prev_conc = NULL;
	}
	prev_conc = (double *)xmalloc(nr_of_species * sizeof(double));
	//if(timeStep < 1e-12) this->set_time_step(timeStep); else this->set_time_step(0.5); // temporary solution
}

Reaction::~Reaction()
{

	if(prev_conc != NULL){
		free(prev_conc);
		prev_conc = NULL;
	}

	//release_reaction_matrix();
}

double **Reaction::compute_reaction(double **concentrations, int loc_el) //multiplication of concentrations array by reaction matrix
{
    cout << "double **Reaction::compute_reaction(double **concentrations, int loc_el) needs to be re-implemented in ancestors." << endl;
	return concentrations;
}

void Reaction::compute_one_step(void)
{
	cout << "Reaction::compute_one_step() needs to be re-implemented in ancestors." << endl;
	 return;
}

void Reaction::set_nr_of_species(int n_substances)
{
	this->nr_of_species = n_substances;
	return;
}

void Reaction::set_nr_of_elements(int nrOfElements)
{
	this->nr_of_elements = nrOfElements;
	return;
}

void Reaction::set_concentration_matrix(double ***ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc)
{
	concentration_matrix = ConcentrationMatrix;
	distribution = conc_distr;
	return;
}

void Reaction::set_time_step(double new_timestep){
	time_step = new_timestep;
	return;
}

void Reaction::set_time_step(Input::Record in_rec)
{
	time_step = in_rec.val<double>("time_step");
	return;
}

//void Reaction::set_mesh_(Mesh *mesh_in){mesh = mesh_in; return;}

void Reaction::set_dual_porosity(Input::Record in_rec)//obsolete function
{
	this->dual_porosity_on = in_rec.val<bool>("dual_porosity"); //OptGetBool("Transport", "Dual_porosity", "no");
	return;
}

double Reaction::get_time_step(void)
{
	return time_step;
}

int Reaction::faktorial(int k)
{
	int faktor = 1;

	if(k < 0)
	{
		//an error message should be placed here
		return 0;
	}

	while(k > 1)
	{
		faktor *= k;
		k--;
	}
	//xprintf(Msg,"\n Koeficient has a value %d.\n",faktor);
	return faktor;
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

/*double **Linear_reaction::modify_reaction_matrix(Input::Record in_rec) //prepare the matrix, which describes reactions
{
	return 0;
}*/
