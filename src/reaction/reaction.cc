#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include "reaction/reaction.hh"
#include "system/system.hh"
#include "materials.hh"
#include "transport/transport.h"
#include "system/par_distribution.hh"
#include "mesh/mesh.h"


using namespace std;

Reaction::Reaction(TimeMarks &marks, Mesh &init_mesh, MaterialDatabase &material_database)//(double timeStep, Mesh * mesh, int nrOfSpecies, bool dualPorosity) //(double timestep, int nrOfElements, double ***ConvectionMatrix)
	: EquationBase(marks, init_mesh, material_database), dual_porosity_on(false), prev_conc(NULL)
{
	nr_of_decays = OptGetInt("Reaction_module","Nr_of_decay_chains","0");
	nr_of_FoR = OptGetInt("Reaction_module","Nr_of_FoR","0");
	nr_of_species = OptGetInt("Transport", "N_substances", "0");
	nom_pol_deg = OptGetInt("Reaction_module","Nom_pol_deg","0");
	den_pol_deg = OptGetInt("Reaction_module","Den_pol_deg","0");
	//set_time_step(); //temporary solution, it reads Global Save_step
	set_dual_porosity();
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
	/*int cols, rows, both;

	if((nr_of_decays > 0) || (nr_of_FoR > 0)){
		for(cols = 0; cols < nr_of_species; cols++){
		prev_conc[cols] = concentrations[cols][loc_el];
		//xprintf(Msg,"\n%d. of %d substances concentration is %f\n", cols,nr_of_species, concentrations[cols][loc_el]); //prev_conc[cols]); //commented to speed the computation up
		concentrations[cols][loc_el] = 0.0;
		}
        for(rows = 0; rows <nr_of_species; rows++){
            for(cols = 0; cols <nr_of_species; cols++){
                concentrations[rows][loc_el] += prev_conc[cols] * reaction_matrix[cols][rows];
            }
            //xprintf(Msg,"\n%d. of %d substances concentration after reaction is %f\n", rows,nr_of_species, concentrations[rows][loc_el]); //commented to speed the computation up
        }
	}*/
	return concentrations;
}

void Reaction::compute_one_step(void)
{
    /*if (reaction_matrix == NULL)   return;

    START_TIMER("decay_step");
	 //for (int loc_el = 0; loc_el < distribution->lsize(distribution->myp()); loc_el++)
	for (int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
	 {
	 	this->compute_reaction(concentration_matrix[MOBILE], loc_el);
	    if (dual_porosity_on == true) {
	     this->compute_reaction(concentration_matrix[IMMOBILE], loc_el);
	    }

	 }
    END_TIMER("decay_step");*/
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
	/*if((nr_of_decays > 0) || (nr_of_FoR > 0)){
		release_reaction_matrix();
		allocate_reaction_matrix();
		if(matrix_exp_on == false)
		{
			modify_reaction_matrix_repeatedly();
		}
	}*/
	return;
}

void Reaction::set_time_step(void)
{
	time_step = OptGetDbl("Global","Save_step","1.0");
	return;
}

//void Reaction::set_mesh_(Mesh *mesh_in){mesh = mesh_in; return;}

void Reaction::set_dual_porosity()
{
	this->dual_porosity_on = OptGetBool("Transport", "Dual_porosity", "no");
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
