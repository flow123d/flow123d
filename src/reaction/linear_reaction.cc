#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>

#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "system/system.hh"
#include "materials.hh"
#include "transport/transport.h"
#include "la/distribution.hh"
#include "mesh/mesh.h"

Input::Type::Record & Linear_reaction::get_one_decay_substep()
{
	using namespace Input::Type;
	static Record rec("Substep", "Equation for reading information about radioactive decays.");

	if(!rec.is_finished()){
		rec.declare_key("parent", String(), Default::obligatory(),
				"Identifier of an isotope.");
        rec.declare_key("half_life", Double(), Default::optional(),
                "Half life of the parent substance.");
        rec.declare_key("kinetic", Double(), Default::optional(),
                "Kinetic constants describing first order reactions.");
		rec.declare_key("products", Array(String()), Default::obligatory(),
				"Identifies isotopes which decays parental atom to.");
		rec.declare_key("branch_ratios", Array(Double()), Default::optional(),
				"Decay chain branching percentage.");
		rec.finish();
	}
	return rec;
}

Input::Type::Record & Linear_reaction::get_input_type()
{
	using namespace Input::Type;
	static Record rec("LinearReactions", "Information for a decision about the way to simulate radioactive decay.");

	if (!rec.is_finished()) {
	    rec.derive_from( Reaction::get_input_type() );
        rec.declare_key("decays", Array( Linear_reaction::get_one_decay_substep() ), Default::obligatory(),
                "Description of particular decay chain substeps.");


		/*rec.declare_key("substances", Array(String()), Default::obligatory(),
								"Names of transported isotopes.");
		rec.declare_key("half_lives", Array(Double()), Default::optional(),
				"Half lives of transported isotopes.");
		rec.declare_key("kinetic_constants", Array(Double()), Default::optional(),
				"Kinetic constants describing first order reactions.");
		rec.declare_key("kinetics", Array(Kinetics()), Default::optional(),
				"Description of particular first order reactions.");*/
		rec.declare_key("matrix_exp_on", Bool(), Default("false"),
				"Enables to use Pade approximant of matrix exponential.");
		/*rec.declare_key("nom_pol_deg", Integer(), Default("2"),
				"Polynomial degree of the nominator of Pade approximant.");
		rec.declare_key("den_pol_deg", Integer(), Default("2"),
				"Polynomial degree of the nominator of Pade approximant");*/
		rec.finish();
	}
	return rec;
}

using namespace std;

Linear_reaction::Linear_reaction(TimeMarks &marks, Mesh &init_mesh, MaterialDatabase &material_database, Input::Record in_rec)//(double timeStep, Mesh * mesh, int nrOfSpecies, bool dualPorosity, Input::Record in_rec) //(double timestep, int nrOfElements, double ***ConvectionMatrix)
	: Reaction(marks, init_mesh, material_database, in_rec), half_lives(NULL), substance_ids(NULL), reaction_matrix(NULL), bifurcation_on(false), prev_conc(NULL), matrix_exp_on(false)
{
	set_indices(in_rec);	//It needs to be called separetelly, earlier.
	set_half_lives(in_rec);
	set_bifurcation(in_rec);
	Input::Array names_array = in_rec.val<Input::Array>("substances");
	nr_of_species = names_array.size();
	Input::Array dec_array = in_rec.val<Input::Array>("decays");
	nr_of_decays = dec_array.size();
	//nr_of_isotopes = OptGetInt("Reaction_module","Nr_of_isotopes","0");
	allocate_reaction_matrix();
	set_time_step(0.5);
}

Linear_reaction::~Linear_reaction()
{
	int i, rows, n_subst;

	//n_subst = sizeof(*reaction_matrix)/sizeof(double *);
	if(half_lives != NULL){
		free(half_lives);
		half_lives = NULL;
	}

	if(substance_ids != NULL){
		free(substance_ids);
		substance_ids = NULL;
	}

	if(prev_conc != NULL){
		free(prev_conc);
		prev_conc = NULL;
	}

	release_reaction_matrix();
}

double **Linear_reaction::allocate_reaction_matrix(void) //reaction matrix initialization
{
	int index, rows, cols, dec_nr, prev_index;
	char dec_name[30];

	cout << "We are going to allocate reaction matrix" << endl;
	if(reaction_matrix == NULL)reaction_matrix = (double **)xmalloc(nr_of_species * sizeof(double*));//allocation section
	for(rows = 0; rows < nr_of_species; rows++){
		reaction_matrix[rows] = (double *)xmalloc(nr_of_species * sizeof(double));
	}
	for(rows = 0; rows < nr_of_species;rows++){
	 for(cols = 0; cols < nr_of_species; cols++){
		 if(rows == cols){
			 reaction_matrix[rows][cols] = 1.0;
		 }
	 }
	}
	//print_reaction_matrix();
	return reaction_matrix;
}

double **Linear_reaction::modify_reaction_matrix(Input::Record in_rec) //prepare the matrix, which describes reactions
{
	int rows,cols, index_par, bif_id;
	int *index_child;
	double rel_step, prev_rel_step;

	half_lives = (double *)xmalloc(nr_of_decays * sizeof(double));
	bifurcation.resize(nr_of_decays);

	if(reaction_matrix == NULL){
		xprintf(Msg,"\nReaction matrix pointer is NULL.\n");
		return NULL;
	}

	set_indices(in_rec);
	set_half_lives(in_rec);
	set_bifurcation(in_rec);

	if(nr_of_decays > 0){
		//pole rozpaduu
		Input::Array decay_array = in_rec.val<Input::Array>("decays");
		int dec_nr = -1;
		//iterator do cyklu
		for(Input::Iterator<Input::Record> dec_it = decay_array.begin<Input::Record>(); dec_it != decay_array.end(); ++dec_it, ++dec_nr)
		{
			index_par = substance_ids[dec_nr][0]; // pole indexu parentu, because indecees in input file run from one whereas indeces in C++ run from ZERO
			//int first_index = substance_ids[0];

			if(cols < (nr_of_isotopes - 1)){
				rel_step = time_step/half_lives[dec_nr];
			}
			reaction_matrix[index_par][index_par] = pow(0.5,rel_step);

			Input::Array prod_array = dec_it->val<Input::Array>("products");

			int i = 0;
			for(Input::Iterator<Input::Array> prod_it = prod_array.begin<Input::Array>(); prod_it != prod_array.end(); ++prod_it, ++i)
			{
				//bif_id = cols;
					reaction_matrix[index_par][substance_ids[dec_nr][i]] += (1 - pow(0.5,rel_step)) * bifurcation[dec_nr][substance_ids[dec_nr][i]];
			}
		}
	}
	print_reaction_matrix();//just for control print
	return reaction_matrix;
}

double **Linear_reaction::modify_reaction_matrix(void) //All the parameters are supposed to be known
{
	int rows,cols, index_par, bif_id;
	int *index_child;
	double rel_step, prev_rel_step;

	half_lives = (double *)xmalloc(nr_of_decays * sizeof(double));
	bifurcation.resize(nr_of_decays);

	if(reaction_matrix == NULL){
		xprintf(Msg,"\nReaction matrix pointer is NULL.\n");
		return NULL;
	}

	if(nr_of_decays > 0){
		//pole rozpadu nemuze byt pouzito
		for(int dec_nr = 0; dec_nr < nr_of_decays; dec_nr++)
		{
			index_par = substance_ids[dec_nr][0];

			if(cols < (nr_of_isotopes - 1)){
				rel_step = time_step/half_lives[dec_nr];
			}
			reaction_matrix[index_par][index_par] = pow(0.5,rel_step);

			int nr_of_indices = sizeof(*(substance_ids[dec_nr]))/sizeof(double);
			for(int i = 0; i < nr_of_indices; ++i)
			{
					reaction_matrix[index_par][substance_ids[dec_nr][i]] += (1 - pow(0.5,rel_step)) * bifurcation[dec_nr][substance_ids[dec_nr][i]];
			}
		}
	}
	print_reaction_matrix();//just for control print
	return reaction_matrix;
}

double **Linear_reaction::compute_reaction(double **concentrations, int loc_el) //multiplication of concentrations array by reaction matrix
{
    int cols, rows, both;

	if(nr_of_decays > 0){
		for(cols = 0; cols < nr_of_species; cols++){
		prev_conc[cols] = concentrations[cols][loc_el];
		concentrations[cols][loc_el] = 0.0;
		}
        for(rows = 0; rows <nr_of_species; rows++){
            for(cols = 0; cols <nr_of_species; cols++){
                concentrations[rows][loc_el] += prev_conc[cols] * reaction_matrix[cols][rows];
            }
        }
	}
	return concentrations;
}

double *Linear_reaction::set_half_lives(Input::Record in_rec)
{
	char  buffer[1024];
	char *pom_buf;
	int j = 0;
	int i=0;

	Input::Array decay_array = in_rec.val<Input::Array>("decays");
	int size = decay_array.size();

	if(half_lives != NULL){
			free(half_lives);
			half_lives = NULL;
	}
	if(half_lives == NULL){
		half_lives = (double *)xmalloc(size * sizeof(double));
	}

	for (Input::Iterator<Input::Record> it = decay_array.begin<Input::Record>(); it != decay_array.end(); ++it)
	{
	  Input::Iterator<double> it_hl = it->find<double>("half_life");
	  if (it_hl) {
           half_lives[i] = *it_hl;
	   } else {
	       it_hl = it->find<double>("kinetic");
	       if (it_hl) {
	    	   half_lives[i] = log(2)/(*it_hl);
	       } else {
	         xprintf(Msg, "You did not specify either the half life nor kinetic konstant for %d substep of decay chain.\n", i);
	         exit(1);
	       }
	   }
	  i++;
	}
	return half_lives;
}

void Linear_reaction::print_half_lives(int nr_of_substances)
{
	int i;

	if(half_lives == NULL)
	{
		xprintf(Msg,"\nHalf-lives are not defined.");
	}else{
		xprintf(Msg,"\nHalf-lives are defined as:");
		for(i=0; i < (nr_of_substances - 1) ; i++)
		{
			if(i < (nr_of_substances  - 2)) //cout << " " << half_lives[i] <<",";
				xprintf(Msg," %f", half_lives[i]);
			if(i == (nr_of_substances  - 2)) //cout << " " << half_lives[i] <<"\n";
				xprintf(Msg," %f\n", this->half_lives[i]);
		}
	}
	return;
}

int **Linear_reaction::set_indices(Input::Record in_rec) //(int index, int nr_of_substances)
{
	char  buffer[1024];
	char *pom_buf;
	int i,j;
	const char *separators = " ,\t";

	if(substance_ids != NULL){
		free(substance_ids);
		substance_ids = NULL;
	}
	if(substance_ids == NULL){
		substance_ids = (int **)xmalloc(nr_of_decays*sizeof(int*));
	}

	Input::Array names_array = in_rec.val<Input::Array>("substances");
	Input::Array decay_array = in_rec.val<Input::Array>("decays");

	int dec_nr = 0;
	for(Input::Iterator<Input::Record> dec_it = decay_array.begin<Input::Record>(); dec_it != decay_array.end(); ++dec_it)
	{
		string parent_name = dec_it->val<string>("parent");
		Input::Array bif_array = dec_it->val<Input::Array>("products");

		if(bif_array.size() > 0)
		{
			substance_ids[dec_nr] = (int*)xmalloc((bif_array.size() + 1)*sizeof(int)); //products Ids + parental atom Id
		}else{
			xprintf(Msg,"For some reason you forget to dedine products of %d-th reaction.\n", dec_nr);
			exit(1);
		}

		/*int pos = -1;
		i = 0;
		for(Input::Iterator<string> name_it = names_array.begin<string>(); name_it != names_array.end() && (pos == -1); ++name_it, ++i)
		{
			//if (strcmp(name_it.c_str(), parent_name.c_str()) == 0)
			if(parent_name.compare(*name_it) == 0)
			{
			        pos = i;
			}
		}*/
		int pos = find_index(names_array, parent_name);
		if(pos > -1)
		{
			substance_ids[dec_nr][0] = pos;
		}else{
			xprintf(Msg,"In %d decay chain substep you used undefined parental atom.", dec_nr);
			exit(1);
		}

		int prod_pos = 0;
		for(Input::Iterator<string> bif_it = bif_array.begin<string>(); bif_it != bif_array.end(); ++bif_it, ++prod_pos)
		{
			/*int pos = -1;
			i = 0;

			for(Input::Iterator<string> name_it = names_array.begin<string>(); name_it != names_array.end() && (pos == -1); ++name_it, ++i)
			{
				if((*bif_it).compare(*name_it) == 0) //if (strcmp(*name_it, *child_name) == 0)
				{
			        pos = i;
				}
			}*/
			int pos = find_index(names_array, *bif_it);
			if(pos > -1)
			{
				substance_ids[dec_nr][prod_pos] = pos;
			}else{
				xprintf(Msg,"In %d decay chain substep you used undefined %d-th product atom.", dec_nr, prod_pos);
				exit(1);
			}
		}

		dec_nr++;
	}
	 return substance_ids;
}

void Linear_reaction::print_indeces(int nr_of_substances)
{
	int i;

	if(substance_ids == NULL)
	{
		xprintf(Msg,"\nReaction/decay has not been defined.");
	}else{
		xprintf(Msg,"\nOrder of substences is defined by %d indeces:", nr_of_isotopes);
		for(i = 0; i < nr_of_substances ; i++)
		{
			if(i < (nr_of_substances  - 1)) xprintf(Msg," %d,",substance_ids[i]);
			if(i == (nr_of_substances  - 1)) xprintf(Msg," %d\n",substance_ids[i]);
		}
	}
	return;
}

void Linear_reaction::print_reaction_matrix(void)
{
	int cols,rows;

	if(reaction_matrix != NULL){
		xprintf(Msg,"\ntime_step %f,Reaction matrix looks as follows:\n",time_step);
		for(rows = 0; rows < nr_of_species; rows++){
			for(cols = 0; cols < nr_of_species; cols++){
				if(cols == (nr_of_species - 1)){
					xprintf(Msg,"%f\n",reaction_matrix[rows][cols]);
				}else{
					xprintf(Msg,"%f\t",reaction_matrix[rows][cols]);
				}
			}
		}
	}else{
		xprintf(Msg,"\nReaction matrix needs to be allocated.\n");
	}
	return;
}

void Linear_reaction::set_bifurcation(Input::Record in_rec) // (int index, Input::Record in_rec)
{

	Input::Array decay_array = in_rec.val<Input::Array>("decays");

		int dec_nr = 0;
		for (Input::Iterator<Input::Record> dec_it = decay_array.begin<Input::Record>(); dec_it != decay_array.end(); ++dec_it)
		{
			Input::Array bif_array = dec_it->val<Input::Array>("branch_ratios");
			int nr_of_prod = bif_array.size();
			bifurcation[dec_nr].resize(nr_of_prod + 1);
			//kdyz bude 1, tak bysem mel dát do radku bifurkaci jednicku, viz par radku nize
	    	if(bif_array.size() > 1)
	    	{
	    		bifurcation[dec_nr].resize(bif_array.size());
		    	bif_array.copy_to(bifurcation[dec_nr]); //atof(pom_buf); //control_sum += bifurcation[dec_nr][j];
			}else{
	    		bifurcation[dec_nr].resize(1);
	    		bifurcation[dec_nr][0]= 1.0;
	    	}
	    	dec_nr++;
		}
	return;
}

void Linear_reaction::set_time_step(double new_timestep, Input::Record in_rec)
{
	time_step = new_timestep;
	release_reaction_matrix();
	allocate_reaction_matrix();
	modify_reaction_matrix(in_rec); //_repeatedly(in_rec);
	return;
}

void Linear_reaction::set_time_step(Input::Record in_rec)
{
	time_step = in_rec.val<double>("time_step"); //OptGetDbl("Global","Save_step","1.0");
	release_reaction_matrix();
	allocate_reaction_matrix();
	modify_reaction_matrix(in_rec); //_repeatedly(in_rec);
	return;
}

void Linear_reaction::set_time_step(double new_timestep)
{
	time_step = new_timestep;
	release_reaction_matrix();
	allocate_reaction_matrix();
	modify_reaction_matrix(); //_repeatedly(in_rec);
	return;
}

void Linear_reaction::compute_one_step(void)
{
    if (reaction_matrix == NULL)   return;

    START_TIMER("decay_step");
	for (int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
	 {
	 	this->compute_reaction(concentration_matrix[MOBILE], loc_el);
	    if (dual_porosity_on == true) {
	     this->compute_reaction(concentration_matrix[IMMOBILE], loc_el);
	    }

	 }
    END_TIMER("decay_step");
	 return;
}

void Linear_reaction::release_reaction_matrix(void)
{
	int i;
	if(reaction_matrix != NULL)
	{
		for(i = 0; i < nr_of_isotopes; i++)
		{
			if(reaction_matrix[i] != NULL)
			{
				free(reaction_matrix[i]);
				reaction_matrix[i] = NULL;
			}
		}
		free(reaction_matrix);
		reaction_matrix = NULL;
	}
}

void Linear_reaction::set_nr_of_isotopes(int Nr_of_isotopes)
{
	nr_of_isotopes = Nr_of_isotopes;
	return;
}
