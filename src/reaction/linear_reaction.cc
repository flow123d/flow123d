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
		rec.declare_key("branching_ratios", Array(Double()), Default::optional(),
				"Decay chain branching percentage.");
		rec.finish();
	}
	return rec;
}

Input::Type::Record & Linear_reaction::get_input_type()
{
	using namespace Input::Type;
	static Record rec("Linear_reactions", "Information for a decision about the way to simulate radioactive decay.");

	if (!rec.is_finished()) {
	    rec.derive_from( Reaction::get_input_type() );
        rec.declare_key("decays", Array( get_one_decay_substep() ), Default::obligatory(),
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
	//nr_of_isotopes = OptGetInt("Reaction_module","Nr_of_isotopes","0");
	allocate_reaction_matrix();
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

double **Linear_reaction::modify_reaction_matrix(void) //prepare the matrix, which describes reactions
{
	int rows,cols, index, prev_index;
	double rel_step, prev_rel_step;

	if(reaction_matrix == NULL){
		xprintf(Msg,"\nReaction matrix pointer is NULL.\n");
		return NULL;
	}
	if((nr_of_decays > 0) || (nr_of_FoR > 0)){
		for(cols = 0; cols < nr_of_isotopes; cols++){
			index = substance_ids[cols] - 1; // because indecees in input file run from one whereas indeces in C++ run from ZERO
			if(cols < (nr_of_isotopes - 1)){
				rel_step = time_step/half_lives[cols];
			}
			if(cols > 0){
				reaction_matrix[prev_index][prev_index] = pow(0.5,prev_rel_step);
				reaction_matrix[prev_index][index] += (1 - pow(0.5,prev_rel_step));
			}
			prev_rel_step = rel_step;
			prev_index = index;
		}
	}
	print_reaction_matrix();//just for control print
	return reaction_matrix;
}

double **Linear_reaction::modify_reaction_matrix(int dec_nr) //prepare the matrix, which describes reactions, takes bifurcation in acount
{
	int rows,cols, index, first_index, bif_id;
	double rel_step, prev_rel_step;

	if(reaction_matrix == NULL){
		xprintf(Msg,"\nReaction matrix pointer is NULL.\n");
		return NULL;
	}

	first_index = substance_ids[0]-1;
	for(cols = 0; cols < nr_of_isotopes; cols++){
		index = substance_ids[cols] - 1; // because indecees in input file run from one whereas indeces in C++ run from ZERO
		if(cols < (nr_of_isotopes -1)){
			rel_step = time_step/half_lives[cols];
			xprintf(Msg,"time_step %f\n", time_step);
		}
		if(cols > 0){
			bif_id = cols -1;
			reaction_matrix[first_index][first_index] = pow(0.5,prev_rel_step); //bifurcation[dec_nr][bif_id] * pow(0.5,prev_rel_step);
			reaction_matrix[first_index][index] += (1 - pow(0.5,prev_rel_step)) * bifurcation[dec_nr][bif_id];
		}
		prev_rel_step = rel_step;
	}
	print_reaction_matrix();//just for control print
	return reaction_matrix;
}

double **Linear_reaction::modify_reaction_matrix_repeatedly(void)
{
	char dec_name[30];
	int rows, cols, dec_nr, dec_name_nr = 1, index, prev_index;

	if(nr_of_decays > 0){
		xprintf(Msg,"\nNumber of decays is %d\n",nr_of_decays);
		if(half_lives != NULL){
					free(half_lives);
					half_lives = NULL;
		}
		half_lives = (double *)xmalloc(nr_of_decays * sizeof(double));
		bifurcation.resize(nr_of_decays);
		for(dec_nr = 0; dec_nr < nr_of_decays; dec_nr++){
			sprintf(dec_name,"Decay_%d", dec_name_nr);
			nr_of_isotopes = OptGetInt(dec_name,"Nr_of_isotopes","0");
			set_half_lives(dec_name);
			set_indeces(dec_name, nr_of_isotopes);
			print_indeces(nr_of_isotopes); //just a control
			print_half_lives(nr_of_isotopes); //just a control
			bifurcation_on = OptGetBool(dec_name,"Bifurcation_on","no");
			if(bifurcation_on == true){
				set_bifurcation(dec_name, dec_nr);
				modify_reaction_matrix(dec_nr);
			}else{
				modify_reaction_matrix();
			}
			dec_name_nr++;
		}
	}
	if(nr_of_FoR > 0){
		xprintf(Msg,"\nNumber of first order reactions is %d\n",nr_of_FoR);
		//half_lives.resize(nr_of_FoR); //does not function at all
		if(half_lives != NULL){
			free(half_lives);
			half_lives = NULL;
		}
		half_lives = (double *)xmalloc(nr_of_FoR * sizeof(double));
		for(dec_nr = 0; dec_nr < nr_of_FoR; dec_nr++){
			sprintf(dec_name,"FoReact_%d", dec_name_nr);
			set_nr_of_isotopes(2);
			set_indeces(dec_name, 2);
			set_kinetic_constants(dec_name, dec_nr);//instead of this line, here should be palced computation of halflives using kinetic constants
			print_indeces(nr_of_isotopes); //just a control
			print_half_lives(2); //just a control
			//modify_reaction_matrix(2);
			modify_reaction_matrix();
			dec_name_nr++;
		}
	}
	return reaction_matrix;
}

double **Linear_reaction::compute_reaction(double **concentrations, int loc_el) //multiplication of concentrations array by reaction matrix
{


    int cols, rows, both;

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
	}
	return concentrations;
}

double *Linear_reaction::set_half_lives(char *decname)//(Input::Record in_rec)
{
	char  buffer[1024];
	char *pom_buf;
	int i,j = 0;
	//const char *separators = " ,\t";

	if(half_lives != NULL){
			free(half_lives);
			half_lives = NULL;
	}
	//(nr_of_isotopes - 1) je velikost pole zadanych rozpadu
	if(half_lives == NULL){
		half_lives = (double *)xmalloc((nr_of_isotopes - 1) * sizeof(double));
	}

	 //strcpy(buffer,OptGetStr(section,"Half_lives",NULL));
	 //pom_buf = strtok( buffer, separators );

	//nutno projit pole rozpadu a ptat se prubezne na polocasy
	//Input::Array decay_array = in_rec.val<Input::Array>("decays");
	// int size = decay_array.size();
	// int i=0;
	 //for (Input::Iterator<Input::Record> it = decay_array.begin<Input::Record>(); it != dacay_array.end(); ++it, ++i) {
	 //  half_lives[i] = it->val<double>("half_life")  /// pouzit find
	 /// Input::Iterator<double> it_hl = it->find<double>("half_life");
	///  if (it_hl) {
    //	        half_lives[i] = *it_hl;
	//   } else {
	//        it_hl = it->find<double>("kinetic");
	//        if (it_hl) {
	//        } else {
	//           error
	//        }


		//if ( pom_buf == NULL )
		{
			xprintf(Msg,"\nHalf-life of %d-th isotope is missing.", j+1);
		}
	    //half_lives[j] = atof(pom_buf);
		//half_lives[j] =
	    xprintf(Msg,"\n %d-th isotopes half-live is %f",j,half_lives[j]);
	    //pom_buf = strtok( NULL, separators );
	 //}
	 /*if ( pom_buf != NULL )
	 {
	    xprintf(Msg,"\nMore parameters then (isotopes -1) has been given. %d", 0);
	 }*/
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

int *Linear_reaction::set_indeces(char *section, int nr_of_substances)
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
		substance_ids = (int *)xmalloc(nr_of_substances*sizeof(int));
	}

	strcpy(buffer,OptGetStr(section,"Substance_ids",NULL));
	pom_buf = strtok( buffer, separators );
	for (j=0; j< nr_of_substances; j++)
	{
	  if ( pom_buf == NULL )
	  {
	    xprintf(Msg,"\nIndex for %d-th substance in %s is missing.", j+1, section);
	  }
	    substance_ids[j] = atoi(pom_buf);
	    pom_buf = strtok( NULL, separators );
	 }
	 if ( pom_buf != NULL )
	 {
	    xprintf(Msg,"\nMore parameters then substances has been given in %s.", section);
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

void Linear_reaction::set_bifurcation(char *section, int dec_nr)
{
	char  buffer[1024];
	char *pom_buf;
	int j;
	const char *separators = " ,\t";
	double control_sum = 0.0;

	if(bifurcation_on == true)
	{
		bifurcation[dec_nr].resize(nr_of_isotopes - 1);
		strcpy(buffer,OptGetStr(section,"Bifurcation",NULL));
		if(buffer == NULL) return;
		pom_buf = strtok( buffer, separators );
		for (j=0; j< (nr_of_isotopes - 1); j++)
		{
			if ( pom_buf == NULL )
			{
				xprintf(Msg,"\nBifurcation parameter of %d-th isotope is missing.", j+1);
			}
	    	bifurcation[dec_nr][j] = atof(pom_buf);
	    	xprintf(Msg,"\n %d-th isotopes bifurcation percentage is %f",j,bifurcation[dec_nr][j]);
	    	pom_buf = strtok( NULL, separators );
	    	if(j > 0)control_sum += bifurcation[dec_nr][j];
	 	 }
	}else{
		bifurcation[dec_nr].resize(1);
		bifurcation[dec_nr][0]= 1.0;
	}
	if(control_sum != 1.0) xprintf(Msg,"\nSum of bifurcation parameters should be 1.0 but it is %f, because of mass conservation law.\n", control_sum);
	if( pom_buf != NULL )
	 {
	    xprintf(Msg,"\nMore parameters then (isotopes -1) has been given. %d", 0);
	 }
	return;
}

void Linear_reaction::set_kinetic_constants(char *section, int react_nr)
{
	char  buffer[1024];
	char *pom_buf;
	//int j;
	const char *separators = " ,\t";

	kinetic_constant.resize(nr_of_FoR);
	strcpy(buffer,OptGetStr(section,"Kinetic_constant",NULL));
	if(buffer == NULL) return;
	pom_buf = strtok( buffer, separators );
	//for (j=0; j< (nr_of_FoR); j++){
		if ( pom_buf == NULL )
		{
			xprintf(Msg,"\nKinetic constant belonging to %d-th reactions is missing.", react_nr+1);
		}
    	kinetic_constant[react_nr] = atof(pom_buf);
    	xprintf(Msg,"\nKinetic constant for %d-th reaction is %f",react_nr,kinetic_constant[react_nr]);
    	pom_buf = strtok( NULL, separators );
    	half_lives[react_nr] = log(2) / kinetic_constant[react_nr];
 	 //}
    return;
}

void Linear_reaction::set_time_step(double new_timestep){
	time_step = new_timestep;
	release_reaction_matrix();
	allocate_reaction_matrix();
	modify_reaction_matrix_repeatedly();
	return;
}

void Linear_reaction::set_time_step(void)
{
	time_step = OptGetDbl("Global","Save_step","1.0");
	release_reaction_matrix();
	allocate_reaction_matrix();
	modify_reaction_matrix_repeatedly();
	return;
}

void Linear_reaction::compute_one_step(void)
{
    if (reaction_matrix == NULL)   return;

    START_TIMER("decay_step");
	 //for (int loc_el = 0; loc_el < distribution->lsize(distribution->myp()); loc_el++)
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
