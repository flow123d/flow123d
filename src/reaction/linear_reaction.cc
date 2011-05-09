#include "linear_reaction.hh"
#include "system/system.hh"
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>

using namespace std;

Linear_reaction::Linear_reaction(int n_subst, double time_step)
	: half_lives(NULL), substance_ids(NULL), reaction_matrix(NULL), bifurcation_on(false)
{
	//decay_on = OptGetBool("Reactions_module","Compute_decay","no");
	//FoR_on = OptGetBool("Reactions_module","Compute_reactions","no");
	//if(decay_on == true)
	nr_of_decays = OptGetInt("Reaction_module","Nr_of_decay_chains","0");
	//if(FoR_on == true)
	nr_of_FoR = OptGetInt("Reaction_module","Nr_of_FoR","0");
	cout << "number of FoR is"<< nr_of_FoR << endl;
	cout << "number of decays is" << nr_of_decays << endl;
	if((nr_of_decays > 0) || (nr_of_FoR > 0)){
		allocate_reaction_matrix(n_subst);
		modify_reaction_matrix_repeatedly(n_subst, time_step);
	}
}

Linear_reaction::~Linear_reaction()
{
	int i, rows, n_subst;

	n_subst = sizeof(*reaction_matrix)/sizeof(double *);
	xprintf(Msg,"\nDestructor is running.");
	free(half_lives);
	half_lives = NULL;

	free(substance_ids);
	this->substance_ids = NULL;

	free(reaction_matrix);
	reaction_matrix = NULL;

	for(rows = 0; rows < n_subst;rows++){
		free(reaction_matrix[rows]);
		reaction_matrix[rows] = NULL;
	}
	free(reaction_matrix);
	reaction_matrix = NULL;
}

double **Linear_reaction::allocate_reaction_matrix(int n_subst) //reaction matrix initialization
{
	int index, rows, cols, dec_nr, dec_name_nr;
	char dec_name[30];

	cout << "We are going to allocate reaction matrix" << endl;
	reaction_matrix = (double **)xmalloc(n_subst * sizeof(double*));//allocation section
	for(rows = 0; rows < n_subst; rows++){
		reaction_matrix[rows] = (double *)xmalloc(n_subst * sizeof(double));
	}

	for(rows = 0; rows < n_subst;rows++){
		for(cols = 0; cols < n_subst; cols++)
			if(cols == rows){
				reaction_matrix[rows][cols] = 1.0;
			}else{
				reaction_matrix[rows][cols] = 0.0;
			}
	}

	dec_name_nr = 1;
	return reaction_matrix;
}

double **Linear_reaction::modify_reaction_matrix(int n_subst, int nr_of_participants, double time_step) //prepare the matrix, which describes reactions
{
	int rows,cols, index, prev_index;
	double rel_step, prev_rel_step;

	if(reaction_matrix == NULL){
		xprintf(Msg,"\nReaction matrix pointer is NULL.\n");
		return NULL;
	}
	if((nr_of_decays > 0) || (nr_of_FoR > 0)){
		for(cols = 0; cols < nr_of_participants; cols++){
			index = substance_ids[cols] - 1; // because indecees in input file run from one whereas indeces in C++ run from ZERO
			if(cols < (nr_of_participants - 1)){
				rel_step = time_step/half_lives[cols];
			}
			if(cols > 0){
				reaction_matrix[prev_index][prev_index] -= pow(0.5,prev_rel_step);
				reaction_matrix[prev_index][index] += pow(0.5,prev_rel_step);
			}
			prev_rel_step = rel_step;
			prev_index = index;
		}
	}
	print_reaction_matrix(n_subst);//just for control print
	return reaction_matrix;
}

double **Linear_reaction::modify_reaction_matrix(int n_subst, double time_step, int dec_nr) //prepare the matrix, which describes reactions, takes bifurcation in acount
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
		}
		if(cols > 0){
			bif_id = cols -1;
			reaction_matrix[first_index][first_index] -= bifurcation[dec_nr][bif_id] * pow(0.5,prev_rel_step);
			reaction_matrix[first_index][index] += bifurcation[dec_nr][bif_id] * pow(0.5,prev_rel_step);
		}
		prev_rel_step = rel_step;
	}
	print_reaction_matrix(n_subst);//just for control print
	return reaction_matrix;
}

double **Linear_reaction::modify_reaction_matrix_repeatedly(int n_subst, double time_step)
{
	char dec_name[30];
	int rows, cols, dec_nr, dec_name_nr = 1;

	if(nr_of_decays > 0){
		xprintf(Msg,"\nNumber of decays is %d\n",nr_of_decays);
		bifurcation.resize(nr_of_decays);
		for(dec_nr = 0; dec_nr < nr_of_decays; dec_nr++){
			sprintf(dec_name,"Decay_%d", dec_name_nr);
			nr_of_isotopes = OptGetInt(dec_name,"Nr_of_isotopes","0");
			set_half_lives(dec_name, nr_of_isotopes);
			set_indeces(dec_name, nr_of_isotopes);
			print_indeces(nr_of_decays); //just a control
			print_half_lives(nr_of_isotopes); //just a control
			bifurcation_on = OptGetBool(dec_name,"Compute_decay","no");
			if(bifurcation_on == true){
				set_bifurcation(dec_name, dec_nr);
				modify_reaction_matrix(n_subst, time_step, dec_nr);
			}else{
				modify_reaction_matrix(n_subst, nr_of_isotopes, time_step);
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
			set_indeces(dec_name, 2);
			set_kinetic_constants(dec_name, dec_nr);//instead of this line, here should be palced computation of halflives using kinetic constants
			print_indeces(nr_of_FoR); //just a control
			print_half_lives(2); //just a control
			modify_reaction_matrix(n_subst, 2, time_step);
			dec_name_nr++;
		}
	}
	return reaction_matrix;
}

double **Linear_reaction::compute_reaction(double **concentrations, int n_subst, int loc_el) //multiplication of concentrations array by reaction matrix
{
	int cols, rows, both;
	double *prev_conc = (double *)xmalloc(n_subst * sizeof(double));

	if(reaction_matrix == NULL){
		xprintf(Msg,"\nReaction matrix pointer is NULL.\n");
		return NULL;
	}

	if((nr_of_decays > 0) || (nr_of_FoR > 0)){
		for(cols = 0; cols < n_subst; cols++){
		prev_conc[cols] = concentrations[cols][loc_el];
		xprintf(Msg,"\n%d. of %d substances concentration is %f\n", cols, n_subst, concentrations[cols][loc_el]); //prev_conc[cols]);
		concentrations[cols][loc_el] = 0.0;
	}
	for(rows = 0; rows < n_subst; rows++){
		for(cols = 0; cols < n_subst; cols++){
			concentrations[rows][loc_el] += prev_conc[cols] * reaction_matrix[cols][rows];
		}
		xprintf(Msg,"\n%d. of %d substances concentration after reaction is %f\n", rows, n_subst, concentrations[rows][loc_el]);
	}
	free(prev_conc);
	prev_conc = NULL;
	}
	return concentrations;
}

double *Linear_reaction::set_half_lives(char *section, int nr_of_substances)
{
	char  buffer[1024];
	char *pom_buf;
	int i,j;
	const char *separators = " ,\t";

	if(half_lives != NULL){
			free(half_lives);
			half_lives = NULL;
	}
	if(half_lives == NULL){
		//xprintf(Msg,"\nAllocation is permited, nr of isotopes %d", nr_of_isotopes);
		half_lives = (double *)xmalloc((nr_of_substances - 1) * sizeof(double));
	}
	 strcpy(buffer,OptGetStr(section,"Half_lives",NULL));
	 pom_buf = strtok( buffer, separators );
	 for (j=0; j< (nr_of_isotopes-1); j++){
		if ( pom_buf == NULL )
		{
			xprintf(Msg,"\nHalf-life of %d-th isotope is missing.", j+1);
		}
	    half_lives[j] = atof(pom_buf);
	    xprintf(Msg,"\n %d-th isotopes half-live is %f",j,half_lives[j]);
	    pom_buf = strtok( NULL, separators );
	 }
	 if ( pom_buf != NULL )
	 {
	    xprintf(Msg,"\nMore parameters then (isotopes -1) has been given. %d", 0);
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
}

void Linear_reaction::print_reaction_matrix(int n_subst)
{
	int cols,rows;

	xprintf(Msg,"\nReaction matrix looks as follows:\n");
	for(rows = 0; rows < n_subst; rows++){
		for(cols = 0; cols < n_subst; cols++){
			if(cols == (n_subst - 1)){
				xprintf(Msg,"%f\n",reaction_matrix[rows][cols]);
			}else{
				xprintf(Msg,"%f\t",reaction_matrix[rows][cols]);
			}
		}
	}
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
}

int Linear_reaction::get_nr_of_decays(void){return nr_of_decays;} // two simple inlinefunction returning private variables
int Linear_reaction::get_nr_of_FoR(void){return nr_of_FoR;}
