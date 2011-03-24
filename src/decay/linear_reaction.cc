#include "linear_reaction.hh"
#include "system.hh"
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>

using namespace std;

Linear_reaction::Linear_reaction()
	: half_lives(NULL), substance_ids(NULL), reaction_matrix(NULL), decaying_isotopes(NULL), bifurcation(NULL), bifurcation_on(false)
{
	/*react_type = decay;
	reaction_matrix = NULL;
	half_lives = NULL;
	substance_ids = NULL;*/
}

Linear_reaction::Linear_reaction(REACTION_TYPE type, int n_subst, char *section, double time_step)
	: half_lives(NULL), substance_ids(NULL), reaction_matrix(NULL), decaying_isotopes(NULL), bifurcation(NULL), bifurcation_on(false)
{
	//react_type = type;
	Set_nr_of_decays(section);
	Prepare_reaction_matrix(n_subst);
	Modify_reaction_matrix_repeatedly(n_subst, section, time_step);
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

	free(decaying_isotopes);
	decaying_isotopes = NULL;

	free(bifurcation);
	bifurcation = NULL;
}

double **Linear_reaction::Prepare_reaction_matrix(int n_subst) //reaction matrix initialization
{
	int index, rows, cols, dec_nr, dec_name_nr;
	char dec_name[30];

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

//double **Linear_reaction::Modify_reaction_matrix(double **reaction_matrix, int n_subst, double time_step) //prepare the matrix, which describes reactions
double **Linear_reaction::Modify_reaction_matrix(int n_subst, double time_step) //prepare the matrix, which describes reactions
{
	int rows,cols, index, prev_index;
	double rel_step, prev_rel_step, decrease;

	if(reaction_matrix == NULL){
		xprintf(Msg,"\nReaction matrix pointer is NULL.\n");
		return NULL;
	}

		for(cols = 0; cols < nr_of_isotopes; cols++){
			rel_step = time_step/half_lives[cols];
			index = substance_ids[cols] - 1; // because indecees in input file run from one whereas indeces in C++ run from ZERO
			decrease = pow(0.5,rel_step);
			if(cols > 0){
				reaction_matrix[prev_index][prev_index] -= pow(0.5,prev_rel_step);
				reaction_matrix[prev_index][index] += pow(0.5,prev_rel_step);
			}
			prev_rel_step = rel_step;
			prev_index = index;
		}
	Print_reaction_matrix(n_subst);//just for control print
	return reaction_matrix;
}

double **Linear_reaction::Modify_reaction_matrix(int n_subst, double time_step, int meaningless) //prepare the matrix, which describes reactions, takes bifurcation in acount
{
	int rows,cols, index, prev_index;
	double rel_step, decrease; //, prev_rel_step

	if(reaction_matrix == NULL){
		xprintf(Msg,"\nReaction matrix pointer is NULL.\n");
		return NULL;
	}

		for(cols = 0; cols < nr_of_isotopes; cols++){
			rel_step = time_step/half_lives[cols];
			index = substance_ids[cols] - 1; // because indecees in input file run from one whereas indeces in C++ run from ZERO
			decrease = pow(0.5,rel_step);
			if(cols > 0){
				reaction_matrix[prev_index][prev_index] -= bifurcation[prev_index] * pow(0.5,rel_step);
				reaction_matrix[prev_index][index] += bifurcation[prev_index] * pow(0.5,rel_step);
			}
			//prev_rel_step = rel_step;
			prev_index = index;
		}
	Print_reaction_matrix(n_subst);//just for control print
	return reaction_matrix;
}

double **Linear_reaction::Modify_reaction_matrix_repeatedly(int n_subst, char *section, double time_step)
{
	char dec_name[30];
	int rows, cols, dec_nr, dec_name_nr = 1;

	if(decaying_isotopes == NULL) Prepare_decaying_isotopes_ids(n_subst);
	xprintf(Msg,"\nNumber of decays is %d\n",nr_of_decays);
		for(dec_nr = 0; dec_nr < nr_of_decays; dec_nr++){
			sprintf(dec_name,"Decay_%d", dec_name_nr);
			Set_nr_of_isotopes(dec_name); //(section);
			Set_half_lives(dec_name); //(section);
			Set_indeces(dec_name); //(section);
			Get_indeces(); //just a control
			Get_half_lives(); //just a control
			Modify_decaying_isotopes_ids();//first place appearence of an isotope in [Decay_i]
			Set_bifurcation_on(dec_name);
			if(bifurcation_on == true){
				Set_bifurcation(dec_name);
				//place for a different implementation of Modify_reaction_matrix
				Modify_reaction_matrix(n_subst, time_step, 1);
			}else{
				Modify_reaction_matrix(n_subst, time_step);
			}
			//free(dec_name);
			dec_name_nr++;
		}
	//}
	return reaction_matrix;
}

double **Linear_reaction::Compute_reaction(double **concentrations, int n_subst, int loc_el) //multiplication of concentrations array by reaction matrix
{
	int cols, rows, both;
	double *prev_conc = (double *)xmalloc(n_subst * sizeof(double));

	if(reaction_matrix == NULL){
		xprintf(Msg,"\nReaction matrix pointer is NULL.\n");
		return NULL;
	}

	for(cols = 0; cols < n_subst; cols++){
		//xprintf(Msg,"\n%d. of %d substances concentration is %f\n", cols, n_subst, concentrations[cols][loc_el]);
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
	return concentrations;
}

REACTION_TYPE Linear_reaction::Set_reaction_type(REACTION_TYPE type) //change of reaction type should be conditionated by generation of new reaction matrix
{
	react_type = decay;
	return react_type;
}

REACTION_TYPE Linear_reaction::Get_reaction_type()
{
	xprintf(Msg,"Type of reaction is: %s.\n",react_type);
	return react_type;
}

int Linear_reaction::Set_nr_of_isotopes(char *section)
{
	nr_of_isotopes = OptGetInt(section,"Nr_of_isotopes","0");
	return nr_of_isotopes;
}

int Linear_reaction::Get_nr_of_isotopes()
{
	return nr_of_isotopes;
}

double *Linear_reaction::Set_half_lives(char *section)
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
		xprintf(Msg,"\nAllocation is permited, nr of isotopes %d", this->nr_of_isotopes);
		this->half_lives = (double *)xmalloc(this->nr_of_isotopes*sizeof(double));
		//this->half_lives = new double[this->nr_of_isotopes];
	}

	strcpy(buffer,OptGetStr(section,"Half_lives",NULL));
	pom_buf = strtok( buffer, separators );
	for (j=0; j< this->nr_of_isotopes; j++)
	{
		if ( pom_buf == NULL )
		{
			xprintf(Msg,"\nHalf-life of %d-th isotope is missing.", j+1);
		}
	    this->half_lives[j] = atof(pom_buf);
	    xprintf(Msg,"\n %d-th isotopes half-live is %f",j,this->half_lives[j]);
	    pom_buf = strtok( NULL, separators );
	 }
	 if ( pom_buf != NULL )
	 {
	    xprintf(Msg,"\nMore parameters then isotopes has been given. %d", 0);
	 }
	return this->half_lives;
}

double *Linear_reaction::Get_half_lives()
{
	int i;

	if(this->half_lives == NULL)
	{
		//std::cout << "\nHalf lives are not defined.\n";
		xprintf(Msg,"\nHalf-lives are not defined.");
	}else{
		//cout << "\nHalf lives are defined as:";
		xprintf(Msg,"\nHalf-lives are defined as:");
		for(i=0; i < this->nr_of_isotopes ; i++)
		{
			if(i < (this->nr_of_isotopes  - 1)) //cout << " " << half_lives[i] <<",";
				xprintf(Msg," %f", half_lives[i]);
			if(i == (this->nr_of_isotopes  - 1)) //cout << " " << half_lives[i] <<"\n";
				xprintf(Msg," %f\n", this->half_lives[i]);
		}
	}
	return half_lives;
}

int *Linear_reaction::Set_indeces(char *section)
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
		substance_ids = (int *)xmalloc(nr_of_isotopes*sizeof(int));
	}

	strcpy(buffer,OptGetStr(section,"Substance_ids",NULL));
	pom_buf = strtok( buffer, separators );
	for (j=0; j< nr_of_isotopes; j++)
	{
	  if ( pom_buf == NULL )
	  {
	    xprintf(Msg,"\nIndex for %d-th isotope is missing.", j+1);
	  }
	    substance_ids[j] = atoi(pom_buf);
	    //xprintf(Msg,"\n P_lat[%d].dGf %Lf\n",j,P_lat[j].dGf);
	    pom_buf = strtok( NULL, separators );
	 }
	 if ( pom_buf != NULL )
	 {
	    xprintf(Msg,"\nMore parameters then isotopes has been given.");
	 }

	 return substance_ids;
}

int *Linear_reaction::Get_indeces()
{
	int i;

	if(substance_ids == NULL)
	{
		xprintf(Msg,"\nDecay chain substances order is not defined.");
	}else{
		xprintf(Msg,"\nDecay chain substences order is defined by %d indeces:", nr_of_isotopes);
		for(i = 0; i < nr_of_isotopes ; i++)
		{
			if(i < (nr_of_isotopes  - 1)) xprintf(Msg," %d,",substance_ids[i]);
			if(i == (nr_of_isotopes  - 1)) xprintf(Msg," %d\n",substance_ids[i]);
		}
	}
	return substance_ids;
}

void Linear_reaction::Print_reaction_matrix(int n_subst)
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
	//return;
}

int Linear_reaction::Set_nr_of_decays(char *section)
{
	nr_of_decays = OptGetInt("Decay_module","Nr_of_decay_chains","1");
	return nr_of_decays;
}

void Linear_reaction::Prepare_decaying_isotopes_ids(int n_subst)
{
	int cols;
	decaying_isotopes = (int *)xmalloc(n_subst * sizeof(int));
	for(cols = 0; cols < n_subst; cols++){
		decaying_isotopes[cols] = 0;
	}
}

void Linear_reaction::Modify_decaying_isotopes_ids(void)
{
	int rows, cols, index, n_subst;

	//for(cols = 0; cols < (n_subst - 1); cols++){//upper limit depends on, if the last one substance in decay chain is stable or not, we think it is
		index = substance_ids[0] - 1;
		decaying_isotopes[index] += 1; //substance appears in some section [Decay_i] at the first position
		n_subst = sizeof(*decaying_isotopes)/sizeof(int);
		xprintf(Msg,"\nThese isotopes are members of decay chains: ");
		for(cols = 1; cols < n_subst; cols++){
			if(cols < (n_subst -1)){
				xprintf(Msg,"%d\t", decaying_isotopes[cols]);
			}else{
				xprintf(Msg,"%d\n", decaying_isotopes[cols]);
			}
		}
	//}
}

void Linear_reaction::Set_bifurcation(char *section)
{
	char  buffer[1024];
	char *pom_buf;
	int j;
	const char *separators = " ,\t";
	double control_sum = 0.0;

	if(bifurcation != NULL){
			free(bifurcation);
			bifurcation = NULL;
	}
	if(bifurcation == NULL){
		bifurcation = (double *)xmalloc(nr_of_isotopes*sizeof(double));
	}

	strcpy(buffer,OptGetStr(section,"Bifurcation",NULL));
	if(buffer == NULL) return;
	pom_buf = strtok( buffer, separators );
	for (j=0; j< nr_of_isotopes; j++)
	{
		if ( pom_buf == NULL )
		{
			xprintf(Msg,"\nBifurcation parameter of %d-th isotope is missing.", j+1);
		}
	    bifurcation[j] = atof(pom_buf);
	    xprintf(Msg,"\n %d-th isotopes bifurcation percentage is %f",j,bifurcation[j]);
	    pom_buf = strtok( NULL, separators );
	    if(j > 0)control_sum += bifurcation[j];
	 }
	if(control_sum != 1.0) xprintf(Msg,"\nSum of bifurcation parameters should be 1.0 but it is %f, because of mass conservation law.\n", control_sum);
	if( pom_buf != NULL )
	 {
	    xprintf(Msg,"\nMore parameters then isotopes has been given. %d", 0);
	 }
}

void Linear_reaction::Set_bifurcation_on(char *section)
{
	bifurcation_on = OptGetBool(section,"Bifurcation_on","no");
}
