#include "linear_reaction.hh"
#include "system.hh"
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>

using namespace std;

Linear_reaction::Linear_reaction()
{
	react_type = decay;
	reaction_matrix = NULL;
	half_lives = NULL;
	substance_ids = NULL;
}

Linear_reaction::Linear_reaction(REACTION_TYPE type, int n_subst, char *section)
	: half_lives(NULL), substance_ids(NULL), reaction_matrix(NULL)
{
	react_type = type;
	reaction_matrix = reaction_matrix;
	nr_of_isotopes = Set_nr_of_isotopes(section);
	half_lives = Set_half_lives(section);
	substance_ids = Set_indeces(section);
	reaction_matrix = Prepare_reaction_matrix(n_subst);
}

Linear_reaction::~Linear_reaction()
{
	int i;

	xprintf(Msg,"\nDestructor is running.");
	free(half_lives);
	half_lives = NULL;

	free(substance_ids);
	this->substance_ids = NULL;

	free(reaction_matrix);
	reaction_matrix = NULL;
}

double **Linear_reaction::Prepare_reaction_matrix(int n_subst) //reaction matrix initialization
{
	int index, rows, cols;

	this->reaction_matrix = (double **)xmalloc(n_subst * sizeof(double*));//allocation section
	for(rows = 0; rows < n_subst; rows++){
		reaction_matrix[rows] = (double *)xmalloc(n_subst * sizeof(double));
	}

	for(rows = 0; rows < n_subst;rows++){
		for(cols = 0; cols < n_subst; cols++) reaction_matrix[rows][cols] = 0.0;
	}

	for(rows = 0; rows < nr_of_isotopes;rows++){
		index = substance_ids[rows] - 1; // because indecees in input file run from one whereas indeces in C++ run from ZERO
		reaction_matrix[index][index] = 1.0;
	}
	return reaction_matrix;
}

//double **Linear_reaction::Modify_reaction_matrix(double **reaction_matrix, int n_subst, double time_step) //prepare the matrix, which describes reactions
double **Linear_reaction::Modify_reaction_matrix(int n_subst, double time_step) //prepare the matrix, which describes reactions
{
	int rows,cols, index, prev_index;
	double rel_step, prev_rel_step;

	if(reaction_matrix == NULL){
		xprintf(Msg,"\nReaction matrix pointer is NULL.\n");
		return NULL;
	}
		for(cols = 0; cols < nr_of_isotopes; cols++){
			rel_step = time_step/half_lives[cols];
			index = substance_ids[cols] - 1; // because indecees in input file run from one whereas indeces in C++ run from ZERO
			if(cols == 0){
				reaction_matrix[index][index] -= pow(0.5,rel_step);
			}else{
				reaction_matrix[index][index] -= pow(0.5,rel_step);
				reaction_matrix[prev_index][index] += pow(0.5,prev_rel_step);
			}
			prev_rel_step = rel_step;
			prev_index = index;
		}
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
		xprintf(Msg,"\n%d. of %d substances concentration is %f\n", cols, n_subst, prev_conc[cols]);
		concentrations[cols][loc_el] = 0.0;
	}
	for(rows = 0; rows < n_subst; rows++){
		for(cols = 0; cols < n_subst; cols++){
			concentrations[rows][loc_el] += prev_conc[cols] * reaction_matrix[cols][rows];
		}
	}
	free(prev_conc);
	prev_conc = NULL;
	return concentrations;
}

REACTION_TYPE Linear_reaction::Set_reaction_type(REACTION_TYPE type) //change of reaction type should be conditionated by generation of new reaction matrix
{
	this->react_type = decay;
	return this->react_type;
}

REACTION_TYPE Linear_reaction::Get_reaction_type()
{
	//std::cout << "\nType of reaction is: " << react_type <<"\n";
	xprintf(Msg,"Type of reaction is: %s.\n",this->react_type);
	return this->react_type;
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

	if(this->half_lives == NULL){
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

	if(this->substance_ids == NULL) this->substance_ids = (int *)xmalloc(this->nr_of_isotopes*sizeof(int));

	strcpy(buffer,OptGetStr(section,"Substance_ids",NULL));
	pom_buf = strtok( buffer, separators );
	for (j=0; j< this->nr_of_isotopes; j++)
	{
	  if ( pom_buf == NULL )
	  {
	    xprintf(Msg,"\nIndex for %d-th isotope is missing.", j+1);
	  }
	    this->substance_ids[j] = atoi(pom_buf);
	    //xprintf(Msg,"\n P_lat[%d].dGf %Lf\n",j,P_lat[j].dGf);
	    pom_buf = strtok( NULL, separators );
	 }
	 if ( pom_buf != NULL )
	 {
	    xprintf(Msg,"\nMore parameters then isotopes has been given.");
	 }

	 return this->substance_ids;
}

int *Linear_reaction::Get_indeces()
{
	int i;

	if(substance_ids == NULL)
	{
		//std::cout << "\nDecay chain substances order is not defined.\n";
		xprintf(Msg,"\nDecay chain substances order is not defined.");
	}else{
		//std::cout << "\nDecay chain substences order is defined by " << nr_of_isotopes << ", " << sizeof(substance_ids) <<", "<< sizeof(*substance_ids) << " indeces:";
		xprintf(Msg,"\nDecay chain substences order is defined by %d indeces:", this->nr_of_isotopes);
		for(i = 0; i < this->nr_of_isotopes ; i++)
		{
			if(i < (this->nr_of_isotopes  - 1)) xprintf(Msg," %d,",this->substance_ids[i]); //std::cout << " " << substance_ids[i] <<",";
			if(i == (this->nr_of_isotopes  - 1)) xprintf(Msg," %d\n",this->substance_ids[i]); //std::cout << " " << substance_ids[i] <<"\n";
		}
	}
	return this->substance_ids;
}

