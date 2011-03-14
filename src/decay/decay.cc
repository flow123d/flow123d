#include "decay.h"
#include "system.hh"
#include <iostream>
#include <cstring>
#include <stdlib.h>

using namespace std;

Linear_reaction::Linear_reaction()
{
	react_type = decay;
	reaction_matrix = NULL;
	half_lives = NULL;
	substance_ids = NULL;
}

Linear_reaction::Linear_reaction(REACTION_TYPE type, double *half_lives, int n_subst)
{
	react_type = type;
	reaction_matrix = Create_reaction_matrix(half_lives, n_subst);
}

Linear_reaction::~Linear_reaction()
{
	int i;

	xprintf(Msg,"\nDestructor is running.");
	free(half_lives);
	half_lives = NULL;

	free(substance_ids);
	substance_ids = NULL;

	free(reaction_matrix);
	reaction_matrix = NULL;
}

double **Linear_reaction::Create_reaction_matrix(double *half_lives, int n_subst) //prepare the matrix, which describes reactions
{
	int i,j;
	double **reaction_matrix = (double**)xmalloc(n_subst*sizeof(double*));
	for(i=0; i<n_subst;i++){
		reaction_matrix[i] = (double*)xmalloc(n_subst*sizeof(double));
		for(j=0; j<n_subst; j++){
			reaction_matrix[i][j] = 0.0; // = half_lives[..]
		}
	}
	return reaction_matrix;
}

double **Linear_reaction::Compute_reaction(double **concentrations) //multiplication of concentrations array by reaction matrix
{
	return concentrations;
}

REACTION_TYPE Linear_reaction::Set_reaction_type(REACTION_TYPE type) //change of reaction type should be conditionated by generation of new reaction matrix
{
	return type;
}

REACTION_TYPE Linear_reaction::Get_reaction_type()
{
	//std::cout << "\nType of reaction is: " << react_type <<"\n";
	xprintf(Msg,"Type of reaction is: %s.\n",react_type);
	return react_type;
}

int Linear_reaction::Set_nr_of_isotopes(char *section)
{
	nr_of_isotopes = OptGetInt(section,"nr_of_isotopes","0");
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

	if(half_lives == NULL) half_lives = (double *)xmalloc(nr_of_isotopes*sizeof(int));

	strcpy(buffer,OptGetStr(section,"Half_lives",NULL));
	pom_buf = strtok( buffer, separators );
	for (j=0; j< nr_of_isotopes; j++)
	{
	  if ( pom_buf == NULL )
	  {
	    xprintf(Msg,"\nHalf-life of %d-th isotope is missing.", j+1);
	  }
	    half_lives[j] = atof(pom_buf);
	    //xprintf(Msg,"\n P_lat[%d].dGf %Lf\n",j,P_lat[j].dGf);
	    pom_buf = strtok( NULL, separators );
	 }
	 if ( pom_buf != NULL )
	 {
	    xprintf(Msg,"\nMore parameters then isotopes has been given. %d", 0);
	 }
	return half_lives;
}

double *Linear_reaction::Get_half_lives()
{
	int i;

	if(half_lives == NULL)
	{
		//std::cout << "\nHalf lives are not defined.\n";
		xprintf(Msg,"\nHalf-lives are not defined.");
	}else{
		//cout << "\nHalf lives are defined as:";
		xprintf(Msg,"\nHalf-lives are defined as:");
		for(i=0; i<nr_of_isotopes ; i++)
		{
			if(i<(nr_of_isotopes  - 1)) //cout << " " << half_lives[i] <<",";
				xprintf(Msg," %f", half_lives[i]);
			if(i == (nr_of_isotopes  - 1)) //cout << " " << half_lives[i] <<"\n";
				xprintf(Msg," %f\n", half_lives[i]);
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

	if(substance_ids == NULL) substance_ids = (int *)xmalloc(nr_of_isotopes*sizeof(int));

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
		//std::cout << "\nDecay chain substances order is not defined.\n";
		xprintf(Msg,"\nDecay chain substances order is not defined.");
	}else{
		//std::cout << "\nDecay chain substences order is defined by " << nr_of_isotopes << ", " << sizeof(substance_ids) <<", "<< sizeof(*substance_ids) << " indeces:";
		xprintf(Msg,"\nDecay chain substences order is defined by %d indeces:", nr_of_isotopes);
		for(i=0; i<nr_of_isotopes ; i++)
		{
			if(i<(nr_of_isotopes  - 1)) xprintf(Msg," %d,",substance_ids[i]); //std::cout << " " << substance_ids[i] <<",";
			if(i == (nr_of_isotopes  - 1)) xprintf(Msg," %d\n",substance_ids[i]); //std::cout << " " << substance_ids[i] <<"\n";
		}
	}
	return substance_ids;
}

