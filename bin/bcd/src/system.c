#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "struct.h"
#include "bcd.h"


//=============================================================================
//		SKIP TO
//=============================================================================
void skip_to( FILE *in, char *sekce)
{
	char line[ MAXBUFF ];
	char string[ MAXBUFF ];
	while( fgets( line, MAXBUFF - 2, in ) != NULL ) {
		sscanf( line, "%s", string );
                if( strcmp( string, sekce ) == 0 )
                	return;
	}
}
//=============================================================================
//		OPEN FILE
//=============================================================================
FILE *open_file(char *file,char *open_param){

        FILE *out;
        out = fopen( file, open_param );
        if(out == NULL){
        	printf("Can't open file %s",file);
        	getchar();
        	exit(EXIT_FAILURE);
        }
        return out;
}
//=============================================================================
// MAKE BRAND NEW COPY OF STRING
//=============================================================================
char *xstrcpy( char *src )
{
	char *rc;
	int length;
	length = strlen( src ) + 1;
	rc = (char*) malloc( length * sizeof( char ) );
	strcpy( rc, src );
	return rc;
}
//=============================================================================
//		GET STRING PARAMETER
//=============================================================================
char* get_string(char *file,char *section,char *key,char* def_value)
{
	FILE *INI;
	char line[MAXBUFF],*string;

	INI = open_file(file,"rt");
	skip_to(INI,section);
	while(fgets(line,MAXBUFF-2,INI) != NULL){
		if(line[0] == '$') break;
		string = strtok(line," =\t\n");
		if(string == NULL) continue;
		if(strcmp(string,key) == 0){
			string = strtok(NULL," =\t\n\r");
			fclose(INI);
			return xstrcpy(string);
		}
	}
	fclose(INI);
	return def_value;
}
//=============================================================================
//		GET INTEGER PARAMETER
//=============================================================================
int get_int(char *file,char *section,char *key,int def_value)
{
	FILE *INI;
	char line[MAXBUFF],*string;
	int n;

	INI = open_file(file,"rt");
	skip_to(INI,section);
	while(fgets(line,MAXBUFF-2,INI) != NULL){
		if(line[0] == '$') break;
		string = strtok(line," =\t\n");
		if(string == NULL) continue;
		if(strcmp(string,key) == 0){
			n = atoi(strtok(NULL," =\t\n\r"));
			fclose(INI);
			return n;
		}
	}
	fclose(INI);
	return def_value;
}
//=============================================================================
//		GET FLOAT PARAMETER
//=============================================================================
double get_float(char *file,char *section,char *key,double def_value)
{
	FILE *INI;
	char line[MAXBUFF],*string;
	double n;

	INI = open_file(file,"rt");
	skip_to(INI,section);
	while(fgets(line,MAXBUFF-2,INI) != NULL){
		if(line[0] == '$') break;
		string = strtok(line," =\t\n");
		if(string == NULL) continue;
		if(strcmp(string,key) == 0){
			n = atof(strtok(NULL," =\t\n\r"));
			fclose(INI);
			return n;
		}
	}
	fclose(INI);
	return def_value;
}
//=============================================================================
//		GET BOOLEAN PARAMETER
//=============================================================================
int get_bool(char *file,char *section,char *key,double def_value)
{
	FILE *INI;
	char line[MAXBUFF],*string,*string1;
	int n;

	INI = open_file(file,"rt");
	skip_to(INI,section);
	while(fgets(line,MAXBUFF-2,INI) != NULL){
		if(line[0] == '$') break;
		string = strtok(line," =\t\n");
		if(string == NULL) continue;
		if(strcmp(string,key) == 0){
			string1 = strtok(NULL," =\t\n\r");
			if(strcmp(string1,"YES") == 0)
				n = 1;
			else
				if(strcmp(string1,"NO") == 0)
					n = 0;
				else{
					printf("%s%s%s%s","\tKey: ",string,", Invalid parameter: ",string1);
					getchar();
					exit(EXIT_FAILURE);
				}
			fclose(INI);
			return n;
		}
	}
	fclose(INI);
	return def_value;
}
//=============================================================================
//		GET CONDITION'S DATA
//=============================================================================
void get_cond(struct Problem *problem,char *secstart,char *secstop,int type)
{
	FILE *INI;
	char line[MAXBUFF],line1[MAXBUFF];//,*string;//,end[MAXBUFF],*string;
	int i;
	struct Condition *cd;

	cd = NULL;

	INI = open_file(problem->ini_file,"rt");
	skip_to(INI,secstart);

	/*if(strcmp(secstart,"$FBC") == 0)
		type = 1;
	if(strcmp(secstart,"$TBC") == 0)
		type = 2;
	if(type == 0)
		exit(EXIT_FAILURE);*/

	fgets(line,MAXBUFF-2,INI);
	sscanf(line,"%s" ,line1);
	while(strcmp(line1,secstop) != 0){
	switch(type){
		case 1:
			cd = new_BC(&problem->fbc,cd,type);
			break;
		case 2:
			cd = new_BC(&problem->tbc,cd,type);
			break;
		case 3:
			cd = new_BC(&problem->tic,cd,type);
			break;
		default:
			exit(EXIT_FAILURE);
			break;
	}


		cd->seting_type = atoi(strtok(line," \t"));
		switch(cd->seting_type){
			case 0:
				break;
			case 1:   		// Element list
				cd->n_elm = atoi(strtok(NULL," \t"));
				cd->elm_list = (int*)malloc(cd->n_elm*sizeof(int));
				for(i=0;i<cd->n_elm;i++)
					cd->elm_list[i] = atoi(strtok(NULL," \t"));
				break;
			case 2: 		// Equation
				cd->boundary_equation = (double*)malloc(4*sizeof(double));
				for(i=0;i<4;i++)
					cd->boundary_equation[i] = atof(strtok(NULL," \t"));
				break;
			case 3:  		// Material
				cd->mid = atoi(strtok(NULL," \t"));
				break;
			case 4:			// Equation + Material
				cd->boundary_equation = (double*)malloc(4*sizeof(double));
					for(i=0;i<4;i++)
				cd->boundary_equation[i] = atof(strtok(NULL," \t"));
				cd->mid = atoi(strtok(NULL," \t"));
				break;
			default:
				printf("\n%s",line);
				printf("%s%d","Unknown condition seting type #",cd->seting_type);
				getchar();
				exit(EXIT_FAILURE);
				break;
		}

		if(cd->condition_type == 1){	//FBC
			cd->bc_type = atoi(strtok(NULL," \t"));
			if(supp_bc_type(cd->bc_type) == 1){
				if(cd->bc_type == 3)
					cd->bc_param = atof(strtok(NULL," \t"));
			}
			else{
				printf("%s%d","Unknown BC type #",cd->bc_type);
				getchar();
				exit(EXIT_FAILURE);
			}
		}

		cd->boundary_rule = (double*)malloc(4*sizeof(double));
		for(i=0;i<4;i++)
			cd->boundary_rule[i] = atof(strtok(NULL," \t"));
		cd->n_tag = atoi(strtok(NULL," \t"));
		cd->tag = (int*)malloc(cd->n_tag*sizeof(int));
		for(i=0;i < cd->n_tag;i++)
			cd->tag[i] = atoi(strtok(NULL," \t"));

		fgets(line,MAXBUFF-2,INI);
		sscanf(line,"%s" ,line1);
	}
}
//=============================================================================
//		NEW CONDITION
//=============================================================================
struct Condition *new_BC(struct Condition **cond,struct Condition *prev_bc,int type)
{
	struct Condition *cd;

	cd =(struct Condition*)malloc(sizeof(struct Condition));
	cd->condition_type = type;
	if(*cond == NULL){
		*cond = cd;
		cd->prev = NULL;
	}
	else{
		prev_bc->next = cd;
		cd->prev = prev_bc;
		cd->next = NULL;
	}
	return cd;
}
