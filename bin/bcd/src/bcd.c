/*
 ============================================================================
 Name        : BCD.c
 Author      : JK
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bcd.h"
#include "struct.h"


int main(int argc, char **argv)
{
	struct Problem *problem;
	struct Condition *cd;

	printf("Boundary condition generator, ver. 1.00 stable");

    problem = (struct Problem*)malloc(sizeof(struct Problem));

    if(argv[1] != NULL)
    	problem->ini_file = xstrcpy(argv[1]);
    else
    	return EXIT_FAILURE;

    printf("\n\nRead ini_file...");
    read_ini(problem);
    printf("OK");


    if(problem->mesh_file != NULL){
    	printf("\nRead msh_file...");
    	read_mesh(problem);
    	printf("OK");
    	create_element_sides(problem);
    	calc_element_center(problem);
    }

    if(problem->ngh_file != NULL){
    	printf("\nRead ngh_file...");
    	read_and_extrude_neigbours(problem);
    	printf("OK");
    	create_boundary_sides_list(problem);
    	init_cond_data(problem);
    }

	if((problem->fbc_file != NULL) && (problem->mesh_file != NULL) && (problem->ngh_file != NULL)){
		printf("\nRead FBC data...");
		get_cond(problem,"$FBC","$EndFBC",1);
		printf("OK");
		FOR_FBC(cd)
			calculate_BC(problem,cd);
		extrude_BC(problem);
		printf("\nWriting fbc_file...");
		print_FBC_file(problem);
		printf("OK");
	}

	if((problem->tbc_file != NULL) && (problem->mesh_file != NULL) && (problem->ngh_file != NULL) && (problem->fbc_file != NULL)){
		printf("\nRead TBC data...");
		get_cond(problem,"$TBC","$EndTBC",2);
		printf("OK");
		FOR_TBC(cd)
			calculate_BC(problem,cd);
		printf("\nWriting tbc_file...");
		print_TBC_file(problem);
		printf("OK");
	}

	if((problem->tic_file != NULL) && (problem->mesh_file != NULL) && (problem->ngh_file != NULL)){
		printf("\nRead TIC data...");
		get_cond(problem,"$TIC","$EndTIC",3);
		printf("OK");
		FOR_TIC(cd)
			calculate_TIC(problem,cd);
		printf("\nWriting tic_file...");
		print_TIC_file(problem);
		printf("OK");
	}

	if((problem->tso_file != NULL) && (problem->mesh_file != NULL) && (problem->ngh_file != NULL)){
			printf("\nRead TSO data...");
			//get_cond(problem,"$TSO","$EndTSO",3);
			printf("OK");
			//FOR_TIC(cd)
			//	calculate_TSO(problem,cd);
			printf("\nWriting tso_file...");
			print_TSO_file(problem);
			printf("OK");
		}


	if((problem->mtr_file != NULL) && (problem->mesh_file != NULL)){
		printf("\nWriting mtr_file...");
		print_MTR_file(problem);
		printf("OK");
	}



	printf("\n"); // \nPress any key to continue...
//	getchar();
    return EXIT_SUCCESS;
}
//=============================================================================
//		READ INI FILE
//=============================================================================
void read_ini(struct Problem *problem)
{
	problem->mesh_file = get_string(problem->ini_file,"$Input","Mesh_file",NULL);
	problem->ngh_file = get_string(problem->ini_file,"$Input","Ngh_file",NULL);
	problem->n_subst = get_int(problem->ini_file,"$Input","N_subst",1);

	problem->output_digits = get_int(problem->ini_file,"$Output","Output_digits",12);
	problem->write_all_BC = get_bool(problem->ini_file,"$Output","Write_all_BC",0);
	problem->accuracy = get_float(problem->ini_file,"$Output","Accuracy",1e-10);
	problem->fbc_file = get_string(problem->ini_file,"$Output","FBC_file",NULL);
	problem->tbc_file = get_string(problem->ini_file,"$Output","TBC_file",NULL);
	problem->tic_file = get_string(problem->ini_file,"$Output","TIC_file",NULL);
	problem->tso_file = get_string(problem->ini_file,"$Output","TSO_file",NULL);
	problem->mtr_file = get_string(problem->ini_file,"$Output","MTR_file",NULL);
	/*
	printf("\n%s",problem->mesh_file);
	printf("\n%s",problem->ngh_file);
	printf("\n%d",problem->n_subst);
	printf("\n%d",problem->output_digits);
	printf("\n%E",problem->accuracy);
	printf("\n%s",problem->fbc_file);
	printf("\n%s",problem->tbc_file);
	printf("\n%s",problem->tic_file);
	printf("\n%s",problem->mtr_file); */
}
