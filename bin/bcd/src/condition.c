#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "struct.h"
#include "bcd.h"

//=============================================================================
//		EXTRUDE HOMOGENEOUS NEUMANN BOUNDARY CONDITION'S
//=============================================================================
void extrude_BC(struct Problem *problem)
{
	struct Side *sde;

	problem->n_wboundaries = problem->n_boundaries;

	if(problem->write_all_BC == NO)
		FOR_BOUNDARY_SIDES(sde)
			if((sde->fbc_type == 2) && (sde->fbc[0] == 0.0)){
				sde->write = NO;
				problem->n_wboundaries--;
			}
}
//=============================================================================
//		CALCULATE BOUNDARY CONDITION'S
//=============================================================================
void calculate_BC(struct Problem *problem,struct Condition *cd)
{
	struct Side *sde;
	int i;

	FOR_BOUNDARY_SIDES(sde)
		switch(cd->seting_type){
			case 0:
				assign_values(sde,NULL,cd);
				break;
			case 1:
				for(i=0;i<cd->n_elm;i++)
					if(sde->element->id == cd->elm_list[i])
						assign_values(sde,NULL,cd);
				break;
			case 2:
				if(equation_side(problem,sde,cd->boundary_equation) == 1)
					assign_values(sde,NULL,cd);
				break;
			case 3:
				if(sde->element->mid == cd->mid)
					assign_values(sde,NULL,cd);
				break;
			case 4:
				if((sde->element->mid == cd->mid) && (equation_side(problem,sde,cd->boundary_equation)))
					assign_values(sde,NULL,cd);
				break;
			default:
				break;
		}
}
//=============================================================================
//		CALCULATE TRANSPORT INITIAL CONDITION'S
//=============================================================================
void calculate_TIC(struct Problem *problem,struct Condition *cd)
{

	struct Element *elm;
	int i;

	FOR_ELEMENTS(elm)
		switch(cd->seting_type){
			case 0:
				assign_values(NULL,elm,cd);
				break;
			case 1:
				for(i=0;i<cd->n_elm;i++)
					if(elm->id == cd->elm_list[i])
						assign_values(NULL,elm,cd);
				break;
			case 2:
				if(equation_element(problem,elm,cd->boundary_equation) == 1)
					assign_values(NULL,elm,cd);
				break;
			case 3:
				if(elm->mid == cd->mid)
					assign_values(NULL,elm,cd);
				break;
			case 4:
				if((elm->mid == cd->mid) && (equation_element(problem,elm,cd->boundary_equation) == 1 ))
					assign_values(NULL,elm,cd);
				break;
			default:
				break;
		}
}
//=============================================================================
//		ASSIGN VALUES
//=============================================================================
void assign_values(struct Side *sde,struct Element *elm,struct Condition *cd)
{
	int i;

	switch(cd->condition_type){
		case 1: 	// FBC
			sde->fbc_type = cd->bc_type;
			sde->fbc[0] = point_value(sde->center,cd->boundary_rule);
			if(sde->fbc_type == 3)
				sde->fbc[1] = cd->bc_param;
			sde->fbc_tag = cd->tag[0];
			break;
		case 2:		// TBC
			for(i=0;i<cd->n_tag;i++)
				sde->tbc[cd->tag[i]] = point_value(sde->center,cd->boundary_rule);
			break;
		case 3:		//TIC
			for(i=0;i<cd->n_tag;i++)
				elm->tic[cd->tag[i]] = point_value(elm->center,cd->boundary_rule);
			break;
	}

}
//=============================================================================
//		IDENTIFY EQUATION SIDE
//=============================================================================
int equation_side(struct Problem *problem,struct Side *sde,double *equation)
{
	int i;
	int bool = 1;

	FOR_SIDE_NODES(sde,i)
		bool *= (int) (fabs(point_value(sde->node[i]->coor,equation)) < problem->accuracy);
	return bool;
}
//=============================================================================
//		IDENTIFY HALF PLANE ELEMENT
//=============================================================================
int equation_element(struct Problem *problem,struct Element *elm,double *equation)
{
	int bool = 1;

	bool *= (int) ((point_value(elm->center,equation)) >= 0);


	return bool;
}
//=============================================================================
//		CALCULATE SIDE VALUE
//=============================================================================
double point_value(double *point,double *equation)
{
	double value;

		value = equation[0] * point[0] +
				equation[1] * point[1] +
				equation[2] * point[2] +
				equation[3];

	return value;
}
//=============================================================================
//		INITIALIZE CONDITION DATA
//=============================================================================
void init_cond_data(struct Problem *problem)
{
	struct Side *sde;
	struct Element *elm;
	int i;

	FOR_BOUNDARY_SIDES(sde){

		if(problem->fbc_file != NULL){
			sde->fbc_type = 2;
			sde->fbc[0] = 0.0;
			sde->fbc[1] = 0.0;
			sde->fbc_tag = 0;
		}

		if(problem->tbc_file != NULL) {
			sde->tbc = (double*)malloc(problem->n_subst * sizeof(double));
			for(i=0;i < problem->n_subst;i++)
				sde->tbc[i] = 0;
		}
	}

	if(problem->tic_file != NULL)
		FOR_ELEMENTS(elm){
			elm->tic = (double*)malloc(problem->n_subst*sizeof(double));
			for(i=0;i < problem->n_subst;i++)
				elm->tic[i] = 0;
	}
}
