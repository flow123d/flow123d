#include <stdlib.h>
#include "struct.h"


//=============================================================================
//		SIDE NODES COUNT
//=============================================================================
int side_nodes_count(struct Side *sde)
{
	switch(sde->element->type){
		case 1:
			return 1;
		case 2:
			return 2;
		case 4:
			return 3;
		default:
			return -1;
	}
}
//=============================================================================
//		ELEMENT NODES COUNT
//=============================================================================
int element_nodes_count(struct Element *elm)
{
	switch(elm->type){
		case 1:
			return 2;
		case 2:
			return 3;
		case 4:
			return 4;
		default:
			return -1;
	}
}
//=============================================================================
//		ELEMENT DIMENSION
//=============================================================================
int element_dim(struct Element *elm)
{
	switch(elm->type){
		case 1:
			return 1;
		case 2:
			return 2;
		case 4:
			return 3;
		default:
			return -1;
	}
}
//=============================================================================
//		ELEMENT SIDES COUNT
//=============================================================================
int element_sides_count(struct Element *elm)
{
	switch(elm->type){
		case 1:
			return 2;
		case 2:
			return 3;
		case 4:
			return 4;
		default:
			return -1;
	}
}
//=============================================================================
//		ASSIGN NODES TO SIDES
//=============================================================================
void assign_nodes_to_sides(struct Problem *problem)
{
	struct Element *elm;
	int i;

	FOR_ELEMENTS(elm)
		FOR_ELEMENT_SIDES(elm,i)
			switch(elm->type){
				case 1:
					elm->side[i]->node[0] = elm->node[i];
					break;
				case 2:
					switch(i){
						case 0:
							elm->side[i]->node[0] = elm->node[0];
							elm->side[i]->node[1] = elm->node[1];
							break;
						case 1:
							elm->side[i]->node[0] = elm->node[1];
							elm->side[i]->node[1] = elm->node[2];
							break;
						case 2:
							elm->side[i]->node[0] = elm->node[2];
							elm->side[i]->node[1] = elm->node[0];
							break;
					}
					break;
				case 4:
					switch(i){
						case 0:
							elm->side[i]->node[0] = elm->node[1];
							elm->side[i]->node[1] = elm->node[2];
							elm->side[i]->node[2] = elm->node[3];
							break;
						case 1:
							elm->side[i]->node[0] = elm->node[0];
							elm->side[i]->node[1] = elm->node[2];
							elm->side[i]->node[2] = elm->node[3];
							break;
						case 2:
							elm->side[i]->node[0] = elm->node[0];
							elm->side[i]->node[1] = elm->node[1];
							elm->side[i]->node[2] = elm->node[3];
							break;
						case 3:
							elm->side[i]->node[0] = elm->node[0];
							elm->side[i]->node[1] = elm->node[1];
							elm->side[i]->node[2] = elm->node[2];
							break;
					}
					break;
				default:
					break;
			}
}
//=============================================================================
//		CALCULATE SIDES CENTER
//=============================================================================
void calc_sides_center(struct Problem *problem)
{
	struct Element *elm;
	int i,j,k;

	FOR_ELEMENTS(elm)
		FOR_ELEMENT_SIDES(elm,i)
			FOR_SIDE_NODES(elm->side[i],j)
				for(k=0;k<3;k++)
					elm->side[i]->center[k] += elm->side[i]->node[j]->coor[k] / elm->side[i]->n_nodes;
}
//=============================================================================
//		SUPPORTED BC TYPE
//=============================================================================
int supp_bc_type(int bc)
{
	switch(bc){
		case 1:
		case 2:
		case 3:
			return 1;
		default:
			return 0;
	}
}
