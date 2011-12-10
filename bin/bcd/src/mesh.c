#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "struct.h"
#include "bcd.h"

//=============================================================================
//		CALCULATE ELEMENT CENTER
//=============================================================================
void calc_element_center(struct Problem *problem)
{
	struct Element *elm;
	int i,j;

	FOR_ELEMENTS(elm)
		FOR_ELEMENT_NODES(elm,i)
			for(j=0;j<3;j++)
				elm->center[j] += elm->node[i]->coor[j] / elm->n_nodes;
}
//=============================================================================
//		CREATE BOUNDARY SIDES LIST
//=============================================================================
void create_boundary_sides_list(struct Problem *problem)
{
	struct Element *elm;
	struct Side *prev_side;
	int i,bool;

		bool = 1;
		prev_side = NULL;
		FOR_ELEMENTS(elm)
			FOR_ELEMENT_SIDES(elm,i)
				if(elm->side[i]->boundary == YES){
					if(bool==1){
						bool = 0;
						problem->boundary_side = elm->side[i];
					}
					else prev_side->next = elm->side[i];
					elm->side[i]->prev = prev_side;
					elm->side[i]->next = NULL;
					prev_side = elm->side[i];
				}
	/*	FOR_BOUNDARY_SIDES(prev_side){
			printf("\n%d\t%d",prev_side->element->id,prev_side->id);
			getchar();
		} */
}
//=============================================================================
//		CREATE ELEMENT SIDES
//=============================================================================
void create_element_sides(struct Problem *problem)
{
	struct Element *elm;
	int i,j;

	elm = problem->element;

	FOR_ELEMENTS(elm)
		{
		for(i=0;i<3;i++)
			elm->center[i] = 0.0;
		elm->n_sides = element_sides_count(elm);
		elm->side = (struct Side**)malloc(elm->n_sides * sizeof(struct Side*));
		for(i=0;i<elm->n_sides;i++)
			{
			problem->n_sides++;
			elm->side[i] = (struct Side*)malloc(sizeof(struct Side));
			elm->side[i]->id = i;
			elm->side[i]->boundary = YES;
			elm->side[i]->element = elm;
			elm->side[i]->write = YES;
			elm->side[i]->n_nodes = side_nodes_count(elm->side[i]);
			for(j=0;j<3;j++)
				elm->side[i]->center[j] = 0;
			elm->side[i]->node = (struct Node**)malloc(elm->side[i]->n_nodes * sizeof(struct Node*));
			for(j=0;j<elm->side[i]->n_nodes;j++)
				elm->side[i]->node[j] = (struct Node*)malloc(sizeof(struct Node));
			}
		}
	assign_nodes_to_sides(problem);
	calc_sides_center(problem);
}
//=============================================================================
//		READ MESH DATA
//=============================================================================
void read_mesh(struct Problem *problem)
{
     char line[MAXBUFF];
     int i,j,n;
     struct Node *node = NULL;
     struct Element *elm = NULL;
     FILE     *MSH;

     init_problem(problem);
     MSH = open_file(problem->mesh_file,"rt");
     skip_to(MSH,"$Nodes");
     fgets( line, MAXBUFF - 2, MSH );
     problem->n_nodes = atoi( strtok( line, " \t" ) );  // number of nodes
     for(i=0;i<problem->n_nodes;i++){    // Read all nodes
    	 node = new_node(&problem->node,node);
    	 fgets( line, MAXBUFF - 2, MSH );
    	 node->id = atoi( strtok( line, " \t" ));
    	 if (problem->max_node_id < node->id) problem->max_node_id = node->id;	// max_node_id
    	 for(j = 0; j < 3; j++)
    		 node->coor[j] = atof( strtok( NULL, " \t" ));
     }
     create_node_hash(problem);
   /*
   FOR_NODES(node){
          printf("%d\t%f\t%f\t%f\n",node->id, node->coor[0], node->coor[1], node->coor[2]);
          getchar();}  */
    /*
     for(i=0;i < problem->max_node_id + 1;i++)
         	 if(problem->node_hash[i] != NULL)
         	 	printf("%d\t%f\t%f\t%f\n",i, problem->node_hash[i]->coor[0], problem->node_hash[i]->coor[1], problem->node_hash[i]->coor[2]);
         	 else continue;
     getchar();
     */


     skip_to(MSH,"$Elements");
     fgets( line, MAXBUFF - 2, MSH );
     problem->n_elm = atoi( strtok( line, " \t" ) );  // number of elements
     for(i = 0; i < problem->n_elm; i++){    // Read all elements
    	 elm = new_element(&problem->element,elm);
         fgets( line, MAXBUFF - 2, MSH );
         elm->id = atoi( strtok( line, " \t" ));
         if (problem->max_elm_id < elm->id) problem->max_elm_id = elm->id;	    // max_elm_id
         elm->type = atoi( strtok( NULL, " \t" ));
         n = atoi( strtok( NULL, " \t" ));
         elm->mid = atoi( strtok( NULL, " \t" ));
         if(elm->mid > problem->max_mat_id)
        	 problem->max_mat_id = elm->mid;
         for(j = 0; j < n-1; j++)
        	 atoi( strtok( NULL, " \t" ));  // unused
         elm->n_nodes = element_nodes_count(elm);
         elm->node = (struct Node**)malloc(elm->n_nodes * sizeof(struct Node*));
         for(j = 0; j < elm->n_nodes; j++)
         		elm->node[j] = problem->node_hash[atoi( strtok( NULL, " \t" ))];
     }
     fclose(MSH);
     create_element_hash(problem);
     create_material_hash(problem);
     create_material_list(problem);
    /*
     FOR_ELEMENTS(elm)

                   	 {
                   	 printf("%d\t%d\t%d\t%d\n",elm->id, elm->n_nodes, elm->mid, elm->n_nodes);
                   	 for(i=0;i < elm->n_nodes;i++)
                   	 printf("\t\t\t%d\t%f\t%f\t%f\n",elm->node[i]->id, elm->node[i]->coor[0], elm->node[i]->coor[1], elm->node[i]->coor[2]);
                   	 getchar();
                   	 }

                    getchar(); */
   /*  for(i=0;i < problem->max_elm_id + 1;i++)
              	 if(problem->element_hash[i] != NULL)
              	 	printf("%d\t%d\t%d\t\n",i, problem->element_hash[i]->id,problem->element_hash[i]->n_nodes);
              	 else continue;
          getchar(); */

}
//=============================================================================
//		CREATE NODE HASH
//=============================================================================
void create_node_hash(struct Problem *problem)
{
	struct Node *node;
	int i;
	problem->node_hash = (struct Node**)malloc((problem->max_node_id + 1) * sizeof(struct Node*)); // create node hash
	for(i=0;i < problem->max_node_id + 1;i++)
		problem->node_hash[i] = NULL;
	FOR_NODES(node)
		problem->node_hash[node->id] = node;

}
//=============================================================================
//		CREATE NODE HASH
//=============================================================================
void create_element_hash(struct Problem *problem)
{
	struct Element *elm;
	int i;
	problem->element_hash = (struct Element**)malloc((problem->max_elm_id + 1) * sizeof(struct Element*));
	for(i=0;i < problem->max_elm_id + 1;i++)
		problem->element_hash[i] = NULL;
	FOR_ELEMENTS(elm)
		problem->element_hash[elm->id]= elm;

}
//=============================================================================
//		CREATE MATERIAL HASH
//=============================================================================
void create_material_hash(struct Problem *problem)
{
	struct Material *mat;
	struct Element *elm;
	int i;
	problem->material_hash = (struct Material**)malloc((problem->max_mat_id + 1) * sizeof(struct Material*));
	for(i=0;i < problem->max_mat_id + 1;i++)
		problem->material_hash[i] = NULL;
	FOR_ELEMENTS(elm)
	{
		if(problem->material_hash[elm->mid] == NULL){
			mat = (struct Material*)malloc(sizeof(struct Material));
			mat->id = elm->mid;
			mat->dim = element_dim(elm);
			problem->material_hash[elm->mid] = mat;
			problem->n_mat++;
		}
	}
}
//=============================================================================
//		CREATE MATERIAL LIST
//=============================================================================
void create_material_list(struct Problem *problem)
{
	struct Material *mat,*prev_mat;
	int i;
	int bool = 1;

	prev_mat = NULL;
	for(i=0;i < problem->max_mat_id + 1;i++){
		if(problem->material_hash[i] == NULL) continue;
		else{
			mat = problem->material_hash[i];
			if(bool == 1){
				bool = 0;
				problem->material = mat;
			}
			else prev_mat->next = mat;
		mat->next = NULL;
		mat->prev = prev_mat;
		prev_mat = mat;
		}
	}
}
//=============================================================================
//		READ AND EXTRUDE NEIGBOURS
//=============================================================================
void read_and_extrude_neigbours(struct Problem *problem)
{
	char line[MAXBUFF];
	int i,j,n,n_bound,id,tag1,tag2,tag3,type;
	FILE     *NGH;

	n = 0;

	problem->n_boundaries = problem->n_sides;
	NGH = open_file(problem->ngh_file,"rt");
	skip_to(NGH,"$Neighbours");
	fgets( line, MAXBUFF - 2, NGH );
	n_bound = atoi( strtok( line, " \t" ) );
	for(i=0;i<n_bound;i++){
		fgets( line, MAXBUFF - 2, NGH );
		id = atoi( strtok( line, " \t" ) );
		type = atoi( strtok( NULL, " \t" ) );
		switch(type){
			case 11:	// [side to side]
				n = atoi( strtok( NULL, " \t" ) );
				problem->n_boundaries -= n;
				for(j=0;j<n;j++){
					tag1 = atoi( strtok( NULL, " \t" ) );
					tag2 = atoi( strtok( NULL, " \t" ) );
					problem->element_hash[tag1]->side[tag2]->boundary = NO;
				}
				break;
			case 20: //compatible  [element to side]
				problem->n_boundaries--;
				tag3 = atoi( strtok( NULL, " \t" ) );	// lower dim. elm.
				tag1 = atoi( strtok( NULL, " \t" ) );
				tag2 = atoi( strtok( NULL, " \t" ) );
				problem->element_hash[tag1]->side[tag2]->boundary = NO;
				break;
			default:
				break;
		}
	}
	fclose(NGH);
}
//=============================================================================
//		INITIALIZE PROBLEM
//=============================================================================
void init_problem(struct Problem *problem)
{
	char buff[MAXBUFF];

	problem->max_elm_id = 0;
	problem->n_elm = 0;
	problem->n_sides = 0; // sides count
	problem->n_boundaries = 0; // boundaries count
	problem->n_wboundaries = 0;
	problem->max_node_id = 0;
	problem->max_mat_id = 0;
	problem->n_nodes = 0;
	problem->n_mat = 0;

	sprintf(buff,"%s%d%s","% .",problem->output_digits,"f");
	problem->dbl_fmt = xstrcpy(buff);

	problem->boundary_side = NULL;
	problem->element = NULL;
	problem->element_hash = NULL;
	problem->node  = NULL;
	problem->node_hash = NULL;
    problem->material = NULL;
    problem->material_hash = NULL;
	problem->fbc = NULL;
	problem->tbc = NULL;
	problem->tic = NULL;
}
//=============================================================================
//		CREATE NEW ELEMENT
//=============================================================================
struct Element* new_element(struct Element **ele,struct Element *p_elm)
{
	struct Element *elm;

	elm = (struct Element*)malloc(sizeof(struct Element));
	if(*ele == NULL)
		*ele = elm;
	else
		p_elm->next = elm;
	elm->next = NULL;
	elm->prev = p_elm;
	return elm;
}
//=============================================================================
//		CREATE NEW NODE
//=============================================================================
struct Node* new_node(struct Node **node,struct Node *p_node)
{
	struct Node *n_node;

	n_node = (struct Node*)malloc(sizeof(struct Node));
	if(*node == NULL)
		*node = n_node;
	else
		p_node->next = n_node;
	n_node->next = NULL;
	n_node->prev = p_node;
	return n_node;
}
