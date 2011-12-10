#ifndef STRUCT_H_
#define STRUCT_H_



#define YES   1
#define NO    0
#define OK   1
#define NOK   0
//#define FBC	 1
//#define TBC	 2
//#define TIC	 3
#define MAXBUFF    10240
#define FOR_ELEMENTS(i)  for((i)=problem->element;(i)!=NULL;(i)=(i)->next)
#define FOR_FBC(i)  for((i)=problem->fbc;(i)!=NULL;(i)=(i)->next)
#define FOR_TBC(i)  for((i)=problem->tbc;(i)!=NULL;(i)=(i)->next)
#define FOR_TIC(i)  for((i)=problem->tic;(i)!=NULL;(i)=(i)->next)
#define FOR_BOUNDARY_SIDES(i)  for((i)=problem->boundary_side;(i)!=NULL;(i)=(i)->next)
#define FOR_NODES(i)  for((i)=problem->node;(i)!=NULL;(i)=(i)->next)
#define FOR_MATERIALS(i)  for((i)=problem->material;(i)!=NULL;(i)=(i)->next)
#define FOR_ELEMENT_SIDES(i,j)  for((j)=0;(j)<(i)->n_sides;(j)++)
#define FOR_ELEMENT_NODES(i,j)  for((j)=0;(j)<(i)->n_nodes;(j)++)
#define FOR_SIDE_NODES(i,j)  for((j)=0;(j)<(i)->n_nodes;(j)++)


struct Problem
{
//	[Internal data]
	int      max_elm_id;
    int      n_elm;
    int      n_sides; // sides count
    int      n_boundaries; // boundaries count
    int		 n_wboundaries; // boundaries count - boundary without homogeneous Neumann's BC
	int      max_node_id;
	int		 max_mat_id;
	int		 n_mat;
	int      n_nodes;
	char	 *dbl_fmt;
    struct Side 	*boundary_side;
	struct Element 	*element;
	struct Element	**element_hash;
	struct Node 	*node;
    struct Node	 	**node_hash;
    struct Material *material;
    struct Material **material_hash;
    struct Condition *fbc;
    struct Condition *tbc;
    struct Condition *tic;
//	[External data]
	int      n_subst;
	int		 output_digits;
	int		 write_all_BC;
	double   accuracy;
    char     *mesh_file;
    char     *ngh_file;
    char     *fbc_file;
    char     *tbc_file;
    char     *tic_file;
    char     *mtr_file;
    char     *tso_file;
    char     *ini_file;
};


struct Element
{
	int      id;
	int      type;
	int     *nid;
	int      mid;
	int 	 n_nodes;
	int 	 n_sides;
	double   center[3];
	struct Element *prev;
	struct Element *next;
	struct Node	   **node;
	struct Side	   **side;
	double 	*tic;
};


struct Node
{
	int      id;
    double  coor[3];
    struct Node *next;
    struct Node *prev;
};

struct Side
{
	int id;
	double	center[3];
    int     n_nodes;
    int 	boundary;
    int		write;
    struct Element *element;
    struct Node **node;
    struct Side *next;
    struct Side *prev;
    int fbc_type;
    int fbc_tag;
    double fbc[2];
    double *tbc;
};
struct Condition
{
	int id;
	int condition_type; // 1 - FBC  2 - TBC  //////3 - TIC
	int seting_type; // 1 - ELM, 2 - EQN, 3 - MAT, 4 - MAT + EQN
	double	*boundary_equation;
    int  	mid;
    int 	*elm_list;
    int		bc_type; 	//only for flow
    double  bc_param;
    int		n_elm;
    int     n_tag;
    int		*tag;
    double	*boundary_rule;
    struct Condition *next;
    struct Condition *prev;

};
struct Material
{
	int id;
	int dim;
    struct Material *next;
    struct Material *prev;
};

#endif /*STRUCT_H_*/
