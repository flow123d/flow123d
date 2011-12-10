#ifndef BCD_H_
#define BCD_H_

struct Problem;
struct Element;
struct Side;
struct Node;
struct Material;
struct Condition;

// --- bcd.c ---
void read_ini(struct Problem *problem);


// --- output.c ---
void print_FBC_file(struct Problem *problem);
void print_TBC_file(struct Problem *problem);
void print_TIC_file(struct Problem *problem);
void print_TSO_file(struct Problem *problem);
void print_MTR_file(struct Problem *problem);

// --- system.c ---
void skip_to( FILE *in, char *sekce);
FILE *open_file(char *file,char *open_param);
char *xstrcpy( char *src );
int get_bool(char *file,char *section,char *key,double def_value);
char* get_string(char *file,char *section,char *key,char* def_value);
int get_int(char *file,char *section,char *key,int def_value);
double get_float(char *file,char *section,char *key,double def_value);
void get_cond(struct Problem *problem,char *secstart,char *secstop,int type);
struct Condition *new_BC(struct Condition **cond,struct Condition *prev_bc,int type);

// --- mesh.c ---
void read_mesh(struct Problem *problem);
void create_element_sides(struct Problem *problem);
void create_node_hash(struct Problem *problem);
void create_element_hash(struct Problem *problem);
void create_material_hash(struct Problem *problem);
void create_material_list(struct Problem *problem);
void read_and_extrude_neigbours(struct Problem *problem);
void create_boundary_sides_list(struct Problem *problem);
void init_problem(struct Problem *problem);
void calc_element_center(struct Problem *problem);
struct Element* new_element(struct Element **ele,struct Element *p_elm);
struct Node* new_node(struct Node **node,struct Node *p_node);

// --- topology_data.c ---
int element_nodes_count(struct Element *elm);
int element_sides_count(struct Element *elm);
int element_dim(struct Element *elm);
int side_nodes_count(struct Side *sde);
void assign_nodes_to_sides(struct Problem *problem);
void calc_sides_center(struct Problem *problem);
int supp_bc_type(int bc);

// --- condition.c ---
void extrude_BC(struct Problem *problem);
void calculate_BC(struct Problem *problem,struct Condition *cd);
void calculate_TIC(struct Problem *problem,struct Condition *cd);
int equation_element(struct Problem *problem,struct Element *elm,double *equation);
int equation_side(struct Problem *problem,struct Side *sde,double *equation);
double point_value(double *point,double *equation);
void assign_values(struct Side *sde,struct Element *elm,struct Condition *cd);
void init_cond_data(struct Problem *problem);

#endif /*BCD_H_*/
