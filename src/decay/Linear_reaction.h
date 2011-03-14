enum REACTION_TYPE { decay = 1, kinetics};

class Linear_reaction
{
	public:
		Linear_reaction();
		Linear_reaction(REACTION_TYPE type, double *order, int n_subst, char* section);
		~Linear_reaction();
		double **Create_reaction_matrix(double *half_lives, int n_subst, int *substance_ids); //prepare the matrix, which describes reactions
		int Waste_reaction_matrix(double **matrix, int n_subst); //is here to dealocate memory and to prevent memory leaks
		int *Get_indeces();
		int Get_nr_of_isotopes();
		double *Get_half_lives();
		REACTION_TYPE Get_reaction_type();
	private:
		int *Set_indeces(char *section);
		int Set_nr_of_isotopes(char *section);
		double *Set_half_lives(char *section);
		double **Compute_reaction(double **concentrations); //multiplication of concentrations array by reaction matrix
		REACTION_TYPE Set_reaction_type(REACTION_TYPE type); //change of reaction type should be conditionated by generation of new reaction matrix
		REACTION_TYPE react_type;
		double **reaction_matrix;
		double *half_lives;
		int *substance_ids;
		int nr_of_isotopes;
};
