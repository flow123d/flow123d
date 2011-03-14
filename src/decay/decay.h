enum REACTION_TYPE { decay = 1, kinetics};

class Linear_reaction
{
	public:
		Linear_reaction();
		Linear_reaction(REACTION_TYPE type, double *order, int n_subst);
		~Linear_reaction();
		double **Create_reaction_matrix(double *order, int n_subst); //prepare the matrix, which describes reactions
		double **Compute_reaction(double **concentrations); //multiplication of concentrations array by reaction matrix
		REACTION_TYPE Set_reaction_type(REACTION_TYPE type); //change of reaction type should be conditionated by generation of new reaction matrix
		REACTION_TYPE Get_reaction_type();
		int *Set_indeces(char *section);
		int *Get_indeces();
		int Set_nr_of_isotopes(char *section);
		int Get_nr_of_isotopes();
		double *Set_half_lives(char *section);
		double *Get_half_lives();
	private:
		REACTION_TYPE react_type;
		double **reaction_matrix;
		double *half_lives;
		int *substance_ids;
		int nr_of_isotopes;
};
