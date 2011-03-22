enum REACTION_TYPE { decay = 1, kinetics};

class Linear_reaction
{
	public:
		Linear_reaction();
		Linear_reaction(REACTION_TYPE type, int n_subst, char* section);
		//Linear_reaction(REACTION_TYPE type, int n_subst, double **reaction_matrix);
		~Linear_reaction();
		//double **Modify_reaction_matrix(double **reaction_matrix, int n_subst, double time_step); //prepare the matrix, which describes reactions
		double **Modify_reaction_matrix(int n_subst, double time_step); //prepare the matrix, which describes reactions
		double **Prepare_reaction_matrix(int n_subst); //reaction matrix initialization
		int Waste_reaction_matrix(double **matrix, int n_subst); //is here to dealocate memory and to prevent memory leaks
		int *Get_indeces();
		int Get_nr_of_isotopes();
		double *Get_half_lives();
		double **Compute_reaction(double **concentrations, int n_subst, int loc_el); //multiplication of concentrations array by reaction matrix
		REACTION_TYPE Get_reaction_type();
		int *Set_indeces(char *section);
		int Set_nr_of_isotopes(char *section);
		double *Set_half_lives(char *section);
	private:
		REACTION_TYPE Set_reaction_type(REACTION_TYPE type); //change of reaction type should be conditionated by generation of new reaction matrix
		REACTION_TYPE react_type;
		double **reaction_matrix;
		double *half_lives;
		int *substance_ids;
		int nr_of_isotopes;
};
