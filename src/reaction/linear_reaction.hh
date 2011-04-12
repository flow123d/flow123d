//enum REACTION_TYPE { decay = 1, kinetics};
#include<vector>

class Linear_reaction
{
	public:
		Linear_reaction(int n_subst, char* section, double time_step);
		~Linear_reaction();
		double **Compute_reaction(double **concentrations, int n_subst, int loc_el); //multiplication of concentrations array by reaction matrix
	private:
		Linear_reaction();
		int *Set_indeces(char *section);
		int Set_nr_of_isotopes(char *section);
		int Set_nr_of_decays(char *section);
		double *Set_half_lives(char *section);
		void Prepare_decaying_isotopes_ids(int n_subst);
		void Modify_decaying_isotopes_ids(void);
		void Set_bifurcation(char *section, int dec_nr);
		void Set_bifurcation_on(char *section);
		double **Allocate_reaction_matrix(int n_subst); //reaction matrix initialization
		double **Modify_reaction_matrix(int n_subst, double time_step); //prepare the matrix, which describes reactions
		double **Modify_reaction_matrix(int n_subst, double time_step, int bifurcation);
		double **Modify_reaction_matrix_repeatedly(int n_subst, char *section, double time_step);
		void Print_reaction_matrix(int n_subst);
		int Waste_reaction_matrix(double **matrix, int n_subst); //is here to dealocate memory and to prevent memory leaks
		int *Get_indeces();
		int Get_nr_of_isotopes();
		double *Get_half_lives();
		double **reaction_matrix;
		double *half_lives;
		int *substance_ids;
		int nr_of_isotopes;
		int nr_of_decays;
		int *decaying_isotopes;
		std::vector<std::vector<double> > bifurcation;
		bool bifurcation_on;
};
