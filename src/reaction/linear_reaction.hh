/** @brief class Linear_reaction is used to enable simulation of simple chemical reactions
 *
 * Class in this file makes it possible to realize  simulation of reaction of the first order by simple matrix multiplication.
 * One step of the linear reaction is represented as a product of a matrix containing concentrations of observed speciesin elements in rows multiplied by so called
 * reaction_matrix.
 * Matrix containing concentrations has a dimension Nxn, where N is a number of elements in mesh and n denotes a number of transported chemical species.
 * The reaction_matrix is a square matrix and it has a dimension nxn.
 *
 */
#include<vector> ///< included to enable saving bifurcation

class Linear_reaction
{
	public:
		Linear_reaction(int n_subst, double time_step); ///< constructor with parameters for initialization of  a new declared class member
		~Linear_reaction(); ///< desctructor
		double **Compute_reaction(double **concentrations, int n_subst, int loc_el); ///< method for multiplication of concentrations array by reaction matrix
	private:
		Linear_reaction(); ///< suppresses a use of constructor without parameters
		int *Set_indeces(char *section); ///< reads a sequence of numbers defining an order of substances in decay chain, from ini-file
		int Set_nr_of_isotopes(char *section); ///< reads an information about number of decay chain members, from ini-file
		int Set_nr_of_decays(void); ///< reads an information about number of decays under consideration, from ini-file
		double *Set_half_lives(char *section); ///< reads a sequence of halflives belonging to separate decay chain steps, from ini-file
		//void Prepare_decaying_isotopes_ids(int n_subst); ///< probably superfluous
		//void Modify_decaying_isotopes_ids(void); ///< probably superfluous
		void Set_bifurcation(char *section, int dec_nr); ///< reads an information for construction of a matrix describing bifurcation of every single decay chain on one row of the matrix, from ini-file
		void Set_bifurcation_on(char *section); ///< reads an information if the bifurcation for a current decay chain is switched on, from ini-file, initialy it is switched off,
		double **Allocate_reaction_matrix(int n_subst); ///< allocates memory for nxn square reaction matrix
		double **Modify_reaction_matrix(int n_subst, double time_step); ///< it is used for reaction matrix modification in cases when a bifurcation for a current decay chain is switched off
		double **Modify_reaction_matrix(int n_subst, double time_step, int bifurcation); ///< it is used for reaction matrix modification in cases when a bifurcation for a current decay chain is switched on
		double **Modify_reaction_matrix_repeatedly(int n_subst, double time_step); ///< calls Modify_reaction_matrix(..) so many times as many nr_of_decays are defined
		void Print_reaction_matrix(int n_subst); ///< it is just control of reaction matrix modifications
		int Waste_reaction_matrix(double **matrix, int n_subst); ///< dealocates reaction_matrix memory to prevent memory leaks
		int *Get_indeces(); ///< prints a sequence of indeces of chemical species contained in decay chain and returns a pointer to an array
		int Get_nr_of_isotopes(); ///< returns number of isotopes without printing
		double *Get_half_lives(); ///< prints a sequence of half-lives belonging to every single step of a current decay chain and returns a pointer to an array
		double **reaction_matrix; ///<two-dimensional array describing together all radioactive decay chains under consideration
		double *half_lives; ///< one-dimensional array containing half-lives belonging to every single step of a current decay chain
		int *substance_ids; ///< one-dimensional array with indeces of transported substances which are taking patr in a current decay chain
		int nr_of_isotopes; ///< number of isotopes in a current decay chain
		int nr_of_decays; ///< how many decay chains are under consideration, it means number of [Decay_i] sections in ini-file
		//int *decaying_isotopes; ///< probably superfluous
		std::vector<std::vector<double> > bifurcation; ///< two dimensional array containing mass percentage of every single decay bifurcation on every single row
		bool bifurcation_on; ///< bifurcation is initialy switched off
};
