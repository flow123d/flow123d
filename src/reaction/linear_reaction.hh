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
		double **compute_reaction(double **concentrations, int n_subst, int loc_el); ///< method for multiplication of concentrations array by reaction matrix
	private:
		Linear_reaction(); ///< suppresses a use of constructor without parameters
		int *set_indeces(char *section, int n_subst); ///< reads a sequence of numbers defining an order of substances in decay chain, from ini-file
		void set_nr_of_isotopes(char *section); ///< reads an information about number of decay chain members, from ini-file
		void set_nr_of_decays(void); ///< reads an information about number of decays under consideration, from ini-file
		void set_nr_of_FoR(void); ///< reads an information about a number of first order reactions under consideration, from ini-file
		double *set_half_lives(char *section, int n_subst); ///< reads a sequence of halflives belonging to separate decay chain steps, from ini-file
		//void Prepare_decaying_isotopes_ids(int n_subst); ///< probably superfluous
		//void Modify_decaying_isotopes_ids(void); ///< probably superfluous
		void set_bifurcation(char *section, int dec_nr); ///< reads an information for construction of a matrix describing bifurcation of every single decay chain on one row of the matrix, from ini-file
		void set_bifurcation_on(char *section); ///< reads an information if the bifurcation for a current decay chain is switched on, from ini-file, initialy it is switched off,
		void set_For_on(void); ///< reads an information if first order reactions simulation is switched on, duplicit to the function in problem.cc
		void set_decay_on(void); ///< reads an information if a decay simulation is switched on, duplicit to the function in problem.cc
		void set_kinetic_constants(char *section); ///< reads an information for construction of a vector containing kinetic constants of every single first order reactions, from those constants half-lives belonging to first order reactions are computed, from ini-file
		double **allocate_reaction_matrix(int n_subst); ///< allocates memory for nxn square reaction matrix
		double **modify_reaction_matrix(int n_subst, double time_step); ///< it is used for reaction matrix modification in cases when a bifurcation for a current decay chain is switched off
		double **modify_reaction_matrix(int n_subst, double time_step, int bifurcation); ///< it is used for reaction matrix modification in cases when a bifurcation for a current decay chain is switched on
		double **modify_reaction_matrix_repeatedly(int n_subst, double time_step); ///< calls Modify_reaction_matrix(..) so many times as many nr_of_decays are defined
		void print_reaction_matrix(int n_subst); ///< it is just control of reaction matrix modifications
		void print_indeces(int n_subst); ///< prints a sequence of indeces of chemical species contained in decay chain and returns a pointer to an array
		void print_half_lives(int n_subst); ///< prints a sequence of half-lives belonging to every single step of a current decay chain and returns a pointer to an array
		double **reaction_matrix; ///<two-dimensional array describing together all radioactive decay chains under consideration
		double *half_lives; ///< one-dimensional array containing half-lives belonging to every single step of a current decay chain
		int *substance_ids; ///< one-dimensional array with indeces of transported substances which are taking patr in a current decay chain or in a first order reaction
		int nr_of_isotopes; ///< number of isotopes in a current decay chain
		int nr_of_decays; ///< how many decay chains are under consideration, it means number of [Decay_i] sections in ini-file
		int nr_of_FoR; ///< how many firts order reactions of the type A -> B are under consideration, it is a number of [FoReaction_i]
		//int *decaying_isotopes; ///< probably superfluous
		std::vector<std::vector<double> > bifurcation; ///< two dimensional array containing mass percentage of every single decay bifurcation on every single row
		std::vector<double> kinetic_constant; ///< one dimensional array of kinetic constants belonging to considered reactions
		bool bifurcation_on; ///< bifurcation is initialy switched off
		bool FoR_on;
		bool decay_on;
};
