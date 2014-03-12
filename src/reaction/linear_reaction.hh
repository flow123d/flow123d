/** @brief class Linear_reaction is used to enable simulation of simple chemical reactions
 *
 * Class in this file makes it possible to realize  simulation of reaction of the first order by simple matrix multiplication.
 * One step of the linear reaction is represented as a product of a matrix containing concentrations of observed speciesin elements in rows multiplied by so called
 * reaction_matrix. Through this way radioactive decay can bee also realized and that was exactly what we did at the begining of journey. :-)
 * Matrix containing concentrations has a dimension Nxn, where N is a number of elements in mesh and n denotes a number of transported chemical species.
 * The reaction_matrix is a square matrix and it has a dimension nxn.
 *
 */
#ifndef LINREACT
#define LINREACT

#include <vector>
#include <input/input_type.hh>

class Mesh;
class Distribution;
class Reaction;

class Linear_reaction: public Reaction
{
	public:
		/*
	 	* Static variable for new input data types input
		*/
		static Input::Type::Record input_type;
		/*
	 	* Static variable gets information about particular decay step
		*/
		static Input::Type::Record input_type_one_decay_substep;
        /**
         *  Constructor with parameter for initialization of a new declared class member
         *  TODO: parameter description
         */
		//Linear_reaction(TimeMarks &marks, Mesh &init_mesh, MaterialDatabase &material_database, Input::Record in_rec, vector<string> &names);
		Linear_reaction(Mesh &init_mesh, Input::Record in_rec, vector<string> &names);
		/**
		*	Destructor.
		*/
		~Linear_reaction(void);

                /**
                *       Fuctions holds together setting of isotopes, bifurcations and substance indices.
                */
                virtual void init_from_input(Input::Record in_rec) override;
                
		/**
		*	For simulation of chemical reaction in just one element either inside of MOBILE or IMMOBILE pores.
		*/
		virtual double **compute_reaction(double **concentrations, int loc_el) override;
		/**
		*	Prepared to compute simple chemical reactions inside all of considered elements. It calls compute_reaction(...) for all the elements controled by concrete processor, when the computation is paralelized.
		*/
		void update_solution(void) override;
		/**
		*	This method modificates reaction matrix as described in ini-file a single section [Decay_i] or [FoReact_i]. It is used when bifurcation is switched off.
		*/
		virtual double **modify_reaction_matrix(void);
             
	protected:

        double **allocate_reaction_matrix(void);

		/**
		*	This method disables to use constructor without parameters.
		*/
		Linear_reaction();
		/**
		*	For control printing of a matrix describing simple chemical raections.
		*/
		void print_reaction_matrix(void);
		/**
		*	For printing nr_of_isotopes identifies of isotopes in a current decay chain.
		*/
		void print_indices(int dec_nr, int n_subst);
		/**
		* Following method releases reaction matrix to make it possible to set a new time step for chemistry.
		*/
		void release_reaction_matrix();
		/**
		*	For printing (nr_of_isotopes - 1) doubles containing half-lives belonging to particular isotopes on screen.
		*/
		void print_half_lives(int n_subst);
                
                /**
                *       Finds a position of a string in specified array.
                */
                unsigned int find_subst_name(const std::string &name);
		/**
		* 	Boolean which indicates the use of Pade approximant of the matrix exponential.
		*/
		//bool matrix_exp_on;
		/**
		*	Small (nr_of_species x nr_of_species) square matrix for realization of radioactive decay and first order reactions simulation.
		*/
		double **reaction_matrix;
                /**
                *       Pointer to reference previous concentration array used in compute_reaction().
                */
                double *prev_conc;
		/**
		*	Sequence of (nr_of_isotopes - 1) doubles containing half-lives belonging to particular isotopes.
		*/
		vector<double> half_lives;
		/**
		*	Sequence of integers describing an order of isotopes in decay chain or first order reaction.
		*/
		vector< vector <unsigned int> >substance_ids;
		/**
		*	Two dimensional array contains mass percentage of every single decay bifurcation on every single row.
		*/
		std::vector<std::vector<double> > bifurcation;
};

#endif
