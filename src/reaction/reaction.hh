/** @brief class Linear_reaction is used to enable simulation of simple chemical reactions
 *
 * Class in this file makes it possible to realize  simulation of reaction of the first order by simple matrix multiplication.
 * One step of the linear reaction is represented as a product of a matrix containing concentrations of observed speciesin elements in rows multiplied by so called
 * reaction_matrix. Through this way radioactive decay can bee also realized and that was exactly what we did at the begining of journey. :-)
 * Matrix containing concentrations has a dimension Nxn, where N is a number of elements in mesh and n denotes a number of transported chemical species.
 * The reaction_matrix is a square matrix and it has a dimension nxn.
 *
 */
#ifndef REACT
#define REACT

#include "input/accessors.hh"
#include "coupling/equation.hh"
//class Mesh;
class Distribution;

enum Reaction_type {No_reaction, Linear_react, Linear_react_Pade, General_react_Semch};

class Reaction: public EquationBase
{
	public:
		/*
		 * Static method for new input data types input
		 */
		static Input::Type::AbstractRecord & get_input_type();
		/*
		 * Static method for new input data types input
		*/
		static Input::Type::Record & get_one_decay_type();
        /**
         *  Constructor with parameter for initialization of a new declared class member
         *  TODO: parameter description
         */
		Reaction(TimeMarks &marks, Mesh &init_mesh, MaterialDatabase &material_database, Input::Record in_rec, const std::vector<string> &names); //(double timeStep, Mesh * mesh, int nrOfSpecies, bool dualPorosity); //(double time_step, int nrOfElements, double ***ConcentrationMatrix);
		/**
		*	Destructor.
		*/
		~Reaction(void);
		/**
		*	For simulation of chemical raection in just one element either inside of MOBILE or IMMOBILE pores.
		*/
		virtual double **compute_reaction(double **concentrations, int loc_el);
		/**
		*	Prepared to compute simple chemical reactions inside all of considered elements. It calls compute_reaction(...) for all the elements controled by concrete processor, when the computation is paralelized.
		*/
		virtual void compute_one_step(void);

		/**
		 * Returns number of substances involved in reactions. This should be same as number of substances in transport.
		 */
        inline unsigned int n_substances()
        { return names_.size(); }


		/**
		* 	It returns current time step used for first order reactions.
		*/
		double get_time_step(void);
		/**
		*	This method enables to change a data source the program is working with, during simulation.
		*/
		void set_concentration_matrix(double ***ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc);
		/**
		*	This method enables to change the timestep for computation of simple chemical reactions. Such a change is conected together with creating of a new reaction matrix necessity.
		*/
		virtual void set_time_step(double new_timestep);
		/**
		* Folowing method enabels the timestep for chemistry to have the value written in ini-file.
		*/
		virtual void set_time_step(Input::Record in_rec);
		//
		void update_solution(void);
		void choose_next_time(void);
		void set_time_step_constrain(double dt);
		void get_parallel_solution_vector(Vec &vc);
		void get_solution_vector(double* &vector, unsigned int &size);
		/**
		* Function for setting dual porosity.
		*/
		void set_dual_porosity(bool dual_porosity_on);//(Input::Record in_rec);
	protected:
		/**
		*	This method disables to use constructor without parameters.
		*/
		Reaction();
		/**
		*	Enables to compute factorial k!.
		*/
		int faktorial(int k);
		/**
		*	Finds a position of a string in specified array.
		*/
		unsigned int find_subst_name(const std::string &name);
		/**
		*	Boolean which enables to compute reactions also in immobile pores.
		*/
		bool dual_porosity_on;
		/**
		*	Holds the double describing time step for radioactive decay or first order reactions simulations.
		*/
		double time_step;
		/**
		*	Informs how many firts order reactions of the type A -> B are under consideration. It is a number of [FoReaction_i] in ini-file.
		*/
		//int nr_of_FoR; //Obsolete variable.
		/**
		*	Pointer to threedimensional array[mobile/immobile][species][elements] containing concentrations.
		*/
		double ***concentration_matrix;
		/**
		*	Pointer to reference to distribution of elements between processors.
		*/
		Distribution *distribution;
		/**
		*	Pointer to reference previous concentration array used in compute_reaction().
		*/
		double *prev_conc;
		/**
		* Names belonging to substances. Should be same as in the transport.
		*/
		vector<string> names_;
};

#endif
