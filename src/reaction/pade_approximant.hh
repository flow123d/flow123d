/*
 * pade_approximant.h
 *
 *  Created on: Apr 2, 2012
 *      Author: lukas
 */

#ifndef PADE_APPROXIMANT_H_
#define PADE_APPROXIMANT_H_

#include <vector>
#include <input/input_type.hh>
#include <input/accessors.hh>

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"


class Mesh;
class Distribution;
class Reaction;
class Linear_reaction;

class Pade_approximant: public Linear_reaction
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
		Pade_approximant(Mesh &mesh, Input::Record in_rec, vector<string> &names);

		/**
		*	Destructor.
		*/
		~Pade_approximant(void);

		/**
		*	For simulation of chemical reaction in just one element either inside of MOBILE or IMMOBILE pores.
		*/
		double **compute_reaction(double **concentrations, int loc_el) override;

		/**
		*	Evaluates Pade approximant from Reaction_matrix.
		*/
		double **modify_reaction_matrix(void) override;
		/**
		* It enables to evaluate matrix nominator and denominator present in Pade approximant.
		*/
		void evaluate_matrix_polynomial(Mat *Polynomial, Mat *Reaction_matrix, PetscScalar *coef);
		
	protected:
                /**
                *       Enables to compute factorial k!.
                */
                int factorial(int k);
                
                /**
                 * Initialize from input interface.
                 */
                void init_from_input(Input::Record in_rec) override;
                
		/**
		*	This method reads a sequence of numbers defining an order of substances in decay chain. The string section defines where too look for indices inside of ini-file, whereas n_subst is a number of isotopes in described decay chain.
		*/
		//int *set_indeces(char *section, int n_subst);
		/**
		*	This method reads an information about a number of isotopes in a decay chain described inside of ini-file in section given as an argument. This method is used for radioactive decay simulation.
		*/
		//void set_nr_of_isotopes(char* section);
		/**
		*	This method sets number of isotopes for the case of first order reaction. The value should be always 2.
		*/
		//void set_nr_of_isotopes(int Nr_of_isotopes);
		/**
		*	This method reads a sequence of (nr_of_isotopes - 1) halflives belonging to separate decay chain step. This information is placed in ini-file in a block starting with a string section.
		*/
		//double *set_half_lives(char *section);
		/**
		*	This method reads form ini-file an information for construction of a matrix describing bifurcation of every single decay chain on one row of the reaction matrix. Informations about bifurcation are placed in a block starting with a string section. dec_nr identifies which one decay chain is handled and which row of twodimensional bifurcation matrix (double **array)should be affected.
		*/
		//void set_bifurcation(char *section, int dec_nr);
		/**
		*	This method reads from ini-file an information if the bifurcation for a current decay chain is switched on in a block starting with a string section. Initialy bifurcation is switched of.
		*/
		//void set_bifurcation_on(char *section);
		/**
		*	This method reads from ini-file an information if first order reactions simulation is switched on.
		*/
		//void set_For_on(void);
		/**
		*	This method reads from ini-file an information if a radioactive decay simulation is switched on.
		*/
		//void set_decay_on(void);
		/**
		*	This method reads from ini-file an information and prepares a vector (onedimensional double *array) containing kinetic constants of every single first order reactions. Those informations are placed in a block with a string section at the beginning. From those constants half-lives belonging to first order reactions are computed.
		*/
		//void set_kinetic_constants(char *section, int reaction_nr);
		/**
		*
		*/
		double **allocate_reaction_matrix(void);
		/**
		*	This method modificates reaction matrix as described in ini-file a single section [Decay_i] or [FoReact_i]. It is used when bifurcation is switched off.
		*/
		//double **modify_reaction_matrix(void);
		/**
		*	For control printing of a matrix describing simple chemical reactions.
		*/
		//void print_reaction_matrix(void);
		/**
		*	For printing nr_of_isotopes identifies of isotopes in a current decay chain.
		*/
		//void print_indeces(int n_subst);
		/**
		* Following method releases reaction matrix to make it possible to set a new time step for chemistry.
		*/
		//void release_reaction_matrix();
		/**
		*	For printing (nr_of_isotopes - 1) doubles containing half-lives belonging to particular isotopes on screen.
		*/
		//void print_half_lives(int n_subst);
		/**
		*	Small (nr_of_species x nr_of_species) square matrix for realization of radioactive decay and first order reactions simulation.
		*/
		//double **reaction_matrix;
		//std::vector<double> half_lives; ///< alternative to following row
		/**
		*	Sequence of (nr_of_isotopes - 1) doubles containing half-lives belonging to particular isotopes.
		*/
		//double *half_lives;
		/**
		*	Sequence of integers describing an order of isotopes in decay chain or first order reaction.
		*/
		//int *substance_ids;
		/**
		*	Informs about the number of isotopes in a current decay chain.
		*/
		//int nr_of_isotopes;
		/**
		*	Two dimensional array contains mass percentage of every single decay bifurcation on every single row.
		*/
		//std::vector<std::vector<double> > bifurcation;
		/**
		*	One dimensional array of kinetic constants belonging to considered reactions.
		*/
		//std::vector<double> kinetic_constant;
		/**
		*	Boolean which enables to turn on branching of considered decay chain.
		*/
		//bool bifurcation_on;
		/**
		*	Integer which informs about the order of a polynomial term in nominator of Pade approximant rational term.
		*/
		int nom_pol_deg;
		/**
		*	Integer which informs about the order of a polynomial term in denominator of Pade approximant rational term.
		*/
		int den_pol_deg;
		/**
		*
		*/
		//void modify_reaction_matrix(void);
		/**
		*
		*/
		//double **allocate_reaction_matrix(void);
		/**
		*
		*/
		//double **modify_reaction_matrix(int dec_nr);
		/**
		* PETSC format of a matrix describing linear chemical reaction.
		*/
		Mat Reaction_matrix;
};

#endif /* PADE_APPROXIMANT_H_ */
