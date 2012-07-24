/*
 * pade_approximant.h
 *
 *  Created on: Apr 2, 2012
 *      Author: lukas
 */

#ifndef PADE_APPROXIMANT_H_
#define PADE_APPROXIMANT_H_

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"
#include "reaction/linear_reaction.hh"

class Pade_approximant: public Linear_reaction
{
	public:
        /**
         *  Constructor with parameter for initialization of a new declared class member
         *  TODO: parameter description
         */
		Pade_approximant(TimeMarks &marks, Mesh &mesh, MaterialDatabase &mat_base); //(double timeStep, Mesh * mesh, int nrOfSpecies, bool dualPorosity); //(double time_step, int nrOfElements, double ***ConcentrationMatrix);
		/**
		*	Destructor.
		*/
		~Pade_approximant(void);
		/**
		*
		*/
		void set_time_step(double new_timestep);
	protected:
		/**
		*
		*/
		void modify_reaction_matrix(void);
		/**
		*
		*/
		double **allocate_reaction_matrix(void);
		/**
		*
		*/
		double **modify_reaction_matrix(int dec_nr);
		/**
		*	Evaluates Pade approximant from Reaction_matrix.
		*/
		double **modify_reaction_matrix_repeatedly(void);
		/**
		* It enables to evaluate matrix nominator and denominator present in Pade approximant.
		*/
		void evaluate_matrix_polynomial(Mat *Polynomial, Mat *Reaction_matrix, PetscScalar *coef);
		/**
		* PETSC format of a matrix describing linear chemical reaction.
		*/
		Mat Reaction_matrix;
};

#endif /* PADE_APPROXIMANT_H_ */
