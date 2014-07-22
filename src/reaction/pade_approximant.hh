/*
 * pade_approximant.h
 *
 *  Created on: Apr 2, 2012
 *      Author: lukas
 */

#ifndef PADE_APPROXIMANT_H_
#define PADE_APPROXIMANT_H_

#include <input/input_type.hh>

// #include "petscvec.h"
// #include "petscmat.h"
// #include "petscksp.h"

class Mesh;
class Distribution;
class LinearReaction;

class PadeApproximant: public LinearReaction
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
   
    /// Constructor.
	PadeApproximant(Mesh &mesh, Input::Record in_rec);

	/// Destructor.
	~PadeApproximant(void);

    void initialize() override;
	void zero_time_step() override;

protected:
    /**
    *   Evaluates Pade approximant from Reaction_matrix.
    */
    void modify_reaction_matrix(void) override;
    void modify_reaction_matrix2(void);
    
    // Evaluate nominator and denominator coeficients of PadeApproximant for exponencial function.
    void compute_exp_coefs(unsigned int nominator_degree, unsigned int denominator_degree,
                           std::vector<double> &nominator_coefs, std::vector<double> &denominator_coefs);
    
    /**
    * It enables to evaluate matrix nominator and denominator present in Pade approximant.
    */
    //void evaluate_matrix_polynomial(Mat *Polynomial, Mat *Reaction_matrix, PetscScalar *coef);
    
    /// Computes factorial of @p k.
    int factorial(int k);
            
	/**
	*	Integer which informs about the order of a polynomial term in nominator of Pade approximant rational term.
	*/
	int nom_pol_deg;
	/**
	*	Integer which informs about the order of a polynomial term in denominator of Pade approximant rational term.
	*/
	int den_pol_deg;

	/**
	* PETSC format of a matrix describing linear chemical reaction.
	*/
// 	//Mat Reaction_matrix;
};

#endif // PADE_APPROXIMANT_H_
