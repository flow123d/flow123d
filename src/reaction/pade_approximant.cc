#include <stdlib.h>
#include <math.h>
#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "system/system.hh"

#include "la/distribution.hh"
#include "mesh/mesh.h"

using namespace std;
using namespace Input::Type;

Record PadeApproximant::input_type_one_decay_substep
	= Record("Substep", "Equation for reading information about radioactive decays.")
	.declare_key("parent", String(), Default::obligatory(),
				"Identifier of an isotope.")
    .declare_key("half_life", Double(), Default::optional(),
                "Half life of the parent substance.")
    .declare_key("kinetic", Double(), Default::optional(),
                "Kinetic constants describing first order reactions.")
	.declare_key("products", Array(String()), Default::obligatory(),
				"Identifies isotopes which decays parental atom to.")
	.declare_key("branch_ratios", Array(Double()), Default("1.0"),  //default is one product, with ratio = 1.0
				"Decay chain branching percentage.");


Record PadeApproximant::input_type
	= Record("PadeApproximant", "Abstract record with an information about pade approximant parameters.")
	.derive_from( ReactionTerm::input_type )
    .declare_key("decays", Array( PadeApproximant::input_type_one_decay_substep ), Default::obligatory(),
                "Description of particular decay chain substeps.")
	.declare_key("nom_pol_deg", Integer(), Default("2"),
				"Polynomial degree of the nominator of Pade approximant.")
	.declare_key("den_pol_deg", Integer(), Default("2"),
				"Polynomial degree of the nominator of Pade approximant");


PadeApproximant::PadeApproximant(Mesh &init_mesh, Input::Record in_rec)
      : LinearReaction(init_mesh, in_rec)
{
}

PadeApproximant::~PadeApproximant()
{
}

void PadeApproximant::initialize()
{
    LinearReaction::initialize();
    
    // init from input
    nom_pol_deg = input_record_.val<int>("nom_pol_deg");
    den_pol_deg = input_record_.val<int>("den_pol_deg");
    if((nom_pol_deg + den_pol_deg) < 0){
        xprintf(UsrErr, "You did not specify Pade approximant required polynomial degrees.");
    }
}

void PadeApproximant::zero_time_step()
{
    LinearReaction::zero_time_step();
}


void PadeApproximant::modify_reaction_matrix(void)
{
	Mat Denominator;
	Mat Nominator;
	Mat Pade_approximant;
	//MatFactorInfo matfact;
	PC Precond;
	//IS rperm, cperm;
	Vec tmp1; //contains the information about concentrations of all the species in one particular element
	Vec tmp2; //the same as tmp1
	//PetscInt n, m = 2;
	PetscScalar nominator_coef[nom_pol_deg];
	PetscScalar denominator_coef[den_pol_deg];
	PetscScalar Hlp_mat[1];
	PetscScalar *Array_hlp;
	//const PetscScalar *Reaction_matrix_row;
	//char dec_name[30];
	int rows, cols, i, j; //int dec_nr, dec_name_nr = 1, index, prev_index;

	//create the matrix Reaction_matrix
	MatCreate(PETSC_COMM_SELF, &Reaction_matrix);
    
    //should be probably multiplied by 2 (which is the value of m)
	MatSetSizes(Reaction_matrix, PETSC_DECIDE, PETSC_DECIDE, n_substances_, n_substances_); 
	MatSetType(Reaction_matrix, MATAIJ);
	MatSetUp(Reaction_matrix);


	//It is necessery to initialize reaction matrix here
	int index_par;
	int index_child;
	PetscScalar rel_step;
	PetscScalar extent;
    for (unsigned int i_decay = 0; i_decay < half_lives_.size(); i_decay++) {
        index_par = substance_ids_[i_decay][0];
        rel_step = time_->dt() / half_lives_[i_decay];
        extent = -log(2)*rel_step; //pow(0.5, rel_step);
        DBGMSG("parental index: %d, extent: %d", index_par, extent);
        MatSetValue(Reaction_matrix, index_par, index_par, extent,INSERT_VALUES);
        for (unsigned int i_product = 1; i_product < substance_ids_[i_decay].size(); ++i_product){
            //reaction_matrix[index_par][ substance_ids_[i_decay][i_product] ]
            extent = log(2)*rel_step* bifurcation_[i_decay][i_product-1];
            index_child = substance_ids_[i_decay][i_product];
        	MatSetValue(Reaction_matrix, index_par, index_child,extent,INSERT_VALUES);
        }
    }

	MatAssemblyBegin(Reaction_matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Reaction_matrix, MAT_FINAL_ASSEMBLY);

	//create the matrix N
    MatDuplicate(Reaction_matrix, MAT_DO_NOT_COPY_VALUES, &Nominator);

    //create the matrix D
    MatDuplicate(Reaction_matrix, MAT_DO_NOT_COPY_VALUES, &Denominator);


	//Computation of nominator in pade approximant follows
	MatZeroEntries(Nominator);
	//MatAssemblyBegin(Nominator, MAT_FINAL_ASSEMBLY);
	//MatAssemblyEnd(Nominator, MAT_FINAL_ASSEMBLY);
	for(j = nom_pol_deg; j >= 0; j--)
	{
		nominator_coef[j] = (PetscScalar) (factorial(nom_pol_deg + den_pol_deg - j) * factorial(nom_pol_deg)) 
                        / (factorial(nom_pol_deg + den_pol_deg) * factorial(j) * factorial(nom_pol_deg - j));
	}
	evaluate_matrix_polynomial(&Nominator, &Reaction_matrix, nominator_coef);
	//MatView(Nominator,PETSC_VIEWER_STDOUT_WORLD);

	//Computation of denominator in pade approximant follows
	MatZeroEntries(Denominator);
	//MatAssemblyBegin(Denominator, MAT_FINAL_ASSEMBLY);
	//MatAssemblyEnd(Denominator, MAT_FINAL_ASSEMBLY);
	for(i = den_pol_deg; i >= 0; i--)
	{
		denominator_coef[i] = (PetscScalar) pow(-1.0,i) * factorial(nom_pol_deg + den_pol_deg - i) 
                              * factorial(den_pol_deg) / (factorial(nom_pol_deg + den_pol_deg) 
                              * factorial(i) * factorial(den_pol_deg - i));
	}
	evaluate_matrix_polynomial(&Denominator, &Reaction_matrix, denominator_coef);
	//MatView(Denominator, PETSC_VIEWER_STDOUT_WORLD);



	PCCreate(PETSC_COMM_WORLD, &Precond);
	PCSetType(Precond, PCLU);
	PCSetOperators(Precond, Denominator, Denominator, DIFFERENT_NONZERO_PATTERN);
	//PCFactorSetMatOrderingType(Precond, MATORDERINGNATURAL);
	PCFactorSetMatOrderingType(Precond, MATORDERINGRCM);
	PCSetUp(Precond);

	VecCreate(PETSC_COMM_WORLD, &tmp1);
	VecSetSizes(tmp1, PETSC_DECIDE, n_substances_);
	VecSetFromOptions(tmp1);
	VecDuplicate(tmp1, &tmp2);


    //create the matrix pade
    MatCreate(PETSC_COMM_SELF, &Pade_approximant);
    
    //should be probably multiplied by 2 (which is the value of m)
    MatSetSizes(Pade_approximant, PETSC_DECIDE, PETSC_DECIDE, n_substances_, n_substances_); 
    MatSetType(Pade_approximant, MATAIJ);
    MatSetUp(Pade_approximant);

	for(rows = 0; rows < (int)( n_substances_ ); rows++){
		MatGetColumnVector(Nominator, tmp1, rows);
		//VecView(tmp1, PETSC_VIEWER_STDOUT_SELF);
		PCApply(Precond, tmp1, tmp2);
		PCView(Precond, PETSC_VIEWER_STDOUT_WORLD);
		//VecView(tmp2, PETSC_VIEWER_STDOUT_SELF);
		VecGetArray(tmp2, &Array_hlp);
		for(cols = 0; cols < (int)( n_substances_ ); cols++)
		{
			MatSetValue(Pade_approximant, rows, cols, Array_hlp[cols], ADD_VALUES);
		}
	}
	MatAssemblyBegin(Pade_approximant, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Pade_approximant, MAT_FINAL_ASSEMBLY);

	//pade assembled to reaction_matrix
	for(rows = 0; rows < (int)( n_substances_ ); rows++)
		{
			for(cols = 0; cols < (int)( n_substances_ ); cols++)
			{
				reaction_matrix_[rows][cols] = 0.0;
			}
		}
	for(rows = 0; rows < (int)( n_substances_ ); rows++)
	{
		for(cols = 0; cols < (int)( n_substances_ ); cols++)
		{
			MatGetValues(Pade_approximant, 1, &rows, 1, &cols, Hlp_mat); //&Hlp_mat[n_substances_*rows + cols]);
			reaction_matrix_[rows][cols] = (double) (Hlp_mat[0]);
		}
	}

	print_reaction_matrix(); //for visual control of equality of reaction_matrix in comparison with pade aproximant*/

	VecDestroy(&tmp1);
	VecDestroy(&tmp2);
	PCDestroy(&Precond);
	MatDestroy(&Denominator);
	MatDestroy(&Nominator);
	MatDestroy(&Pade_approximant);
}

void PadeApproximant::evaluate_matrix_polynomial(Mat *Polynomial, Mat *Reaction_matrix, PetscScalar *coef)
{
	Mat Identity;

	//create Identity matrix
	MatCreate(PETSC_COMM_SELF, &Identity);
	MatSetSizes(Identity, PETSC_DECIDE, PETSC_DECIDE, n_substances_, n_substances_); //should be probably multiplied by 2 (which is the value of m)
	MatSetType(Identity, MATAIJ);
	MatSetUp(Identity);

	MatAssemblyBegin(Identity, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Identity, MAT_FINAL_ASSEMBLY);
	MatShift(Identity, 1.0);

	for(int i = den_pol_deg; i >= 0; i--)
		{
			MatMatMult(*Polynomial, *Reaction_matrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, Polynomial);
			MatAXPY(*Polynomial, coef[i], Identity, DIFFERENT_NONZERO_PATTERN);
		}

	MatDestroy(&Identity);
}


double **PadeApproximant::compute_reaction(double **concentrations, int loc_el) //multiplication of concentrations array by reaction matrix
{
    unsigned int cols, rows;

	for(rows = 0; rows < n_substances_; rows++){
        
        prev_conc_[rows] = concentrations[rows][loc_el];
        concentrations[rows][loc_el] = 0.0;
        
        for(cols = 0; cols < n_substances_; cols++){
            concentrations[rows][loc_el] += reaction_matrix_[rows][cols]*prev_conc_[cols];
        }
    }

	return concentrations;
}



int PadeApproximant::factorial(int k)
{
    int faktor = 1;

    if(k < 0)
    {
            //an error message should be placed here
            return 0;
    }

    while(k > 1)
    {
            faktor *= k;
            k--;
    }
    //xprintf(Msg,"\n Koeficient has a value %d.\n",faktor);
    return faktor;
}
