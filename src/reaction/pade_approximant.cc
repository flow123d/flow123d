#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "system/system.hh"
//#include "transport/transport.h"
//#include "system/par_distribution.hh"
#include "la/distribution.hh"
#include "mesh/mesh.h"

#define MOBILE 0
#define IMMOBILE 1

using namespace Input::Type;

Record Pade_approximant::input_type_one_decay_substep
	= Record("Substep", "Equation for reading information about radioactive decays.")
	.declare_key("parent", String(), Default::obligatory(),
				"Identifier of an isotope.")
    .declare_key("half_life", Double(), Default::optional(),
                "Half life of the parent substance.")
    .declare_key("kinetic", Double(), Default::optional(),
                "Kinetic constants describing first order reactions.")
	.declare_key("products", Array(String()), Default::obligatory(),
				"Identifies isotopes which decays parental atom to.")
	.declare_key("branch_ratios", Array(Double()), Default("1.0"),   // default is one product, with ratio == 1.0
				"Decay chain branching percentage.");


Record Pade_approximant::input_type
	= Record("PadeApproximant", "Abstract record with an information about pade approximant parameters.")
	.derive_from( Reaction::input_type )
    .declare_key("decays", Array( Pade_approximant::input_type_one_decay_substep ), Default::obligatory(),
                "Description of particular decay chain substeps.")
	.declare_key("nom_pol_deg", Integer(), Default("2"),
				"Polynomial degree of the nominator of Pade approximant.")
	.declare_key("den_pol_deg", Integer(), Default("2"),
				"Polynomial degree of the nominator of Pade approximant");


using namespace std;


Pade_approximant::Pade_approximant(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)
      : Linear_reaction(init_mesh, in_rec, names)
{
	init_from_input(in_rec);
	DBGMSG("Pade_approximant constructor is running.\n");
	allocate_reaction_matrix();
}

Pade_approximant::~Pade_approximant()
{
	/*int i, rows, n_subst;

	if(half_lives != NULL){
		free(half_lives);
		half_lives = NULL;
	}

	if(substance_ids != NULL){
		free(substance_ids);
		substance_ids = NULL;
	}

	release_reaction_matrix();*/

	DBGMSG("Pade approximant destructor is running.\n");
}

double **Pade_approximant::allocate_reaction_matrix(void) //reaction matrix initialization
{
	unsigned int rows, cols;

	cout << "We are going to allocate reaction matrix" << endl;
	if (reaction_matrix == NULL) reaction_matrix = (double **)xmalloc(n_substances() * sizeof(double*));//allocation section
	for(rows = 0; rows < n_substances(); rows++){
		reaction_matrix[rows] = (double *)xmalloc(n_substances() * sizeof(double));
	}
	for(rows = 0; rows < n_substances();rows++){
	 for(cols = 0; cols < n_substances(); cols++){
		 reaction_matrix[rows][cols] = 0.0;
	 }
	}
	print_reaction_matrix();
	return reaction_matrix;
}

double **Pade_approximant::modify_reaction_matrix(void)
{
	Mat Denominator;
	Mat Nominator;
	//Mat Reaction_matrix;
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
	MatSetSizes(Reaction_matrix, PETSC_DECIDE, PETSC_DECIDE, n_substances(), n_substances()); //should be probably multiplied by 2 (which is the value of m)
	MatSetType(Reaction_matrix, MATAIJ);
	MatSetUp(Reaction_matrix);


	//It is necessery to initialize reaction matrix here
	int index_par;
	int index_child;
	PetscScalar rel_step;
	PetscScalar extent;
    for (unsigned int i_decay = 0; i_decay < half_lives.size(); i_decay++) {
        index_par = substance_ids[i_decay][0];
        rel_step = time_->dt() / half_lives[i_decay];
        extent = -log(2)*rel_step; //pow(0.5, rel_step);
        cout<<"parental index" << index_par << ", extent "<< extent << endl;
        MatSetValue(Reaction_matrix, index_par, index_par, extent,INSERT_VALUES);
        for (unsigned int i_product = 1; i_product < substance_ids[i_decay].size(); ++i_product){
            //reaction_matrix[index_par][ substance_ids[i_decay][i_product] ]
            extent = log(2)*rel_step* bifurcation[i_decay][i_product-1];
            index_child = substance_ids[i_decay][i_product];
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
		nominator_coef[j] = (PetscScalar) (factorial(nom_pol_deg + den_pol_deg - j) * factorial(nom_pol_deg)) / (factorial(nom_pol_deg + den_pol_deg) * factorial(j) * factorial(nom_pol_deg - j));
	}
	evaluate_matrix_polynomial(&Nominator, &Reaction_matrix, nominator_coef);
	//MatView(Nominator,PETSC_VIEWER_STDOUT_WORLD);

	//Computation of denominator in pade approximant follows
	MatZeroEntries(Denominator);
	//MatAssemblyBegin(Denominator, MAT_FINAL_ASSEMBLY);
	//MatAssemblyEnd(Denominator, MAT_FINAL_ASSEMBLY);
	for(i = den_pol_deg; i >= 0; i--)
	{
		denominator_coef[i] = (PetscScalar) pow(-1.0,i) * factorial(nom_pol_deg + den_pol_deg - i) * factorial(den_pol_deg) / (factorial(nom_pol_deg + den_pol_deg) * factorial(i) * factorial(den_pol_deg - i));
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
	VecSetSizes(tmp1, PETSC_DECIDE, n_substances());
	VecSetFromOptions(tmp1);
	VecDuplicate(tmp1, &tmp2);


    //create the matrix pade
    MatCreate(PETSC_COMM_SELF, &Pade_approximant);
    MatSetSizes(Pade_approximant, PETSC_DECIDE, PETSC_DECIDE, n_substances(), n_substances()); //should be probably multiplied by 2 (which is the value of m)
    MatSetType(Pade_approximant, MATAIJ);
    MatSetUp(Pade_approximant);

	for(rows = 0; rows < (int)( n_substances() ); rows++){
		MatGetColumnVector(Nominator, tmp1, rows);
		//VecView(tmp1, PETSC_VIEWER_STDOUT_SELF);
		PCApply(Precond, tmp1, tmp2);
		PCView(Precond, PETSC_VIEWER_STDOUT_WORLD);
		//VecView(tmp2, PETSC_VIEWER_STDOUT_SELF);
		VecGetArray(tmp2, &Array_hlp);
		for(cols = 0; cols < (int)( n_substances() ); cols++)
		{
			MatSetValue(Pade_approximant, rows, cols, Array_hlp[cols], ADD_VALUES);
		}
	}
	MatAssemblyBegin(Pade_approximant, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Pade_approximant, MAT_FINAL_ASSEMBLY);

	//pade assembled to reaction_matrix
	for(rows = 0; rows < (int)( n_substances() ); rows++)
		{
			for(cols = 0; cols < (int)( n_substances() ); cols++)
			{
				reaction_matrix[rows][cols] = 0.0;
			}
		}
	for(rows = 0; rows < (int)( n_substances() ); rows++)
	{
		for(cols = 0; cols < (int)( n_substances() ); cols++)
		{
			MatGetValues(Pade_approximant, 1, &rows, 1, &cols, Hlp_mat); //&Hlp_mat[n_substances()*rows + cols]);
			reaction_matrix[rows][cols] = (double) (Hlp_mat[0]);
		}
	}

	print_reaction_matrix(); //for visual control of equality of reaction_matrix in comparison with pade aproximant*/

	VecDestroy(&tmp1);
	VecDestroy(&tmp2);
	PCDestroy(&Precond);
	MatDestroy(&Denominator);
	MatDestroy(&Nominator);
	//MatDestroy(Reaction_matrix);
	MatDestroy(&Pade_approximant);

	return reaction_matrix;
}

void Pade_approximant::evaluate_matrix_polynomial(Mat *Polynomial, Mat *Reaction_matrix, PetscScalar *coef)
{
	Mat Identity;

	//create Identity matrix
	MatCreate(PETSC_COMM_SELF, &Identity);
	MatSetSizes(Identity, PETSC_DECIDE, PETSC_DECIDE, n_substances(), n_substances()); //should be probably multiplied by 2 (which is the value of m)
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

	return;
}


double **Pade_approximant::compute_reaction(double **concentrations, int loc_el) //multiplication of concentrations array by reaction matrix
{
    unsigned int cols, rows;

    if (reaction_matrix == NULL) return concentrations;

	for(cols = 0; cols < n_substances(); cols++){
		prev_conc[cols] = concentrations[cols][loc_el];
		concentrations[cols][loc_el] = 0.0;
	}

	for(cols = 0; cols < n_substances(); cols++){
		for(rows = 0; rows < n_substances(); rows++){
            concentrations[cols][loc_el] += reaction_matrix[cols][rows]*prev_conc[rows];
        }
    }

	return concentrations;
}


void Pade_approximant::init_from_input(Input::Record in_rec)
{
	nom_pol_deg = in_rec.val<int>("nom_pol_deg");
	den_pol_deg = in_rec.val<int>("den_pol_deg");
	if((nom_pol_deg + den_pol_deg) < 0){
		cout << "You did not specify Pade approximant required polynomial degrees." << endl;
		//TODO: This occasion should cause an error.
		//break;
	}
}


int Pade_approximant::factorial(int k)
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
