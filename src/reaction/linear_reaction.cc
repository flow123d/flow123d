#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include "reaction/linear_reaction.hh"
#include "system/system.hh"
#include "materials.hh"
#include "transport/transport.h"
#include "system/par_distribution.hh"
#include "mesh/mesh.h"


using namespace std;

Linear_reaction::Linear_reaction(double timeStep, Mesh * mesh, int nrOfSpecies, bool dualPorosity) //(double timestep, int nrOfElements, double ***ConvectionMatrix)
	: half_lives(NULL), substance_ids(NULL), reaction_matrix(NULL), bifurcation_on(false), dual_porosity_on(false), prev_conc(NULL), matrix_exp_on(false)
{
	nr_of_decays = OptGetInt("Reaction_module","Nr_of_decay_chains","0");
	nr_of_FoR = OptGetInt("Reaction_module","Nr_of_FoR","0");
	nr_of_species = OptGetInt("Transport", "N_substances", "0");
	matrix_exp_on = OptGetBool("Reaction_module","Matrix_exp_on","No");
	if(matrix_exp_on == true)
	{
		nom_pol_deg = OptGetInt("Reaction_module","Nom_pol_deg","2");
		den_pol_deg = OptGetInt("Reaction_module","Den_pol_deg","2");
	}
	//dual_porosity_on = dualPorosity;
	set_dual_porosity();
	set_mesh_(mesh);
	set_nr_of_elements(mesh->n_elements());
	cout << "number of FoR is "<< nr_of_FoR << endl;
	cout << "number of decays is " << nr_of_decays << endl;
	cout << "number of species is " << nr_of_species << endl;
	if(prev_conc != NULL){
		free(prev_conc);
		prev_conc = NULL;
	}
	prev_conc = (double *)xmalloc(nr_of_species * sizeof(double));
	if(timeStep > 1e-12) this->set_time_step(timeStep); else this->set_time_step(0.5);
	if((nr_of_decays > 0) || (nr_of_FoR > 0)){
		allocate_reaction_matrix();
		if(matrix_exp_on == true)
		{
			modify_reaction_matrix_using_pade();
		}else{
			modify_reaction_matrix_repeatedly();
		}
	}
}

Linear_reaction::~Linear_reaction()
{
	int i, rows, n_subst;

	//n_subst = sizeof(*reaction_matrix)/sizeof(double *);
	if(half_lives != NULL){
		free(half_lives);
		half_lives = NULL;
	}

	if(substance_ids != NULL){
		free(substance_ids);
		substance_ids = NULL;
	}

	if(prev_conc != NULL){
		free(prev_conc);
		prev_conc = NULL;
	}

	release_reaction_matrix();
}

double **Linear_reaction::allocate_reaction_matrix(void) //reaction matrix initialization
{
	int index, rows, cols, dec_nr, prev_index;
	char dec_name[30];

	cout << "We are going to allocate reaction matrix" << endl;
	reaction_matrix = (double **)xmalloc(nr_of_species * sizeof(double*));//allocation section
	for(rows = 0; rows < nr_of_species; rows++){
		reaction_matrix[rows] = (double *)xmalloc(nr_of_species * sizeof(double));
	}
	for(rows = 0; rows < nr_of_species;rows++){
	 for(cols = 0; cols < nr_of_species; cols++){
		 if(rows == cols){
			 if(matrix_exp_on == false) reaction_matrix[rows][cols] = 1.0;
		 }else{
			reaction_matrix[rows][cols] = 0.0;
		 }
	 }
	}
	//print_reaction_matrix();
	return reaction_matrix;
}

double **Linear_reaction::modify_reaction_matrix(void) //prepare the matrix, which describes reactions
{
	int rows,cols, index, prev_index;
	double rel_step, prev_rel_step;

	if(reaction_matrix == NULL){
		xprintf(Msg,"\nReaction matrix pointer is NULL.\n");
		return NULL;
	}
	if((nr_of_decays > 0) || (nr_of_FoR > 0)){
		for(cols = 0; cols < nr_of_isotopes; cols++){
			index = substance_ids[cols] - 1; // because indecees in input file run from one whereas indeces in C++ run from ZERO
			if(cols < (nr_of_isotopes - 1)){
				rel_step = time_step/half_lives[cols];
			}
			if(cols > 0){
				reaction_matrix[prev_index][prev_index] = pow(0.5,prev_rel_step);
				reaction_matrix[prev_index][index] += (1 - pow(0.5,prev_rel_step));
			}
			prev_rel_step = rel_step;
			prev_index = index;
		}
	}
	print_reaction_matrix();//just for control print
	return reaction_matrix;
}

double **Linear_reaction::modify_reaction_matrix(int dec_nr) //prepare the matrix, which describes reactions, takes bifurcation in acount
{
	int rows,cols, index, first_index, bif_id;
	double rel_step, prev_rel_step;

	if(reaction_matrix == NULL){
		xprintf(Msg,"\nReaction matrix pointer is NULL.\n");
		return NULL;
	}

	first_index = substance_ids[0]-1;
	for(cols = 0; cols < nr_of_isotopes; cols++){
		index = substance_ids[cols] - 1; // because indecees in input file run from one whereas indeces in C++ run from ZERO
		if(cols < (nr_of_isotopes -1)){
			rel_step = time_step/half_lives[cols];
			xprintf(Msg,"time_step %f\n", time_step);
		}
		if(cols > 0){
			bif_id = cols -1;
			reaction_matrix[first_index][first_index] = pow(0.5,prev_rel_step); //bifurcation[dec_nr][bif_id] * pow(0.5,prev_rel_step);
			reaction_matrix[first_index][index] += (1 - pow(0.5,prev_rel_step)) * bifurcation[dec_nr][bif_id];
		}
		prev_rel_step = rel_step;
	}
	print_reaction_matrix();//just for control print
	return reaction_matrix;
}

double **Linear_reaction::modify_reaction_matrix_repeatedly(void)
{
	char dec_name[30];
	int rows, cols, dec_nr, dec_name_nr = 1, index, prev_index;

	if(nr_of_decays > 0){
		xprintf(Msg,"\nNumber of decays is %d\n",nr_of_decays);
		if(half_lives != NULL){
					free(half_lives);
					half_lives = NULL;
		}
		half_lives = (double *)xmalloc(nr_of_decays * sizeof(double));
		bifurcation.resize(nr_of_decays);
		for(dec_nr = 0; dec_nr < nr_of_decays; dec_nr++){
			sprintf(dec_name,"Decay_%d", dec_name_nr);
			nr_of_isotopes = OptGetInt(dec_name,"Nr_of_isotopes","0");
			set_half_lives(dec_name);
			set_indeces(dec_name, nr_of_isotopes);
			print_indeces(nr_of_isotopes); //just a control
			print_half_lives(nr_of_isotopes); //just a control
			bifurcation_on = OptGetBool(dec_name,"Bifurcation_on","no");
			if(bifurcation_on == true){
				set_bifurcation(dec_name, dec_nr);
				modify_reaction_matrix(dec_nr);
			}else{
				modify_reaction_matrix();
			}
			dec_name_nr++;
		}
	}
	if(nr_of_FoR > 0){
		xprintf(Msg,"\nNumber of first order reactions is %d\n",nr_of_FoR);
		//half_lives.resize(nr_of_FoR); //does not function at all
		if(half_lives != NULL){
			free(half_lives);
			half_lives = NULL;
		}
		half_lives = (double *)xmalloc(nr_of_FoR * sizeof(double));
		for(dec_nr = 0; dec_nr < nr_of_FoR; dec_nr++){
			sprintf(dec_name,"FoReact_%d", dec_name_nr);
			set_nr_of_isotopes(2);
			set_indeces(dec_name, 2);
			set_kinetic_constants(dec_name, dec_nr);//instead of this line, here should be palced computation of halflives using kinetic constants
			print_indeces(nr_of_isotopes); //just a control
			print_half_lives(2); //just a control
			//modify_reaction_matrix(2);
			modify_reaction_matrix();
			dec_name_nr++;
		}
	}
	return reaction_matrix;
}

void Linear_reaction::modify_reaction_matrix(Mat *R) //prepare the matrix, which describes reactions
{
	int rows,cols, index, prev_index;
	double rel_step, prev_rel_step, Hlp_kin;

	if(R == NULL){
		xprintf(Msg,"\nPointer to the reaction matrix R is NULL.\n");
		return;
	}
	if((nr_of_decays > 0) || (nr_of_FoR > 0)){
		for(cols = 0; cols < nr_of_isotopes; cols++){
			index = substance_ids[cols] - 1; // because indecees in input file run from one whereas indeces in C++ run from ZERO
			if(cols > 0){
				Hlp_kin = time_step*log(2)/half_lives[cols-1];
				//xprintf(Msg, "\ncols identifier is %d\n",cols);
				//xprintf(Msg, "\nHlp_kin coef has the value %f an time_step is equal to %f, half-live is %f\n",Hlp_kin, time_step, half_lives[cols-1]);
				MatSetValue(*R, prev_index, prev_index, Hlp_kin, INSERT_VALUES);
				//MatSetValues(*R, 1, &prev_index, 1, &prev_index, &Hlp_kin, INSERT_VALUES);
				//reaction_matrix[prev_index][prev_index] = pow(0.5,prev_rel_step);
				MatSetValue(*R, index, prev_index, Hlp_kin, INSERT_VALUES);
				//reaction_matrix[prev_index][index] += (1 - pow(0.5,prev_rel_step));
			}
			prev_index = index;
		}
	}
	return;
}

double **Linear_reaction::modify_reaction_matrix(Mat *R, int dec_nr) //prepare the matrix, which describes reactions, takes bifurcation in acount
{
	int rows,cols, index, first_index, bif_id;
	double rel_step, prev_rel_step;
	PetscScalar Hlp_kin;

	if(reaction_matrix == NULL){
		xprintf(Msg,"\nReaction matrix pointer is NULL.\n");
		return NULL;
	}

	first_index = substance_ids[0]-1;
	Hlp_kin = time_step*log(2)/half_lives[0];
	MatSetValue(*R,first_index,first_index,Hlp_kin,INSERT_VALUES);
	for(cols = 0; cols < nr_of_isotopes; cols++){
		index = substance_ids[cols] - 1; // because indecees in input file run from one whereas indeces in C++ run from ZERO
		if(cols > 0){
			bif_id = cols -1;
			Hlp_kin = time_step*log(2)/half_lives[cols-1] * bifurcation[dec_nr][bif_id];
			MatSetValue(*R, index, first_index, Hlp_kin, INSERT_VALUES);
		}
	}
	return reaction_matrix;
}

double **Linear_reaction::modify_reaction_matrix_using_pade(void)
{
	Mat D;
	Mat N;
	Mat R;
	Mat pade;
	Mat Hlp;
	Mat Hlp2;
	Mat B;
	PetscInt n, m = 2;
	PetscScalar koef_hlp;
	char dec_name[30];
	int rows, cols, dec_nr, dec_name_nr = 1, index, prev_index;
	PetscScalar *Hlp_mat;
	IS rperm, cperm;
	MatFactorInfo matfact;

	//create the matrix D
	MatCreate(PETSC_COMM_WORLD, &D);
	MatSetSizes(D, PETSC_DECIDE, PETSC_DECIDE, nr_of_species, nr_of_species); //nr_of_species should be probably multiplied by 2 (which is the value of m), but I do not know why
	MatSetType(D, MATAIJ);
	MatZeroEntries(D);

	//create the matrix N
	//MatDuplicate(D, MAT_COPY_VALUES, &N);
	MatCreate(PETSC_COMM_WORLD, &N);
	MatSetSizes(N, PETSC_DECIDE, PETSC_DECIDE, nr_of_species, nr_of_species); //should be probably multiplied by 2 (which is the value of m)
	MatSetType(N, MATAIJ);
	MatZeroEntries(N);

	//create the matrix R
	//MatDuplicate(D, MAT_COPY_VALUES, &R);
	MatCreate(PETSC_COMM_WORLD, &R);
	MatSetSizes(R, PETSC_DECIDE, PETSC_DECIDE, nr_of_species, nr_of_species); //should be probably multiplied by 2 (which is the value of m)
	MatSetType(R, MATAIJ);
	MatZeroEntries(R);

	//create the matrix pade
	//MatDuplicate(D, MAT_COPY_VALUES, &pade);
	MatCreate(PETSC_COMM_WORLD, &pade);
	MatSetSizes(pade, PETSC_DECIDE, PETSC_DECIDE, nr_of_species, nr_of_species); //should be probably multiplied by 2 (which is the value of m)
	MatSetType(pade, MATAIJ);
	MatZeroEntries(pade);

	//create the matrix Hlp
	//MatDuplicate(D, MAT_COPY_VALUES, &Hlp);
	MatCreate(PETSC_COMM_WORLD, &Hlp);
	MatSetSizes(Hlp, PETSC_DECIDE, PETSC_DECIDE, nr_of_species, nr_of_species); //should be probably multiplied by 2 (which is the value of m)
	MatSetType(Hlp, MATAIJ);
	MatZeroEntries(Hlp);

	//create the matrix Hlp2
	//MatDuplicate(D, MAT_COPY_VALUES, &Hlp2);
	MatCreate(PETSC_COMM_WORLD, &Hlp2);
	MatSetSizes(Hlp2, PETSC_DECIDE, PETSC_DECIDE, nr_of_species, nr_of_species); //should be probably multiplied by 2 (which is the value of m)
	MatSetType(Hlp2, MATAIJ);
	MatZeroEntries(Hlp2);

	if(nr_of_decays > 0){
		xprintf(Msg,"\nNumber of decays is %d\n",nr_of_decays);
		if(half_lives != NULL){
			free(half_lives);
			half_lives = NULL;
		}
		half_lives = (double *)xmalloc(nr_of_decays * sizeof(double));
		bifurcation.resize(nr_of_decays);
		for(dec_nr = 0; dec_nr < nr_of_decays; dec_nr++){
			sprintf(dec_name,"Decay_%d", dec_name_nr);
			nr_of_isotopes = OptGetInt(dec_name,"Nr_of_isotopes","0");
			set_half_lives(dec_name);
			set_indeces(dec_name, nr_of_isotopes);
			print_indeces(nr_of_isotopes); //just a control
			print_half_lives(nr_of_isotopes); //just a control
			bifurcation_on = OptGetBool(dec_name,"Bifurcation_on","no");
			if(bifurcation_on == true){
				set_bifurcation(dec_name, dec_nr);
				modify_reaction_matrix(&R, dec_nr);
			}else{
				if(&R != NULL)
				{
					modify_reaction_matrix(&R);
					xprintf(Msg,"Reaction matrix R is has been allocated. The addres is %d.\n", R);
				}else{
					xprintf(Msg,"Reaction matrix R is has not been allocated.\n"); //cout << "Reaction matrix R is has not been allocated." << endl;
				}
			}
			dec_name_nr++;
		}
		//MatAssemblyBegin(R, MAT_FINAL_ASSEMBLY);
		//MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY);
	}
	if(nr_of_FoR > 0){
		xprintf(Msg,"\nNumber of first order reactions is %d\n",nr_of_FoR);
		//half_lives.resize(nr_of_FoR); //does not function at all
		if(half_lives != NULL){
			free(half_lives);
			half_lives = NULL;
		}
		half_lives = (double *)xmalloc(nr_of_FoR * sizeof(double));
		for(dec_nr = 0; dec_nr < nr_of_FoR; dec_nr++){
			sprintf(dec_name,"FoReact_%d", dec_name_nr);
			set_nr_of_isotopes(2);
			set_indeces(dec_name, 2);
			set_kinetic_constants(dec_name, dec_nr);//instead of this line, here should be palced computation of halflives using kinetic constants
			print_indeces(nr_of_isotopes); //just a control
			print_half_lives(2); //just a control
			//modify_reaction_matrix(2);
			modify_reaction_matrix(&R, 2);
			dec_name_nr++;
		}
		//MatAssemblyBegin(R, MAT_FINAL_ASSEMBLY);
		//MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY);
	}
	MatAssemblyBegin(R, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY);
	//MatView(R,PETSC_VIEWER_STDOUT_WORLD);

	//Computation of nominator in pade approximant follows
	MatZeroEntries(Hlp);
	MatZeroEntries(Hlp2);
	MatZeroEntries(N);
	MatAssemblyBegin(Hlp,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Hlp,MAT_FINAL_ASSEMBLY);
	MatShift(Hlp, 1.0); //identity matrix
	//MatCopy(R, Hlp, DIFFERENT_NONZERO_PATTERN);
	MatAssemblyBegin(N, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(N, MAT_FINAL_ASSEMBLY);
	for(rows = 0; rows <= nom_pol_deg; rows++)
	{
		koef_hlp = (PetscScalar) (faktorial(nom_pol_deg + den_pol_deg - rows) * faktorial(nom_pol_deg)) / (faktorial(nom_pol_deg + den_pol_deg) * faktorial(rows) * faktorial(nom_pol_deg - rows));
		xprintf(Msg,"\n koeficient has a value %f\n", koef_hlp);
		//if(rows > 0)
		MatAXPY(N, koef_hlp, Hlp, DIFFERENT_NONZERO_PATTERN);
		//else MatAXPY(N, koef_hlp, Hlp, SUBSET_NONZERO_PATTERN);
		MatMatMult(R,Hlp,MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Hlp2);
		MatZeroEntries(Hlp);
		MatCopy(Hlp2, Hlp, DIFFERENT_NONZERO_PATTERN);
		MatZeroEntries(Hlp2);
	}
	//MatView(N,PETSC_VIEWER_STDOUT_WORLD);

	//Computation of denominator in pade approximant follows
	MatZeroEntries(Hlp);
	MatZeroEntries(Hlp2);
	MatShift(Hlp, 1.0); //identity matrix
	//MatCopy(R, Hlp, DIFFERENT_NONZERO_PATTERN);
	MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(D, MAT_FINAL_ASSEMBLY);
	for(rows = 0; rows <= den_pol_deg; rows++)
	{
		koef_hlp = (-1.0) * faktorial(nom_pol_deg + den_pol_deg - rows) * faktorial(den_pol_deg) / (faktorial(nom_pol_deg + den_pol_deg) * faktorial(rows) * faktorial(den_pol_deg - rows));
		//if(rows > 0)
		MatAXPY(D, koef_hlp, Hlp, DIFFERENT_NONZERO_PATTERN);
		//else MatAXPY(D, koef_hlp, Hlp, SUBSET_NONZERO_PATTERN);
		MatMatMult(R,Hlp,MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Hlp2);
		MatZeroEntries(Hlp);
		MatCopy(Hlp2, Hlp, DIFFERENT_NONZERO_PATTERN);
		MatZeroEntries(Hlp2);
	}
	//MatView(D,PETSC_VIEWER_STDOUT_WORLD);

	//create the matrix B
	MatCreate(PETSC_COMM_WORLD, &B); //MatCreateSeqDense
	MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, nr_of_species, nr_of_species); //nr_of_species should be probably multiplied by 2 (which is the value of m), but I do not know why
	MatSetType(B, MATSEQAIJ);
	MatZeroEntries(B);

	//pade approximant
	MatZeroEntries(B);//MatZeroEntries(Hlp);
	MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
	MatShift(B,1.0); //MatShift(Hlp,1.0); //It should result in identity matrix.
	MatZeroEntries(Hlp2);
	MatFactorInfoInitialize(&matfact);
	MatGetOrdering(D, MATORDERINGNATURAL, &rperm, &cperm);
	MatLUFactor(D, rperm, cperm, &matfact);
	MatMatSolve(D, B, Hlp2);//MatMatSolve(D, Hlp, Hlp2); //D^{-1} into Hlp2
	MatZeroEntries(Hlp);
	MatMatMult(Hlp2, N , MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Hlp);
	MatCreateTranspose(Hlp, &pade);

	MatAssemblyBegin(pade, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(pade, MAT_FINAL_ASSEMBLY);

	//pade assembled to reaction_matrix
	//MatGetArray(pade, &Hlp_mat);
	for(rows = 0; rows < nr_of_species; rows++)
	{
		for(cols = 0; cols < nr_of_species; cols++)
		{
			//MatZeroEntries(Hlp_mat);
			//reaction_matrix[rows][cols] = Hlp_mat[rows * nr_of_species + cols]; //values from pade approximant are assembled to reaction_matrix, here
			MatGetValues(pade, 1, &rows, 1, &cols, Hlp_mat); //&Hlp_mat[nr_of_species*rows + cols]);
			reaction_matrix[rows][cols] = (double) (*Hlp_mat);
			/*if(cols == (nr_of_species - 1))
			{
				PetscPrintf(PETSC_COMM_WORLD,"%f\n", Hlp_mat[nr_of_species*rows + cols]);
			}else{
				PetscPrintf(PETSC_COMM_WORLD,"%f", Hlp_mat[nr_of_species*rows + cols]);
			}*/
		}
	}

	PetscPrintf(PETSC_COMM_WORLD,"pade matrix looks as follows:\n");
	MatView(pade,PETSC_VIEWER_STDOUT_WORLD);

	print_reaction_matrix(); //for visual control of equality of reaction_matrix in comparison with pade aproximant

	MatDestroy(&B);
	MatDestroy(&D);
	MatDestroy(&N);
	MatDestroy(&R);
	MatDestroy(&pade);
	MatDestroy(&Hlp);
	MatDestroy(&Hlp2);

	return reaction_matrix;
}

double **Linear_reaction::compute_reaction(double **concentrations, int loc_el) //multiplication of concentrations array by reaction matrix
{


    int cols, rows, both;

	if((nr_of_decays > 0) || (nr_of_FoR > 0)){
		for(cols = 0; cols < nr_of_species; cols++){
		prev_conc[cols] = concentrations[cols][loc_el];
		//xprintf(Msg,"\n%d. of %d substances concentration is %f\n", cols,nr_of_species, concentrations[cols][loc_el]); //prev_conc[cols]); //commented to speed the computation up
		concentrations[cols][loc_el] = 0.0;
		}
        for(rows = 0; rows <nr_of_species; rows++){
            for(cols = 0; cols <nr_of_species; cols++){
                concentrations[rows][loc_el] += prev_conc[cols] * reaction_matrix[cols][rows];
            }
            //xprintf(Msg,"\n%d. of %d substances concentration after reaction is %f\n", rows,nr_of_species, concentrations[rows][loc_el]); //commented to speed the computation up
        }
	}
	return concentrations;
}

double *Linear_reaction::set_half_lives(char *section)
{
	char  buffer[1024];
	char *pom_buf;
	int i,j;
	const char *separators = " ,\t";

	if(half_lives != NULL){
			free(half_lives);
			half_lives = NULL;
	}
	if(half_lives == NULL){
		//xprintf(Msg,"\nAllocation is permited, nr of isotopes %d", nr_of_isotopes);
		half_lives = (double *)xmalloc((nr_of_isotopes - 1) * sizeof(double));
	}
	 strcpy(buffer,OptGetStr(section,"Half_lives",NULL));
	 pom_buf = strtok( buffer, separators );
	 for (j=0; j< (nr_of_isotopes-1); j++){
		if ( pom_buf == NULL )
		{
			xprintf(Msg,"\nHalf-life of %d-th isotope is missing.", j+1);
		}
	    half_lives[j] = atof(pom_buf);
	    xprintf(Msg,"\n %d-th isotopes half-live is %f",j,half_lives[j]);
	    pom_buf = strtok( NULL, separators );
	 }
	 if ( pom_buf != NULL )
	 {
	    xprintf(Msg,"\nMore parameters then (isotopes -1) has been given. %d", 0);
	 }
	return half_lives;
}

void Linear_reaction::print_half_lives(int nr_of_substances)
{
	int i;

	if(half_lives == NULL)
	{
		xprintf(Msg,"\nHalf-lives are not defined.");
	}else{
		xprintf(Msg,"\nHalf-lives are defined as:");
		for(i=0; i < (nr_of_substances - 1) ; i++)
		{
			if(i < (nr_of_substances  - 2)) //cout << " " << half_lives[i] <<",";
				xprintf(Msg," %f", half_lives[i]);
			if(i == (nr_of_substances  - 2)) //cout << " " << half_lives[i] <<"\n";
				xprintf(Msg," %f\n", this->half_lives[i]);
		}
	}
	return;
}

int *Linear_reaction::set_indeces(char *section, int nr_of_substances)
{
	char  buffer[1024];
	char *pom_buf;
	int i,j;
	const char *separators = " ,\t";

	if(substance_ids != NULL){
		free(substance_ids);
		substance_ids = NULL;
	}
	if(substance_ids == NULL){
		substance_ids = (int *)xmalloc(nr_of_substances*sizeof(int));
	}

	strcpy(buffer,OptGetStr(section,"Substance_ids",NULL));
	pom_buf = strtok( buffer, separators );
	for (j=0; j< nr_of_substances; j++)
	{
	  if ( pom_buf == NULL )
	  {
	    xprintf(Msg,"\nIndex for %d-th substance in %s is missing.", j+1, section);
	  }
	    substance_ids[j] = atoi(pom_buf);
	    pom_buf = strtok( NULL, separators );
	 }
	 if ( pom_buf != NULL )
	 {
	    xprintf(Msg,"\nMore parameters then substances has been given in %s.", section);
	 }

	 return substance_ids;
}

void Linear_reaction::print_indeces(int nr_of_substances)
{
	int i;

	if(substance_ids == NULL)
	{
		xprintf(Msg,"\nReaction/decay has not been defined.");
	}else{
		xprintf(Msg,"\nOrder of substences is defined by %d indeces:", nr_of_isotopes);
		for(i = 0; i < nr_of_substances ; i++)
		{
			if(i < (nr_of_substances  - 1)) xprintf(Msg," %d,",substance_ids[i]);
			if(i == (nr_of_substances  - 1)) xprintf(Msg," %d\n",substance_ids[i]);
		}
	}
	return;
}

void Linear_reaction::print_reaction_matrix(void)
{
	int cols,rows;

	if(reaction_matrix != NULL){
		xprintf(Msg,"\ntime_step %f,Reaction matrix looks as follows:\n",time_step);
		for(rows = 0; rows < nr_of_species; rows++){
			for(cols = 0; cols < nr_of_species; cols++){
				if(cols == (nr_of_species - 1)){
					xprintf(Msg,"%f\n",reaction_matrix[rows][cols]);
				}else{
					xprintf(Msg,"%f\t",reaction_matrix[rows][cols]);
				}
			}
		}
	}else{
		xprintf(Msg,"\nReaction matrix needs to be allocated.\n");
	}
	return;
}

void Linear_reaction::set_bifurcation(char *section, int dec_nr)
{
	char  buffer[1024];
	char *pom_buf;
	int j;
	const char *separators = " ,\t";
	double control_sum = 0.0;

	if(bifurcation_on == true)
	{
		bifurcation[dec_nr].resize(nr_of_isotopes - 1);
		strcpy(buffer,OptGetStr(section,"Bifurcation",NULL));
		if(buffer == NULL) return;
		pom_buf = strtok( buffer, separators );
		for (j=0; j< (nr_of_isotopes - 1); j++)
		{
			if ( pom_buf == NULL )
			{
				xprintf(Msg,"\nBifurcation parameter of %d-th isotope is missing.", j+1);
			}
	    	bifurcation[dec_nr][j] = atof(pom_buf);
	    	xprintf(Msg,"\n %d-th isotopes bifurcation percentage is %f",j,bifurcation[dec_nr][j]);
	    	pom_buf = strtok( NULL, separators );
	    	if(j > 0)control_sum += bifurcation[dec_nr][j];
	 	 }
	}else{
		bifurcation[dec_nr].resize(1);
		bifurcation[dec_nr][0]= 1.0;
	}
	if(control_sum != 1.0) xprintf(Msg,"\nSum of bifurcation parameters should be 1.0 but it is %f, because of mass conservation law.\n", control_sum);
	if( pom_buf != NULL )
	 {
	    xprintf(Msg,"\nMore parameters then (isotopes -1) has been given. %d", 0);
	 }
	return;
}

void Linear_reaction::set_kinetic_constants(char *section, int react_nr)
{
	char  buffer[1024];
	char *pom_buf;
	//int j;
	const char *separators = " ,\t";

	kinetic_constant.resize(nr_of_FoR);
	strcpy(buffer,OptGetStr(section,"Kinetic_constant",NULL));
	if(buffer == NULL) return;
	pom_buf = strtok( buffer, separators );
	//for (j=0; j< (nr_of_FoR); j++){
		if ( pom_buf == NULL )
		{
			xprintf(Msg,"\nKinetic constant belonging to %d-th reactions is missing.", react_nr+1);
		}
    	kinetic_constant[react_nr] = atof(pom_buf);
    	xprintf(Msg,"\nKinetic constant for %d-th reaction is %f",react_nr,kinetic_constant[react_nr]);
    	pom_buf = strtok( NULL, separators );
    	half_lives[react_nr] = log(2) / kinetic_constant[react_nr];
 	 //}
    return;
}

void Linear_reaction::compute_one_step(void)
{
    if (reaction_matrix == NULL)   return;

    START_TIMER("decay_step");
	 //for (int loc_el = 0; loc_el < distribution->lsize(distribution->myp()); loc_el++)
	for (int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
	 {
	 	this->compute_reaction(concentration_matrix[MOBILE], loc_el);
	    if (dual_porosity_on == true) {
	     this->compute_reaction(concentration_matrix[IMMOBILE], loc_el);
	    }

	 }
    END_TIMER("decay_step");
	 return;
}

void Linear_reaction::set_nr_of_species(int n_substances)
{
	this->nr_of_species = n_substances;
	return;
}

void Linear_reaction::set_nr_of_elements(int nrOfElements)
{
	this->nr_of_elements = nrOfElements;
	return;
}

void Linear_reaction::set_concentration_matrix(double ***ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc)
{
	concentration_matrix = ConcentrationMatrix;
	distribution = conc_distr;
	return;
}

void Linear_reaction::set_time_step(double new_timestep){
	time_step = new_timestep;
	if((nr_of_decays > 0) || (nr_of_FoR > 0)){
		release_reaction_matrix();
		allocate_reaction_matrix();
		if(matrix_exp_on == false)
		{
			modify_reaction_matrix_repeatedly();
		}else{
			modify_reaction_matrix_using_pade();
		}
	}
	return;
}
int Linear_reaction::get_nr_of_decays(void){return nr_of_decays;} // two simple inlinefunction returning private variables
int Linear_reaction::get_nr_of_FoR(void){return nr_of_FoR;}
void Linear_reaction::set_mesh_(Mesh *mesh_in){mesh = mesh_in; return;}

void Linear_reaction::set_dual_porosity()
{
	this->dual_porosity_on = OptGetBool("Transport", "Dual_porosity", "no");
	return;
}

void Linear_reaction::release_reaction_matrix(void)
{
	int i;
	if(reaction_matrix != NULL)
	{
		for(i = 0; i < nr_of_isotopes; i++)
		{
			if(reaction_matrix[i] != NULL)
			{
				free(reaction_matrix[i]);
				reaction_matrix[i] = NULL;
			}
		}
		free(reaction_matrix);
		reaction_matrix = NULL;
	}
}

double Linear_reaction::get_time_step(void)
{
	return time_step;
}

void Linear_reaction::set_nr_of_isotopes(int Nr_of_isotopes)
{
	nr_of_isotopes = Nr_of_isotopes;
	return;
}

void Linear_reaction::set_nr_of_decays(void)
{
	nr_of_decays = OptGetInt("Reaction_module","Nr_of_decay_chains","0");
	return;
}

void Linear_reaction::set_nr_of_FoR(void)
{
	nr_of_FoR = OptGetInt("Reaction_module","Nr_of_FoR","0");
	return;
}

int Linear_reaction::faktorial(int k)
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
