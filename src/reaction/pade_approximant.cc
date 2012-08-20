#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "system/system.hh"
#include "materials.hh"
#include "transport/transport.h"
//#include "system/par_distribution.hh"
#include "la/distribution.hh"
#include "mesh/mesh.h"

//class Distribution;

Input::Type::Record & Pade_approximant::get_one_decay_substep()
{
	using namespace Input::Type;
	static Record rec("Substep", "Equation for reading information about radioactive decays.");

	if(!rec.is_finished()){
		rec.declare_key("parent", String(), Default::obligatory(),
				"Identifier of an isotope.");
        rec.declare_key("half_life", Double(), Default::optional(),
                "Half life of the parent substance.");
        rec.declare_key("kinetic", Double(), Default::optional(),
                "Kinetic constants describing first order reactions.");
		rec.declare_key("products", Array(String()), Default::obligatory(),
				"Identifies isotopes which decays parental atom to.");
		rec.declare_key("branching_ratios", Array(Double()), Default::optional(),
				"Decay chain branching percentage.");
		rec.finish();
	}
	return rec;
}

Input::Type::Record & Pade_approximant::get_input_type()
{
	using namespace Input::Type;
	static Record rec("Approximant", "Abstract record with an information about pade approximant parameters.");

	if (!rec.is_finished()) {
	    rec.derive_from( Reaction::get_input_type() );
		/*rec.declare_key("substances", Array(String()), Default::obligatory(),
								"Names of transported isotopes.");

		rec.declare_key("half_lives", Array(Double()), Default::optional(),
				"Half lives of transported isotopes.");
		rec.declare_key("decays", Array(Decay()), Default::optional(),
				"Description of particular decay chain substeps.");
		rec.declare_key("kinetic_constants", Array(Double()), Default::optional(),
				"Kinetic constants describing first order reactions.");
		rec.declare_key("kinetics", Array(Kinetics()), Default::optional(),
				"Description of particular first order reactions.");
		rec.declare_key("matrix_exp_on", Bool(), Default("false"),
				"Enables to use Pade approximant of matrix exponential.");*/
		rec.declare_key("nom_pol_deg", Integer(), Default("0"),
				"Polynomial degree of the nominator of Pade approximant.");
		rec.declare_key("den_pol_deg", Integer(), Default("0"),
				"Polynomial degree of the nominator of Pade approximant");
		rec.finish();

		//Linear_reaction::get_input_type();

		//rec.no_more_descendants();
	}
	return rec;
}

using namespace std;

Pade_approximant::Pade_approximant(TimeMarks &marks, Mesh &init_mesh, MaterialDatabase &material_database, Input::Record in_rec) //(double timeStep, Mesh * mesh, int nrOfSpecies, bool dualPorosity) //(double timestep, int nrOfElements, double ***ConvectionMatrix)
			:Reaction(marks, init_mesh, material_database, in_rec)//, Reaction_matrix(NULL)
{
	nom_pol_deg = in_rec.val<int>("nom_pol_deg");
	den_pol_deg = in_rec.val<int>("den_pol_deg");
	if((nom_pol_deg + den_pol_deg) < 0){
		cout << "You did not specify Pade approximant required polynomial degrees." << endl;
		//This occasion should cause an error.
		//break;
	}
	cout << "Pade_approximant constructor is running." << endl;
}

Pade_approximant::~Pade_approximant()
{
	int i, rows, n_subst;

	if(half_lives != NULL){
		free(half_lives);
		half_lives = NULL;
	}

	if(substance_ids != NULL){
		free(substance_ids);
		substance_ids = NULL;
	}

	release_reaction_matrix();

	cout << "Pade approximant destructor is running."  << endl;
}

void Pade_approximant::set_time_step(double new_timestep){
	time_step = new_timestep;
	release_reaction_matrix();
	this->allocate_reaction_matrix();
	this->modify_reaction_matrix_repeatedly();
	return;
}

double **Pade_approximant::allocate_reaction_matrix(void) //reaction matrix initialization
{
	int index, rows, cols, dec_nr, prev_index;
	char dec_name[30];

	cout << "We are going to allocate reaction matrix" << endl;
	if(reaction_matrix == NULL)
	{
		; //Do nothing!!!
	}else{
		release_reaction_matrix();
	}
	reaction_matrix = (double **)xmalloc(nr_of_species * sizeof(double*));//allocation section
	for(rows = 0; rows < nr_of_species; rows++){
		reaction_matrix[rows] = (double *)xmalloc(nr_of_species * sizeof(double));
	}
	for(rows = 0; rows < nr_of_species;rows++){
	 for(cols = 0; cols < nr_of_species; cols++){
		 /*if(rows == cols){
			 ;//if((nom_pol_deg + den_pol_deg) == 0) reaction_matrix[rows][cols] = 1.0; //this row is different in comparison to the case of Linear_reaction
		 }else*/{
			reaction_matrix[rows][cols] = 0.0;
		 }
	 }
	}
	//print_reaction_matrix();
	return reaction_matrix;
}

void Pade_approximant::modify_reaction_matrix(void) //prepare the matrix, which describes reactions
{
	int rows,cols;
	double rel_step, prev_rel_step;
	PetscScalar Hlp_kin, index, prev_index;

	if(Reaction_matrix == NULL){
		xprintf(Msg,"\nPointer to the Reaction matrix is NULL.\n");
		return;
	}
	if(nr_of_decays > 0){
		for(cols = 0; cols < nr_of_isotopes; cols++){
			index = substance_ids[cols] - 1; // because indecees in input file run from one whereas indeces in C++ run from ZERO
			if(cols > 0){
				Hlp_kin = (PetscScalar) time_step*log(2)/half_lives[cols-1];
				MatSetValue(Reaction_matrix, prev_index, prev_index, ((-1.0) * Hlp_kin), INSERT_VALUES);
				MatSetValue(Reaction_matrix, index, prev_index, Hlp_kin, INSERT_VALUES);
			}
			prev_index = index;
		}
	}
	return;
}

double **Pade_approximant::modify_reaction_matrix(int dec_nr) //prepare the matrix, which describes reactions, takes bifurcation in acount
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
	MatSetValue(Reaction_matrix,first_index,first_index,Hlp_kin,ADD_VALUES);
	for(cols = 0; cols < nr_of_isotopes; cols++){
		index = substance_ids[cols] - 1; // because indecees in input file run from one whereas indeces in C++ run from ZERO
		if(cols > 0)
		{
			bif_id = cols -1;
			Hlp_kin = time_step*log(2)/half_lives[cols-1] * bifurcation[dec_nr][bif_id];
			MatSetValue(Reaction_matrix, index, first_index, Hlp_kin, INSERT_VALUES);
		}
	}
	return reaction_matrix;
}

double **Pade_approximant::modify_reaction_matrix_repeatedly(void)
{
	Mat Denominator;
	Mat Nominator;
	//Mat Reaction_matrix;
	Mat Pade_approximant;
	MatFactorInfo matfact;
	PC Precond;
	IS rperm, cperm;
	Vec tmp1; //contains the information about concentrations of all the species in one particular element
	Vec tmp2; //the same as tmp1
	PetscInt n, m = 2;
	PetscScalar nominator_coef[nom_pol_deg];
	PetscScalar denominator_coef[den_pol_deg];
	PetscScalar Hlp_mat[1];
	PetscScalar *Array_hlp;
	const PetscScalar *Reaction_matrix_row;
	char dec_name[30];
	int rows, cols, dec_nr, dec_name_nr = 1, index, prev_index, i, j;

	//create the matrix Reaction_matrix
	MatCreate(PETSC_COMM_SELF, &Reaction_matrix);
	MatSetSizes(Reaction_matrix, PETSC_DECIDE, PETSC_DECIDE, nr_of_species, nr_of_species); //should be probably multiplied by 2 (which is the value of m)
	MatSetType(Reaction_matrix, MATAIJ);
	MatAssemblyBegin(Reaction_matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Reaction_matrix, MAT_FINAL_ASSEMBLY);

	//create the matrix N
	MatDuplicate(Reaction_matrix, MAT_COPY_VALUES, &Nominator);

	//create the matrix D
	MatDuplicate(Reaction_matrix, MAT_COPY_VALUES, &Denominator);

	//create the matrix pade
	MatDuplicate(Reaction_matrix, MAT_COPY_VALUES, &Pade_approximant);
	MatAssemblyBegin(Pade_approximant, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Pade_approximant, MAT_FINAL_ASSEMBLY);

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
			nr_of_isotopes = OptGetInt(dec_name,"Nr_of_isotopes","0"); //sizeof(substances) or sum(sign(half_lives + 1))
			set_half_lives(dec_name);
			set_indeces(dec_name, nr_of_isotopes);
			print_indeces(nr_of_isotopes); //just a control
			print_half_lives(nr_of_isotopes); //just a control
			bifurcation_on = OptGetBool(dec_name,"Bifurcation_on","no"); //
			if(bifurcation_on == true){
				set_bifurcation(dec_name, dec_nr);
				modify_reaction_matrix(dec_nr);
			}else{
				if(&Reaction_matrix != NULL)
				{
					modify_reaction_matrix();
					xprintf(Msg,"Reaction matrix R is has been allocated. The addres is %d.\n", &Reaction_matrix);
				}else{
					xprintf(Msg,"Reaction matrix R has not been allocated.\n"); //cout << "Reaction matrix R is has not been allocated." << endl;
				}
			}
			dec_name_nr++;
		}
	}
	/*if(nr_of_FoR > 0){
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
			modify_reaction_matrix(2);
			dec_name_nr++;
		}
	}*/
	MatAssemblyBegin(Reaction_matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Reaction_matrix, MAT_FINAL_ASSEMBLY);
	MatView(Reaction_matrix,PETSC_VIEWER_STDOUT_SELF);

	//Computation of nominator in pade approximant follows
	MatZeroEntries(Nominator);
	MatAssemblyBegin(Nominator, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Nominator, MAT_FINAL_ASSEMBLY);
	for(j = nom_pol_deg; j >= 0; j--)
	{
		nominator_coef[j] = (PetscScalar) (faktorial(nom_pol_deg + den_pol_deg - j) * faktorial(nom_pol_deg)) / (faktorial(nom_pol_deg + den_pol_deg) * faktorial(j) * faktorial(nom_pol_deg - j));
	}
	evaluate_matrix_polynomial(&Nominator, &Reaction_matrix, nominator_coef);
	MatView(Nominator,PETSC_VIEWER_STDOUT_WORLD);

	//Computation of denominator in pade approximant follows
	MatZeroEntries(Denominator);
	MatAssemblyBegin(Denominator, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Denominator, MAT_FINAL_ASSEMBLY);
	for(i = den_pol_deg; i >= 0; i--)
	{
		denominator_coef[i] = (PetscScalar) pow(-1.0,i) * faktorial(nom_pol_deg + den_pol_deg - i) * faktorial(den_pol_deg) / (faktorial(nom_pol_deg + den_pol_deg) * faktorial(i) * faktorial(den_pol_deg - i));
	}
	evaluate_matrix_polynomial(&Denominator, &Reaction_matrix, denominator_coef);
	MatView(Denominator, PETSC_VIEWER_STDOUT_WORLD);

	PCCreate(PETSC_COMM_WORLD, &Precond);
	PCSetType(Precond, PCLU);
	PCSetOperators(Precond, Denominator, Denominator, DIFFERENT_NONZERO_PATTERN);
	//PCFactorSetMatOrderingType(Precond, MATORDERINGNATURAL);
	PCFactorSetMatOrderingType(Precond, MATORDERINGRCM);
	PCSetUp(Precond);

	VecCreate(PETSC_COMM_WORLD, &tmp1);
	VecSetSizes(tmp1, PETSC_DECIDE, nr_of_species);
	VecSetFromOptions(tmp1);
	VecDuplicate(tmp1, &tmp2);

	for(rows = 0; rows < nr_of_species; rows++){
		MatGetColumnVector(Nominator, tmp1, rows);
		//VecView(tmp1, PETSC_VIEWER_STDOUT_SELF);
		PCApply(Precond, tmp1, tmp2);
		PCView(Precond, PETSC_VIEWER_STDOUT_WORLD);
		//VecView(tmp2, PETSC_VIEWER_STDOUT_SELF);
		VecGetArray(tmp2, &Array_hlp);
		for(cols = 0; cols < nr_of_species; cols++)
		{
			MatSetValue(Pade_approximant, rows, cols, Array_hlp[cols], ADD_VALUES);
		}
	}
	MatAssemblyBegin(Pade_approximant, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Pade_approximant, MAT_FINAL_ASSEMBLY);

	//pade assembled to reaction_matrix
	for(rows = 0; rows < nr_of_species; rows++)
	{
		for(cols = 0; cols < nr_of_species; cols++)
		{
			MatGetValues(Pade_approximant, 1, &rows, 1, &cols, Hlp_mat); //&Hlp_mat[nr_of_species*rows + cols]);
			reaction_matrix[rows][cols] = (double) (Hlp_mat[0]);
		}
	}

	/*PetscPrintf(PETSC_COMM_SELF,"pade matrix looks as follows:\n");
	MatView(Pade_approximant,PETSC_VIEWER_STDOUT_WORLD);*/

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
	MatSetSizes(Identity, PETSC_DECIDE, PETSC_DECIDE, nr_of_species, nr_of_species); //should be probably multiplied by 2 (which is the value of m)
	MatSetType(Identity, MATAIJ);
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


    int cols, rows, both;

	if(nr_of_decays > 0){
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

double *Pade_approximant::set_half_lives(char *section)
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

void Pade_approximant::print_half_lives(int nr_of_substances)
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

int *Pade_approximant::set_indeces(char *section, int nr_of_substances)
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

void Pade_approximant::print_indeces(int nr_of_substances)
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

void Pade_approximant::print_reaction_matrix(void)
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

void Pade_approximant::set_bifurcation(char *section, int dec_nr)
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

/*void Pade_approximant::set_kinetic_constants(char *section, int react_nr)
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
}*/

void Pade_approximant::set_time_step(void)
{
	time_step = OptGetDbl("Global","Save_step","1.0");
	release_reaction_matrix();
	allocate_reaction_matrix();
	modify_reaction_matrix_repeatedly();
	return;
}

void Pade_approximant::compute_one_step(void)
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

void Pade_approximant::release_reaction_matrix(void)
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

void Pade_approximant::set_nr_of_isotopes(int Nr_of_isotopes)
{
	nr_of_isotopes = Nr_of_isotopes;
	return;
}
