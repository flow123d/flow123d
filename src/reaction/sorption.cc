#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>

#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "reaction/sorption.hh"
#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "la/distribution.hh"
#include "mesh/mesh.h"

using namespace Input::Type;

Record Sorption::input_type_isotherm // input_type_one_decay_substep
	= Record("Isotherm", "Equation for reading information about limmited solubility affected sorption.")
	.declare_key("specie", String(), Default::obligatory(),
				"Identifier of a sorbing isotope.")
	.declare_key("molar_mass", Double(), Default("1.0"),
				"Molar mass.")
    .declare_key("region", Integer(), Default::optional(), // not considered at first
                "Arrea where considered sorption appears.")
    .declare_key("type", String(), Default::obligatory(), // this item is obsolete may be
                "Declaration of an isotherm type either lang or lin.")
    .declare_key("direction", Double(), Default("0.0"), // parameter of a linear isotherm
				"Defines a slope of a linear isotherm.")
	.declare_key("omega", Double(), Default("1.0"), // multiplicative parameter of a langmuir isotherm
							"Defines a multiplicative parameter omega of a langmuir isotherm c_s = omega * (alpha*c_a)/(1- alpha*c_a).")
	.declare_key("alpha", Double(), Default("0.0"), // parameter of a langmuir isotherm
				"Defines a parameter alpha of a langmuir isotherm c_s = omega * (alpha*c_a)/(1- alpha*c_a).")
	.declare_key("solvable", Double(), Default("1.0"),   // concentration limit for a solubility of the specie under concideration
				"Solubility limit.");


Record Sorption::input_type
	= Record("Sorptions", "Information about all the limited solubility affected sorptions.")
	.derive_from( Reaction::input_type )
    .declare_key("sorptions", Array( Sorption::input_type_isotherm ), Default::obligatory(),
                "Description of particular sorption cases under consideration.");


using namespace std;

Sorption::Sorption(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)//(double timeStep, Mesh * mesh, int nrOfSpecies, bool dualPorosity, Input::Record in_rec) //(double timestep, int nrOfElements, double ***ConvectionMatrix)
      : Reaction(init_mesh, in_rec, names)
{
	prepare_inputs(in_rec);
	//determine_crossection()
}

/*Sorption::~Sorption()
{
	//int i, rows, n_subst;
	//
	//release_reaction_matrix();
}*/

double **Sorption::compute_reaction(double **concentrations, int loc_el) // Sorptions are realized just for one element.
{
    int cols, rows;

    //  Measurements [c_a,c_s] will be rotated
    //  Rotated measurements must be projected on rotated isotherm
    //  Projections need to be transformed back to original CS

    /*if (reaction_matrix == NULL) return concentrations;

	for(cols = 0; cols < n_substances(); cols++){
		prev_conc[cols] = concentrations[cols][loc_el];
		concentrations[cols][loc_el] = 0.0;
	}

	for(rows = 0; rows < n_substances(); rows++){
        for(cols = 0; cols < n_substances(); cols++){
            concentrations[rows][loc_el] += prev_conc[cols] * reaction_matrix[cols][rows];
        }
    }*/

	return concentrations;
}

void Sorption::compute_one_step(void) // Computes sorption simulation over all the elements.
{
    //DBGMSG("decay step\n");
    //if (reaction_matrix == NULL)   return;

    /*START_TIMER("sorption_step");
	for (int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
	 {
	 	this->compute_reaction(concentration_matrix[MOBILE], loc_el);
	    if (dual_porosity_on == true) {
	     this->compute_reaction(concentration_matrix[IMMOBILE], loc_el);
	    }

	 }
    END_TIMER("sorption_step");*/
	 return;
}


void Sorption::print_sorption_parameters(int nr_of_substances) {
    int i;

    xprintf(Msg, "\nSorption parameters are defined as:");
    for (i = 0; i < (nr_of_substances - 1); i++) {
        if (i < (nr_of_substances - 2)) ; //cout << " " << half_lives[i] <<",";
        	// describing table header
        	// name(specie)
        	// molar mass
        	// isotherm type
        	// region
        	// parameter(s)
        	// solubility limit
            //xprintf(Msg, " %f", half_lives[i]);
        if (i == (nr_of_substances - 2)) ; //cout << " " << half_lives[i] <<"\n";
            // xprintf(Msg, " %f\n", this->half_lives[i]);
    }
}

// TODO: check duplicity of parents
//       raise warning if sum of ratios is not one
void Sorption::prepare_inputs(Input::Record in_rec)
{
    unsigned int idx;

	Input::Array sorption_array = in_rec.val<Input::Array>("sorptions");

	substance_ids.resize( sorption_array.size() );
	//half_lives.resize( sorption_array.size() );
	//bifurcation.resize( sorption_array.size() );

	int i_sorption=0;
	for (Input::Iterator<Input::Record> sorp_it = sorption_array.begin<Input::Record>(); sorp_it != sorption_array.end(); ++sorp_it, ++i_sorption)
	{
		//region determining part
		/*Input::Iterator<double> it_hl = sorp_it->find<double>("half_life");
		if (it_hl) {
		   half_lives[i_sorption] = *it_hl;
		} else {
		   it_hl = sorp_it->find<double>("kinetic");
		   if (it_hl) {
			   half_lives[i_sorption] = log(2)/(*it_hl);
		   } else {
		    xprintf(UsrErr, "Missing half-life or kinetic in the %d-th reaction.\n", i_sorption);
		  }
		}*/

		//isotherm type determining part
		/*string parent_name = sorp_it->val<string>("type");
		Input::Array product_array = sorp_it->val<Input::Array>("products");
		Input::Array ratio_array = sorp_it->val<Input::Array>("branch_ratios"); // has default value [ 1.0 ]*/

		// linear isotherm direction (slope)
		/*if (product_array.size() > 0)   substance_ids[i_sorption].resize( product_array.size()+1 );
		else			xprintf(UsrErr,"Empty array of products in the %d-th reaction.\n", i_sorption);*/


		// multiplicative coefficient omega for Langmuir isotherm, c_s = omega*(alpha*c_a)/(1 - alpha*c_a)
		/*idx = find_subst_name(parent_name);
		if (idx < n_substances())	substance_ids[i_sorption][0] = idx;
		else                		xprintf(UsrErr,"Wrong name of parent substance in the %d-th reaction.\n", i_sorption);*/

		// alpha parameter for Langmuir isotherm, c_s = omega*(alpha*c_a)/(1 - alpha*c_a)
		/*unsigned int i_product = 1;
		for(Input::Iterator<string> product_it = product_array.begin<string>(); product_it != product_array.end(); ++product_it, i_product++)
		{
			idx = find_subst_name(*product_it);
			if (idx < n_substances())   substance_ids[i_sorption][i_product] = idx;
			else                    	xprintf(Msg,"Wrong name of %d-th product in the %d-th reaction.\n", i_product-1 , i_sorption);
		}*/

		//Critical concentrations solvable in water.
        /*if (ratio_array.size() == product_array.size() )   ratio_array.copy_to( bifurcation[i_sorption] );
        else            xprintf(UsrErr,"Number of branches %d has to match number of products %d in the %d-th reaction.\n",
                                       ratio_array.size(), product_array.size(), i_sorption);*/

		// Molar mass of particular adsorbent.

	}
}

void Sorption::determine_crossections(int k_points)
{
	;
}

void rotate_points(double angle, double **points)
{
	;
}

void interpolate_datapoints(void)
{
	;
}
