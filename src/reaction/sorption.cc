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

namespace it=Input::Type;

Sorption::EqData::EqData(const std::string &name)
: EqDataBase(name)
{
    ADD_FIELD(nr_of_points, "Number of crossections allong isotherm", Input::Type::Default("10"));
    ADD_FIELD(region_ident, "Rock matrix area identifier.", Input::Type::Default("0"));
    ADD_FIELD(rock_density, "Rock matrix density.", Input::Type::Default("0.0"));

    ADD_FIELD(sorption_type,"Considered adsorption is described by selected isotherm.", it::Default("none") );
              sorption_type.set_selection(&sorption_type_selection);

    ADD_FIELD(slopes,"Directions of linear isotherm.",it::Default("0.0"));
    std::vector<FieldEnum> list; list.push_back(none); list.push_back(Langmuir);
    slopes.disable_where(& sorption_type, list );

    ADD_FIELD(omegas,"Langmuir isotherm multiplication parameters in c_s = omega * (alpha*c_a)/(1- alpha*c_a).");
    list.clear(); list.push_back(none); list.push_back(Linear);
    omegas.disable_where(& sorption_type, list );

    ADD_FIELD(alphas,"Langmuir isotherm alpha parameters in c_s = omega * (alpha*c_a)/(1- alpha*c_a).");
    list.clear(); list.push_back(none); list.push_back(Linear);
    alphas.disable_where(& sorption_type, list );
}

/*RegionSet Sorption::EqData::read_bulk_list_item(Input::Record rec) {
    RegionSet domain=EqDataBase::read_bulk_list_item(rec);
    Input::AbstractRecord field_a_rec;
    if (rec.opt_val("init_piezo_head", field_a_rec)) {
                init_pressure.set_field(domain, boost::make_shared< FieldAddPotential<3, FieldValue<3>::Scalar > >( this->gravity_, field_a_rec) );
    }
    return domain;
}*/

using namespace Input::Type;

Record Sorption::input_type_isotherm // region independent parameters
	= Record("Isotherm", "Equation for reading information about limmited solubility affected sorption.")
	.declare_key("specie", String(), Default::obligatory(),
				"Identifier of a sorbing isotope.")
	.declare_key("molar_mass", Double(), Default("1.0"),
				"Molar mass.")
	.declare_key("solvable", Double(), Default("1.0"),   // concentration limit for a solubility of the specie under concideration
				"Solubility limit.")
	.declare_key("substeps", Integer(), Default("10"),
				"Number of equidistant substeps, molar mass and isotherm intersections");


Record Sorption::input_type
	= Record("Sorptions", "Information about all the limited solubility affected sorptions.")
	.derive_from( Reaction::input_type )
    .declare_key("sorptions", Array( Sorption::input_type_isotherm ), Default::obligatory(),
                "Description of particular sorption cases under consideration.");

using namespace std;

Sorption::Sorption(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)
	: Reaction(init_mesh, in_rec, names)
{
	nr_of_regions = 1; // temporarirly suppresed //material_database.size();
	nr_of_substances = names.size();

	//both following operations connected together in compute_isotherm(..) method
	//prepare_inputs(in_rec);
	//determine_crossection()
	//compute_isotherms(in_rec);
}

// TODO: check duplicity of parents
//       raise warning if sum of ratios is not one
void Sorption::compute_isotherms(Input::Record in_rec)
{
    unsigned int idx;

	Input::Array sorption_array = in_rec.val<Input::Array>("sorptions");
	/*n_points = in_rec.val<int>("substeps");

	isotherms.resize(nr_of_regions*nr_of_substances);
	for(int i = 0; i < nr_of_regions; i++)
	{
		for(int j = 0; j < nr_of_substances; j++)
		{
			isotherms[i][j] = 0; //initialization
		}
	}

	int i_sorption=0;
	for (Input::Iterator<Input::Record> sorp_it = sorption_array.begin<Input::Record>(); sorp_it != sorption_array.end(); ++sorp_it, ++i_sorption)
	{
		//isotherm determining part
		int region_id = sorp_it->find<int>("region");
		int substance_id = sorp_it->find<int>("specie");
		isotherms[region_id][specie_id] = new Isotherm(init_mesh, *sorp_it, names);
		(&isotherms[region_id][specie_id])->set_parameters(*sorp_it);
		determine_crossections(n_points);//
		rotate_points(0.70711, isotherms[region_id][substance_id]);*/

		//Sorption_type sorption_type = sorp_it->find<enum Sorption_type>("type");
		/*double slope = sorp_it->find<double>("direction");
		if (slope > 0.0) {
		   ; //intersections of an isotherm with mass balance lines can be found exactly
		} else {
		   double alpha = sorp_it->find<double>("alpha");
		   double omega = sorp_it->find<double>("omega");
		   if ((alpha == -1.0) || (omega == -1.0)) {
			   xprintf(UsrErr, "Some of parameters for isotherm are either missing or incorect.\n");
		   } else {
			   int substance_id = sorp_it->find<int>("specie");
			   int region_id = sorp_it->find<int>("region");
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

	//}
}

double **Sorption::compute_reaction(double **concentrations, int loc_el) // Sorptions are realized just for one element.
{
    int cols, rows;

    //  If intersections of isotherm with mass balance lines are known, then interpolate.
    	//  Measurements [c_a,c_s] will be rotated
    	//  Rotated measurements must be projected on rotated isotherm, interpolate_datapoints()
    	//  Projections need to be transformed back to original CS
    //	If intersections are not known then solve the problem analytically (toms748_solve).

    //*if (reaction_matrix == NULL) return concentrations;

	/*for(cols = 0; cols < n_substances(); cols++){
		prev_conc[cols] = concentrations[cols][loc_el];
		concentrations[cols][loc_el] = 0.0;
	}*/

	/*for(rows = 0; rows < n_substances(); rows++){
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


void Sorption::print_sorption_parameters(void)
{

    xprintf(Msg, "\nSorption parameters are defined as:");
    /*for (int i = 0; i < (nr_of_substances - 1); i++) {
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
    }*/
}

void Sorption::determine_crossections(void)
{
	;
}

void Sorption::rotate_point(double angle, std::vector<std::vector<double> > points)
{
	;
}

void Sorption::interpolate_datapoints(void)
{
	;
}
