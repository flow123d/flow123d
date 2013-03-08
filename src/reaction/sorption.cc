#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <boost/foreach.hpp>

#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "reaction/sorption.hh"
#include "reaction/sorption_impl.hh"
#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "la/distribution.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "transport/transport.h" //because of definition of constants MOBILE, IMMOBILE,

const double pi = 3.1415;
namespace it=Input::Type;

Sorption::EqData::EqData(const std::string &name)
: EqDataBase(name)
{
    //ADD_FIELD(nr_of_points, "Number of crossections allong isotherm", Input::Type::Default("10"));
    //ADD_FIELD(region_ident, "Rock matrix area identifier.", Input::Type::Default("0"));
    ADD_FIELD(rock_density, "Rock matrix density.", Input::Type::Default("0.0"));
    //ADD_FIELD(solvent_density, "Solvent density.", Input::Type::Default("1.0"));

    ADD_FIELD(sorption_types,"Considered adsorption is described by selected isotherm.", it::Default("none") );
              sorption_types.set_selection(&sorption_type_selection);

    ADD_FIELD(mult_coefs,"Multiplication parameters (k, omega) in either Langmuir c_s = omega * (alpha*c_a)/(1- alpha*c_a) or in linear c_s = k * c_a isothermal description.");
    //std::vector<FieldEnum> list; list.push_back(none); list.push_back(Langmuir); list.push_back(Freundlich);
    //slopes.disable_where(& sorption_type, list );

    //ADD_FIELD(omegas,"Langmuir isotherm multiplication parameters in c_s = omega * (alpha*c_a)/(1- alpha*c_a).");
    //list.clear(); list.push_back(none); list.push_back(Linear); list.push_back(Freundlich);
    //omegas.disable_where(& sorption_type, list );

    ADD_FIELD(alphas,"Second parameters (alpha, ...) defining isotherm  c_s = omega * (alpha*c_a)/(1- alpha*c_a).");
    //list.clear(); list.push_back(none); list.push_back(Linear); list.push_back(Freundlich);
    //alphas.disable_where(& sorption_type, list );
    ADD_FIELD(mob_porosity,"Mobile porosity of the rock matrix.");
    ADD_FIELD(immob_porosity,"Immobile porosity of the rock matrix.", Input::Type::Default("0.0"));
}

using namespace Input::Type;

Record Sorption::input_type
	= Record("Sorptions", "Information about all the limited solubility affected sorptions.")
	.derive_from( Reaction::input_type )
	.declare_key("solv_dens", Double(), Default("1.0"),
				"Density of the solvent.")
	.declare_key("substeps", Integer(), Default("10"),
				"Number of equidistant substeps, molar mass and isotherm intersections")
	.declare_key("species", Array(String()), Default::obligatory(),
							"Names of all the sorbing species")
	.declare_key("molar_mass", Array(Double()), Default::obligatory(),
							"Specifies molar masses of all the sorbing species")
	.declare_key("solubility", Array(Double()), Default::obligatory(),
							"Specifies solubility limits of all the sorbing species")
    /*.declare_key("sorptions", Array( Sorption::input_type_isotherm ), Default::obligatory(),
                "Description of particular sorption cases under consideration.")*/
    .declare_key("bulk_data", Array(Sorption::EqData().bulk_input_type()), Default::obligatory(),
                   	   	   "Containes region specific data necessery to construct isotherms.");

using namespace std;

Sorption::Sorption(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)
	: Reaction(init_mesh, in_rec, names)
{
	nr_of_regions = init_mesh.n_materials;
	nr_of_substances = names.size();
	nr_of_points = in_rec.val<int>("substeps");
	solvent_dens = in_rec.val<int>("solv_dens");
	//Simple vectors holding  common informations.
	region_ids.resize( nr_of_regions ); // ( nr_of_regions );
	substance_ids.resize(nr_of_substances); // ( nr_of_substances );
	molar_masses.resize( nr_of_substances );
	c_aq_max.resize( nr_of_substances );
	//isotherms array resized bellow
	isotherms_mob.resize(nr_of_regions*nr_of_substances);
	if(dual_porosity_on) isotherms_immob.resize(nr_of_regions*nr_of_substances);
	//Multidimensional array isotherm
	for(int i_reg = 0; i_reg < nr_of_regions; i_reg++)
	{
		//isotherm[i_mob].resize(nr_of_regions);
		for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
		{

			//Isotherm iso_hlp();
			Isotherm (*isotherms_mob[i_reg][i_subst])();
			//----------------------------------------
			//Isotherm *iso_hlp = new Isotherm();
			//isotherms_mob[i_reg][i_subst] = iso_hlp;
			//----------------------------------------
			//isotherms[i_reg][i_subst].reinit(rock_dens[i_subst], solvent_dens, double porosity, double molar_mass, double c_aqua_limit);
			//isotherms_mob[i_reg][i_subst].intersections.resize(n_points); // obsolete nonsense
			if(dual_porosity_on) Isotherm (*isotherms_immob[i_reg][i_subst])(); // obsolete nonsense
		}
	}

	TimeGovernor tg(0.0, 1.0);

    //data.init_conc.set_n_comp(4);        // set number of substances posibly read from elsewhere
    //data.bc_conc.set_n_comp(4);

    data_.set_mesh(&init_mesh);
    data_.init_from_input( in_rec.val<Input::Array>("bulk_data"),Input::Array() );
    data_.set_time(tg);

    //cycle over regions and species to create particular isotherms_mob and isotherms_immob objects of Isotherm type
    	//Isotherm(double rock_density,double solvent_density, double porosity, double molar_mass, double step_length);
}

void Sorption::prepare_inputs(Input::Record in_rec)
{
    unsigned int idx;

	Input::Array sorption_array = in_rec.val<Input::Array>("sorptions");

	//Multidimensional array isotherm, initialization
	/*for(int i_reg = 0; i_reg < nr_of_regions; i_reg++)
	{
		//isotherm[i_mob].resize(nr_of_regions);
		for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
		{
			//isotherms[i_reg][i_subst].set_step_length(c_)
			for(int i_point = 0; i_point < nr_of_points; i_point++)
			{
				isotherms_mob[i_reg][i_subst].intersections[i_point] = 0;
				if(dual_porosity_on)isotherms_mob[i_reg][i_subst].intersections[i_point] = 0;
			}
		}
	}*/

	int i_sorp = 0;
	for (Input::Iterator<Input::Record> sorp_it = sorption_array.begin<Input::Record>(); sorp_it != sorption_array.end(); ++sorp_it, ++i_sorp)
	{
		int idx;
		double mol_mass, solub;

		//indices determining part
		string specie_name = sorp_it->val<string>("specie");
		idx = find_subst_name(specie_name);
		if (idx < n_substances())
		{
			substance_ids[i_sorp] = idx;
		}else{
			//xprintf(UsrErr,"Unknown name %s of substance undergoing the %d-th sorption.\n", specie_name, i_sorp);
			xprintf(UsrErr,"Unknown name (identifier) of the substance undergoing the %d-th sorption.\n", i_sorp);
		}

		// solubulity limit
		solub = sorp_it->val<double>("solubility");
		if (solub > 0.0) {
		   c_aq_max[idx] = solub;
		} else {
			//xprintf(UsrErr, "Unknown solubility limit of substance %s undergoing the %d-th sorption.\n", specie_name, i_sorp);
			xprintf(UsrErr, "Unknown solubility limit of the substance undergoing the %d-th sorption.\n", i_sorp);
		}

		//molar masses determining part
		mol_mass = sorp_it->val<double>("molar_mass");
		if (mol_mass) {
		   molar_masses[idx] = mol_mass;
		} else {
			//xprintf(UsrErr, "Unknown molar mass of substance %s undergoing the %d-th sorption.\n", specie_name, i_sorp);
			xprintf(UsrErr, "Unknown molar mass of the substance undergoing the %d-th sorption.\n", i_sorp);
		}

		//region determining part
		idx = sorp_it->val<int>("region");
		if(idx)
		{
			region_ids[i_sorp] = idx;
		}else{
			xprintf(UsrErr, "Undefined region identifier where the %d-th sorption takes place.\n", i_sorp);
		}

	}

}

/*void Sorption::precompute_isotherm_tables() {
    BOOST_FOREACH(const Region & reg, this->mesh_->region_db().get_region_set("BULK")) {
        arma::Col<unsigned int> sorption_type_vec;
        arma::Col<double> scale_vec, alpha_vec;
        if (data_.sorption_types.get_const_value(reg, sorption_type_vec))
        if (data_.mult_coefs.get_const_value(reg, scale_vec))
        if (data_.alphas.get_const_value(reg, alpha_vec)) {
            // precompute isotherm
        }
        // else leave isotherm empty for this region
    }
}*/

// TODO: check duplicity of parents
//       raise warning if sum of ratios is not one

/*void Sorption::handle_datapoints(double rock_density, double porosity, std::vector<double> &prev_conc, std::vector<double> isotherm, int reg_id_nr, int i_subst)
{
	double coef_hlp, conc_hlp;
	int iso_ind_floor, iso_ind_ceil;
	std::vector<double> rot_point;
	rot_point.resize(2);

	//compute_rot_coefs(porosity, rock_density, i_subst); // computes rotation matrix entries
	rot_coefs[0] = porosity*solvent_dens;
	rot_coefs[1] = molar_masses[i_subst]*(1-porosity)*rock_density;
	//rot_point = rotate_point(previous_conc); //counterclockwise rotation to mass balancing coordination system
	rot_point[0] = rot_coefs[0]*prev_conc[0] + rot_coefs[1]*prev_conc[1]; //coordinate x^R or c_a^R or m'_x(aqueous == dissolved)
	//rot_point[1] = interpolate_datapoint(rot_point, reg_id_nr, i_subst); // interpolation in mass ballancing coordination system
	iso_ind_floor = (int)(rot_point[0]/(step_length)); iso_ind_ceil = iso_ind_floor + 1;
	rot_point[1] = isotherm[iso_ind_floor] + (rot_point[0] - isotherm[iso_ind_floor])*(isotherm[iso_ind_ceil] - isotherm[iso_ind_floor])/step_length;
	//switch_rot_coefs();
	coef_hlp = rot_coefs[0];
	rot_coefs[0] = (-1.0)*rot_coefs[1];
	rot_coefs[1] = coef_hlp;

	//=========================================================
	//		Alternative to backward rotation and scaling
	//=========================================================
	prev_conc[0] = (rot_point[0] + rot_point[1])/(2*rot_coef[0]); // coordinate x^R or c_a^R (aqueous == dissolved)
	prev_conc[1] = (rot_point[0] - rot_point[1])/(2*rot_coef[1]); // coordinate y^R or c_s^R (sorbed == solid)

	return;
}*/

double **Sorption::compute_reaction(double **concentrations, int loc_el) // Sorptions are realized just for one element.
{
    ElementFullIter elem = mesh_->element(el_4_loc[loc_el]);
    double mob_porosity, immob_porosity; // = data_.mob_porosity.value(elem->centre(),elem->element_accessor());
    double rock_density; // = data_.rock_density.value(elem->centre(),elem->element_accessor());
    double k_rep;
    Region region = elem->region();
    int reg_id_nr = region.idx(); //->region_->reg_id;

    //	Identify loc_el region.
    //  If intersections of isotherm with mass balance lines are known, then interpolate.
    	//  Measurements [c_a,c_s] will be rotated
    	//  Rotated measurements must be projected on rotated isotherm, interpolate_datapoints()
    	//  Projections need to be transformed back to original CS
    //	If intersections are not known then solve the problem analytically (toms748_solve).

	if( data_.rock_density.get_const_value(region, rock_density)) // constant value of rock density over whole the region
	{
		if(data_.mob_porosity.get_const_value(region, mob_porosity)) // constant values of porosity over whole the region
		{
			for(int i_subst = 0; i_subst < n_substances(); i_subst++)
			{
				isotherms_mob[reg_id_nr][i_subst].compute_projection(concentration_matrix[IMMOBILE][i_subst][loc_el], concentration_matrix[IMMOBILE_SORB][i_subst][loc_el]);
			}
		}else{
			cout << "It is not possible in this time to compute sorption if the porosity is not constant over whole the region " << reg_id_nr << endl;
		}
		if(dual_porosity_on && data_.immob_porosity.get_const_value(region, immob_porosity))
		{
			for(int i_subst = 0; i_subst < n_substances(); i_subst++)
			{
				isotherms_immob[reg_id_nr][i_subst].compute_projection(concentration_matrix[IMMOBILE][i_subst][loc_el], concentration_matrix[IMMOBILE_SORB][i_subst][loc_el]);
			}
		}else{
			cout << "It is not possible in this time to compute sorption in immobile pores if the dual porosity  is not constant over whole the region " << reg_id_nr << endl;
		}
	}else{
		cout << "It is not possible in this time to compute sorption in if the rock density is not constant over whole the region " << reg_id_nr << endl;
	}

	return concentrations;
}

void Sorption::compute_one_step(void) // Computes sorption simulation over all the elements.
{
    START_TIMER("sorption_step");
	for (int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
	 {
	 	this->compute_reaction(concentration_matrix[0], loc_el); //MOBILE and IMMOBILE 	computed

	 }
    END_TIMER("sorption_step");
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

/*Isotherm::Isotherm(double rock_density,double solvent_density, double porosity, double molar_mass, double step_length)
{
	k_W = porosity*solvent_density;
	k_H = molar_mass*(1-porosity)*rock_density;
}

void Isotherm::compute_projection(double &c_aq, double &c_sorb)
{
	int iso_ind_floor, iso_ind_ceil;
	double conc_hlp = c_aq;

	//coordinates are transformed
	conc_hlp = c_aq;
	c_aq = k_W*c_aq + k_H*c_sorb; //coordinate x^R or c_a^R (aqueous == dissolved)

	//interpolation is realized
	iso_ind_floor = (int)(c_aq/(step_length)); iso_ind_ceil = iso_ind_floor + 1;
	c_sorb = intersections[iso_ind_floor] + (c_aq - iso_ind_floor*step_length)*(intersections[iso_ind_ceil] - intersections[iso_ind_floor])/step_length;

	//tarnsformation back to original coordination system
	conc_hlp = c_aq;
	c_aq = (c_aq + c_sorb)/(2*k_W); // coordinate x or c_a (aqueous == dissolved)
	c_sorb = (conc_hlp - c_sorb)/(2*k_H); // coordinate y or c_s (sorbed == solid)

	return;
}

void Isotherm::set_step_length(double c_aq_max, int n_steps)
{
	double c_sorb_max; //= some value

	step_length = (k_W*c_aq_max + k_H*c_sorb_max)/n_steps;
}*/
