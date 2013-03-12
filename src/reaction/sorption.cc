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
	.declare_key("solvent_dens", Double(), Default("1.0"),
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
	TimeGovernor tg(0.0, 1.0);

    //data.init_conc.set_n_comp(4);        // set number of substances posibly read from elsewhere
    //data.bc_conc.set_n_comp(4);

    data_.set_mesh(&init_mesh);
    //Input::Array isotherms =
    		data_.init_from_input( in_rec.val<Input::Array>("bulk_data"),Input::Array() );
    data_.set_time(tg);

    Input::Array sorptions_array = in_rec.val<Input::Array>("bulk_data"); // no idea how to get infos from data_.init_from_input( in_rec.val<Input::Array>("bulk_data"),Input::Array() );

	nr_of_regions = init_mesh.n_materials;
	nr_of_substances = names.size();
	nr_of_points = in_rec.val<int>("substeps");

	//Simple vectors holding  common informations.
	region_ids.resize( nr_of_regions ); // ( nr_of_regions );
	substance_ids.resize(nr_of_substances); // ( nr_of_substances );
	molar_masses.resize( nr_of_substances );
	c_aq_max.resize( nr_of_substances );

	//common data for all the isotherms loaded bellow
	solvent_dens = in_rec.val<double>("solvent_dens");

	Input::Array molar_mass_array = in_rec.val<Input::Array>("molar_masses");
	//molar_masses = in_rec.val<Array(Double())>("molar_masses");
	if (molar_mass_array.size() == molar_masses.size() )   molar_mass_array.copy_to( molar_masses );
	  else  xprintf(UsrErr,"Number of molar masses %d has to match number of sorbing species %d.\n",
	                                       molar_mass_array.size(), molar_masses.size());
	//c_aq_max = in_rec.val<Array(Double())>("solubility");

	Input::Array solub_limit_array = in_rec.val<Input::Array>("solubility");
	if (solub_limit_array.size() == c_aq_max.size() )   solub_limit_array.copy_to( c_aq_max );
	  else  xprintf(UsrErr,"Number of given solubility limits %d has to match number of sorbing species %d.\n",
	                                       solub_limit_array.size(), c_aq_max.size());

	Input::Array species_array = in_rec.val<Input::Array>("species");
	unsigned int idx, i_spec = 1;
	for(Input::Iterator<string> spec_iter = species_array.begin<string>(); spec_iter != species_array.end(); ++spec_iter, i_spec++)
	{
		idx = find_subst_name(*spec_iter);
		if (idx < n_substances())   substance_ids[i_spec] = idx;
		else	xprintf(Msg,"Wrong name of %d-th sorbing specie.\n", i_spec);
	}

	//isotherms array resized bellow
	//isotherms_mob.resize(nr_of_regions*nr_of_substances);
	isotherms_mob.resize(nr_of_regions);
	for(int i_reg = 0; i_reg < nr_of_regions; i_reg++)
	{
		isotherms_mob[i_reg].resize(nr_of_substances);
	}
	if(dual_porosity_on)
	{
		//isotherms_immob.resize(nr_of_regions*nr_of_substances);
		isotherms_immob.resize(nr_of_regions);
			for(int i_reg = 0; i_reg < nr_of_regions; i_reg++)
			{
				isotherms_immob[i_reg].resize(nr_of_substances);
			}
	}

	// list of types of isotherms in particular regions
	std::vector<SorptionType> iso_type; iso_type.resize(nr_of_substances);
	//Multidimensional array
	int i_reg = 0;
	for(Input::Iterator<Input::Record> reg_iter = sorptions_array.begin<Input::Record>(); reg_iter != sorptions_array.end(); ++reg_iter, i_reg++)
	{
		// list of types of isotherms in particular regions, initialization
		Input::Array sorption_types_array = reg_iter->val<Input::Array>("sorption_types");
		if (sorption_types_array.size() == iso_type.size() )   sorption_types_array.copy_to( iso_type );
		  else  xprintf(UsrErr,"Number of sorption types %d has to match number of sorbing species %d.\n",
		                                       sorption_types_array.size(), iso_type.size());
		//isotherm[i_mob].resize(nr_of_regions);
		double rock_density = reg_iter->val<double>("rock_density");
		double mobile_porosity = reg_iter->val<double>("mob_porosity");
		double immobile_porosity = reg_iter->val<double>("immob_porosity");

		for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
		{
			// reinit isotherm, what about to define a type of isotherm in reinit
			isotherms_mob[i_reg][i_subst].reinit(iso_type[i_subst],rock_density,solvent_dens,mobile_porosity, molar_masses[i_subst], c_aq_max[i_subst]);
			//precompute necessary multiplication coefficient
			double k = isotherms_mob[i_reg][i_subst].get_scale_sorbed()/isotherms_mob[i_reg][i_subst].get_scale_aqua();
			//create temporary class Linear object
			//Linear lin_iso(k);
			switch (iso_type[i_subst])
			{
			 case 1: //Linear:
				//call appropriate function for isotherm descritption
				//isotherms_mob[i_reg][i_subst].make_table(lin_iso, nr_of_points);
			 break;
			 case 2: //Langmuir:
			 	 ;
			 break;
			 case 3: //Freundlich:
				 ;
			 break;
			 default:
				xprintf(UsrErr, "Unknown type of sorption of %d-th specie in %d-th region.\n", i_subst, i_reg);
			}
			if(dual_porosity_on)
			{
				//Isotherm (*isotherms_immob[i_reg][i_subst])();
				isotherms_immob[i_reg][i_subst].reinit(iso_type[i_subst],rock_density,solvent_dens,immobile_porosity, molar_masses[i_subst], c_aq_max[i_subst]);
			}
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
