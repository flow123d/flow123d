#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <boost/foreach.hpp>

#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "reaction/isotherm.hh"
#include "reaction/sorption.hh"
#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "la/distribution.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "transport/transport.h" //because of definition of constants MOBILE, IMMOBILE,
#include "input/type_selection.hh"

const double pi = 3.1415;
namespace it=Input::Type;

it::Selection Sorption::EqData::sorption_type_selection = it::Selection("SorptionType")
	.add_value(none,"none","No sorption considered")
	.add_value(linear,"linear","Linear isotherm described sorption considered.")
	.add_value(langmuir,"langmuir","Langmuir isotherm described sorption considered")
	.add_value(freundlich,"freundlich","Freundlich isotherm described sorption considered");
//.finish();

Sorption::EqData::EqData()
: EqDataBase("Sorption")
{
    ADD_FIELD(rock_density, "Rock matrix density.", Input::Type::Default("0.0"));

    ADD_FIELD(sorption_types,"Considered adsorption is described by selected isotherm." );
              sorption_types.set_selection(&sorption_type_selection);
    //ADD_FIELD(sorption_types,"Considered adsorption is described by selected isotherm.", it::Default(0) );

    ADD_FIELD(mult_coefs,"Multiplication parameters (k, omega) in either Langmuir c_s = omega * (alpha*c_a)/(1- alpha*c_a) or in linear c_s = k * c_a isothermal description.", Input::Type::Default("1.0"));
    //std::vector<FieldEnum> list; list.push_back(none); //SorptionType
    //mult_coefs.disable_where(&sorption_types, list ); //function disable where requires different parameters

    ADD_FIELD(second_params,"Second parameters (alpha, ...) defining isotherm  c_s = omega * (alpha*c_a)/(1- alpha*c_a).", Input::Type::Default("1.0"));
    //list.clear(); list.push_back(none); list.push_back(linear);
    //alphas.disable_where(&sorption_types, list );

    ADD_FIELD(mob_porosity,"Mobile porosity of the rock matrix.", Input::Type::Default("0.0000001"));
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
	.declare_key("molar_masses", Array(Double()), Default::obligatory(),
							"Specifies molar masses of all the sorbing species")
	.declare_key("solubility", Array(Double()), Default::obligatory(),
							"Specifies solubility limits of all the sorbing species")
    .declare_key("bulk_data", Array(Sorption::EqData().bulk_input_type()), Default::obligatory(), //
                   	   	   "Containes region specific data necessery to construct isotherms.");

using namespace std;

Sorption::Sorption(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)
	: Reaction(init_mesh, in_rec, names)
{
	cout << "Sorption constructor is running." << endl;
	TimeGovernor tg(0.0, 1.0);
    nr_of_regions = init_mesh.region_db().bulk_size();
    nr_of_substances = in_rec.val<Input::Array>("species").size();
    nr_of_points = in_rec.val<int>("substeps");


    data_.sorption_types.set_n_comp(nr_of_substances);        // set number of substances posibly read from elsewhere
    data_.mult_coefs.set_n_comp(nr_of_substances);
    data_.second_params.set_n_comp(nr_of_substances);
    data_.set_mesh(&init_mesh);
    data_.init_from_input( in_rec.val<Input::Array>("bulk_data"), Input::Array());
    data_.set_time(tg);

	//Simple vectors holding  common informations.
	region_ids.resize( nr_of_regions ); // ( nr_of_regions );
	substance_ids.resize(nr_of_substances); // ( nr_of_substances );
	molar_masses.resize( nr_of_substances );
	c_aq_max.resize( nr_of_substances );

	//isotherms array resized bellow
	//isotherms_mob.resize(nr_of_regions*nr_of_substances);
	//free(isotherms_mob);
	isotherms_mob.resize(nr_of_regions);
	//isotherms_mob = (std::vector *) malloc(nr_of_regions*sizeof(std::vector *));
	for(int i_reg = 0; i_reg < nr_of_regions; i_reg++)
	{
		//check the size of isotherms_mob
			//if it is smaller then i_reg then add i-th isotherm
			//else use reinit().
		Isotherm iso_mob;
		isotherms_mob[i_reg].push_back(iso_mob);
		isotherms_mob[i_reg].resize(nr_of_substances);
		/*for(int i_spec = 0; i_spec < nr_of_substances; i_spec++)
		{
			//isotherms_mob[i_reg][i_spec] = *(new Isotherm);
		}*/
}
	if(dual_porosity_on)
	{
		//isotherms_immob.resize(nr_of_regions*nr_of_substances);
		//free(isotherms_immob);
		isotherms_immob.resize(nr_of_regions);
			for(int i_reg = 0; i_reg < nr_of_regions; i_reg++)
			{
				Isotherm iso_immob;
				isotherms_immob[i_reg].push_back(iso_immob);
				isotherms_immob[i_reg].resize(nr_of_substances);
				/*for(int i_spec = 0; i_spec < nr_of_substances; i_spec++)
				{
					//isotherms_immob[i_reg][i_spec] = *(new Isotherm);
				}*/
			}
	}
	prepare_inputs(in_rec);
}

Sorption::~Sorption(void)
{
}

void Sorption::init_from_input(Input::Array bulk_list)
{
	//Not sure what to write here.
	return;
}

void Sorption::prepare_inputs(Input::Record in_rec)
{
    Input::Array sorptions_array = in_rec.val<Input::Array>("bulk_data");
	//common data for all the isotherms loaded bellow
	solvent_dens = in_rec.val<double>("solvent_dens");

	Input::Array molar_mass_array = in_rec.val<Input::Array>("molar_masses");
	//molar_masses = in_rec.val<Array(Double())>("molar_masses");
	if (molar_mass_array.size() == molar_masses.size() )   molar_mass_array.copy_to( molar_masses );
	  else  xprintf(UsrErr,"Number of molar masses %d has to match number of sorbing species %d.\n", molar_mass_array.size(), molar_masses.size());

	Input::Array solub_limit_array = in_rec.val<Input::Array>("solubility");
	if (solub_limit_array.size() == c_aq_max.size() )   solub_limit_array.copy_to( c_aq_max );
	  else  xprintf(UsrErr,"Number of given solubility limits %d has to match number of sorbing species %d.\n", solub_limit_array.size(), c_aq_max.size());

	Input::Array species_array = in_rec.val<Input::Array>("species");
	unsigned int idx, i_spec = 0;
	for(Input::Iterator<string> spec_iter = species_array.begin<string>(); spec_iter != species_array.end(); ++spec_iter, i_spec++)
	{
		idx = find_subst_name(*spec_iter);
		if ((idx < n_substances()) && (idx >= 0))   substance_ids[i_spec] = idx;
		else	xprintf(Msg,"Wrong name of %d-th sorbing specie.\n", i_spec);
	}

	// list of types of isotherms in particular regions
	//FieldValue<3>::EnumVector::return_type iso_type;
	arma::Col<unsigned int> iso_type;
	//FieldValue<3>::Vector::return_type iso_type;
	cout << "there are " << nr_of_substances <<" substances under concideration." << endl;
	iso_type.resize(nr_of_substances); //arma::Col<unsigned int> je ten typ ze začátku řádku, std::vector<SorptionType>
	//std::vector<FieldEnum> iso_type; iso_type.resize(nr_of_substances);
	// list of sorption parameters
	FieldValue<3>::Vector::return_type mult_param;
	//arma::Col<double> mult_param;
	mult_param.resize(nr_of_substances);
	FieldValue<3>::Vector::return_type second_coef;
	//arma::Col<double> second_coef;
	second_coef.resize(nr_of_substances);
	double rock_density, mobile_porosity, immobile_porosity;
	//Multidimensional array
	int i_reg = 0;
	//std::map<SorptionType, std::string>;
	//for(Input::Iterator<Input::Record> reg_iter = sorptions_array.begin<Input::Record>(); reg_iter != sorptions_array.end(); ++reg_iter, i_reg++)
	BOOST_FOREACH(const Region &reg_iter, this->mesh_->region_db().get_region_set("BULK") )
	{
		// list of types of isotherms in particular regions, initialization
		if(data_.sorption_types.get_const_value(reg_iter, iso_type))
		{
			; //for(int index_latky = 0; index_latky < nr_of_substances; index_latky++) cout << "Type of isotherm of " << index_latky << " specie is " << iso_type[index_latky] << endl;
			//xprintf(Msg,"Type of isotherm of %d-th specie is %d.\n", index_latky, iso_type[index_latky]);
		}else  xprintf(UsrErr,"Type of isotherm must be the same all over the %d-th region, but it is not.", reg_iter.id());

		// multiplication coefficient parameter follows
		if(data_.mult_coefs.get_const_value(reg_iter, mult_param)) ;
		  else  xprintf(UsrErr,"All the sorption parameters must be constant all over the %d-th region, but it is not.", reg_iter.id());

		// second sorption parameter for isothermal description
		if(data_.second_params.get_const_value(reg_iter, second_coef)) ;
		  else  xprintf(UsrErr,"All the sorption parameters must be constant all over the %d-th region, but it is not.", reg_iter.id());

		// rock material characteristics are read bellow
		if(data_.rock_density.get_const_value(reg_iter, rock_density)) ;
		  else xprintf(UsrErr,"Rock density is suppossed to be constant all over the %d-th region, but it is not.", reg_iter.id());
		if(data_.mob_porosity.get_const_value(reg_iter, mobile_porosity)) ;
		  else xprintf(UsrErr,"Mobile porosity is suppossed to be constant all over the %d-th region, but it is not.", reg_iter.id());
		if(data_.immob_porosity.get_const_value(reg_iter, immobile_porosity)) ;
		  else xprintf(UsrErr,"Immobile porosity is suppossed to be constant all over the %d-th region, but it is not.", reg_iter.id());

		for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
		{
			// reinit isotherm, what about to define a type of isotherm in reinit
			SorptionType hlp_iso_type =  SorptionType(iso_type[i_subst]);
			//int hlp_iso_type =  int(iso_type[i_subst]);
			// did not function //SorptionType hlp_iso_type =  (SorptionType) iso_type[i_subst];
			//xprintf(Msg,"Sorption type of %d-th substance is %d.\n",i_subst, iso_type[i_subst]);

			int reg_idx=reg_iter.bulk_idx();
			isotherms_mob[reg_idx][i_subst].reinit(hlp_iso_type,rock_density,solvent_dens,mobile_porosity, molar_masses[i_subst], c_aq_max[i_subst]);
			//cout << "This message should indicate fault." << endl;
			if(dual_porosity_on)
			{
				isotherms_immob[reg_idx][i_subst].reinit(hlp_iso_type,rock_density,solvent_dens,immobile_porosity, molar_masses[i_subst], c_aq_max[i_subst]);
			}
			switch(hlp_iso_type)
			{
			 case none: // 0: //
			 {
				 //xprintf(Msg,"No sorption is considered for %d-th specie in %d-th region.\n", i_subst, reg_idx);
				 int one_point = 1;
				 isotherms_mob[reg_idx][i_subst].make_one_point_table();
				 if(dual_porosity_on)
				 {
					const Linear obj_isotherm_immob(mult_param[i_subst]);
					isotherms_immob[reg_idx][i_subst].make_one_point_table();
				 }

			 }
			 break;
			 case  linear: // 1: //
			 {
				Linear obj_isotherm(mult_param[i_subst]);
				isotherms_mob[reg_idx][i_subst].make_table(obj_isotherm, nr_of_points);
				if(dual_porosity_on)
				{
					const Linear obj_isotherm_immob(mult_param[i_subst]);
					isotherms_immob[reg_idx][i_subst].make_table(obj_isotherm_immob, nr_of_points);
				}
			 }
			 break;
			 case langmuir: // 2: //
			 {
			 	Langmuir obj_isotherm(mult_param[i_subst], second_coef[i_subst]);
				isotherms_mob[reg_idx][i_subst].make_table(obj_isotherm, nr_of_points);
			 	if(dual_porosity_on)
			 	{
				 	Langmuir obj_isotherm_immob(mult_param[i_subst], second_coef[i_subst]);
					isotherms_mob[reg_idx][i_subst].make_table(obj_isotherm_immob, nr_of_points);
			 	}
			 }
			 break;
			 case freundlich: // 3: //
			 {
				 xprintf(Msg,"Freundlich isotherm is not implemented yet.");
			 }
			 break;
			 default:
			 {
				 xprintf(Msg,"Sorption of %d-th specie in %d-th region has type %s.", i_subst, i_reg, hlp_iso_type);
			 }
			 break;
			}
		}
	}
}

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
    int reg_id_nr = region.bulk_idx(); //->region_->reg_id;

    //	Identify loc_el region.
    //  If intersections of isotherm with mass balance lines are known, then interpolate.
    	//  Measurements [c_a,c_s] will be rotated
    	//  Rotated measurements must be projected on rotated isotherm, interpolate_datapoints()
    	//  Projections need to be transformed back to original CS
    //	If intersections are not known then solve the problem analytically (toms748_solve).

	if( data_.rock_density.get_const_value(region, rock_density)) // constant value of rock density over the whole region
	{
		if(data_.mob_porosity.get_const_value(region, mob_porosity)) // constant values of porosity over the whole region
		{
			for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
			{
			    if(isotherms_mob[reg_id_nr][i_subst].get_interpolation_table_size() >= 2)
			    {
					int subst_id =substance_ids[i_subst];
					isotherms_mob[reg_id_nr][i_subst].compute_projection(concentration_matrix[IMMOBILE][subst_id][loc_el], concentration_matrix[IMMOBILE_SORB][subst_id][loc_el]);
			    }else{
			    	cout << "The isotherm is either not specified or it is defined as 'none'" << endl;
			    }
			}
		}else{
			cout << "It is not possible in this time to compute sorption if the porosity is not constant over the whole region " << reg_id_nr << endl;
		}
		if(dual_porosity_on && data_.immob_porosity.get_const_value(region, immob_porosity))
		{
			for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
			{
			    if(isotherms_immob[reg_id_nr][i_subst].get_interpolation_table_size() >= 2)
			    {
					int subst_id =substance_ids[i_subst];
					isotherms_immob[reg_id_nr][i_subst].compute_projection(concentration_matrix[IMMOBILE][subst_id][loc_el], concentration_matrix[IMMOBILE_SORB][subst_id][loc_el]);
			    }else{
			    	cout << "The isotherm is either not specified or it is defined as 'none'" << endl;
			    }
			}
		}else{
			cout << "It is not possible in this time to compute sorption in immobile pores if the dual porosity  is not constant over the whole region " << reg_id_nr << endl;
		}
	}else{
		cout << "It is not possible in this time to compute sorption in if the rock density is not constant over the whole region " << reg_id_nr << endl;
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

/**
* Meaningless inherited methods.
*/
void Sorption::update_solution(void)
{
	cout << "Meaningless inherited method." << endl;
	return;
}
void Sorption::choose_next_time(void)
{
	cout << "Meaningless inherited method." << endl;
	return;
}

void Sorption::set_time_step_constrain(double dt)
{
	cout << "Meaningless inherited method." << endl;
	return;
}

void Sorption::get_parallel_solution_vector(Vec &vc)
{
	cout << "Meaningless inherited method." << endl;
	return;
}

void Sorption::get_solution_vector(double* &vector, unsigned int &size)
{
	cout << "Meaningless inherited method." << endl;
	return;
}

void Sorption::set_time_step(double new_timestep)
{
	//cout << "This method is obsolete for equilibrial sorptions and reactions, but it must be implemented." << endl;
	return;
}

void Sorption::set_time_step(Input::Record in_rec)
{
	cout << "This method is obsolete for equilibrial sorptions and reactions, but it must be implemented." << endl;
	return;
}
