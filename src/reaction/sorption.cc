#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <boost/foreach.hpp>

#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
//#include "reaction/isotherm.hh"
#include "reaction/sorption.hh"
#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "la/distribution.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "mesh/region.hh"
//#include "transport/transport.h" //because of definition of constants MOBILE, IMMOBILE,
#include "input/type_selection.hh"

const double pi = 3.1415;
namespace it=Input::Type;

it::Selection Sorption::EqData::sorption_type_selection = it::Selection("SorptionType")
	.add_value(none,"none","No sorption considered")
	.add_value(linear,"linear","Linear isotherm described sorption considered.")
	.add_value(langmuir,"langmuir","Langmuir isotherm described sorption considered")
	.add_value(freundlich,"freundlich","Freundlich isotherm described sorption considered");

using namespace Input::Type;

Record Sorption::input_type
	= Record("Sorptions", "Information about all the limited solubility affected adsorptions.")
	.derive_from( Reaction::input_type )
	.declare_key("solvent_dens", Double(), Default("1.0"),
				"Density of the solvent.")
	.declare_key("substeps", Integer(), Default("100"),
				"Number of equidistant substeps, molar mass and isotherm intersections")
	.declare_key("species", Array(String()), Default::obligatory(),
							"Names of all the sorbing species")
	.declare_key("molar_masses", Array(Double()), Default::obligatory(),
							"Specifies molar masses of all the sorbing species")
	.declare_key("solubility", Array(Double()), Default::obligatory(),
							"Specifies solubility limits of all the sorbing species")
    .declare_key("bulk_data", Array(Sorption::EqData().bulk_input_type()), Default::obligatory(), //
                   	   	   "Containes region specific data necessery to construct isotherms.");

Sorption::EqData::EqData()
: EqDataBase("Sorption")
{
    ADD_FIELD(rock_density, "Rock matrix density.", Input::Type::Default("0.0"));

    ADD_FIELD(sorption_types,"Considered adsorption is described by selected isotherm."); //, Input::Type::Default("none"));
              sorption_types.set_selection(&sorption_type_selection);

    ADD_FIELD(mult_coefs,"Multiplication parameters (k, omega) in either Langmuir c_s = omega * (alpha*c_a)/(1- alpha*c_a) or in linear c_s = k * c_a isothermal description.", Input::Type::Default("1.0"));

    ADD_FIELD(second_params,"Second parameters (alpha, ...) defining isotherm  c_s = omega * (alpha*c_a)/(1- alpha*c_a).", Input::Type::Default("1.0"));

    ADD_FIELD(alpha, "Diffusion coefficient of non-equilibrium linear exchange between mobile and immobile zone (dual porosity)."
            " Vector, one value for every substance.", IT::Default("0"));
}

using namespace std;

Sorption::Sorption(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)//, pScalar mob_porosity, pScalar immob_porosity)
	: Reaction(init_mesh, in_rec, names)
{
	cout << "Sorption constructor is running." << endl;
	TimeGovernor tg(0.0, 1.0);
    nr_of_regions = init_mesh.region_db().bulk_size();
    nr_of_substances = in_rec.val<Input::Array>("species").size();
    nr_of_points = in_rec.val<int>("substeps");

    data_.sorption_types.set_n_comp(nr_of_substances);
    data_.mult_coefs.set_n_comp(nr_of_substances);
    data_.second_params.set_n_comp(nr_of_substances);
    data_.set_mesh(&init_mesh);
    data_.init_from_input( in_rec.val<Input::Array>("bulk_data"), Input::Array());
    data_.set_time(tg);

	//Simple vectors holding  common informations.
	substance_ids.resize(nr_of_substances);
	molar_masses.resize( nr_of_substances );
	c_aq_max.resize( nr_of_substances );

	//isotherms array resized bellow
	isotherms.resize(nr_of_regions);
	for(int i_reg = 0; i_reg < nr_of_regions; i_reg++)
	{
		for(int i_spec = 0; i_spec < nr_of_substances; i_spec++)
		{
			Isotherm iso_mob;
			isotherms[i_reg].push_back(iso_mob);
		}
	}
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

    // Common data for all the isotherms loaded bellow
	solvent_dens = in_rec.val<double>("solvent_dens");

	Input::Array molar_mass_array = in_rec.val<Input::Array>("molar_masses");

	// Molar_masses = in_rec.val<Array(Double())>("molar_masses");
	if (molar_mass_array.size() == molar_masses.size() )   molar_mass_array.copy_to( molar_masses );
	  else  xprintf(UsrErr,"Number of molar masses %d has to match number of adsorbing species %d.\n", molar_mass_array.size(), molar_masses.size());

	Input::Array solub_limit_array = in_rec.val<Input::Array>("solubility");
	if (solub_limit_array.size() == c_aq_max.size() )   solub_limit_array.copy_to( c_aq_max );
	  else  xprintf(UsrErr,"Number of given solubility limits %d has to match number of adsorbing species %d.\n", solub_limit_array.size(), c_aq_max.size());

	Input::Array species_array = in_rec.val<Input::Array>("species");
	unsigned int idx, i_spec = 0;
	for(Input::Iterator<string> spec_iter = species_array.begin<string>(); spec_iter != species_array.end(); ++spec_iter, i_spec++)
	{
		idx = find_subst_name(*spec_iter);
		if ((idx < n_substances()) && (idx >= 0))   substance_ids[i_spec] = idx;
		else	xprintf(Msg,"Wrong name of %d-th adsorbing specie.\n", i_spec);
	}

	// List of types of isotherms in particular regions
	arma::Col<unsigned int> iso_type;
	cout << "there are " << nr_of_substances <<" substances under concideration." << endl;
	iso_type.resize(nr_of_substances);

	// List of sorption parameters
	FieldValue<3>::Vector::return_type mult_param;
	mult_param.resize(nr_of_substances);
	FieldValue<3>::Vector::return_type second_coef;
	second_coef.resize(nr_of_substances);
	double rock_density;

	BOOST_FOREACH(const Region &reg_iter, this->mesh_->region_db().get_region_set("BULK") )
	{
		double por_m, por_imm, phi;

		int reg_idx=reg_iter.bulk_idx();

		// List of types of isotherms in particular regions, initialization
		if(data_.sorption_types.get_const_value(reg_iter, iso_type))
		{
		}else  xprintf(UsrErr,"Type of isotherm must be the same all over the %d-th region, but it is not.", reg_iter.id());

		// Multiplication coefficient parameter follows, initialization
		if(data_.mult_coefs.get_const_value(reg_iter, mult_param)) ;
		  else  xprintf(UsrErr,"All the sorption parameters must be constant all over the %d-th region, but it is not.", reg_iter.id());

		// Second sorption parameter for isothermal description, initialization
		if(data_.second_params.get_const_value(reg_iter, second_coef)) ;
		  else  xprintf(UsrErr,"All the sorption parameters must be constant all over the %d-th region, but it is not.", reg_iter.id());

		for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
		{
			SorptionType hlp_iso_type =  SorptionType(iso_type[i_subst]);

			// Initializes isotherm parameters
			switch(hlp_iso_type)
			{
				case 0: // none:
				{
			 		;
				}
				break;
				case 1: // linear
				{
					isotherms[reg_idx][i_subst].set_mult_coef_(mult_param[i_subst]);
				}
				case 2:  // freundlich
				case 3: // langmuir
				{
					isotherms[reg_idx][i_subst].set_mult_coef_(mult_param[i_subst]);
					isotherms[reg_idx][i_subst].set_second_coef_(second_coef[i_subst]);
				}
				break;
			 	default:
			 	{
			 		xprintf(UsrErr,"1) Sorption of %d-th specie in %d-th region has unknown type %d.", i_subst, reg_idx, hlp_iso_type);
			 	}
			 	break;
			}

			// Creates interpolation tables in the case of constant rock matrix parameters
			if((data_.rock_density.get_const_value(reg_iter, rock_density)) && (this->porosity_->get_const_value(reg_iter, por_m)) && (this->immob_porosity_->get_const_value(reg_iter, por_imm)) && (this->phi_->get_const_value(reg_iter, phi)))
			{
				// here should be some kind of set_scales(..) method instead of following two lines
				double scale_aqua; // = por_m;
				double scale_sorbed; //= phi * (1 - por_m - por_imm) * rock_density * molar_mass;
				this->set_scales(scale_aqua, scale_sorbed, por_m, por_imm, phi, rock_density, molar_masses[i_subst]);

				isotherms[reg_idx][i_subst].reinit(hlp_iso_type, rock_density, solvent_dens, scale_aqua, scale_sorbed, molar_masses[i_subst], c_aq_max[i_subst]);
				switch(hlp_iso_type)
				{
				case 0: // none:
				{
				 	isotherms[reg_idx][i_subst].make_one_point_table();
			 	 }
				 break;
			 	 case 1: //  linear:
			 	 {
				 	Linear obj_isotherm(mult_param[i_subst]);
					isotherms[reg_idx][i_subst].make_table(obj_isotherm, nr_of_points);
			 	 }
			 	 break;
			 	 case 2: // freundlich:
			 	 {
				 	Freundlich obj_isotherm(mult_param[i_subst], second_coef[i_subst]);
					isotherms[reg_idx][i_subst].make_table(obj_isotherm, nr_of_points);
			 	 }
			 	 break;
			 	 case 3: // langmuir:
			 	 {
				 	Langmuir obj_isotherm(mult_param[i_subst], second_coef[i_subst]);
					isotherms[reg_idx][i_subst].make_table(obj_isotherm, nr_of_points);
			 	 }
			 	 break;
			 	 default:
			 	 {
				 	 xprintf(UsrErr,"2) Sorption of %d-th specie in %d-th region has unknown type %d.", i_subst, reg_idx, hlp_iso_type);
			 	 }
			 	 break;
				}
			}else{
				isotherms[reg_idx][i_subst].set_scale_aqua(-1.0);
			}
			isotherms[reg_idx][i_subst].set_sorption_type(hlp_iso_type);
		}
	}
}

// TODO: check duplicity of parents
//       raise warning if sum of ratios is not one

double **Sorption::compute_reaction(double **concentrations, int loc_el) // Sorption simulations are realized just for one element.
{
    ElementFullIter elem = mesh_->element(el_4_loc[loc_el]);
    double porosity;
    double rock_density;
    Region region = elem->region();
    int reg_id_nr = region.bulk_idx();
    int variabl_int = 0;

    if(reg_id_nr != 0) cout << "region id is " << reg_id_nr << endl;

    //  If intersections of isotherm with mass balance lines are known, then interpolate.
    	//  Measurements [c_a,c_s] will be rotated
    	//  Rotated measurements must be projected on rotated isotherm, interpolate_datapoints()
    	//  Projections need to be transformed back to original CS
    //	If intersections are not known then solve the problem analytically (toms748_solve).

    // Constant value of rock density and mobile porosity over the whole region
    if( (data_.rock_density.get_const_value(region, rock_density)) &&  (this->porosity_->get_const_value(region, porosity)) )
	{
		START_TIMER("new-sorption interpolation");
		{
			for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
			{
				int subst_id = substance_ids[i_subst];
			    if( (this->isotherms[reg_id_nr][i_subst].get_sorption_type() > 0) )
			    {
					if((isotherms[reg_id_nr][i_subst].compute_projection((concentration_matrix[subst_id][loc_el]), sorbed_conc_array[i_subst][loc_el]) ) == false)
					{
						cout << "Sorption computed using interpolation failed for substance " << subst_id << " == " << i_subst << endl;
					}
			    }else{
					 ;
			    }
			}
		}
		END_TIMER("new-sorption interpolation");
	}else{
		START_TIMER("new-sorption toms748_solve values-readed");
		if( !(data_.rock_density.get_const_value(region, rock_density)) ) rock_density = data_.rock_density.value(elem->centre(),elem->element_accessor());
		double porosity;
		if( !(this->porosity_->get_const_value(region, porosity)) )
			porosity = this->porosity_->value(elem->centre(),elem->element_accessor());
		variabl_int = 1;
		END_TIMER("new-sorption toms748_solve values-readed");
	}

	if(variabl_int == 1)
	{
		START_TIMER("new-sorption toms748_solve");
		double phi = this->phi_->value(elem->centre(),elem->element_accessor());
		double por_m = this->porosity_->value(elem->centre(),elem->element_accessor());
		double por_imm = this->immob_porosity_->value(elem->centre(),elem->element_accessor());
		double scale_aqua, scale_sorbed;

		for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
		{
		    if(this->isotherms[reg_id_nr][i_subst].get_sorption_type() > 0)
		    {
		    	SorptionType elem_sorp_type = this->isotherms[reg_id_nr][i_subst].get_sorption_type();

				this->set_scales(scale_aqua, scale_sorbed, por_m, por_imm, phi, rock_density, molar_masses[i_subst]);

		    	{
		    		this->isotherms[reg_id_nr][i_subst].reinit(elem_sorp_type, rock_density, solvent_dens, scale_aqua, scale_sorbed, molar_masses[i_subst], c_aq_max[i_subst]);
		    	}
				int subst_id = substance_ids[i_subst];
				switch(elem_sorp_type)
				{
				 case 0:
				 {
					 ;
				 }
				 break;
				 case 1: //  linear:
				 {
					Linear obj_isotherm(isotherms[reg_id_nr][i_subst].get_mult_coef_());
					isotherms[reg_id_nr][i_subst].solve_conc((concentration_matrix[subst_id][loc_el]), sorbed_conc_array[i_subst][loc_el], obj_isotherm);
				 }
				 break;
				 case 2: // freundlich
				 {
					Freundlich obj_isotherm(isotherms[reg_id_nr][i_subst].get_mult_coef_(), isotherms[reg_id_nr][i_subst].get_second_coef_());
					isotherms[reg_id_nr][i_subst].solve_conc((concentration_matrix[subst_id][loc_el]), sorbed_conc_array[i_subst][loc_el], obj_isotherm);
				 }
				 break;
				 case 3:  // langmuir:
				 {
					Langmuir obj_isotherm(isotherms[reg_id_nr][i_subst].get_mult_coef_(), isotherms[reg_id_nr][i_subst].get_second_coef_());
					isotherms[reg_id_nr][i_subst].solve_conc((concentration_matrix[subst_id][loc_el]), sorbed_conc_array[i_subst][loc_el], obj_isotherm);
				 }
				 break;
				 default:
				 {
					 xprintf(UsrErr,"4) Sorption of %d-th specie in %d-th region has unknown type %d.", i_subst, reg_id_nr, isotherms[reg_id_nr][i_subst].get_sorption_type());
				 }
				 break;
				}
		    }
		}
		END_TIMER("new-sorption toms748_solve");
	}

	return concentrations;
}

// Computes sorption simulation over all the elements.
void Sorption::compute_one_step(void)
{
    START_TIMER("new_sorp_step");
	for (int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
	 {
	 	this->compute_reaction(concentration_matrix, loc_el);
	 }
    END_TIMER("new_sorp_step");
	 return;
}


void Sorption::print_sorption_parameters(void)
{

    xprintf(Msg, "\nSorption parameters are defined as follows:\n");
}

void Sorption::set_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc_)
{
	concentration_matrix = ConcentrationMatrix;
	distribution = conc_distr;
	el_4_loc = el_4_loc_;
	return;
}

void Sorption::set_sorb_conc_array(double** sorb_conc_array)
{
	sorbed_conc_array = sorb_conc_array;
	return;
}

void Sorption::set_sorb_conc_array(unsigned int nr_of_local_elm)
{
	this->sorbed_conc_array = new double * [nr_of_substances];
    for (unsigned int sbi = 0; sbi < nr_of_substances; sbi++)
    {
      sorbed_conc_array[sbi] = new double[ nr_of_local_elm ];
      for (unsigned int i = 0; i < nr_of_local_elm; i++)
      {
        sorbed_conc_array[sbi][i] = 0.0;
      }
    }
}

void Sorption::set_porosity(pScalar porosity)
{
	this->porosity_ = porosity;
	return;
}



void Sorption::set_porosity(pScalar porosity, pScalar immob_porosity)
{
	this->porosity_ = porosity;
	this->immob_porosity_ = immob_porosity;
	return;
}

void Sorption::set_phi(pScalar phi)
{
	this->phi_ = phi;
	return;
}

pScalar Sorption::get_phi(void)
{
	return phi_;
}

void Sorption::set_scales(double &scale_aqua, double &scale_sorbed, double por_m, double por_imm, double phi, double rock_density, double molar_mass)
{
	scale_aqua = por_m;

	if(dual_porosity_on) scale_sorbed = phi * (1 - por_m - por_imm) * rock_density * molar_mass;
	  else scale_sorbed = (1 - por_m) * rock_density * molar_mass;

	return;
}

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
	cout << "Meaningless inherited method." << endl;
	return;
}

void Sorption::set_time_step(Input::Record in_rec)
{
	cout << "This method is obsolete for equilibrial sorptions and reactions, but it must be implemented." << endl;
	return;
}
