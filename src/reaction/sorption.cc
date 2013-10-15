#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <boost/foreach.hpp>

#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "reaction/sorption.hh"
#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "la/distribution.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "mesh/region.hh"
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

    ADD_FIELD(sorption_types,"Considered adsorption is described by selected isotherm."); //
              sorption_types.set_selection(&sorption_type_selection);

    ADD_FIELD(mult_coefs,"Multiplication parameters (k, omega) in either Langmuir c_s = omega * (alpha*c_a)/(1- alpha*c_a) or in linear c_s = k * c_a isothermal description.", Input::Type::Default("1.0"));

    ADD_FIELD(second_params,"Second parameters (alpha, ...) defining isotherm  c_s = omega * (alpha*c_a)/(1- alpha*c_a).", Input::Type::Default("1.0"));

    ADD_FIELD(alphas, "Diffusion coefficient of non-equilibrium linear exchange between mobile and immobile zone (dual porosity)."
            " Vector, one value for every substance.", Input::Type::Default("0"));
}

using namespace std;

Sorption::Sorption(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)//
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
    int nr_transp_subst = names.size();
    data_.alphas.set_n_comp(nr_transp_subst);
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

void Sorption::prepare_inputs(Input::Record in_rec, int porosity_type)
{

    // Common data for all the isotherms loaded bellow
	solvent_dens = in_rec.val<double>("solvent_dens");

	Input::Array molar_mass_array = in_rec.val<Input::Array>("molar_masses");

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
	//cout << "there are " << nr_of_substances <<" substances under concideration." << endl;
	iso_type.resize(nr_of_substances);

	/*arma::Col<double> mult_coef;
	mult_coef.resize(nr_of_substances);

	arma::Col<double> second_coef;
	second_coef.resize(nr_of_substances);*/

	// List of sorption parameters
	FieldValue<3>::Vector::return_type mult_param;
	mult_param.resize(nr_of_substances);
	FieldValue<3>::Vector::return_type second_coef;
	second_coef.resize(nr_of_substances);
	// Mass transfer coeffs between mobile and immobile pores
	FieldValue<3>::Vector::return_type mass_transfer_coeffs;
	mass_transfer_coeffs.resize(nr_of_substances);

	double rock_density;

	BOOST_FOREACH(const Region &reg_iter, this->mesh_->region_db().get_region_set("BULK") )
	{
		double por_m, por_imm, phi; //, mult_coef, second_coef;

		int reg_idx = reg_iter.bulk_idx();

		ElementAccessor<3> elm;

		// List of types of isotherms in particular regions, initialization
		if(data_.sorption_types.get_const_accessor(reg_iter, elm)) ;
		  else  xprintf(UsrErr,"Type of isotherm must be the same all over the %d-th region, but it is not.", reg_iter.id());


		// Creates interpolation tables in the case of constant rock matrix parameters
		if((data_.rock_density.get_const_accessor(reg_iter, elm)) &&
				(data_.mult_coefs.get_const_accessor(reg_iter, elm)) &&
				(data_.second_params.get_const_accessor(reg_iter, elm)) &&
				(this->porosity_->get_const_accessor(reg_iter, elm)) &&
				(this->immob_porosity_->get_const_accessor(reg_iter, elm)) &&
				(this->phi_->get_const_accessor(reg_iter, elm)))
		{
		    // TODO: call isotherm_reinit(isotherms[reg_idx], elm) here to get constant values and call reinit for isotherms

			for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
			{
				SorptionType hlp_iso_type =  SorptionType(iso_type[i_subst]);
				Isotherm & isotherm = this->isotherms[reg_idx][i_subst];

				//scales for the case of sorption in mobile pores follows
				double scale_aqua = por_m;
				double scale_sorbed;
				if((scale_sorbed = phi * (1 - por_m - por_imm) * rock_density * molar_masses[i_subst]) == 0.0)
					xprintf(UsrErr, "Sorption::prepare_inputs() failed. Parameter scale_sorbed (phi * (1 - por_m - por_imm) * rock_density * molar_masses[i_subst]) is equal to zero.");

				isotherm.reinit(hlp_iso_type, solvent_dens, scale_aqua, scale_sorbed, c_aq_max[i_subst], mult_param[i_subst], second_coef[i_subst]); // hlp_iso_type, rock_density, solvent_dens, por_m, por_imm, phi, molar_masses[i_subst], c_aq_max[i_subst]);
				isotherm.make_table(nr_of_points);
			}
		}else{
			for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
			{
				SorptionType hlp_iso_type =  SorptionType(iso_type[i_subst]);
				Isotherm & isotherm = this->isotherms[reg_idx][i_subst];

				// TODO: should do nothing herre

				// We need to get coeficients belonging to particular isotherm. Solution realized bellow is just temporary. If it would be final, then the initialization could preceed if-condition above.
				//data_.mult_coefs.get_const_accessor(reg_iter, mult_param);
				//data_.second_params.get_const_accessor(reg_iter, second_coef);
				xprintf(Msg,"Sorption::prepare_inputs(), mult_coef_ %f, second_coef_ %f\n", mult_param[i_subst], second_coef[i_subst]);
				isotherm.set_iso_params(hlp_iso_type,mult_param[i_subst],second_coef[i_subst]);
			}
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
    if( this->isotherms[reg_id_nr][0].is_precomputed() ) //(data_.rock_density.get_const_accessor(region, rock_density)) &&  (this->porosity_->get_const_accessor(region, porosity)) )
	{
		START_TIMER("new-sorption interpolation");
		{
			for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
			{
				Isotherm & isotherm = this->isotherms[reg_id_nr][i_subst];
				int subst_id = substance_ids[i_subst];
			    isotherm.compute_projection((concentration_matrix[subst_id][loc_el]), sorbed_conc_array[i_subst][loc_el]);
			}
		}
		END_TIMER("new-sorption interpolation");
	}else{
		START_TIMER("new-sorption toms748_solve values-readed");
		//if( !(data_.rock_density.get_const_accessor(region, rock_density)) )
			rock_density = data_.rock_density.value(elem->centre(),elem->element_accessor());
		double porosity;
		//if( !(this->porosity_->get_const_accessor(region, porosity)) )
			porosity = this->porosity_->value(elem->centre(),elem->element_accessor());

		double phi = this->phi_->value(elem->centre(),elem->element_accessor());
		double por_m = this->porosity_->value(elem->centre(),elem->element_accessor());
		double por_imm = this->immob_porosity_->value(elem->centre(),elem->element_accessor());
		END_TIMER("new-sorption toms748_solve values-readed");

		double scale_aqua = por_m;

		for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
		{
			double scale_sorbed;
			// following condition is valid for adsorption in mobile pores
			if((scale_sorbed = phi * (1 - por_m - por_imm) * rock_density * molar_masses[i_subst]) == 0.0)
								xprintf(UsrErr, "Sorption::compute_reaction() failed. Parameter scale_sorbed (phi * (1 - por_m - por_imm) * rock_density * molar_masses[i_subst]) is equal to zero.");

			Isotherm & isotherm = this->isotherms[reg_id_nr][i_subst];
		    SorptionType elem_sorp_type = isotherm.adsorption_type_; //this->isotherms[reg_id_nr][i_subst].get_sorption_type();

		    // Coeficients describing isoterms must be constant at the moment I think.
		    double mult_coef = isotherm.mult_coef_;
		    double second_coef = isotherm.second_coef_;
		    //xprintf(Msg,"coeficients %f %f\n", mult_coef, second_coef);

			isotherm.reinit(elem_sorp_type, solvent_dens, scale_aqua, scale_sorbed, c_aq_max[i_subst], mult_coef, second_coef);
		    int subst_id = substance_ids[i_subst];

		    ConcPair conc(concentration_matrix[subst_id][loc_el], sorbed_conc_array[i_subst][loc_el]);
			START_TIMER("new-sorption toms748_solve");
		    conc = isotherm.solve_conc(conc);
			END_TIMER("new-sorption toms748_solve");
		    concentration_matrix[subst_id][loc_el] = conc.first;
		    sorbed_conc_array[i_subst][loc_el] = conc.second;
		}
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

void Sorption::set_sorb_conc_array(unsigned int nr_of_local_elm) // could be transposed to optimize computation speed
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

void Sorption::set_immob_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc_)
{
	immob_concentration_matrix = ConcentrationMatrix;
	distribution = conc_distr;
	el_4_loc = el_4_loc_;
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
	phi_ = phi;
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
