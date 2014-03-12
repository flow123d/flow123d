#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <boost/foreach.hpp>

#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
//#include "reaction/isotherm.hh"
#include "reaction/sorption_immob.hh"
#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "la/distribution.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "mesh/region.hh"
#include "input/type_selection.hh"

#include "coupling/time_governor.hh"

const double pi = 3.1415;
/*namespace it=Input::Type;

/*it::Selection SorptionImmob::EqData::sorption_type_selection = it::Selection("SorptionType")
	.add_value(Isotherm::none,"none", "No adsorption considered")
	.add_value(Isotherm::linear, "linear",
			"Linear isotherm described adsorption considered.")
	.add_value(Isotherm::langmuir, "langmuir",
			"Langmuir isotherm described adsorption considered")
	.add_value(Isotherm::freundlich, "freundlich",
			"Freundlich isotherm described adsorption considered");/**/

//using namespace Input::Type;
/*Record SorptionImmob::input_type
	= Record("SorptionsImmob", "Information about all the limited solubility affected adsorptions.")
	.derive_from( Reaction::input_type )
	.declare_key("solvent_dens", Double(), Default("1.0"),
				"Density of the solvent.")
	.declare_key("substeps", Integer(), Default("1000"),
				"Number of equidistant substeps, molar mass and isotherm intersections")
	.declare_key("species", Array(String()), Default::obligatory(),
							"Names of all the adsorbing species")
	.declare_key("molar_masses", Array(Double()), Default::obligatory(),
							"Specifies molar masses of all the sorbing species")
	.declare_key("solubility", Array(Double(0.0)), Default::optional(), //("-1.0"), //
							"Specifies solubility limits of all the sorbing species")
	.declare_key("table_limits", Array(Double(0.0)), Default::optional(), //("-1.0"), //
							"Specifies highest aqueous concentration in interpolation table.")
    .declare_key("bulk_data", Array(Sorption::EqData().bulk_input_type()), Default::obligatory(), //
                   	   	   "Containes region specific data necessery to construct isotherms.")//;
	.declare_key("time", Double(), Default("1.0"),
			"Key called time required by TimeGovernor in Sorption constructor.");/**/

/*SorptionImmob::EqData::EqData()
: EqDataBase("SorptionImmob")
{
    /*ADD_FIELD(rock_density, "Rock matrix density.", Input::Type::Default("0.0"));

    ADD_FIELD(sorption_types,"Considered adsorption is described by selected isotherm."); //
              sorption_types.set_selection(&sorption_type_selection);

    ADD_FIELD(mult_coefs,"Multiplication parameters (k, omega) in either Langmuir c_s = omega * (alpha*c_a)/(1- alpha*c_a) or in linear c_s = k * c_a isothermal description.", Input::Type::Default("1.0"));

    ADD_FIELD(second_params,"Second parameters (alpha, ...) defining isotherm  c_s = omega * (alpha*c_a)/(1- alpha*c_a).", Input::Type::Default("1.0"));

    //ADD_FIELD(alphas, "Diffusion coefficient of non-equilibrium linear exchange between mobile and immobile zone (dual porosity)."
            //" Vector, one value for every substance.", Input::Type::Default("0"));
}/**/

using namespace std;

SorptionImmob::SorptionImmob(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)//
	: SorptionDual(init_mesh, in_rec, names)
{
	cout << "SorptionImmob constructor is running." << endl;
}

SorptionImmob::~SorptionImmob(void)
{
}

/*void Sorption::init_from_input(Input::Array bulk_list)
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

	Input::Iterator<Input::Array> solub_iter = in_rec.find<Input::Array>("solubility");
	if( solub_iter )
	{
		solub_iter->copy_to(solubility_vec_);
		if (solubility_vec_.size() != nr_of_substances)
		{
			xprintf(UsrErr,"Number of given solubility limits %d has to match number of adsorbing species %d.\n", solubility_vec_.size(), nr_of_substances);
		}
	}else{
		// fill solubility_vec_ with zeros or resize it at least
		solubility_vec_.resize(nr_of_substances);
	}

	Input::Iterator<Input::Array> interp_table_limits = in_rec.find<Input::Array>("table_limits");
	if( interp_table_limits )
	{
		interp_table_limits->copy_to(table_limit_);
		if (table_limit_.size() != nr_of_substances)
		{
			xprintf(UsrErr,"Number of given table limits %d has to match number of adsorbing species %d.\n", table_limit_.size(), nr_of_substances);
		}
	}else{
		// fill table_limit_ with zeros or resize it at least
		table_limit_.resize(nr_of_substances);
	}

	Input::Array species_array = in_rec.val<Input::Array>("species");
	unsigned int idx, i_spec = 0;
	for(Input::Iterator<string> spec_iter = species_array.begin<string>(); spec_iter != species_array.end(); ++spec_iter, i_spec++)
	{
		idx = find_subst_name(*spec_iter);
		if ((idx < n_substances()) && (idx >= 0)) substance_ids[i_spec] = idx;
		else	xprintf(UsrErr,"Wrong name of %d-th adsorbing specie.\n", i_spec);
	}

	make_tables();
}/**/

void SorptionImmob::make_tables(void)
{
	/*ElementAccessor<3> elm;

	BOOST_FOREACH(const Region &reg_iter, this->mesh_->region_db().get_region_set("BULK") )
	{
		int reg_idx = reg_iter.bulk_idx();

		// Creates interpolation tables in the case of constant rock matrix parameters
		if((data_.rock_density.get_const_accessor(reg_iter, elm)) &&
				(data_.mult_coefs.get_const_accessor(reg_iter, elm)) &&
				(data_.second_params.get_const_accessor(reg_iter, elm)) &&
				(this->porosity_->get_const_accessor(reg_iter, elm)) &&
				(this->immob_porosity_->get_const_accessor(reg_iter, elm)) &&
				(this->phi_->get_const_accessor(reg_iter, elm)))
		{
			isotherm_reinit(isotherms[reg_idx],elm);
			xprintf(Msg,"parameters are constant\n");
			for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
			{
				isotherms[reg_idx][i_subst].make_table(nr_of_points);
			}
		}
	}/**/
	return;
}/**/

void SorptionImmob::isotherm_reinit(std::vector<Isotherm> &isotherms_vec, const ElementAccessor<3> &elem)
{
	START_TIMER("SorptionImmob::isotherm_reinit");

	const double &rock_density = data_.rock_density.value(elem.centre(),elem);
	//double porosity = this->porosity_->value(elem.centre(),elem);

	double phi = this->phi_->value(elem.centre(),elem);
	double por_m = this->porosity_->value(elem.centre(),elem);
	double por_imm = this->immob_porosity_->value(elem.centre(),elem);

	// List of types of isotherms in particular regions
	arma::uvec adsorption_type = data_.sorption_types.value(elem.centre(),elem);
	arma::Col<double> mult_coef_vec = data_.mult_coefs.value(elem.centre(),elem);
	arma::Col<double> second_coef_vec = data_.second_params.value(elem.centre(),elem);

	for(int i_subst = 0; i_subst < n_substances_; i_subst++)
	{
		double mult_coef = mult_coef_vec[i_subst];
		double second_coef = second_coef_vec[i_subst];
		Isotherm & isotherm = isotherms_vec[i_subst];

		//scales are different for the case of sorption in mobile and immobile pores
		double scale_aqua, scale_sorbed;
		scale_aqua = por_imm;
		scale_sorbed = (1 - phi) * (1 - por_m - por_imm) * rock_density * molar_masses[i_subst];
		if(scale_sorbed == 0.0)
			xprintf(UsrErr, "Sorption::prepare_inputs() failed. Parameter scale_sorbed ((1 - phi) * (1 - por_m - por_imm) * rock_density * molar_masses[i_subst]) is equal to zero.");
		bool limited_solubility_on;
		double table_limit;
		if (solubility_vec_[i_subst] <= 0.0) {
			limited_solubility_on = false;
			table_limit=table_limit_[i_subst];

		} else {
			limited_solubility_on = true;
			table_limit=solubility_vec_[i_subst];
		}
		isotherm.reinit(Isotherm::SorptionType(adsorption_type[i_subst]), limited_solubility_on,
					solvent_dens, scale_aqua, scale_sorbed, table_limit, mult_coef, second_coef);

	}

	END_TIMER("SorptionImmob::isotherm_reinit");

	return;
}

// TODO: check duplicity of parents
//       raise warning if sum of ratios is not one

/*double **Sorption::compute_reaction(double **concentrations, int loc_el) // Sorption simulations are realized just for one element.
{
    ElementFullIter elem = mesh_->element(el_4_loc[loc_el]);
    double porosity;
    double rock_density;
    Region region = elem->region();
    int reg_id_nr = region.bulk_idx();
    int variabl_int = 0;

	std::vector<Isotherm> & isotherms_vec = isotherms[reg_id_nr];

    if(reg_id_nr != 0) cout << "region id is " << reg_id_nr << endl;

    // Constant value of rock density and mobile porosity over the whole region => interpolation_table is precomputed
    if (isotherms_vec[0].is_precomputed()) {
    	for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
    	{
    		Isotherm & isotherm = this->isotherms[reg_id_nr][i_subst];
    		int subst_id = substance_ids[i_subst];
            isotherm.interpolate((concentration_matrix[subst_id][loc_el]), sorbed_conc_array[i_subst][loc_el]);
    	}
    } else {
		isotherm_reinit(isotherms_vec, elem->element_accessor());
    	for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
    	{
    		Isotherm & isotherm = this->isotherms[reg_id_nr][i_subst];
    		int subst_id = substance_ids[i_subst];
            isotherm.compute((concentration_matrix[subst_id][loc_el]), sorbed_conc_array[i_subst][loc_el]);
    	}
    }

	return concentrations;
}/**/

// Computes adsorption simulation over all the elements.
void SorptionImmob::update_solution(void)
{
    data_.set_time(*time_); // set to the last computed time
    //if parameters changed during last time step, reinit isotherms and eventualy update interpolation tables in the case of constant rock matrix parameters
	if((data_.rock_density.changed_during_set_time) &&
		(data_.mult_coefs.changed_during_set_time) &&
		(data_.second_params.changed_during_set_time) &&
		(this->porosity_->changed_during_set_time) &&
		(this->immob_porosity_->changed_during_set_time) &&
		(this->phi_->changed_during_set_time))
	{
		make_tables();
	}

    START_TIMER("Computes reaction");
	for (int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
	 {
	 	this->compute_reaction(concentration_matrix, loc_el);
	 }
    END_TIMER("Computes reaction");

	return;
}


/*void Sorption::print_sorption_parameters(void)
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
}/**/

void SorptionImmob::set_immob_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc_)
{
	immob_concentration_matrix = ConcentrationMatrix;
	distribution = conc_distr;
	el_4_loc = el_4_loc_;
	return;
}

/*void Sorption::set_porosity(pScalar porosity, pScalar immob_porosity)
{
	this->porosity_ = porosity;
	this->immob_porosity_ = immob_porosity;
	return;
}
void Sorption::set_phi(pScalar phi)
{
	phi_ = phi;
	return;
}/**/

/*void Sorption::update_solution(void)
{
	//cout << "1) Meaningless inherited method." << endl;
	return;
}
void Sorption::choose_next_time(void)
{
	//cout << "2) Meaningless inherited method." << endl;
	return;
}

void Sorption::set_time_step_constrain(double dt)
{
	//cout << "3) Meaningless inherited method." << endl;
	return;
}

void Sorption::get_parallel_solution_vector(Vec &vc)
{
	//cout << "4) Meaningless inherited method." << endl;
	return;
}

void Sorption::get_solution_vector(double* &vector, unsigned int &size)
{
	//cout << "5) Meaningless inherited method." << endl;
	return;
}

void Sorption::set_time_step(double new_timestep)
{
	//cout << "6) Meaningless inherited method." << endl;
	return;
}

/*void Sorption::set_time_step(Input::Record in_rec)
{
	//cout << "This method is obsolete for equilibrial sorptions and reactions, but it must be implemented." << endl;
	return;
}/**/
