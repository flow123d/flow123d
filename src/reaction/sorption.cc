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
#include "mesh/region.hh"
#include "input/type_selection.hh"

#include "coupling/time_governor.hh"

const double pi = 3.1415;


using namespace std;

SorptionSimple::SorptionSimple(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)//
	: SorptionBase(init_mesh, in_rec, names)
{
	/*cout << "Sorption constructor is running." << endl;
	TimeGovernor tg(0.0, 1.0);
    nr_of_regions = init_mesh.region_db().bulk_size();
    nr_of_substances = in_rec.val<Input::Array>("species").size();
    nr_of_points = in_rec.val<int>("substeps");

    data_.sorption_types.set_n_comp(nr_of_substances);
    data_.mult_coefs.set_n_comp(nr_of_substances);
    data_.second_params.set_n_comp(nr_of_substances);
    int nr_transp_subst = names.size();
    //data_.alphas.set_n_comp(nr_transp_subst);
    data_.set_mesh(&init_mesh);
    data_.init_from_input( in_rec.val<Input::Array>("bulk_data"), Input::Array());

	//Simple vectors holding  common informations.
	substance_ids.resize(nr_of_substances);
	molar_masses.resize( nr_of_substances );

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

    time_ = new TimeGovernor(in_rec.val<double>("time"), TimeGovernor::marks().type_fixed_time());/**/
}

SorptionSimple::~SorptionSimple(void)
{
}

/* //implemented in base
void SorptionSimple::init_from_input(Input::Record in_rec)
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
}
//*/

void SorptionSimple::make_tables(void)
{
  DBGMSG("make_tables\n");
  ElementAccessor<3> elm;

  BOOST_FOREACH(const Region &reg_iter, this->mesh_->region_db().get_region_set("BULK") )
	{
		int reg_idx = reg_iter.bulk_idx();
                
		// Creates interpolation tables in the case of constant rock matrix parameters
		if((data_.rock_density.get_const_accessor(reg_iter, elm)) &&
				(data_.mult_coefs.get_const_accessor(reg_iter, elm)) &&
				(data_.second_params.get_const_accessor(reg_iter, elm)) &&
				(this->porosity_->get_const_accessor(reg_iter, elm))) /* &&
				(this->immob_porosity_->get_const_accessor(reg_iter, elm)) &&
				(this->phi_->get_const_accessor(reg_iter, elm)))/**/
		{
			isotherm_reinit(isotherms[reg_idx],elm);
			xprintf(MsgDbg,"parameters are constant\n");
			for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
			{
				isotherms[reg_idx][i_subst].make_table(nr_of_points);
			}
		}
	}
}/**/

void SorptionSimple::isotherm_reinit(std::vector<Isotherm> &isotherms_vec, const ElementAccessor<3> &elem)
{
	START_TIMER("SorptionSimple::isotherm_reinit");

	double rock_density = data_.rock_density.value(elem.centre(),elem);
	//double porosity = this->porosity_->value(elem.centre(),elem);

	//double phi = this->phi_->value(elem.centre(),elem);
	double por_m = this->porosity_->value(elem.centre(),elem);
	//double por_imm = this->immob_porosity_->value(elem.centre(),elem);

	// List of types of isotherms in particular regions
	arma::uvec adsorption_type = data_.sorption_types.value(elem.centre(),elem);
	arma::Col<double> mult_coef_vec = data_.mult_coefs.value(elem.centre(),elem);
	arma::Col<double> second_coef_vec = data_.second_params.value(elem.centre(),elem);

	for(int i_subst = 0; i_subst < nr_of_substances; i_subst++)
	{
		double mult_coef = mult_coef_vec[i_subst];
		double second_coef = second_coef_vec[i_subst];
		Isotherm & isotherm = isotherms_vec[i_subst];

		//scales are different for the case of sorption in mobile and immobile pores
		double scale_aqua = por_m, 
                       scale_sorbed = scale_sorbed = (1 - por_m) * rock_density * molar_masses[i_subst];

		if ( scale_sorbed == 0.0)
			xprintf(UsrErr, "Sorption::init_from_input() failed. Parameter scale_sorbed (phi * (1 - por_m - por_imm) * rock_density * molar_masses[i_subst]) is equal to zero.");
		bool limited_solubility_on = false;
		double table_limit;
		if (solubility_vec_[i_subst] <= 0.0) {
			limited_solubility_on = false;
			table_limit=table_limit_[i_subst];
		} else {
			limited_solubility_on = true;
			table_limit=solubility_vec_[i_subst];
		}
		
//                 DBGMSG("limited_solubility_on: %d\n",limited_solubility_on);
//                 DBGMSG("solvent_dens: %f\n",solvent_dens);
//                 DBGMSG("scale_aqua: %f\n",scale_aqua);
//                 DBGMSG("scale_sorbed: %f\n",scale_sorbed);
//                 DBGMSG("table_limit: %f\n",table_limit);
//                 DBGMSG("mult_coef: %f\n",mult_coef);
//                 DBGMSG("second_coef: %f\n",second_coef);
		isotherm.reinit(Isotherm::SorptionType(adsorption_type[i_subst]), limited_solubility_on,
					solvent_dens, scale_aqua, scale_sorbed, table_limit, mult_coef, second_coef);

	}

	END_TIMER("SorptionSimple::isotherm_reinit");

	return;
}

/*/ TODO: check duplicity of parents
//       raise warning if sum of ratios is not one

double **SorptionSimple::compute_reaction(double **concentrations, int loc_el) // Sorption simulations are realized just for one element.
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
void SorptionSimple::update_solution(void)
{
  DBGMSG("update_solution\n");
  data_.set_time(*time_); // set to the last computed time
  
  DBGMSG("update_solution - data-set_time\n");
  DBGMSG("rock_density name: %s\n",data_.rock_density.name().c_str());
  //if parameters changed during last time step, reinit isotherms and eventualy update interpolation tables in the case of constant rock matrix parameters
  if( (data_.rock_density.changed_during_set_time) &&
      (data_.mult_coefs.changed_during_set_time) &&
      (data_.second_params.changed_during_set_time) &&
      (this->porosity_->changed_during_set_time)) /* &&
		(this->immob_porosity_->changed_during_set_time) &&
		(this->phi_->changed_during_set_time))/**/
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


/*void SorptionSimple::print_sorption_parameters(void)
{
    xprintf(Msg, "\nSorption parameters are defined as follows:\n");
}

void SorptionSimple::set_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc_)
{
	concentration_matrix = ConcentrationMatrix;
	distribution = conc_distr;
	el_4_loc = el_4_loc_;
	return;
}

void SorptionSimple::set_sorb_conc_array(double** sorb_conc_array)
{
	sorbed_conc_array = sorb_conc_array;
	return;
}

void SorptionSimple::set_sorb_conc_array(unsigned int nr_of_local_elm) // could be transposed to optimize computation speed
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

void SorptionSimple::set_immob_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc_)
{
	immob_concentration_matrix = ConcentrationMatrix;
	distribution = conc_distr;
	el_4_loc = el_4_loc_;
	return;
}/**/


/*void SorptionSimple::set_phi(pScalar phi)
{
	phi_ = phi;
	return;
}




void Sorption::choose_next_time(void)
{
	//cout << "2) Meaningless inherited method." << endl;
	return;
}

void SorptionSimple::set_time_step_constrain(double dt)
{
	//cout << "3) Meaningless inherited method." << endl;
	return;
}

void SorptionSimple::get_parallel_solution_vector(Vec &vc)
{
	//cout << "4) Meaningless inherited method." << endl;
	return;
}

void SorptionSimple::get_solution_vector(double* &vector, unsigned int &size)
{
	//cout << "5) Meaningless inherited method." << endl;
	return;
}

void SorptionSimple::set_time_step(double new_timestep)
{
	//cout << "6) Meaningless inherited method." << endl;
	return;
}

void SorptionSimple::set_time_step(Input::Record in_rec)
{
	//cout << "This method is obsolete for equilibrial sorptions and reactions, but it must be implemented." << endl;
	return;
<<<<<<< HEAD
}/**/

void SorptionSimple::set_concentration_vector(Vec &vc)
{
        //cout << "7) Meaningless inherited method." << endl;
        return;
}
