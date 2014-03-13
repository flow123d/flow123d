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

#include "fields/field_set.hh"
#include "coupling/time_governor.hh"

using namespace std;
namespace it=Input::Type;

it::Selection SorptionBase::EqData::sorption_type_selection = it::Selection("SorptionType")
	.add_value(Isotherm::none,"none", "No adsorption considered")
	.add_value(Isotherm::linear, "linear",
			"Linear isotherm described adsorption considered.")
	.add_value(Isotherm::langmuir, "langmuir",
			"Langmuir isotherm described adsorption considered")
	.add_value(Isotherm::freundlich, "freundlich",
			"Freundlich isotherm described adsorption considered");

using namespace Input::Type;

Record SorptionBase::input_type
	= Record("Sorption", "Information about all the limited solubility affected adsorptions.")
	.derive_from( Reaction::input_type )
	.declare_key("solvent_dens", Double(), Default("1.0"),
				"Density of the solvent.")
	.declare_key("substeps", Integer(), Default("1000"),
				"Number of equidistant substeps, molar mass and isotherm intersections")
	.declare_key("molar_masses", Array(Double()), Default::obligatory(),
							"Specifies molar masses of all the sorbing species")
	.declare_key("solubility", Array(Double(0.0)), Default::optional(), //("-1.0"), //
							"Specifies solubility limits of all the sorbing species")
	.declare_key("table_limits", Array(Double(0.0)), Default::optional(), //("-1.0"), //
							"Specifies highest aqueous concentration in interpolation table.")
    .declare_key("data", Array(SorptionBase::EqData().make_field_descriptor_type("Sorption")), Default::obligatory(), //
                    "Containes region specific data necessary to construct isotherms.")//;
	.declare_key("time", Double(), Default("1.0"),
			"Key called time required by TimeGovernor in Sorption constructor.");/**/

SorptionBase::EqData::EqData()
{
    ADD_FIELD(rock_density, "Rock matrix density.", "0.0");

    ADD_FIELD(sorption_types,"Considered adsorption is described by selected isotherm."); //
              sorption_types.input_selection(&sorption_type_selection);

    ADD_FIELD(mult_coefs,"Multiplication parameters (k, omega) in either Langmuir c_s = omega * (alpha*c_a)/(1- alpha*c_a) or in linear c_s = k * c_a isothermal description.","1.0");

    ADD_FIELD(second_params,"Second parameters (alpha, ...) defining isotherm  c_s = omega * (alpha*c_a)/(1- alpha*c_a).","1.0");
}


SorptionBase::SorptionBase(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)//
	: Reaction(init_mesh, in_rec, names)
{
  DBGMSG("SorptionBase constructor.\n");
  
  nr_of_regions = init_mesh.region_db().bulk_size();
  nr_of_points = in_rec.val<int>("substeps");

  data_.sorption_types.n_comp(n_substances_);
  data_.mult_coefs.n_comp(n_substances_);
  data_.second_params.n_comp(n_substances_);
    
  data_.set_mesh(init_mesh);
  data_.set_input_list( in_rec.val<Input::Array>("data"));
  
  data_+=(data_.porosity
          .name("porosity")
          .units("0")
         );
  
  data_.set_limit_side(LimitSide::right);
  
  //Simple vectors holding  common informations.
  molar_masses.resize( n_substances_ );

  //isotherms array resized bellow
  isotherms.resize(nr_of_regions);
  for(int i_reg = 0; i_reg < nr_of_regions; i_reg++)
    for(int i_spec = 0; i_spec < n_substances_; i_spec++)
    {
      Isotherm iso_mob;
      isotherms[i_reg].push_back(iso_mob);
    }
    
  init_from_input(in_rec);
}

SorptionBase::~SorptionBase(void)
{
}


void SorptionBase::initialize(void )
{
  ASSERT(distribution != nullptr, "Distribution has not been set yet.\n");
  ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  
    //allocating new array for sorbed concentrations
    unsigned int nr_of_local_elm = distribution->lsize();
    sorbed_conc_array = new double * [n_substances_];
    for (unsigned int sbi = 0; sbi < n_substances_; sbi++)
    {
      sorbed_conc_array[sbi] = new double[ nr_of_local_elm ];
      for (unsigned int i = 0; i < nr_of_local_elm; i++)
      {
        sorbed_conc_array[sbi][i] = 0.0;
      }
    }
    
}


void SorptionBase::init_from_input(Input::Record in_rec)
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
		if (solubility_vec_.size() != n_substances_)
		{
			xprintf(UsrErr,"Number of given solubility limits %d has to match number of adsorbing species %d.\n", solubility_vec_.size(), n_substances_);
		}
	}else{
		// fill solubility_vec_ with zeros or resize it at least
		solubility_vec_.resize(n_substances_);
	}

	Input::Iterator<Input::Array> interp_table_limits = in_rec.find<Input::Array>("table_limits");
	if( interp_table_limits )
	{
		interp_table_limits->copy_to(table_limit_);
		if (table_limit_.size() != n_substances_)
		{
			xprintf(UsrErr,"Number of given table limits %d has to match number of adsorbing species %d.\n", table_limit_.size(), n_substances_);
		}/**/
	}else{
		// fill table_limit_ with zeros or resize it at least
		table_limit_.resize(n_substances_);
	}
}


void SorptionBase::update_solution(void)
{
  DBGMSG("SorptionSimple - update_solution\n");
  data_.set_time(*time_); // set to the last computed time
  
  // if parameters changed during last time step, reinit isotherms and eventualy 
  // update interpolation tables in the case of constant rock matrix parameters
  if(data_.changed())
    make_tables();
    

  START_TIMER("Sorption");
  for (int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
  {
    this->compute_reaction(concentration_matrix, loc_el);
  }
  END_TIMER("Sorption");

  return;
}

void SorptionBase::make_tables(void)
{
  ElementAccessor<3> elm;
  BOOST_FOREACH(const Region &reg_iter, this->mesh_->region_db().get_region_set("BULK") )
  {
    int reg_idx = reg_iter.bulk_idx();

    if(data_.is_constant(reg_iter))
    {
      ElementAccessor<3> elm(this->mesh_, reg_iter); // constant element accessor
      isotherm_reinit(isotherms[reg_idx],elm);
      xprintf(MsgDbg,"parameters are constant\n");
      for(int i_subst = 0; i_subst < n_substances_; i_subst++)
      {
        isotherms[reg_idx][i_subst].make_table(nr_of_points);
      }
    }
  }
}

double **SorptionBase::compute_reaction(double **concentrations, int loc_el) // Sorption simulations are realized just for one element.
{
  //DBGMSG("compute_reaction\n");
    ElementFullIter elem = mesh_->element(el_4_loc[loc_el]);
    double porosity;
    double rock_density;
    Region region = elem->region();
    int reg_id_nr = region.bulk_idx();
    int variabl_int = 0;

	std::vector<Isotherm> & isotherms_vec = isotherms[reg_id_nr];

    //if(reg_id_nr != 0) cout << "region id is " << reg_id_nr << endl;
    
    // Constant value of rock density and mobile porosity over the whole region => interpolation_table is precomputed
    if (isotherms_vec[0].is_precomputed()) {
    	for(int i_subst = 0; i_subst < n_substances_; i_subst++)
    	{
    		Isotherm & isotherm = this->isotherms[reg_id_nr][i_subst];
    		int subst_id = substance_id[i_subst];
            isotherm.interpolate((concentration_matrix[subst_id][loc_el]), sorbed_conc_array[i_subst][loc_el]);
    	}
    } else {
		isotherm_reinit(isotherms_vec, elem->element_accessor());
    	for(int i_subst = 0; i_subst < n_substances_; i_subst++)
    	{
            Isotherm & isotherm = this->isotherms[reg_id_nr][i_subst];
            int subst_id = substance_id[i_subst];
            isotherm.compute((concentration_matrix[subst_id][loc_el]), sorbed_conc_array[i_subst][loc_el]);
    	}
    }

	return concentrations;
}


void SorptionBase::print_sorption_parameters(void)
{
    xprintf(Msg, "\nSorption parameters are defined as follows:\n");
}
