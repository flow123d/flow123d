#include <boost/foreach.hpp>

#include "reaction/reaction.hh"
#include "reaction/isotherm.hh"
#include "reaction/sorption_immob.hh"

#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "mesh/mesh.h"

using namespace std;


IT::Record SorptionImmob::input_type
	= IT::Record("SorptionImmobile", "Information about all the limited solubility affected adsorptions.")
	.derive_from( ReactionTerm::input_type )
	.copy_keys(SorptionBase::input_type)
	.declare_key("output_fields", IT::Array(make_output_selection("conc_immobile_solid", "SorptionImmobile_Output")),
            IT::Default("conc_immobile_solid"), "List of fields to write to output stream.");



SorptionImmob::SorptionImmob(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)
	: SorptionDual(init_mesh, in_rec, names)
{
  //DBGMSG("SorptionImmob constructor.\n");
	data_ = new EqData("conc_immobile_solid");
	output_selection = make_output_selection("conc_immobile_solid", "SorptionImmobile_Output");
}

SorptionImmob::~SorptionImmob(void)
{
}

/*
double SorptionImmob::compute_sorbing_scale(double por_m, double por_imm)
{
  double phi = por_imm / (por_m + por_imm);
  return phi;
}
*/

void SorptionImmob::isotherm_reinit(std::vector<Isotherm> &isotherms_vec, const ElementAccessor<3> &elem)
{
	START_TIMER("SorptionImmob::isotherm_reinit");

	double rock_density = data_->rock_density.value(elem.centre(),elem);
	
	double por_m = data_->porosity.value(elem.centre(),elem);
	double por_imm = immob_porosity_.value(elem.centre(),elem);
        //double phi = por_imm / (por_m + por_imm);     // = 1-phi
        double phi = por_m/(por_m + por_imm);
        
	// List of types of isotherms in particular regions
	arma::uvec adsorption_type = data_->adsorption_type.value(elem.centre(),elem);
	arma::Col<double> mult_coef_vec = data_->isotherm_mult.value(elem.centre(),elem);
	arma::Col<double> second_coef_vec = data_->isotherm_other.value(elem.centre(),elem);

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
			xprintf(UsrErr, "Parameter scale_sorbed ((1 - phi) * (1 - por_m - por_imm) * rock_density * molar_masses[i_subst]) is equal to zero.");
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
					solvent_density, scale_aqua, scale_sorbed, table_limit, mult_coef, second_coef);

	}

	END_TIMER("SorptionImmob::isotherm_reinit");
}



