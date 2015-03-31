#include <vector>
#include <limits>

#include "reaction/isotherm.hh"
#include "reaction/sorption.hh"
#include "system/sys_profiler.hh"
#include "mesh/accessors.hh"
#include "input/factory.hh"

/*********************************                 *********************************************************/
/********************************* SORPTION_SIMPLE *********************************************************/
/*********************************                 *********************************************************/

IT::Record SorptionSimple::input_type = SorptionBase::record_factory(SorptionRecord::simple);

SorptionSimple::SorptionSimple(Mesh &init_mesh, Input::Record in_rec)
  : SorptionBase(init_mesh, in_rec)
{
	data_ = new EqData("conc_solid");
    this->eq_data_ = data_;
	output_selection = make_output_selection("conc_solid", "SorptionSimple_Output");
}

const int SorptionSimple::registrar =
		Input::register_class< SorptionSimple, Mesh &, Input::Record >("Sorption");

SorptionSimple::~SorptionSimple(void)
{}

void SorptionSimple::isotherm_reinit(std::vector<Isotherm> &isotherms_vec, const ElementAccessor<3> &elem)
{
	START_TIMER("SorptionSimple::isotherm_reinit");

	double rock_density = data_->rock_density.value(elem.centre(),elem);
	double por_m = data_->porosity.value(elem.centre(),elem);

	// List of types of isotherms in particular regions
	arma::uvec adsorption_type = data_->sorption_type.value(elem.centre(),elem);
	arma::Col<double> mult_coef_vec = data_->isotherm_mult.value(elem.centre(),elem);
	arma::Col<double> second_coef_vec = data_->isotherm_other.value(elem.centre(),elem);

	for(unsigned int i_subst = 0; i_subst < n_substances_; i_subst++)
	{
		double mult_coef = mult_coef_vec[i_subst];
		double second_coef = second_coef_vec[i_subst];
		Isotherm & isotherm = isotherms_vec[i_subst];

		//scales are different for the case of sorption in mobile and immobile pores
		double scale_aqua = por_m,
               scale_sorbed = (1 - por_m) * rock_density * substances_[substance_global_idx_[i_subst]].molar_mass();

		bool limited_solubility_on = false;
		double table_limit;
		if (solubility_vec_[i_subst] <= 0.0) {
			limited_solubility_on = false;
			table_limit=table_limit_[i_subst];
		} else {
			limited_solubility_on = true;
			table_limit=solubility_vec_[i_subst];
		}
		
		if( (1-por_m) <= std::numeric_limits<double>::epsilon()) //means there is no sorbing surface
        {
            isotherm.reinit(Isotherm::none, false, solvent_density_, scale_aqua, scale_sorbed,table_limit,0,0);
            continue;
        }

        if ( scale_sorbed <= 0.0)
            xprintf(UsrErr, "Scaling parameter in sorption is not positive. Check the input for rock density and molar mass of %d. substance.",i_subst);
        
		isotherm.reinit(Isotherm::SorptionType(adsorption_type[i_subst]), limited_solubility_on,
					solvent_density_, scale_aqua, scale_sorbed, table_limit, mult_coef, second_coef);

	}

	END_TIMER("SorptionSimple::isotherm_reinit");
}


/***********************************               *********************************************************/
/*********************************** SORPTION_DUAL *********************************************************/
/***********************************               *********************************************************/

SorptionDual::SorptionDual(Mesh &init_mesh, Input::Record in_rec,
                           const string &output_conc_name,
                           const string &output_selection_name)
    : SorptionBase(init_mesh, in_rec)
{
    data_ = new EqData(output_conc_name);
    *data_+=immob_porosity_
        .flags_add(FieldFlag::input_copy)
        .name("porosity_immobile");
    this->eq_data_ = data_;
    output_selection = make_output_selection(output_conc_name, output_selection_name);
}

SorptionDual::~SorptionDual(void)
{}

/**********************************                  *******************************************************/
/*********************************** SORPTION_MOBILE *******************************************************/
/**********************************                  *******************************************************/

IT::Record SorptionMob::input_type = SorptionBase::record_factory(SorptionRecord::mobile);


const int SorptionMob::registrar =
		Input::register_class< SorptionMob, Mesh &, Input::Record >("SorptionMobile");


SorptionMob::SorptionMob(Mesh &init_mesh, Input::Record in_rec)
    : SorptionDual(init_mesh, in_rec, "conc_solid", "SorptionMobile_Output")
{}


SorptionMob::~SorptionMob(void)
{}

/*
double SorptionMob::compute_sorbing_scale(double por_m, double por_imm)
{
  double phi = por_m/(por_m + por_imm);
  return phi;
}
*/

void SorptionMob::isotherm_reinit(std::vector<Isotherm> &isotherms_vec, const ElementAccessor<3> &elem)
{
    START_TIMER("SorptionMob::isotherm_reinit");

    double rock_density = data_->rock_density.value(elem.centre(),elem);

    double por_m = data_->porosity.value(elem.centre(),elem);
    double por_imm = immob_porosity_.value(elem.centre(),elem);
    double phi = por_m/(por_m + por_imm);

    // List of types of isotherms in particular regions
    arma::uvec adsorption_type = data_->sorption_type.value(elem.centre(),elem);
    arma::Col<double> mult_coef_vec = data_->isotherm_mult.value(elem.centre(),elem);
    arma::Col<double> second_coef_vec = data_->isotherm_other.value(elem.centre(),elem);

    for(unsigned int i_subst = 0; i_subst < n_substances_; i_subst++)
    {
        double mult_coef = mult_coef_vec[i_subst];
        double second_coef = second_coef_vec[i_subst];
        Isotherm & isotherm = isotherms_vec[i_subst];

        //scales are different for the case of sorption in mobile and immobile pores
        double scale_aqua = por_m,
               scale_sorbed = phi * (1 - por_m - por_imm) * rock_density * substances_[substance_global_idx_[i_subst]].molar_mass();

        bool limited_solubility_on;
        double table_limit;
        if (solubility_vec_[i_subst] <= 0.0) {
                limited_solubility_on = false;
                table_limit=table_limit_[i_subst];

        } else {
                limited_solubility_on = true;
                table_limit=solubility_vec_[i_subst];
        }
        
        if( (1-por_m-por_imm) <= std::numeric_limits<double>::epsilon()) //means there is no sorbing surface
        {
            isotherm.reinit(Isotherm::none, false, solvent_density_, scale_aqua, scale_sorbed,table_limit,0,0);
            continue;
        }
        
        if ( scale_sorbed <= 0.0)
            xprintf(UsrErr, "Scaling parameter in sorption is not positive. Check the input for rock density and molar mass of %d. substance.",i_subst);
        
        isotherm.reinit(Isotherm::SorptionType(adsorption_type[i_subst]), limited_solubility_on,
                        solvent_density_, scale_aqua, scale_sorbed, table_limit, mult_coef, second_coef);

    }

    END_TIMER("SorptionMob::isotherm_reinit");
}


/***********************************                   *****************************************************/
/*********************************** SORPTION_IMMOBILE *****************************************************/
/***********************************                   *****************************************************/

IT::Record SorptionImmob::input_type = SorptionBase::record_factory(SorptionRecord::immobile);

const int SorptionImmob::registrar =
		Input::register_class< SorptionImmob, Mesh &, Input::Record >("SorptionImmobile");

SorptionImmob::SorptionImmob(Mesh &init_mesh, Input::Record in_rec)
: SorptionDual(init_mesh, in_rec, "conc_immobile_solid", "SorptionImmobile_Output")
{}

SorptionImmob::~SorptionImmob(void)
{}

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
    double phi = por_m/(por_m + por_imm);
        
    // List of types of isotherms in particular regions
    arma::uvec adsorption_type = data_->sorption_type.value(elem.centre(),elem);
    arma::Col<double> mult_coef_vec = data_->isotherm_mult.value(elem.centre(),elem);
    arma::Col<double> second_coef_vec = data_->isotherm_other.value(elem.centre(),elem);

    for(unsigned int i_subst = 0; i_subst < n_substances_; i_subst++)
    {
        double mult_coef = mult_coef_vec[i_subst];
        double second_coef = second_coef_vec[i_subst];
        Isotherm & isotherm = isotherms_vec[i_subst];

        //scales are different for the case of sorption in mobile and immobile pores
        double scale_aqua = por_imm,
               scale_sorbed = (1 - phi) * (1 - por_m - por_imm) * rock_density * substances_[substance_global_idx_[i_subst]].molar_mass();

        bool limited_solubility_on;
        double table_limit;
        if (solubility_vec_[i_subst] <= 0.0) {
            limited_solubility_on = false;
            table_limit=table_limit_[i_subst];

        } else {
            limited_solubility_on = true;
            table_limit=solubility_vec_[i_subst];
        }
        
        if( (1-por_m-por_imm) <= std::numeric_limits<double>::epsilon()) //means there is no sorbing surface
        {
            isotherm.reinit(Isotherm::none, false, solvent_density_, scale_aqua, scale_sorbed,0,0,0);
            continue;
        }
        
        if ( scale_sorbed <= 0.0)
            xprintf(UsrErr, "Scaling parameter in sorption is not positive. Check the input for rock density and molar mass of %d. substance.",i_subst);
        
        isotherm.reinit(Isotherm::SorptionType(adsorption_type[i_subst]), limited_solubility_on,
                    solvent_density_, scale_aqua, scale_sorbed, table_limit, mult_coef, second_coef);

    }

    END_TIMER("SorptionImmob::isotherm_reinit");
}
