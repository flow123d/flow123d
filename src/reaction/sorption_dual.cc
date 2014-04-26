#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <boost/foreach.hpp>

//#include "reaction/reaction.hh"
//#include "reaction/isotherm.hh"
#include "reaction/sorption_dual.hh"
//#include "system/system.hh"
//#include "system/sys_profiler.hh"

//#include "la/distribution.hh"
#include "mesh/mesh.h"
//#include "mesh/elements.h"
//#include "mesh/region.hh"

//#include "coupling/time_governor.hh"



using namespace std;

SorptionDual::SorptionDual(Mesh &init_mesh, Input::Record in_rec)//
	: SorptionBase(init_mesh, in_rec)
{
  immob_porosity_.just_copy()
    .name("porosity_immobile");
}

SorptionDual::~SorptionDual(void)
{
}

/*
void SorptionDual::isotherm_reinit(std::vector<Isotherm> &isotherms_vec, const ElementAccessor<3> &elem)
{
        START_TIMER("SorptionMob::isotherm_reinit");

        double rock_density = data_.rock_density.value(elem.centre(),elem);

        double por_m = data_.porosity.value(elem.centre(),elem);
        double por_imm = immob_porosity_.value(elem.centre(),elem);

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
                scale_sorbed = compute_sorbing_scale(por_m,por_imm) * (1 - por_m - por_imm) * rock_density * molar_masses[i_subst];
                if(scale_sorbed == 0.0)
                        xprintf(UsrErr, "Parameter scale_sorbed (phi * (1 - por_m - por_imm) * rock_density * molar_masses[i_subst]) is equal to zero.");
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

        END_TIMER("SorptionMob::isotherm_reinit");

        return;
}
*/
/*
// Computes adsorption simulation over all the elements.
void SorptionDual::update_solution(void)
{
    data_.set_time(*time_); // set to the last computed time
    immob_porosity_.set_time(*time_);
    
    //if parameters changed during last time step, reinit isotherms and eventualy update interpolation tables in the case of constant rock matrix parameters
    if(data_.changed()
      || immob_porosity_.changed()
    )
      make_tables();

    START_TIMER("SorptionDual");
    for (int loc_el = 0; loc_el < distribution->lsize(); loc_el++)
    {
      compute_reaction(concentration_matrix, loc_el);
    }
    END_TIMER("SorptionDual");
}
*/
