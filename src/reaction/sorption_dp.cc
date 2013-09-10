#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <boost/foreach.hpp>

#include "reaction/reaction.hh"
#include "reaction/sorption_dp.hh"
#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "la/distribution.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "mesh/region.hh"
#include "input/type_selection.hh"

namespace it=Input::Type;

using namespace Input::Type;

using namespace std;

Sorption_dp::Sorption_dp(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)
	: Sorption(init_mesh, in_rec, names)
{
	cout << "Sorption_dp constructor is running." << endl;
	TimeGovernor tg(0.0, 1.0);
    nr_of_regions = init_mesh.region_db().bulk_size();
    nr_of_substances = in_rec.val<Input::Array>("species").size();
    nr_of_points = in_rec.val<int>("substeps");

    data_.sorption_types.set_n_comp(nr_of_substances);
    data_.mult_coefs.set_n_comp(nr_of_substances);
    data_.second_params.set_n_comp(nr_of_substances);
    // alpha must be set for all transported substances
    int nr_transp_subst = names.size();
    data_.alpha.set_n_comp(nr_transp_subst);
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
			Isotherm iso_immob;
			isotherms[i_reg].push_back(iso_immob);
		}
	}
}

Sorption_dp::~Sorption_dp(void)
{
}

void Sorption_dp::transport_dual_porosity(void) {

    double conc_avg = 0.0;
    unsigned int loc_el,sbi;
    double cm, pcm, ci, pci, por_m, por_imm, alpha;

    for (sbi = 0; sbi < nr_transp_subst_; sbi++) {
    	for (loc_el = 0; loc_el < distribution->lsize(); loc_el++) {
    		ElementFullIter elem = mesh_->element(el_4_loc[loc_el]);
    		por_m = this->porosity_->value(elem->centre(),elem->element_accessor());
    		por_imm = this->porosity_->value(elem->centre(),elem->element_accessor());
    		alpha = data_.alpha.value(elem->centre(), elem->element_accessor())(sbi);
    		pcm = this->concentration_matrix[sbi][loc_el];
    		pci = this->immob_concentration_matrix[sbi][loc_el];

    		// ---compute average concentration------------------------------------------
    		conc_avg = ((por_m * pcm) + (por_imm * pci)) / (por_m + por_imm);

    		if ((conc_avg != 0.0) && (por_imm != 0.0)) {
        		// ---compute concentration in mobile area-----------------------------------
        		cm = (pcm - conc_avg) * exp(-alpha * ((por_m + por_imm) / (por_m * por_imm)) * time_->dt()) + conc_avg;

        		// ---compute concentration in immobile area---------------------------------
        		ci = (pci - conc_avg) * exp(-alpha * ((por_m + por_imm) / (por_m * por_imm)) * time_->dt()) + conc_avg;
        		// --------------------------------------------------------------------------

        		this->concentration_matrix[sbi][loc_el] = cm;
        		this->concentration_matrix[sbi][loc_el] = ci;
    		}
    	}
    }

}

void Sorption_dp::set_scales(double &scale_aqua, double &scale_sorbed, double por_m, double por_imm, double phi, double rock_density, double molar_mass)
{
	scale_aqua = por_imm;
	scale_sorbed = (1 - phi) * (1 - por_m - por_imm) * rock_density * molar_mass;
	return;
}

void Sorption_dp::set_nr_transp(int nr_transp_subst)
{
	nr_transp_subst_ = nr_transp_subst;
	return;
}

void Sorption_dp::set_immob_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc_)
{
	immob_concentration_matrix = ConcentrationMatrix;
	distribution = conc_distr;
	el_4_loc = el_4_loc_;
	return;
}
