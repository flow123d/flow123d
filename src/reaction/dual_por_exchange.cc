#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <boost/foreach.hpp>

#include "reaction/reaction.hh"
#include "reaction/dual_por_exchange.hh"
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

Dual_por_exchange::EqData::EqData()
{
	ADD_FIELD(alphas, "Diffusion coefficient of non-equilibrium linear exchange between mobile and immobile zone (dual porosity)."
            " Vector, one value for every substance.", Input::Type::Default("0"));
}

Dual_por_exchange::Dual_por_exchange(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)
	: Reaction(init_mesh, in_rec, names)
{
	cout << "Dual_por_exchange constructor is running." << endl;
	TimeGovernor tg(0.0, 1.0);
    nr_of_regions = init_mesh.region_db().bulk_size();

    data_.set_mesh(&init_mesh);
    data_.init_from_input( in_rec.val<Input::Array>("data"));
    data_.set_time(tg);

    nr_of_substances = names.size();
}

Dual_por_exchange::~Dual_por_exchange(void)
{
}

void Dual_por_exchange::compute_one_step(void) {

    double conc_avg = 0.0;
    unsigned int loc_el,sbi;
    double cm, pcm, ci, pci, por_m, por_imm, alpha;

    START_TIMER("dual_por_exchange_step");
    for (sbi = 0; sbi < nr_of_substances; sbi++) {
    	for (loc_el = 0; loc_el < distribution->lsize(); loc_el++) {
    		ElementFullIter elem = mesh_->element(el_4_loc[loc_el]);
    		por_m = this->porosity_->value(elem->centre(),elem->element_accessor());
    		por_imm = this->porosity_->value(elem->centre(),elem->element_accessor());
    		alpha = data_.alphas.value(elem->centre(), elem->element_accessor())(sbi);
    		pcm = this->concentration_matrix[sbi][loc_el];
    		pci = this->immob_concentration_matrix[sbi][loc_el];

    		// ---compute average concentration------------------------------------------
    		conc_avg = ((por_m * pcm) + (por_imm * pci)) / (por_m + por_imm);

    		if ((conc_avg != 0.0) && (por_imm != 0.0)) {
        		// ---compute concentration in mobile area-----------------------------------
        		cm = (pcm - conc_avg) * exp(-alpha * ((por_m + por_imm) / (por_m * por_imm)) * time_step) + conc_avg;

        		// ---compute concentration in immobile area---------------------------------
        		ci = (pci - conc_avg) * exp(-alpha * ((por_m + por_imm) / (por_m * por_imm)) * time_step) + conc_avg;
        		// --------------------------------------------------------------------------

        		this->concentration_matrix[sbi][loc_el] = cm;
        		this->concentration_matrix[sbi][loc_el] = ci;
    		}
    	}
    }
    START_TIMER("dual_por_exchange_step");
}

void Dual_por_exchange::set_immob_concentration_matrix(double **ConcentrationMatrix, Distribution *conc_distr, int *el_4_loc_)
{
	immob_concentration_matrix = ConcentrationMatrix;
	distribution = conc_distr;
	el_4_loc = el_4_loc_;
	return;
}

void Dual_por_exchange::set_porosity(pScalar porosity, pScalar immob_porosity)
{
	this->porosity_ = porosity;
	this->immob_porosity_ = immob_porosity;
	return;
}
