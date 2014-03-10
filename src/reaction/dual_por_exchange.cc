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

Record Dual_por_exchange::input_type
        = Record("TransportOperatorSplitting",
            "Explicit FVM transport (no diffusion)\n"
            "coupled with reaction and adsorption model (ODE per element)\n"
            " via operator splitting.")
    .derive_from(Reaction::input_type)
    .declare_key("alpha", Dual_por_exchange::input_type, Default::obligatory(),
                "Diffusion parameter between mobile and immobile zone.")
    .declare_key("init_conc_immobile", Dual_por_exchange::input_type, Default::optional(),
                        "Initial value of the concentration in the immobile zone.");
    
Dual_por_exchange::EqData::EqData()
: EqDataBase("Exchange_dp")
{
  ADD_FIELD(alpha, "Diffusion coefficient of non-equilibrium linear exchange between mobile and immobile zone (dual porosity)."
            " Vector, one value for every substance.", Input::Type::Default("0"));
  ADD_FIELD(immob_porosity, "Porosity of the immobile zone.", Input::Type::Default("0"));
  ADD_FIELD(init_conc_immobile, "Initial concentration of substances in the immobile zone."
            " Vector, one value for every substance.", Input::Type::Default("0"));
}

Dual_por_exchange::Dual_por_exchange(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)
	: Reaction(init_mesh, in_rec, names)
{
    cout << "Dual_por_exchange constructor is running." << endl;

    data_.set_mesh(&init_mesh);
    data_.init_from_input( in_rec.val<Input::Array>("bulk_data"), Input::Array());
    
    time_ = new TimeGovernor();
    data_.set_time(*time_);

    nr_of_substances = names.size();
}

Dual_por_exchange::~Dual_por_exchange(void)
{
}

void Dual_por_exchange::update_solution(void) {

    data_.set_time(*time_);
    double conc_avg = 0.0;
    unsigned int loc_el,sbi;
    double cm, pcm, ci, pci, por_m, por_imm, alpha;
    
    START_TIMER("dual_por_exchange_step");
    for (sbi = 0; sbi < nr_of_substances; sbi++) {
    	for (loc_el = 0; loc_el < distribution->lsize(); loc_el++) {
    		ElementFullIter elem = mesh_->element(el_4_loc[loc_el]);
    		por_m = data_.porosity->value(elem->centre(),elem->element_accessor());
    		por_imm = data_.porosity->value(elem->centre(),elem->element_accessor());
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
    START_TIMER("dual_por_exchange_step");
}


void Dual_por_exchange::set_porosity(pScalar porosity)
{
  data_.porosity = porosity;
  return;
}

void Dual_por_exchange::init_from_input(Input::Record in_rec)
{
	// Initialize from input interface
}
