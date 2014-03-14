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

#include "reaction/sorption_dual.hh"
#include "reaction/sorption_immob.hh"
#include "reaction/sorption_mob.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "semchem/semchem_interface.hh"

using namespace Input::Type;
using namespace std;


const double Dual_por_exchange::min_dt = 1e-5;

Record Dual_por_exchange::input_type
        = Record("DualPorosity",
            "Dual porosity model in transport problems.\n"
            "Provides computing the concentration of substances in mobile and immobile zone.\n"
            )
    .derive_from(Reaction::input_type)
    .declare_key("data", Array(Dual_por_exchange::EqData().make_field_descriptor_type("DualPorosity")), Default::obligatory(),
                    "Containes region specific data necessary to construct dual porosity model.")
    
    .declare_key("reactions_mob", Reaction::input_type, Default::optional(), "Reaction model in mobile zone.")
    .declare_key("reactions_immob", Reaction::input_type, Default::optional(), "Reaction model in immobile zone.");
    
Dual_por_exchange::EqData::EqData()
{
  ADD_FIELD(alpha, "Diffusion coefficient of non-equilibrium linear exchange between mobile and immobile zone (dual porosity)."
            " Vector, one value for every substance.", "0");
  ADD_FIELD(immob_porosity, "Porosity of the immobile zone.", "0");
  ADD_FIELD(init_conc_immobile, "Initial concentration of substances in the immobile zone."
            " Vector, one value for every substance.", "0");
}

Dual_por_exchange::Dual_por_exchange(Mesh &init_mesh, Input::Record in_rec, vector<string> &names)
	: Reaction(init_mesh, in_rec, names)
{
    //DBGMSG("DualPorosity - constructor\n");

    data_.alpha.n_comp(n_all_substances_);
    data_.init_conc_immobile.n_comp(n_all_substances_);
    
    
    data_.set_input_list(in_rec.val<Input::Array>("data"));
    
    data_+=(data_.porosity
          .name("porosity")
          .units("0")
         );
    
    data_.set_mesh(init_mesh);
    
    data_.set_limit_side(LimitSide::right);
}

Dual_por_exchange::~Dual_por_exchange(void)
{
  if(reaction_mob != nullptr) delete reaction_mob;
  if(reaction_immob != nullptr) delete reaction_immob;
  
  for (unsigned int sbi = 0; sbi < n_all_substances_; sbi++) 
  {
      //no mpi vectors
      xfree(immob_concentration_matrix[sbi]);
  }

  xfree(immob_concentration_matrix);
}


void Dual_por_exchange::init_from_input(Input::Record in_rec)
{  
  DBGMSG("dual_por init_from_input\n");
  Input::Iterator<Input::AbstractRecord> reactions_it = in_rec.find<Input::AbstractRecord>("reactions_mob");
  if ( reactions_it ) 
  {
    if (reactions_it->type() == Linear_reaction::input_type ) {
        reaction_mob =  new Linear_reaction(*mesh_, *reactions_it, names_);
                
    } else
    if (reactions_it->type() == Pade_approximant::input_type) {
        reaction_mob = new Pade_approximant(*mesh_, *reactions_it, names_ );
    } else
    if (reactions_it->type() == SorptionBase::input_type ) {
        reaction_mob =  new SorptionMob(*mesh_, *reactions_it, names_);
                
       static_cast<SorptionMob *> (reaction_mob) -> set_porosity(data_.porosity);
       static_cast<SorptionMob *> (reaction_mob) -> set_porosity_immobile(data_.immob_porosity);
                
    } else
    if (reactions_it->type() == Dual_por_exchange::input_type ) {
        xprintf(UsrErr, "Dual porosity model cannot have another descendant dual porosity model.\n");
    } else
    if (reactions_it->type() == Semchem_interface::input_type ) 
    {
        xprintf(UsrErr, "Semchem chemistry model is not supported at current time.\n");
    } else 
    {
        xprintf(UsrErr, "Wrong reaction type in DualPorosity model.\n");
    }
    
  } else
  {
    reaction_mob = nullptr;
  }
  
  reactions_it = in_rec.find<Input::AbstractRecord>("reactions_immob");
  if ( reactions_it ) 
  {
    if (reactions_it->type() == Linear_reaction::input_type ) {
        reaction_immob =  new Linear_reaction(*mesh_, *reactions_it, names_);
                
    } else
    if (reactions_it->type() == Pade_approximant::input_type) {
        reaction_immob = new Pade_approximant(*mesh_, *reactions_it, names_ );
    } else
    if (reactions_it->type() == SorptionBase::input_type ) {
        reaction_immob =  new SorptionImmob(*mesh_, *reactions_it, names_);
        
       static_cast<SorptionImmob *> (reaction_immob) -> set_porosity(data_.porosity);        
       static_cast<SorptionImmob *> (reaction_immob) -> set_porosity_immobile(data_.immob_porosity);
                
    } else
    if (reactions_it->type() == Dual_por_exchange::input_type ) {
        xprintf(UsrErr, "Dual porosity model cannot have another descendant dual porosity model.\n");
    } else
    if (reactions_it->type() == Semchem_interface::input_type ) 
    {
        xprintf(UsrErr, "Semchem chemistry model is not supported at current time.\n");
    } else 
    {
        xprintf(UsrErr, "Unknown reactions type in DualPorosity model.\n");
    }
    
  } else
  {
    reaction_immob = nullptr;
  }
}


void Dual_por_exchange::initialize(void )
{ 
  ASSERT(distribution != nullptr, "Distribution has not been set yet.\n");
  ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  
  data_.set_time(*time_);
    
  //allocating memory for immobile concentration matrix
  immob_concentration_matrix = (double**) xmalloc(n_all_substances_ * sizeof(double*));
  for (unsigned int sbi = 0; sbi < n_all_substances_; sbi++)
    immob_concentration_matrix[sbi] = (double*) xmalloc(distribution->lsize() * sizeof(double));
   
  //DBGMSG("DualPorosity - init_conc_immobile.\n");
    
  //copied from convection set_initial_condition
  //setting initial condition for immobile concentration matrix
  FOR_ELEMENTS(mesh_, elem)
  {
    if (!distribution->is_local(el_4_loc[elem.index()])) continue;

    unsigned int index = el_4_loc[elem.index()] - distribution->begin();
    ElementAccessor<3> ele_acc = mesh_->element_accessor(elem.index());
    arma::vec value = data_.init_conc_immobile.value(elem->centre(), ele_acc);
        
    for (int sbi=0; sbi < n_all_substances_; sbi++)
    {
      immob_concentration_matrix[sbi][index] = value(sbi);
    }
  }
  
  // creating reactions from input and setting their parameters
  init_from_input(input_record_);
  
  if(reaction_mob != nullptr)
  { 
    reaction_mob->set_time_governor(*time_);
    reaction_mob->set_concentration_matrix(concentration_matrix, distribution, el_4_loc);
    reaction_mob->initialize();
  }
    
  if(reaction_immob != nullptr) 
  {
    reaction_immob->set_time_governor(*time_);
    reaction_immob->set_concentration_matrix(immob_concentration_matrix, distribution, el_4_loc);
    reaction_immob->initialize();
  }
}


void Dual_por_exchange::update_solution(void) 
{
  DBGMSG("DualPorosity - update solution\n");
  data_.set_time(*time_);
 
  START_TIMER("dual_por_exchange_step");
  for (unsigned int loc_el = 0; loc_el < distribution->lsize(); loc_el++) 
  {
    compute_reaction(immob_concentration_matrix, loc_el);
  }
  END_TIMER("dual_por_exchange_step");
  
  if(reaction_mob != nullptr) reaction_mob->update_solution();
  if(reaction_immob != nullptr) reaction_immob->update_solution();
}


double **Dual_por_exchange::compute_reaction(double **concentrations, int loc_el) 
{
  double conc_avg = 0.0;
  unsigned int sbi, sbi_loc;
  double cm, pcm, ci, pci, por_m, por_imm, temp_exp;
   
  if(time_->dt() >= min_dt)
  {
      ElementFullIter ele = mesh_->element(el_4_loc[loc_el]);
      por_m = data_.porosity.value(ele->centre(),ele->element_accessor());
      por_imm = data_.immob_porosity.value(ele->centre(),ele->element_accessor());
      arma::Col<double> alpha_vec = data_.alpha.value(ele->centre(), ele->element_accessor());
      
      //TODO:
      //for (sbi = 0; sbi < n_substances_; sbi++) //over substances involved in dual porosity model
      for (sbi = 0; sbi < n_all_substances_; sbi++) //over all substances
      {
        //sbi_loc = substance_id[sbi];    //mapping to global substance index
                //previous values
                pcm = concentration_matrix[sbi][loc_el];
                pci = immob_concentration_matrix[sbi][loc_el];

                // ---compute average concentration------------------------------------------
                conc_avg = ((por_m * pcm) + (por_imm * pci)) / (por_m + por_imm);

                if ((conc_avg != 0.0) && (por_imm != 0.0)) {
                        temp_exp = exp(-alpha_vec[sbi] * ((por_m + por_imm) / (por_m * por_imm)) * time_->dt());
                        // ---compute concentration in mobile area-----------------------------------
                        cm = (pcm - conc_avg) * temp_exp + conc_avg;

                        // ---compute concentration in immobile area---------------------------------
                        ci = (pci - conc_avg) * temp_exp + conc_avg;
                        // --------------------------------------------------------------------------
//                         DBGMSG("cm: %f  ci: %f  pcm: %f  pci: %f  conc_avg: %f  alpha: %f  por_m: %f  por_imm: %f  time_dt: %f\n",
//                                 cm, ci, pcm, pci, conc_avg, alpha_vec[sbi], por_m, por_imm, time_->dt());
                        concentration_matrix[sbi][loc_el] = cm;
                        immob_concentration_matrix[sbi][loc_el] = ci;
                }
        }

  }
  else{
      ElementFullIter ele = mesh_->element(el_4_loc[loc_el]);
      por_m = data_.porosity.value(ele->centre(),ele->element_accessor());
      por_imm = data_.immob_porosity.value(ele->centre(),ele->element_accessor());
      arma::Col<double> alpha_vec = data_.alpha.value(ele->centre(), ele->element_accessor());
      
      for (sbi = 0; sbi < n_all_substances_; sbi++) {
                //previous values
                pcm = concentration_matrix[sbi][loc_el];
                pci = immob_concentration_matrix[sbi][loc_el];

                if (por_imm != 0.0) {
                        temp_exp = alpha_vec[sbi]*(pci - pcm);
                        // ---compute concentration in mobile area-----------------------------------
                        cm = temp_exp / por_m + pcm;

                        // ---compute concentration in immobile area---------------------------------
                        ci = -temp_exp / por_imm + pci;
                        // --------------------------------------------------------------------------

                        concentration_matrix[sbi][loc_el] = cm;
                        immob_concentration_matrix[sbi][loc_el] = ci;
                }
        }
  }
  
  return immob_concentration_matrix;
}
