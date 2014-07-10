#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "reaction/reaction.hh"
#include "reaction/dual_por_exchange.hh"
#include "system/system.hh"
#include "system/sys_profiler.hh"
#include <petscmat.h>

#include "la/distribution.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "mesh/region.hh"
#include "fields/field_elementwise.hh" 

#include "reaction/sorption.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "semchem/semchem_interface.hh"

using namespace Input::Type;
using namespace std;


Selection DualPorosity::EqData::output_selection
		= EqData().output_fields.make_output_field_selection("DualPorosity_Output")
		.close();

Record DualPorosity::input_type
        = Record("DualPorosity",
            "Dual porosity model in transport problems.\n"
            "Provides computing the concentration of substances in mobile and immobile zone.\n"
            )
    .derive_from(ReactionTerm::input_type)
    .declare_key("input_fields", Array(DualPorosity::EqData().make_field_descriptor_type("DualPorosity")), Default::obligatory(),
                    "Containes region specific data necessary to construct dual porosity model.")
    .declare_key("scheme_tolerance", Double(0.0), Default("1e-3"), 
                 "Tolerance according to which the explicit Euler scheme is used or not."
                 "Set 0.0 to use analytic formula only (can be slower).")
    
    .declare_key("reaction_mobile", ReactionTerm::input_type, Default::optional(), "Reaction model in mobile zone.")
    .declare_key("reaction_immobile", ReactionTerm::input_type, Default::optional(), "Reaction model in immobile zone.")
    
    .declare_key("output_fields", Array(EqData::output_selection),
                Default("conc_immobile"), "List of fields to write to output stream.");
    
DualPorosity::EqData::EqData()
{
  *this += diffusion_rate_immobile
           .name("diffusion_rate_immobile")
           .description("Diffusion coefficient of non-equilibrium linear exchange between mobile and immobile zone.")
           .input_default("0")
           .units("");
  
  *this += porosity_immobile
          .name("porosity_immobile")
          .description("Porosity of the immobile zone.")
          .input_default("0")
          .units("1");

  *this += init_conc_immobile
          .name("init_conc_immobile")
          .description("Initial concentration of substances in the immobile zone.")
          .units("M/L^3");

  //creating field for porosity that is set later from the governing equation (transport)
  *this +=porosity
        .name("porosity")
        .units("1")
        .flags( FieldFlag::input_copy );

  output_fields += *this;
  output_fields += conc_immobile.name("conc_immobile").units("M/L^3");
}

DualPorosity::DualPorosity(Mesh &init_mesh, Input::Record in_rec)
	: ReactionTerm(init_mesh, in_rec)
{
  //set pointer to equation data fieldset
  this->eq_data_ = &data_;
  
  //reads input and creates possibly other reaction terms
  make_reactions();
  //read value from input
  scheme_tolerance_ = input_record_.val<double>("scheme_tolerance");
}

DualPorosity::~DualPorosity(void)
{
  if(reaction_mobile != nullptr) delete reaction_mobile;
  if(reaction_immobile != nullptr) delete reaction_immobile;

  VecScatterDestroy(&(vconc_out_scatter));
  VecDestroy(vconc_immobile);
  VecDestroy(vconc_immobile_out);

  for (unsigned int sbi = 0; sbi < names_.size(); sbi++)
  {
      //no mpi vectors
      xfree(conc_immobile[sbi]);
      xfree(conc_immobile_out[sbi]);
  }

  xfree(conc_immobile);
  xfree(conc_immobile_out);
}


void DualPorosity::make_reactions() {
    Input::Iterator<Input::AbstractRecord> reactions_it = input_record_.find<Input::AbstractRecord>("reaction_mobile");
    if ( reactions_it )
    {
      if (reactions_it->type() == LinearReaction::input_type ) {
          reaction_mobile =  new LinearReaction(*mesh_, *reactions_it);

      } else
      if (reactions_it->type() == PadeApproximant::input_type) {
          reaction_mobile = new PadeApproximant(*mesh_, *reactions_it);
      } else
      if (reactions_it->type() == SorptionMob::input_type ) {
          reaction_mobile =  new SorptionMob(*mesh_, *reactions_it);
      } else
      if (reactions_it->type() == DualPorosity::input_type ) {
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
      reaction_mobile = nullptr;
    }

    reactions_it = input_record_.find<Input::AbstractRecord>("reaction_immobile");
    if ( reactions_it )
    {
      if (reactions_it->type() == LinearReaction::input_type ) {
          reaction_immobile =  new LinearReaction(*mesh_, *reactions_it);

      } else
      if (reactions_it->type() == PadeApproximant::input_type) {
          reaction_immobile = new PadeApproximant(*mesh_, *reactions_it);
      } else
      if (reactions_it->type() == SorptionImmob::input_type ) {
          reaction_immobile =  new SorptionImmob(*mesh_, *reactions_it);
      } else
      if (reactions_it->type() == DualPorosity::input_type ) {
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
      reaction_immobile = nullptr;
    }

}

void DualPorosity::initialize()
{
  //DBGMSG("DualPorosity - initialize.\n");
  ASSERT(distribution_ != nullptr, "Distribution has not been set yet.\n");
  ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  ASSERT(output_stream_,"Null output stream.");
  ASSERT_LESS(0, names_.size());
  
  //allocating memory for immobile concentration matrix
  conc_immobile = (double**) xmalloc(names_.size() * sizeof(double*));
  conc_immobile_out = (double**) xmalloc(names_.size() * sizeof(double*));
  for (unsigned int sbi = 0; sbi < names_.size(); sbi++)
  {
    conc_immobile[sbi] = (double*) xmalloc(distribution_->lsize() * sizeof(double));
    conc_immobile_out[sbi] = (double*) xmalloc(distribution_->size() * sizeof(double));
  }
  allocate_output_mpi();
  
  initialize_fields();

  if(reaction_mobile != nullptr)
  {
    reaction_mobile->names(names_)
                .output_stream(*output_stream_)
                .concentration_matrix(concentration_matrix_, distribution_, el_4_loc_, row_4_el_)
                .set_time_governor(*time_);
    reaction_mobile->initialize();
  }

  if(reaction_immobile != nullptr)
  {
    reaction_immobile->names(names_)
                .output_stream(*output_stream_)
                .concentration_matrix(conc_immobile, distribution_, el_4_loc_, row_4_el_)
                .set_time_governor(*time_);
    reaction_immobile->initialize();
  }

}

void DualPorosity::initialize_fields()
{
  //setting fields in data
  data_.set_n_components(names_.size());

  //setting fields that are set from input file
  input_data_set_+=data_;
  input_data_set_.set_input_list(input_record_.val<Input::Array>("input_fields"));

  data_.set_mesh(*mesh_);
  data_.set_limit_side(LimitSide::right);
  
  //initialization of output
  output_array = input_record_.val<Input::Array>("output_fields");
  
  //initialization of output
  data_.conc_immobile.init(names_);
  data_.conc_immobile.set_mesh(*mesh_);
  data_.output_fields.output_type(OutputTime::ELEM_DATA);

  for (unsigned int sbi=0; sbi<names_.size(); sbi++)
  {
    // create shared pointer to a FieldElementwise and push this Field to output_field on all regions
    std::shared_ptr<FieldElementwise<3, FieldValue<3>::Scalar> > output_field_ptr(
        new FieldElementwise<3, FieldValue<3>::Scalar>(conc_immobile_out[sbi], names_.size(), mesh_->n_elements()));
    data_.conc_immobile[sbi].set_field(mesh_->region_db().get_region_set("ALL"), output_field_ptr, 0);
  }
  data_.output_fields.set_limit_side(LimitSide::right);
  output_stream_->add_admissible_field_names(output_array, data_.output_selection);
}


void DualPorosity::zero_time_step()
{
  //DBGMSG("DualPorosity - zero_time_step.\n");
  ASSERT(distribution_ != nullptr, "Distribution has not been set yet.\n");
  ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  ASSERT(output_stream_,"Null output stream.");
  ASSERT_LESS(0, names_.size());
 
  //coupling - passing fields
  if(reaction_mobile)
  if (typeid(*reaction_mobile) == typeid(SorptionMob))
  {
          reaction_mobile->data().set_field("porosity", data_["porosity"]);
          reaction_mobile->data().set_field("porosity_immobile", data_["porosity_immobile"]);
  }
  if(reaction_immobile)
  if (typeid(*reaction_immobile) == typeid(SorptionImmob))
  {
      reaction_immobile->data().set_field("porosity", data_["porosity"]);
      reaction_immobile->data().set_field("porosity_immobile", data_["porosity_immobile"]);
  }
  
  data_.set_time(*time_);
  set_initial_condition();
  
  // write initial condition
  output_vector_gather();
  data_.output_fields.set_time(*time_);
  data_.output_fields.output(output_stream_);
  
  if(reaction_mobile != nullptr)
    reaction_mobile->zero_time_step();

  if(reaction_immobile != nullptr)
    reaction_immobile->zero_time_step();
}

void DualPorosity::set_initial_condition()
{
  //DBGMSG("DualPorosity - init_conc_immobile.\n");
  //setting initial condition for immobile concentration matrix
  for (unsigned int loc_el = 0; loc_el < distribution_->lsize(); loc_el++)
  {
    unsigned int index = el_4_loc_[loc_el];
    ElementAccessor<3> ele_acc = mesh_->element_accessor(index);
    arma::vec value = data_.init_conc_immobile.value(ele_acc.centre(), ele_acc);
        
    for (unsigned int sbi=0; sbi < names_.size(); sbi++)
    {
      conc_immobile[sbi][loc_el] = value(sbi);
    }
  }
}

void DualPorosity::update_solution(void) 
{
  //DBGMSG("DualPorosity - update solution\n");
  data_.set_time(*time_);
 
  START_TIMER("dual_por_exchange_step");
  for (unsigned int loc_el = 0; loc_el < distribution_->lsize(); loc_el++) 
  {
    compute_reaction(conc_immobile, loc_el);
  }
  END_TIMER("dual_por_exchange_step");
  
  if(reaction_mobile != nullptr) reaction_mobile->update_solution();
  if(reaction_immobile != nullptr) reaction_immobile->update_solution();
}


double **DualPorosity::compute_reaction(double **concentrations, int loc_el) 
{
  unsigned int sbi;
  double conc_average, // weighted (by porosity) average of concentration
         conc_mob, conc_immob,  // new mobile and immobile concentration
         previous_conc_mob, previous_conc_immob, // mobile and immobile concentration in previous time step
         conc_max, //difference between concentration and average concentration
         por_mob, por_immob; // mobile and immobile porosity
   
  // get data from fields
  ElementFullIter ele = mesh_->element(el_4_loc_[loc_el]);
  por_mob = data_.porosity.value(ele->centre(),ele->element_accessor());
  por_immob = data_.porosity_immobile.value(ele->centre(),ele->element_accessor());
  arma::Col<double> diff_vec = data_.diffusion_rate_immobile.value(ele->centre(), ele->element_accessor());
 
// OLD CODE:  
//     for (sbi = 0; sbi < names_.size(); sbi++) //over all substances
//     {
//         //sbi_loc = substance_id[sbi];    //mapping to global substance index
//         //previous values
//         previous_conc_mob = concentration_matrix_[sbi][loc_el];
//         previous_conc_immob = conc_immobile[sbi][loc_el];
// 
//         // ---compute average concentration------------------------------------------
//         conc_average = ((por_mob * previous_conc_mob) + (por_immob * previous_conc_immob)) / (por_mob + por_immob);
// 
//         if ((conc_average != 0.0) && (por_immob != 0.0)) {
//                 temp_exp = exp(-diff_vec[sbi] * ((por_mob + por_immob) / (por_mob * por_immob)) * time_->dt());
//                 // ---compute concentration in mobile area-----------------------------------
//                 conc_mob = (previous_conc_mob - conc_average) * temp_exp + conc_average;
// 
//                 // ---compute concentration in immobile area---------------------------------
//                 conc_immob = (previous_conc_immob - conc_average) * temp_exp + conc_average;
//                 // --------------------------------------------------------------------------
// //                 DBGMSG("conc_mob: %f  conc_immob: %f  previous_conc_mob: %f  previous_conc_immob: %f  conc_average: %f  diff: %f  por_mob: %f  por_immob: %f  time_dt: %f\n",
// //                         conc_mob, conc_immob, previous_conc_mob, previous_conc_immob, conc_average, diff_vec[sbi], por_mob, por_immob, time_->dt());
//                 concentration_matrix_[sbi][loc_el] = conc_mob;
//                 conc_immobile[sbi][loc_el] = conc_immob;
//         }
//     }

    
    // if porosity_immobile == 0 then mobile concentration stays the same 
    // and immobile concentration cannot change
    if (por_immob == 0.0) return conc_immobile;
    
    double exponent,
           temp_exponent = (por_mob + por_immob) / (por_mob * por_immob) * time_->dt();
  
    for (sbi = 0; sbi < names_.size(); sbi++) //over all substances
    {
        exponent = diff_vec[sbi] * temp_exponent;
        //previous values
        previous_conc_mob = concentration_matrix_[sbi][loc_el];
        previous_conc_immob = conc_immobile[sbi][loc_el];
        
        // ---compute average concentration------------------------------------------
        conc_average = ((por_mob * previous_conc_mob) + (por_immob * previous_conc_immob)) 
                       / (por_mob + por_immob);
        
        conc_max = std::max(previous_conc_mob-conc_average, previous_conc_immob-conc_average);
        
        if( conc_max <= (2*scheme_tolerance_/(exponent*exponent)*conc_average) )               // forward euler
        {
            //DBGMSG("forward euler\n");
            double temp = diff_vec[sbi]*(previous_conc_immob - previous_conc_mob) * time_->dt();
            // ---compute concentration in mobile area
            conc_mob = temp / por_mob + previous_conc_mob;

            // ---compute concentration in immobile area
            conc_immob = -temp / por_immob + previous_conc_immob;
        }
        else                                                        //analytic solution
        {
            //DBGMSG("analytic\n");
            double temp = exp(-exponent);
            // ---compute concentration in mobile area
            conc_mob = (previous_conc_mob - conc_average) * temp + conc_average;

            // ---compute concentration in immobile area
            conc_immob = (previous_conc_immob - conc_average) * temp + conc_average;
            
//          DBGMSG("conc_mob: %f  conc_immob: %f  previous_conc_mob: %f  previous_conc_immob: %f  conc_average: %f  diff: %f  por_mob: %f  por_immob: %f  time_dt: %f\n",
//                  conc_mob, conc_immob, previous_conc_mob, previous_conc_immob, conc_average, diff_vec[sbi], por_mob, por_immob, time_->dt()); 
        }
        
        concentration_matrix_[sbi][loc_el] = conc_mob;
        conc_immobile[sbi][loc_el] = conc_immob;
    }
  //*/
  
  return conc_immobile;
}


void DualPorosity::allocate_output_mpi(void )
{
    //DBGMSG("DualPorosity - allocate_output_mpi.\n");
    int sbi, n_subst, ierr;
    n_subst = names_.size();

    vconc_immobile = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vconc_immobile_out = (Vec*) xmalloc(n_subst * (sizeof(Vec))); // extend to all


    for (sbi = 0; sbi < n_subst; sbi++) {
        ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1, distribution_->lsize(), mesh_->n_elements(), conc_immobile[sbi],
                &vconc_immobile[sbi]);
        VecZeroEntries(vconc_immobile[sbi]);

        //  if(rank == 0)
        ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,1, mesh_->n_elements(), conc_immobile_out[sbi], &vconc_immobile_out[sbi]);
        VecZeroEntries(vconc_immobile_out[sbi]);
    }
    
    // create output vector scatter
    IS is;
    ISCreateGeneral(PETSC_COMM_SELF, mesh_->n_elements(), row_4_el_, PETSC_COPY_VALUES, &is); //WithArray
    VecScatterCreate(vconc_immobile[0], is, vconc_immobile_out[0], PETSC_NULL, &vconc_out_scatter);
    ISDestroy(&(is));
}


void DualPorosity::output_vector_gather() 
{
    //DBGMSG("DualPorosity - output_vector_gather.\n");
    unsigned int sbi;

    //PetscViewer inviewer;

    for (sbi = 0; sbi < names_.size(); sbi++) {
        VecScatterBegin(vconc_out_scatter, vconc_immobile[sbi], vconc_immobile_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(vconc_out_scatter, vconc_immobile[sbi], vconc_immobile_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
    }
    //VecView(transport->vconc[0],PETSC_VIEWER_STDOUT_WORLD);
    //VecView(transport->vconc_out[0],PETSC_VIEWER_STDOUT_WORLD);
}


void DualPorosity::output_data(void )
{
    //DBGMSG("DualPorosity output\n");
    output_vector_gather();

    // Register fresh output data
    data_.output_fields.set_time(*time_);
    data_.output_fields.output(output_stream_);
    
    if (reaction_mobile) reaction_mobile->output_data();
    if (reaction_immobile) reaction_immobile->output_data();
}

