#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "reaction/dual_porosity.hh"
#include "reaction/reaction_term.hh"
#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "la/distribution.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "mesh/region.hh"
#include "fields/field_elementwise.hh" 

#include "reaction/sorption.hh"
#include "reaction/first_order_reaction.hh"
#include "reaction/radioactive_decay.hh"
#include "semchem/semchem_interface.hh"
#include "transport/mass_balance.hh"

using namespace Input::Type;


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
           .units( UnitSI().s(-1) );
  
  *this += porosity_immobile
          .name("porosity_immobile")
          .description("Porosity of the immobile zone.")
          .input_default("0")
          .units( UnitSI::dimensionless() );

  *this += init_conc_immobile
          .name("init_conc_immobile")
          .description("Initial concentration of substances in the immobile zone.")
          .units( UnitSI().kg().m(-3) );

  //creating field for porosity that is set later from the governing equation (transport)
  *this +=porosity
        .name("porosity")
        .units( UnitSI::dimensionless() )
        .flags( FieldFlag::input_copy );

	//creating field for cross section that is set later from the governing equation (transport)
	*this +=cross_section
			.name("cross_section")
	        .units("")
	        .flags( FieldFlag::input_copy );


  output_fields += *this;
  output_fields += conc_immobile.name("conc_immobile").units( UnitSI().kg().m(-3) );
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

  for (unsigned int sbi = 0; sbi < substances_.size(); sbi++)
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
      if (reactions_it->type() == FirstOrderReaction::input_type ) {
          reaction_mobile =  new FirstOrderReaction(*mesh_, *reactions_it);

      } else
      if (reactions_it->type() == RadioactiveDecay::input_type) {
          reaction_mobile = new RadioactiveDecay(*mesh_, *reactions_it);
      } else
      if (reactions_it->type() == SorptionMob::input_type ) {
          reaction_mobile =  new SorptionMob(*mesh_, *reactions_it);
      } else
      if (reactions_it->type() == DualPorosity::input_type ) {
        THROW( ReactionTerm::ExcWrongDescendantModel() 
                << ReactionTerm::EI_Model((*reactions_it).type().type_name()) 
                << (*reactions_it).ei_address());
      } else
      if (reactions_it->type() == Semchem_interface::input_type )
      { THROW( ReactionTerm::ExcWrongDescendantModel() 
                << ReactionTerm::EI_Model((*reactions_it).type().type_name())
                << EI_Message("This model is not currently supported!") 
                << (*reactions_it).ei_address());
      } else
      { //This point cannot be reached. The TYPE_selection will throw an error first. 
        THROW( ExcMessage() 
                << EI_Message("Descending model type selection failed (SHOULD NEVER HAPPEN).") 
                << (*reactions_it).ei_address());
      }
    } else
    {
      reaction_mobile = nullptr;
    }

    reactions_it = input_record_.find<Input::AbstractRecord>("reaction_immobile");
    if ( reactions_it )
    {
      if (reactions_it->type() == FirstOrderReaction::input_type ) {
          reaction_immobile =  new FirstOrderReaction(*mesh_, *reactions_it);

      } else
      if (reactions_it->type() == RadioactiveDecay::input_type) {
          reaction_immobile = new RadioactiveDecay(*mesh_, *reactions_it);
      } else
      if (reactions_it->type() == SorptionImmob::input_type ) {
          reaction_immobile =  new SorptionImmob(*mesh_, *reactions_it);
      } else
      if (reactions_it->type() == DualPorosity::input_type ) {
        THROW( ReactionTerm::ExcWrongDescendantModel() 
                << ReactionTerm::EI_Model((*reactions_it).type().type_name()) 
                << (*reactions_it).ei_address());
      } else
      if (reactions_it->type() == Semchem_interface::input_type )
      { THROW( ReactionTerm::ExcWrongDescendantModel() 
                << ReactionTerm::EI_Model((*reactions_it).type().type_name())
                << EI_Message("This model is not currently supported!") 
                << (*reactions_it).ei_address());
      } else
      { //This point cannot be reached. The TYPE_selection will throw an error first. 
        THROW( ExcMessage() 
                << EI_Message("Descending model type selection failed (SHOULD NEVER HAPPEN).") 
                << (*reactions_it).ei_address());
      }
    } else
    {
      reaction_immobile = nullptr;
    }

}

void DualPorosity::initialize()
{
  ASSERT(distribution_ != nullptr, "Distribution has not been set yet.\n");
  ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  ASSERT(output_stream_,"Null output stream.");
  ASSERT_LESS(0, substances_.size());
  
  //allocating memory for immobile concentration matrix
  conc_immobile = (double**) xmalloc(substances_.size() * sizeof(double*));
  conc_immobile_out = (double**) xmalloc(substances_.size() * sizeof(double*));
  for (unsigned int sbi = 0; sbi < substances_.size(); sbi++)
  {
    conc_immobile[sbi] = (double*) xmalloc(distribution_->lsize() * sizeof(double));
    conc_immobile_out[sbi] = (double*) xmalloc(distribution_->size() * sizeof(double));
  }
  allocate_output_mpi();
  
  initialize_fields();

  if(reaction_mobile != nullptr)
  {
    reaction_mobile->substances(substances_)
                .output_stream(*output_stream_)
                .concentration_matrix(concentration_matrix_, distribution_, el_4_loc_, row_4_el_)
                .set_time_governor(*time_);
    reaction_mobile->initialize();
  }

  if(reaction_immobile != nullptr)
  {
    reaction_immobile->substances(substances_)
                .output_stream(*output_stream_)
                .concentration_matrix(conc_immobile, distribution_, el_4_loc_, row_4_el_)
                .set_time_governor(*time_);
    reaction_immobile->initialize();
  }

}

void DualPorosity::initialize_fields()
{
  //setting fields in data
  data_.set_n_components(substances_.size());

  //setting fields that are set from input file
  input_data_set_+=data_;
  input_data_set_.set_input_list(input_record_.val<Input::Array>("input_fields"));

  data_.set_mesh(*mesh_);
  data_.set_limit_side(LimitSide::right);
  
  //initialization of output
  output_array = input_record_.val<Input::Array>("output_fields");
  
  //initialization of output
  data_.conc_immobile.init(substances_.names());
  data_.conc_immobile.set_mesh(*mesh_);
  data_.output_fields.output_type(OutputTime::ELEM_DATA);

  for (unsigned int sbi=0; sbi<substances_.size(); sbi++)
  {
    // create shared pointer to a FieldElementwise and push this Field to output_field on all regions
    std::shared_ptr<FieldElementwise<3, FieldValue<3>::Scalar> > output_field_ptr(
        new FieldElementwise<3, FieldValue<3>::Scalar>(conc_immobile_out[sbi], substances_.size(), mesh_->n_elements()));
    data_.conc_immobile[sbi].set_field(mesh_->region_db().get_region_set("ALL"), output_field_ptr, 0);
  }
  data_.output_fields.set_limit_side(LimitSide::right);
  output_stream_->add_admissible_field_names(output_array, data_.output_selection);

}


void DualPorosity::set_balance_object(boost::shared_ptr<Balance> &balance,
		const std::vector<unsigned int> &subst_idx)
{
	ReactionTerm::set_balance_object(balance, subst_idx);

	subst_idx_mob_ = subst_idx_;
	for (auto name : substances_.names())
		subst_idx_immob_.push_back(balance_->add_quantity(name + "_immobile"));

	if (reaction_mobile)
		reaction_mobile->set_balance_object(balance_, subst_idx_mob_);
	if (reaction_immobile)
		reaction_immobile->set_balance_object(balance_, subst_idx_immob_);

	sources_mob.resize(substances_.size(), vector<double>(mesh_->region_db().bulk_size(), 0));
	sources_in_mob.resize(substances_.size(), vector<double>(mesh_->region_db().bulk_size(), 0));
	sources_out_mob.resize(substances_.size(), vector<double>(mesh_->region_db().bulk_size(), 0));

	sources_immob.resize(substances_.size(), vector<double>(mesh_->region_db().bulk_size(), 0));
	sources_in_immob.resize(substances_.size(), vector<double>(mesh_->region_db().bulk_size(), 0));
	sources_out_immob.resize(substances_.size(), vector<double>(mesh_->region_db().bulk_size(), 0));

	old_mass_coef_mob_.resize(distribution_->lsize(), 0);
	old_mass_coef_immob_.resize(distribution_->lsize(), 0);
}



void DualPorosity::zero_time_step()
{
  ASSERT(distribution_ != nullptr, "Distribution has not been set yet.\n");
  ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  ASSERT(output_stream_,"Null output stream.");
  ASSERT_LESS(0, substances_.size());
 
  //coupling - passing fields
  if(reaction_mobile)
  {
	  if (typeid(*reaction_mobile) == typeid(SorptionMob))
	  {
			  reaction_mobile->data().set_field("porosity", data_["porosity"]);
			  ((SorptionMob *)reaction_mobile)->set_porosity_immobile(data_.porosity_immobile);
			  reaction_mobile->data().set_field("cross_section", data_["cross_section"]);
	  }
	  if (typeid(*reaction_mobile) == typeid(LinearReaction) ||
		  typeid(*reaction_mobile) == typeid(DecayChain))
	  {
		  reaction_mobile->data().set_field("porosity", data_["porosity"]);
		  reaction_mobile->data().set_field("cross_section", data_["cross_section"]);
	  }
  }
  if(reaction_immobile)
  {
	  if (typeid(*reaction_immobile) == typeid(SorptionImmob))
	  {
		  reaction_immobile->data().set_field("porosity", data_["porosity"]);
		  ((SorptionImmob *)reaction_immobile)->set_porosity_immobile(data_.porosity_immobile);
		  reaction_immobile->data().set_field("cross_section", data_["cross_section"]);
	  }
	  if (typeid(*reaction_immobile) == typeid(LinearReaction) ||
		  typeid(*reaction_immobile) == typeid(DecayChain))
	  {
		  reaction_immobile->data().set_field("porosity", data_["porosity_immobile"]);
		  reaction_immobile->data().set_field("cross_section", data_["cross_section"]);
	  }
  }
  
  data_.set_time(*time_);
  set_initial_condition();
  
	if (data_.porosity_immobile.changed())
	{
		if (data_.cross_section.changed())
			assemble_balance_matrix();

		if (balance_ != nullptr)
		{
			for (unsigned int loc_el=0; loc_el<distribution_->lsize(); ++loc_el)
			{
				ElementFullIter ele = mesh_->element(el_4_loc_[loc_el]);
				double csection = data_.cross_section.value(ele->centre(), ele->element_accessor());
				double por_mob = data_.porosity.value(ele->centre(),ele->element_accessor());
				double por_immob = data_.porosity_immobile.value(ele->centre(),ele->element_accessor());
				old_mass_coef_mob_[loc_el] = por_mob*csection;
				old_mass_coef_immob_[loc_el] = por_immob*csection;
			}
		}
	}

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
  //setting initial condition for immobile concentration matrix
  for (unsigned int loc_el = 0; loc_el < distribution_->lsize(); loc_el++)
  {
    unsigned int index = el_4_loc_[loc_el];
    ElementAccessor<3> ele_acc = mesh_->element_accessor(index);
    arma::vec value = data_.init_conc_immobile.value(ele_acc.centre(), ele_acc);
        
    for (unsigned int sbi=0; sbi < substances_.size(); sbi++)
    {
      conc_immobile[sbi][loc_el] = value(sbi);
    }
  }
}

void DualPorosity::assemble_balance_matrix()
{
    if (balance_ != nullptr)
    {
    	balance_->start_mass_assembly(subst_idx_immob_);

    	for (unsigned int loc_el = 0; loc_el < distribution_->lsize(); loc_el++)
    	{
    		ElementFullIter elm = mesh_->element(el_4_loc_[loc_el]);
    		double csection = data_.cross_section.value(elm->centre(), elm->element_accessor());
    		double por_imm = data_.porosity_immobile.value(elm->centre(), elm->element_accessor());

        	for (unsigned int sbi=0; sbi<substances_.size(); ++sbi)
        		balance_->add_mass_matrix_values(subst_idx_immob_[sbi], elm->region().bulk_idx(), {row_4_el_[el_4_loc_[loc_el]]}, {csection*por_imm*elm->measure()} );
        }

    	balance_->finish_mass_assembly(subst_idx_immob_);
    }
}

void DualPorosity::update_solution(void) 
{
  data_.set_time(*time_);

  if (balance_ != nullptr)
  {
	vector<double> tmp(mesh_->region_db().bulk_size(), 0);
	fill_n(sources_mob.begin(), substances_.size(), tmp);
	fill_n(sources_in_mob.begin(), substances_.size(), tmp);
	fill_n(sources_out_mob.begin(), substances_.size(), tmp);
	fill_n(sources_immob.begin(), substances_.size(), tmp);
	fill_n(sources_in_immob.begin(), substances_.size(), tmp);
	fill_n(sources_out_immob.begin(), substances_.size(), tmp);

	if (data_.porosity_immobile.changed() ||
		data_.cross_section.changed())
		assemble_balance_matrix();
  }
 
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
  double csection = data_.cross_section.value(ele->centre(), ele->element_accessor());
  por_mob = data_.porosity.value(ele->centre(),ele->element_accessor());
  por_immob = data_.porosity_immobile.value(ele->centre(),ele->element_accessor());
  arma::Col<double> diff_vec = data_.diffusion_rate_immobile.value(ele->centre(), ele->element_accessor());
 
    // if porosity_immobile == 0 then mobile concentration stays the same 
    // and immobile concentration cannot change
    if (por_immob == 0.0) return conc_immobile;
    
    double exponent,
           temp_exponent = (por_mob + por_immob) / (por_mob * por_immob) * time_->dt();
  
    for (sbi = 0; sbi < substances_.size(); sbi++) //over all substances
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
            double temp = diff_vec[sbi]*(previous_conc_immob - previous_conc_mob) * time_->dt();
            // ---compute concentration in mobile area
            conc_mob = temp / por_mob + previous_conc_mob;

            // ---compute concentration in immobile area
            conc_immob = -temp / por_immob + previous_conc_immob;
        }
        else                                                        //analytic solution
        {
            double temp = exp(-exponent);
            // ---compute concentration in mobile area
            conc_mob = (previous_conc_mob - conc_average) * temp + conc_average;

            // ---compute concentration in immobile area
            conc_immob = (previous_conc_immob - conc_average) * temp + conc_average;
        }
        
        concentration_matrix_[sbi][loc_el] = conc_mob;
        conc_immobile[sbi][loc_el] = conc_immob;

        if (balance_ != nullptr)
        {
        	double source_mob = (conc_mob*csection*por_mob - previous_conc_mob*old_mass_coef_mob_[loc_el])
        			*ele->measure()/time_->dt();
        	double source_immob = (conc_immob*por_immob*csection - previous_conc_immob*old_mass_coef_immob_[loc_el])
        			*ele->measure()/time_->dt();
			sources_mob[sbi][ele->region().bulk_idx()] += source_mob;
			sources_immob[sbi][ele->region().bulk_idx()] += source_immob;
			if (source_mob > 0)
				sources_in_mob[sbi][ele->region().bulk_idx()] += source_mob;
			else
				sources_out_mob[sbi][ele->region().bulk_idx()] += source_mob;
			if (source_immob > 0)
				sources_in_immob[sbi][ele->region().bulk_idx()] += source_immob;
			else
				sources_out_immob[sbi][ele->region().bulk_idx()] += source_immob;

        }
    }
    if (balance_ != nullptr)
    {
    	old_mass_coef_mob_[loc_el] = por_mob*csection;
    	old_mass_coef_immob_[loc_el] = por_immob*csection;
    }
  
  return conc_immobile;
}

void DualPorosity::update_instant_balance()
{
    for (unsigned int sbi = 0; sbi < substances_.size(); sbi++)
    {
    	balance_->add_instant_sources(subst_idx_mob_[sbi], sources_mob[sbi], sources_in_mob[sbi], sources_out_mob[sbi]);
    	balance_->add_instant_sources(subst_idx_immob_[sbi], sources_immob[sbi], sources_in_immob[sbi], sources_out_immob[sbi]);
    	balance_->calculate_mass(subst_idx_immob_[sbi], vconc_immobile[sbi]);
    }

    if (reaction_mobile) reaction_mobile->update_instant_balance();
    if (reaction_immobile) reaction_immobile->update_instant_balance();
}


void DualPorosity::update_cumulative_balance()
{
    for (unsigned int sbi = 0; sbi < substances_.size(); sbi++)
    {
    	balance_->add_cumulative_sources(subst_idx_mob_[sbi], sources_mob[sbi], time_->dt());
    	balance_->add_cumulative_sources(subst_idx_immob_[sbi], sources_immob[sbi], time_->dt());
    }

    if (reaction_mobile) reaction_mobile->update_cumulative_balance();
    if (reaction_immobile) reaction_immobile->update_cumulative_balance();
}



void DualPorosity::allocate_output_mpi(void )
{
    int sbi, n_subst;
    n_subst = substances_.size();

    vconc_immobile = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vconc_immobile_out = (Vec*) xmalloc(n_subst * (sizeof(Vec))); // extend to all


    for (sbi = 0; sbi < n_subst; sbi++) {
        VecCreateMPIWithArray(PETSC_COMM_WORLD,1, distribution_->lsize(), mesh_->n_elements(), conc_immobile[sbi],
                &vconc_immobile[sbi]);
        VecZeroEntries(vconc_immobile[sbi]);

        //  if(rank == 0)
        VecCreateSeqWithArray(PETSC_COMM_SELF,1, mesh_->n_elements(), conc_immobile_out[sbi], &vconc_immobile_out[sbi]);
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
    unsigned int sbi;

    for (sbi = 0; sbi < substances_.size(); sbi++) {
        VecScatterBegin(vconc_out_scatter, vconc_immobile[sbi], vconc_immobile_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(vconc_out_scatter, vconc_immobile[sbi], vconc_immobile_out[sbi], INSERT_VALUES, SCATTER_FORWARD);
    }
}


void DualPorosity::output_data(void )
{
    output_vector_gather();

    // Register fresh output data
    data_.output_fields.set_time(*time_);
    data_.output_fields.output(output_stream_);
    
    if (reaction_mobile) reaction_mobile->output_data();
    if (reaction_immobile) reaction_immobile->output_data();
}

