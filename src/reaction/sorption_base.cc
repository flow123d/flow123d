#include <boost/foreach.hpp>

#include "reaction/sorption_base.hh"
#include "reaction/reaction_term.hh"
#include "reaction/first_order_reaction.hh"
#include "reaction/radioactive_decay.hh"
#include "reaction/dual_porosity.hh"
#include "semchem/semchem_interface.hh"
#include "reaction/isotherm.hh"

#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "la/distribution.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "input/type_selection.hh"

#include "fields/field_set.hh"
#include "fields/field_elementwise.hh" 

using namespace Input::Type;

Selection SorptionBase::EqData::sorption_type_selection = Selection("AdsorptionType")
	.add_value(Isotherm::none,"none", "No adsorption considered.")
	.add_value(Isotherm::linear, "linear",
			"Linear isotherm runs the concentration exchange between liquid and solid.")
	.add_value(Isotherm::langmuir, "langmuir",
			"Langmuir isotherm runs the concentration exchange between liquid and solid.")
	.add_value(Isotherm::freundlich, "freundlich",
			"Freundlich isotherm runs the concentration exchange between liquid and solid.");



Record SorptionBase::input_type
	= Record("Adsorption", "AUXILIARY RECORD. Should not be directly part of the input tree.")
    .declare_key("substances", Array(String(),1), Default::obligatory(),
                 "Names of the substances that take part in the adsorption model.")
	.declare_key("solvent_density", Double(0.0), Default("1.0"),
				"Density of the solvent.")
	.declare_key("substeps", Integer(1), Default("1000"),
				"Number of equidistant substeps, molar mass and isotherm intersections")
	.declare_key("solubility", Array(Double(0.0)), Default::optional(), //("-1.0"), //
							"Specifies solubility limits of all the adsorbing species.")
	.declare_key("table_limits", Array(Double(0.0)), Default::optional(), //("-1.0"), //
							"Specifies highest aqueous concentration in interpolation table.")
    .declare_key("input_fields", Array(EqData("").input_data_set_.make_field_descriptor_type("Sorption")), Default::obligatory(), //
                    "Containes region specific data necessary to construct isotherms.")//;
    .declare_key("reaction_liquid", ReactionTerm::input_type, Default::optional(), "Reaction model following the sorption in the liquid.")
    .declare_key("reaction_solid", ReactionTerm::input_type, Default::optional(), "Reaction model following the sorption in the solid.");
    

SorptionBase::EqData::EqData(const string &output_field_name)
{
    ADD_FIELD(rock_density, "Rock matrix density.", "0.0");
    	rock_density.units( UnitSI().kg().m(-3) );

    ADD_FIELD(sorption_type,"Considered adsorption is described by selected isotherm. If porosity on an element is equal or even higher than 1.0 (meaning no sorbing surface), then type 'none' will be selected automatically."); //
        sorption_type.input_selection(&sorption_type_selection);
        sorption_type.units( UnitSI::dimensionless() );

    ADD_FIELD(isotherm_mult,"Multiplication parameters (k, omega) in either Langmuir c_s = omega * (alpha*c_a)/(1- alpha*c_a) or in linear c_s = k * c_a isothermal description.","1.0");
    	isotherm_mult.units( UnitSI().mol().kg(-1) );

    ADD_FIELD(isotherm_other,"Second parameters (alpha, ...) defining isotherm  c_s = omega * (alpha*c_a)/(1- alpha*c_a).","1.0");
    	isotherm_other.units( UnitSI::dimensionless() );

    ADD_FIELD(init_conc_solid, "Initial solid concentration of substances."
            " Vector, one value for every substance.", "0");
    	init_conc_solid.units( UnitSI().mol().kg(-1) );

    input_data_set_ += *this;

    // porosity field is set from governing equation (transport) later
    // hence we do not add it to the input_data_set_
    *this += porosity
            .name("porosity")
            .units( UnitSI::dimensionless() )
            .flags(FieldFlag::input_copy);
    
    output_fields += *this;
    output_fields += conc_solid.name(output_field_name).units( UnitSI().kg().m(-3) );
}

Record SorptionBase::record_factory(SorptionBase::SorptionRecord::Type fact)
{
    IT::Record rec;
    switch(fact)
    {
        case SorptionRecord::mobile:
            rec = IT::Record("SorptionMobile", "Adsorption model in the mobile zone, following the dual porosity model.")
                .derive_from( ReactionTerm::input_type )
                .copy_keys(SorptionBase::input_type)
                .declare_key("output_fields", IT::Array(make_output_selection("conc_solid", "SorptionMobile_Output")),
                    IT::Default("conc_solid"), "List of fields to write to output stream.");
            break;
        case SorptionRecord::immobile:  
            rec = IT::Record("SorptionImmobile", "Adsorption model in the immobile zone, following the dual porosity model.")
                .derive_from( ReactionTerm::input_type )
                .copy_keys(SorptionBase::input_type)
                .declare_key("output_fields", IT::Array(make_output_selection("conc_immobile_solid", "SorptionImmobile_Output")),
                    IT::Default("conc_immobile_solid"), "List of fields to write to output stream.");
            break;
            
        default:
            rec = IT::Record("Sorption", "Adsorption model in the reaction term of transport.")
                .derive_from( ReactionTerm::input_type )
                .copy_keys(SorptionBase::input_type)
                .declare_key("output_fields", IT::Array(make_output_selection("conc_solid", "Sorption_Output")),
                             IT::Default("conc_solid"), "List of fields to write to output stream.");
            break;
    }
    return rec;
}


SorptionBase::SorptionBase(Mesh &init_mesh, Input::Record in_rec)//
	: ReactionTerm(init_mesh, in_rec),
	  data_(nullptr)
{
  // creating reaction from input and setting their parameters
  make_reactions();
}


SorptionBase::~SorptionBase(void)
{
  if(reaction_liquid != nullptr) delete reaction_liquid;
  if(reaction_solid != nullptr) delete reaction_solid;
  if (data_ != nullptr) delete data_;

  VecScatterDestroy(&(vconc_out_scatter));
  VecDestroy(vconc_solid);

  for (unsigned int sbi = 0; sbi < substances_.size(); sbi++)
  {
    //no mpi vectors
    xfree(conc_solid[sbi]);
  }
  xfree(conc_solid);
}

void SorptionBase::make_reactions()
{
  Input::Iterator<Input::AbstractRecord> reactions_it;
  
  reactions_it = input_record_.find<Input::AbstractRecord>("reaction_liquid");
  if ( reactions_it )
  {
    if (reactions_it->type() == FirstOrderReaction::input_type ) {
        reaction_liquid =  new FirstOrderReaction(*mesh_, *reactions_it);
    } else
    if (reactions_it->type() == RadioactiveDecay::input_type) {
        reaction_liquid = new RadioactiveDecay(*mesh_, *reactions_it);
    } else
    if (reactions_it->type() == SorptionBase::input_type ) {
        THROW( ReactionTerm::ExcWrongDescendantModel() 
                << ReactionTerm::EI_Model((*reactions_it).type().type_name()) 
                << (*reactions_it).ei_address());
    } else
    if (reactions_it->type() == DualPorosity::input_type ) {
        THROW( ReactionTerm::ExcWrongDescendantModel() 
                << ReactionTerm::EI_Model((*reactions_it).type().type_name()) 
                << (*reactions_it).ei_address());
    } else
    if (reactions_it->type() == Semchem_interface::input_type )
    {   THROW( ReactionTerm::ExcWrongDescendantModel() 
                << ReactionTerm::EI_Model((*reactions_it).type().type_name())
                << EI_Message("This model is not currently supported!") 
                << (*reactions_it).ei_address());
    } else
    {   //This point cannot be reached. The TYPE_selection will throw an error first. 
        THROW( ExcMessage() 
                << EI_Message("Descending model type selection failed (SHOULD NEVER HAPPEN).") 
                << (*reactions_it).ei_address());
    }
  } else
  {
    reaction_liquid = nullptr;
  }
  
  reactions_it = input_record_.find<Input::AbstractRecord>("reaction_solid");
  if ( reactions_it )
  {
    if (reactions_it->type() == FirstOrderReaction::input_type ) {
        reaction_solid =  new FirstOrderReaction(*mesh_, *reactions_it);
    } else
    if (reactions_it->type() == RadioactiveDecay::input_type) {
        reaction_solid = new RadioactiveDecay(*mesh_, *reactions_it);
    } else
    if (reactions_it->type() == SorptionBase::input_type ) {
        THROW( ReactionTerm::ExcWrongDescendantModel() 
                << ReactionTerm::EI_Model((*reactions_it).type().type_name()) 
                << (*reactions_it).ei_address());
    } else
    if (reactions_it->type() == DualPorosity::input_type ) {
        THROW( ReactionTerm::ExcWrongDescendantModel() 
                << ReactionTerm::EI_Model((*reactions_it).type().type_name()) 
                << (*reactions_it).ei_address());
    } else
    if (reactions_it->type() == Semchem_interface::input_type )
    {   THROW( ReactionTerm::ExcWrongDescendantModel() 
                << ReactionTerm::EI_Model((*reactions_it).type().type_name())
                << EI_Message("This model is not currently supported!") 
                << (*reactions_it).ei_address());
    } else
    {   //This point cannot be reached. The TYPE_selection will throw an error first. 
        THROW( ExcMessage() 
                << EI_Message("Descending model type selection failed (SHOULD NEVER HAPPEN).") 
                << (*reactions_it).ei_address());
    }
  } else
  {
    reaction_solid = nullptr;
  }
}

void SorptionBase::initialize()
{
  ASSERT(distribution_ != nullptr, "Distribution has not been set yet.\n");
  ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  ASSERT(output_stream_,"Null output stream.");
  ASSERT_LESS(0, substances_.size());
  
  initialize_substance_ids(); //computes present substances and sets indices
  initialize_from_input();          //reads non-field data from input
  
  //isotherms array resized bellow
  unsigned int nr_of_regions = mesh_->region_db().bulk_size();
  isotherms.resize(nr_of_regions);
  for(unsigned int i_reg = 0; i_reg < nr_of_regions; i_reg++)
  {
    isotherms[i_reg].resize(n_substances_);
    for(unsigned int i_subst = 0; i_subst < n_substances_; i_subst++)
    {
      isotherms[i_reg][i_subst] = Isotherm();
    }
  }   
  
  //allocating new array for sorbed concentrations
  conc_solid = (double**) xmalloc(substances_.size() * sizeof(double*));//new double * [n_substances_];
  conc_solid_out.clear();
  conc_solid_out.resize( substances_.size() );
  for (unsigned int sbi = 0; sbi < substances_.size(); sbi++)
  {
    conc_solid[sbi] = (double*) xmalloc(distribution_->lsize() * sizeof(double));//new double[ nr_of_local_elm ];
    conc_solid_out[sbi].resize( distribution_->size() );
    //zero initialization of solid concentration for all substances
    for(unsigned int i=0; i < distribution_->lsize(); i++)
      conc_solid[sbi][i] = 0;
  }

  allocate_output_mpi();
  
  initialize_fields();
  
  if(reaction_liquid != nullptr)
  {
    reaction_liquid->substances(substances_)
      .concentration_matrix(concentration_matrix_, distribution_, el_4_loc_, row_4_el_)
      .set_time_governor(*time_);
    reaction_liquid->initialize();
  }
  if(reaction_solid != nullptr)
  {
    reaction_solid->substances(substances_)
      .concentration_matrix(conc_solid, distribution_, el_4_loc_, row_4_el_)
      .set_time_governor(*time_);
    reaction_solid->initialize();
  }
}


void SorptionBase::initialize_substance_ids()
{
  Input::Array substances_array = input_record_.val<Input::Array>("substances");
  unsigned int k, global_idx, i_subst = 0;
  bool found;
  Input::Iterator<string> spec_iter = substances_array.begin<string>();
  
  for(; spec_iter != substances_array.end(); ++spec_iter, i_subst++)
  {
    //finding the name of a substance in the global array of names
    found = false;
    for(k = 0; k < substances_.size(); k++)
    {
      if (*spec_iter == substances_[k].name())
      {
        global_idx = k;
        found = true;
        break;
      }
    }
    
    if(!found)
        THROW(ReactionTerm::ExcUnknownSubstance() 
                << ReactionTerm::EI_Substance(*spec_iter) 
                << substances_array.ei_address());
    
    //finding the global index of substance in the local array
    found = false;
    for(k = 0; k < substance_global_idx_.size(); k++)
    {
      if(substance_global_idx_[k] == global_idx)
      {
        found = true;
        break;
      }
    }
    
    if(!found)
    {
      substance_global_idx_.push_back(global_idx);
    }

  }  
  n_substances_ = substance_global_idx_.size();
}

void SorptionBase::initialize_from_input()
{
    // read number of interpolation steps - value checked by the record definition
    n_interpolation_steps_ = input_record_.val<int>("substeps");
    
    // read the density of solvent - value checked by the record definition
	solvent_density_ = input_record_.val<double>("solvent_density");

    // read the solubility vector
	Input::Iterator<Input::Array> solub_iter = input_record_.find<Input::Array>("solubility");
	if( solub_iter )
	{
		if (solub_iter->Array::size() != n_substances_)
		{
            THROW(SorptionBase::ExcSubstanceCountMatch() 
                << SorptionBase::EI_ArrayName("solubility") 
                << input_record_.ei_address());
            // there is no way to get ei_address from 'solub_iter', only as a string
		}
		
		else solub_iter->copy_to(solubility_vec_);
	}
	else{
		// fill solubility_vec_ with zeros
        solubility_vec_.clear();
		solubility_vec_.resize(n_substances_,0.0);
	}

	// read the interpolation table limits
	Input::Iterator<Input::Array> interp_table_limits = input_record_.find<Input::Array>("table_limits");
	if( interp_table_limits )
	{
		if (interp_table_limits->Array::size() != n_substances_)
		{
            THROW(SorptionBase::ExcSubstanceCountMatch() 
                << SorptionBase::EI_ArrayName("table_limits") 
                << input_record_.ei_address());
            // there is no way to get ei_address from 'interp_table_limits', only as a string
		}
		
		else interp_table_limits->copy_to(table_limit_);
	}
	else{
		// fill table_limit_ with zeros
        table_limit_.clear();
		table_limit_.resize(n_substances_,0.0);
	}
}

void SorptionBase::initialize_fields()
{
  ASSERT(n_substances_ > 0, "Number of substances is wrong, they might have not been set yet.\n");

  // read fields from input file
  data_->input_data_set_.set_input_list(input_record_.val<Input::Array>("input_fields"));

  data_->set_components( std::vector<std::string>(n_substances_, "") );
  data_->set_mesh(*mesh_);
  data_->set_limit_side(LimitSide::right);

  //initialization of output
  output_array = input_record_.val<Input::Array>("output_fields");
  data_->conc_solid.set_components(substances_.names());
  data_->output_fields.set_mesh(*mesh_);
  data_->output_fields.set_limit_side(LimitSide::right);
  data_->output_fields.output_type(OutputTime::ELEM_DATA);
  data_->conc_solid.set_up_components();
  for (unsigned int sbi=0; sbi<substances_.size(); sbi++)
  {
      // create shared pointer to a FieldElementwise and push this Field to output_field on all regions
	  auto output_field_ptr = conc_solid_out[sbi].create_field<3, FieldValue<3>::Scalar>(substances_.size());
      data_->conc_solid[sbi].set_field(mesh_->region_db().get_region_set("ALL"), output_field_ptr, 0);
  }
  output_stream_->add_admissible_field_names(output_array);
}


void SorptionBase::zero_time_step()
{
  ASSERT(distribution_ != nullptr, "Distribution has not been set yet.\n");
  ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  ASSERT(output_stream_,"Null output stream.");
  ASSERT_LESS(0, substances_.size());
  
  data_->set_time(time_->step());
  set_initial_condition();
  make_tables();
    
  // write initial condition
  output_vector_gather();
  data_->output_fields.set_time(time_->step());
  data_->output_fields.output(output_stream_);
  
  if(reaction_liquid != nullptr) reaction_liquid->zero_time_step();
  if(reaction_solid != nullptr) reaction_solid->zero_time_step();
}

void SorptionBase::set_initial_condition()
{
  for (unsigned int loc_el = 0; loc_el < distribution_->lsize(); loc_el++)
  {
    unsigned int index = el_4_loc_[loc_el];
    ElementAccessor<3> ele_acc = mesh_->element_accessor(index);
    arma::vec value = data_->init_conc_solid.value(ele_acc.centre(),
        ele_acc);

    //setting initial solid concentration for substances involved in adsorption
    for (unsigned int sbi = 0; sbi < n_substances_; sbi++)
    {
      int subst_id = substance_global_idx_[sbi];
      conc_solid[subst_id][loc_el] = value(sbi);
    }
  }
}


void SorptionBase::update_solution(void)
{
  data_->set_time(time_->step()); // set to the last computed time

  // if parameters changed during last time step, reinit isotherms and eventualy 
  // update interpolation tables in the case of constant rock matrix parameters
  if(data_->changed())
    make_tables();
    

  START_TIMER("Sorption");
  for (unsigned int loc_el = 0; loc_el < distribution_->lsize(); loc_el++)
  {
    compute_reaction(concentration_matrix_, loc_el);
  }
  END_TIMER("Sorption");
  
  if(reaction_liquid != nullptr) reaction_liquid->update_solution();
  if(reaction_solid != nullptr) reaction_solid->update_solution();
}

void SorptionBase::make_tables(void)
{
    try
    {
        ElementAccessor<3> elm;
        BOOST_FOREACH(const Region &reg_iter, this->mesh_->region_db().get_region_set("BULK") )
        {
            int reg_idx = reg_iter.bulk_idx();

            if(data_->is_constant(reg_iter))
            {
                ElementAccessor<3> elm(this->mesh_, reg_iter); // constant element accessor
                isotherm_reinit(isotherms[reg_idx],elm);
                for(unsigned int i_subst = 0; i_subst < n_substances_; i_subst++)
                {
                    isotherms[reg_idx][i_subst].make_table(n_interpolation_steps_);
                }
            }
        }
    }
    catch(ExceptionBase const &e)
    {
        e << input_record_.ei_address();
        throw;
    }
}

double **SorptionBase::compute_reaction(double **concentrations, int loc_el)
{
    ElementFullIter elem = mesh_->element(el_4_loc_[loc_el]);
    int reg_idx = elem->region().bulk_idx();
    unsigned int i_subst, subst_id;

    std::vector<Isotherm> & isotherms_vec = isotherms[reg_idx];
    
    try{
        // Constant value of rock density and mobile porosity over the whole region 
        // => interpolation_table is precomputed
        if (isotherms_vec[0].is_precomputed()) 
        {
            for(i_subst = 0; i_subst < n_substances_; i_subst++)
            {
                subst_id = substance_global_idx_[i_subst];
            
                isotherms_vec[i_subst].interpolate(concentration_matrix_[subst_id][loc_el], 
                                                conc_solid[subst_id][loc_el]);
            }
        }
        else 
        {
            isotherm_reinit(isotherms_vec, elem->element_accessor());
        
            for(i_subst = 0; i_subst < n_substances_; i_subst++)
            {
                subst_id = substance_global_idx_[i_subst];
                isotherms_vec[i_subst].compute(concentration_matrix_[subst_id][loc_el], 
                                            conc_solid[subst_id][loc_el]);
            }
        }
    }
    catch(ExceptionBase const &e)
    {
        e << input_record_.ei_address();
        throw;
    }

  return concentrations;
}


/**************************************** OUTPUT ***************************************************/

void SorptionBase::allocate_output_mpi(void )
{
    int sbi, n_subst;
    n_subst = substances_.size();

    vconc_solid = (Vec*) xmalloc(n_subst * (sizeof(Vec)));

    for (sbi = 0; sbi < n_subst; sbi++) {
        VecCreateMPIWithArray(PETSC_COMM_WORLD,1, distribution_->lsize(), mesh_->n_elements(), conc_solid[sbi],
                &vconc_solid[sbi]);
        VecZeroEntries(vconc_solid[sbi]);

        VecZeroEntries(conc_solid_out[sbi].get_data_petsc());
    }
    
    // creating output vector scatter
    IS is;
    ISCreateGeneral(PETSC_COMM_SELF, mesh_->n_elements(), row_4_el_, PETSC_COPY_VALUES, &is); //WithArray
    VecScatterCreate(vconc_solid[0], is, conc_solid_out[0].get_data_petsc(), PETSC_NULL, &vconc_out_scatter);
    ISDestroy(&(is));
}


void SorptionBase::output_vector_gather() 
{
    unsigned int sbi;

    for (sbi = 0; sbi < substances_.size(); sbi++) {
        VecScatterBegin(vconc_out_scatter, vconc_solid[sbi], conc_solid_out[sbi].get_data_petsc(), INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(vconc_out_scatter, vconc_solid[sbi], conc_solid_out[sbi].get_data_petsc(), INSERT_VALUES, SCATTER_FORWARD);
    }
}


void SorptionBase::output_data(void )
{
    output_vector_gather();

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0)
    {
      // Register fresh output data
      data_->output_fields.set_time(time_->step());
      data_->output_fields.output(output_stream_);
    }

    //for synchronization when measuring time by Profiler
    MPI_Barrier(MPI_COMM_WORLD);
}
