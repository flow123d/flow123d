/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    sorption_base.cc
 * @brief   
 */

#include <boost/foreach.hpp>

#include "reaction/sorption_base.hh"
#include "reaction/reaction_term.hh"
#include "reaction/first_order_reaction.hh"
#include "reaction/radioactive_decay.hh"
#include "reaction/dual_porosity.hh"
#include "reaction/isotherm.hh"

#include "system/system.hh"
#include "system/sys_profiler.hh"

#include "la/distribution.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "mesh/accessors.hh"
#include "input/type_selection.hh"

#include "fields/field_set.hh"
#include "fields/field_fe.hh"
#include "fields/fe_value_handler.hh"

using namespace Input::Type;

const Selection & SorptionBase::EqData::get_sorption_type_selection() {
	return Selection("SorptionType")
		.add_value(Isotherm::none,"none", "No sorption considered.")
		.add_value(Isotherm::linear, "linear",
				"Linear isotherm runs the concentration exchange between liquid and solid.")
		.add_value(Isotherm::langmuir, "langmuir",
				"Langmuir isotherm runs the concentration exchange between liquid and solid.")
		.add_value(Isotherm::freundlich, "freundlich",
				"Freundlich isotherm runs the concentration exchange between liquid and solid.")
		.close();
}



const Record & SorptionBase::get_input_type() {
	return Record("SorptionAux", "AUXILIARY RECORD. Should not be directly part of the input tree.")
		.declare_key("substances", Array(String(),1), Default::obligatory(),
					 "Names of the substances that take part in the sorption model.")
		.declare_key("solvent_density", Double(0.0), Default("1.0"),
					"Density of the solvent.")
		.declare_key("substeps", Integer(1), Default("1000"),
					"Number of equidistant substeps, molar mass and isotherm intersections")
		.declare_key("solubility", Array(Double(0.0)), Default::optional(), //("-1.0"), //
								"Specifies solubility limits of all the sorbing species.")
		.declare_key("table_limits", Array(Double(-1.0)), Default::optional(), //("-1.0"), //
                             "Specifies the highest aqueous concentration in the isotherm function interpolation table. "
                             "Use any negative value for an automatic choice according to current maximal concentration (default and recommended). "
                             "Use '0' to always evaluate isotherm function directly (can be very slow). "
                             "Use a positive value to set the interpolation table limit manually "
                             "(if aqueous concentration is higher, then the isotherm function is evaluated directly).")
		.declare_key("input_fields", Array(EqData("","").input_data_set_.make_field_descriptor_type("Sorption")), Default::obligatory(), //
						"Containes region specific data necessary to construct isotherms.")//;
		.declare_key("reaction_liquid", ReactionTerm::it_abstract_reaction(), Default::optional(), "Reaction model following the sorption in the liquid.")
		.declare_key("reaction_solid", ReactionTerm::it_abstract_reaction(), Default::optional(), "Reaction model following the sorption in the solid.")
		.close();
}
    

SorptionBase::EqData::EqData(const string &output_field_name, const string &output_field_desc)
{
    *this += rock_density.name("rock_density")
            .description("Rock matrix density.")
            .input_default("0.0")
            .units( UnitSI().kg().m(-3) );

    *this += sorption_type.name("sorption_type")
            .description("Considered sorption is described by selected isotherm.\n"
                "If porosity on an element is equal to 1.0 (or even higher), meaning no sorbing surface, then type 'none' will be selected automatically.")
            .input_selection(get_sorption_type_selection())
            .units( UnitSI::dimensionless() );

    *this += distribution_coefficient.name("distribution_coefficient")
            .description("Distribution coefficient (( $k_l, k_F, k_L $)) of linear, Freundlich or Langmuir isotherm respectively.")
            .input_default("1.0")
            .units( UnitSI().m(3).kg(-1) );

    *this += isotherm_other.name("isotherm_other")
            .description("Additional parameter (($ \\alpha $)) of nonlinear isotherms.")
            .input_default("1.0")
            .units( UnitSI::dimensionless() );

    *this += init_conc_solid.name("init_conc_solid")
            .description("Initial solid concentration of substances. It is a vector: one value for every substance.")
            .input_default("0")
            .units( UnitSI().dimensionless() );

    input_data_set_ += *this;

    // porosity field is set from governing equation (transport) later
    // hence we do not add it to the input_data_set_
    *this += porosity
            .name("porosity")
            .units( UnitSI::dimensionless() )
            .flags(FieldFlag::input_copy)
			.set_limits(0.0);
    
    output_fields += *this;
    output_fields += conc_solid.name(output_field_name)
                     .description(output_field_desc)
                     .units( UnitSI().dimensionless() )
                     .flags(FieldFlag::equation_result);
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
    if (data_ != nullptr) delete data_;
}

void SorptionBase::make_reactions()
{
  Input::Iterator<Input::AbstractRecord> reactions_it;
  
  reactions_it = input_record_.find<Input::AbstractRecord>("reaction_liquid");
  if ( reactions_it )
  {
    // TODO: allowed instances in this case are only
    // FirstOrderReaction, RadioactiveDecay
    reaction_liquid = (*reactions_it).factory< ReactionTerm, Mesh &, Input::Record >(*mesh_, *reactions_it);
  } else
  {
    reaction_liquid = nullptr;
  }
  
  reactions_it = input_record_.find<Input::AbstractRecord>("reaction_solid");
  if ( reactions_it )
  {
    // TODO: allowed instances in this case are only
    // FirstOrderReaction, RadioactiveDecay
	reaction_solid = (*reactions_it).factory< ReactionTerm, Mesh &, Input::Record >(*mesh_, *reactions_it);
  } else
  {
    reaction_solid = nullptr;
  }
}

void SorptionBase::initialize()
{
  ASSERT_PTR(time_).error("Time governor has not been set yet.\n");
  ASSERT(output_stream_).error("Null output stream.\n");
  ASSERT_LT(0, substances_.size());
  
  initialize_substance_ids(); //computes present substances and sets indices
  initialize_from_input();          //reads non-field data from input
  
  //isotherms array resized bellow
  unsigned int nr_of_regions = mesh_->region_db().bulk_size();
  isotherms.resize(nr_of_regions);
  max_conc.resize(nr_of_regions);
  for(unsigned int i_reg = 0; i_reg < nr_of_regions; i_reg++)
  {
    isotherms[i_reg].resize(n_substances_);
    max_conc[i_reg].resize(n_substances_, 0.0);
    for(unsigned int i_subst = 0; i_subst < n_substances_; i_subst++)
    {
      isotherms[i_reg][i_subst] = Isotherm();
    }
  }   
  
  initialize_fields();
  
  if(reaction_liquid)
  {
    reaction_liquid->substances(substances_)
      .concentration_fields(conc_mobile_fe)
      .set_time_governor(*time_);
    reaction_liquid->initialize();
  }
  if(reaction_solid)
  {
    reaction_solid->substances(substances_)
      .concentration_fields(data_->conc_solid_fe)
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
		// fill table_limit_ with negative values -> automatic limit
        table_limit_.clear();
		table_limit_.resize(n_substances_,-1.0);
	}
}

void SorptionBase::initialize_fields()
{
  ASSERT_GT(n_substances_, 0).error("Number of substances is wrong, they might have not been set yet.\n");

  // create vector of substances that are involved in sorption
  // and initialize data_ with their names
  std::vector<std::string> substances_sorption;
  for (unsigned int i : substance_global_idx_)
    substances_sorption.push_back(substances_[i].name());
  data_->set_components(substances_sorption);
  
  // read fields from input file
  data_->input_data_set_.set_input_list(input_record_.val<Input::Array>("input_fields"), *time_);
  
  data_->set_mesh(*mesh_);

  //initialization of output
  //output_array = input_record_.val<Input::Array>("output_fields");
  data_->conc_solid.set_components(substances_.names());
  data_->output_fields.set_mesh(*mesh_);
  data_->output_fields.output_type(OutputTime::ELEM_DATA);
  data_->conc_solid.setup_components();

  //creating field fe and output multifield for sorbed concentrations
  data_->conc_solid_fe.resize(substances_.size());
  for (unsigned int sbi = 0; sbi < substances_.size(); sbi++)
  {
      // create shared pointer to a FieldFE and push this Field to output_field on all regions
      data_->conc_solid_fe[sbi] = create_field_fe< 3, FieldValue<3>::Scalar >(dof_handler_);
      data_->conc_solid[sbi].set(data_->conc_solid_fe[sbi], 0);
  }
  //output_stream_->add_admissible_field_names(output_array);
  data_->output_fields.initialize(output_stream_, mesh_, input_record_.val<Input::Record>("output"), time());
}


void SorptionBase::zero_time_step()
{
  ASSERT_PTR(time_).error("Time governor has not been set yet.\n");
  ASSERT(output_stream_).error("Null output stream.\n");
  ASSERT_LT(0, substances_.size());
  
  data_->set_time(time_->step(), LimitSide::right);
  std::stringstream ss; // print warning message with table of uninitialized fields
  if ( FieldCommon::print_message_table(ss, "sorption") ) {
      WarningOut() << ss.str();
  }
  set_initial_condition();
  
  update_max_conc();
  make_tables();
    
  // write initial condition
  //data_->output_fields.set_time(time_->step(), LimitSide::right);
  //data_->output_fields.output(output_stream_);
  
  if(reaction_liquid) reaction_liquid->zero_time_step();
  if(reaction_solid) reaction_solid->zero_time_step();

  output_data();
}

void SorptionBase::set_initial_condition()
{
    for ( DHCellAccessor dh_cell : dof_handler_->own_range() ) {
        IntIdx dof_p0 = dh_cell.get_loc_dof_indices()[0];
        const ElementAccessor<3> ele = dh_cell.elm();

        //setting initial solid concentration for substances involved in adsorption
        for (unsigned int sbi = 0; sbi < n_substances_; sbi++)
        {
            int subst_id = substance_global_idx_[sbi];
            data_->conc_solid_fe[subst_id]->vec()[dof_p0] = data_->init_conc_solid[sbi].value(ele.centre(), ele);
        }
    }
}


void SorptionBase::update_solution(void)
{
  data_->set_time(time_->step(), LimitSide::right); // set to the last computed time

  // if parameters changed during last time step, reinit isotherms and eventualy 
  // update interpolation tables in the case of constant rock matrix parameters
  make_tables();
  clear_max_conc();

  START_TIMER("Sorption");
  for ( DHCellAccessor dh_cell : dof_handler_->own_range() )
  {
      compute_reaction(dh_cell);
  }
  END_TIMER("Sorption");
  
  if(reaction_liquid) reaction_liquid->update_solution();
  if(reaction_solid) reaction_solid->update_solution();
}

void SorptionBase::isotherm_reinit(unsigned int i_subst, const ElementAccessor<3> &elem)
{
    START_TIMER("SorptionBase::isotherm_reinit");
    
    double mult_coef = data_->distribution_coefficient[i_subst].value(elem.centre(),elem);
    double second_coef = data_->isotherm_other[i_subst].value(elem.centre(),elem);
    
    int reg_idx = elem.region().bulk_idx();
    Isotherm & isotherm = isotherms[reg_idx][i_subst];
    
    bool limited_solubility_on = solubility_vec_[i_subst] > 0.0;
    
    // in case of no sorbing surface, set isotherm type None
    if( common_ele_data.no_sorbing_surface_cond <= std::numeric_limits<double>::epsilon())
    {
        isotherm.reinit(Isotherm::none, false, solvent_density_,
                        common_ele_data.scale_aqua, common_ele_data.scale_sorbed,
                        0,0,0);
        return;
    }
    
    if ( common_ele_data.scale_sorbed <= 0.0)
        xprintf(UsrErr, "Scaling parameter in sorption is not positive. Check the input for rock density and molar mass of %d. substance.",i_subst);
    
    isotherm.reinit(Isotherm::SorptionType(data_->sorption_type[i_subst].value(elem.centre(),elem)),
                    limited_solubility_on, solvent_density_,
                    common_ele_data.scale_aqua, common_ele_data.scale_sorbed,
                    solubility_vec_[i_subst], mult_coef, second_coef);
}

void SorptionBase::isotherm_reinit_all(const ElementAccessor<3> &elem)
{
    for(unsigned int i_subst = 0; i_subst < n_substances_; i_subst++)
    {
        isotherm_reinit(i_subst, elem);
    }
}

void SorptionBase::clear_max_conc()
{
    unsigned int reg_idx, i_subst;
    
    // clear max concetrations array
    unsigned int nr_of_regions = mesh_->region_db().bulk_size();
    for(reg_idx = 0; reg_idx < nr_of_regions; reg_idx++)
        for(i_subst = 0; i_subst < n_substances_; i_subst++)
            max_conc[reg_idx][i_subst] = 0.0;
}

void SorptionBase::update_max_conc()
{
    unsigned int reg_idx, i_subst, subst_id;
    
    clear_max_conc();
    
    for ( DHCellAccessor dh_cell : dof_handler_->own_range() ) {
        IntIdx dof_p0 = dh_cell.get_loc_dof_indices()[0];
        reg_idx = dh_cell.elm().region().bulk_idx();
        for(i_subst = 0; i_subst < n_substances_; i_subst++){
            subst_id = substance_global_idx_[i_subst];
            max_conc[reg_idx][i_subst] = std::max(max_conc[reg_idx][i_subst], conc_mobile_fe[subst_id]->vec()[dof_p0]);
      }
    }
}

void SorptionBase::make_tables(void)
{
    START_TIMER("SorptionBase::make_tables");
    try
    {
        ElementAccessor<3> elm;
        for(const Region &reg_iter: this->mesh_->region_db().get_region_set("BULK"))
        {
            int reg_idx = reg_iter.bulk_idx();
            // true if data has been changed and are constant on the region
            bool call_reinit = data_->changed() && data_->is_constant(reg_iter);
            
            if(call_reinit)
            {
                ElementAccessor<3> elm(this->mesh_, reg_iter);
//                 DebugOut().fmt("isotherm reinit\n");
                compute_common_ele_data(elm);
                isotherm_reinit_all(elm);
            }
            
            // find table limit and create interpolation table for every substance
            for(unsigned int i_subst = 0; i_subst < n_substances_; i_subst++){
                
                // clear interpolation tables, if not spacially constant OR switched off
                if(! data_->is_constant(reg_iter) || table_limit_[i_subst] == 0.0){
                    isotherms[reg_idx][i_subst].clear_table();
//                     DebugOut().fmt("limit: 0.0 -> clear table\n");
                    continue;
                }
                
                // if true then make_table will be called at the end
                bool call_make_table = call_reinit;
                // initialy try to keep the current table limit (it is zero at zero time step)
                double subst_table_limit = isotherms[reg_idx][i_subst].table_limit();
                
                // if automatic, possibly remake tables with doubled range when table maximum was reached
                if(table_limit_[i_subst] < 0.0)
                {
                    if(subst_table_limit < max_conc[reg_idx][i_subst])
                    {
                        call_make_table = true;
                        subst_table_limit = 2*max_conc[reg_idx][i_subst];
//                         DebugOut().fmt("limit: max conc\n");
                    }
                }
                // if not automatic, set given table limit
                else
                {
                    subst_table_limit = table_limit_[i_subst];
                }
                
                if(call_make_table){
                    isotherms[reg_idx][i_subst].make_table(n_interpolation_steps_, subst_table_limit);
//                     DebugOut().fmt("reg: {} i_subst {}: table_limit = {}\n", reg_idx, i_subst, isotherms[reg_idx][i_subst].table_limit());
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

void SorptionBase::compute_reaction(const DHCellAccessor& dh_cell)
{
    const ElementAccessor<3> ele = dh_cell.elm();
    int reg_idx = ele.region().bulk_idx();
    IntIdx dof_p0 = dh_cell.get_loc_dof_indices()[0];
    unsigned int i_subst, subst_id;
    // for checking, that the common element data are computed once at maximum
    bool is_common_ele_data_valid = false;
    
    try{
        for(i_subst = 0; i_subst < n_substances_; i_subst++)
        {
            subst_id = substance_global_idx_[i_subst];
            Isotherm & isotherm = isotherms[reg_idx][i_subst];
            if (isotherm.is_precomputed()){
//                 DebugOut().fmt("isotherms precomputed - interpolate, subst[{}]\n", i_subst);
                isotherm.interpolate(conc_mobile_fe[subst_id]->vec()[dof_p0],
                                     data_->conc_solid_fe[subst_id]->vec()[dof_p0]);
            }
            else{
//                 DebugOut().fmt("isotherms reinit - compute , subst[{}]\n", i_subst);
                if(! is_common_ele_data_valid){
                    compute_common_ele_data(ele);
                    is_common_ele_data_valid = true;
                }
                
                isotherm_reinit(i_subst, ele);
                isotherm.compute(conc_mobile_fe[subst_id]->vec()[dof_p0],
                                 data_->conc_solid_fe[subst_id]->vec()[dof_p0]);
            }
            
            // update maximal concentration per region (optimization for interpolation)
            if(table_limit_[i_subst] < 0)
                max_conc[reg_idx][i_subst] = std::max(max_conc[reg_idx][i_subst],
                                                      conc_mobile_fe[subst_id]->vec()[dof_p0]);
        }
    }
    catch(ExceptionBase const &e)
    {
        e << input_record_.ei_address();
        throw;
    }
}


/**************************************** OUTPUT ***************************************************/

void SorptionBase::output_data(void )
{
    data_->output_fields.set_time(time().step(), LimitSide::right);
    // Register fresh output data
    data_->output_fields.output(time().step());
}
