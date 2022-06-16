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
 * @file    first_order_reaction_base.cc
 * @brief   
 */

#include "reaction/first_order_reaction_base.hh"
#include "reaction/reaction_term.hh"

#include "reaction/linear_ode_solver.hh"


#include "system/global_defs.h"
#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "la/distribution.hh"
#include "input/accessors.hh"



using namespace Input::Type;


FirstOrderReactionBase::FirstOrderReactionBase(Mesh &init_mesh, Input::Record in_rec)
    : ReactionTerm(init_mesh, in_rec)
{
    linear_ode_solver_ = std::make_shared<LinearODESolver>();
    this->eq_fields_base_ = std::make_shared<EqFields>();
    this->eq_data_base_ = std::make_shared<EqData>();
}

FirstOrderReactionBase::~FirstOrderReactionBase()
{
}

void FirstOrderReactionBase::initialize()
{
	ASSERT_PERMANENT(time_ != nullptr).error("Time governor has not been set yet.\n");
	ASSERT_PERMANENT_LT(0, eq_data_base_->substances_.size()).error("No substances for rection term.\n");
    
    n_substances_ = eq_data_base_->substances_.size();
    initialize_from_input();

    // allocation
    prev_conc_.resize(n_substances_);
    reaction_matrix_.resize(n_substances_, n_substances_);
    molar_matrix_.resize(n_substances_, n_substances_);
    molar_mat_inverse_.resize(n_substances_, n_substances_);

    // initialize diagonal matrices with molar masses
    molar_matrix_.zeros();
    molar_mat_inverse_.zeros();
    for (unsigned int i=0; i<n_substances_; ++i)
    {
    	molar_matrix_(i,i) = eq_data_base_->substances_[i].molar_mass();
    	molar_mat_inverse_(i,i) = 1./eq_data_base_->substances_[i].molar_mass();
    }
}


void FirstOrderReactionBase::zero_time_step()
{
    ASSERT(time_ != nullptr).error("Time governor has not been set yet.\n");
	ASSERT_LT(0, eq_data_base_->substances_.size()).error("No substances for rection term.\n");

    assemble_ode_matrix();
    // make scaling that takes into account different molar masses of substances
    reaction_matrix_ = molar_matrix_ * reaction_matrix_ * molar_mat_inverse_;
    
    linear_ode_solver_->set_system_matrix(reaction_matrix_);
}


void FirstOrderReactionBase::compute_reaction(const DHCellAccessor& dh_cell)
{      
    unsigned int sbi;  // row in the concentration matrix, regards the substance index
    arma::vec new_conc;
    
    IntIdx dof_p0 = dh_cell.get_loc_dof_indices()[0];

    // save previous concentrations to column vector
    for(sbi = 0; sbi < n_substances_; sbi++)
        prev_conc_(sbi) = this->eq_fields_base_->conc_mobile_fe[sbi]->vec().get(dof_p0);
    
    // compute new concetrations R*c
    linear_ode_solver_->update_solution(prev_conc_, new_conc);
    
    // save new concentrations to the concentration matrix
    for(sbi = 0; sbi < n_substances_; sbi++)
        this->eq_fields_base_->conc_mobile_fe[sbi]->vec().set( dof_p0, new_conc(sbi) );
}

void FirstOrderReactionBase::update_solution(void)
{
    //DebugOut() << "FirstOrderReactionBases - update solution\n";
    if(time_->is_changed_dt())
    {
        linear_ode_solver_->set_step(time_->dt());
    }

    START_TIMER("linear reaction step");

    for ( DHCellAccessor dh_cell : eq_data_base_->dof_handler_->own_range() )
    {
        compute_reaction(dh_cell);
    }
    END_TIMER("linear reaction step");
}


unsigned int FirstOrderReactionBase::find_subst_name(const string &name)
{
    unsigned int k=0;
        for(; k < n_substances_; k++)
                if (name == eq_data_base_->substances_[k].name()) return k;

        return k;
}
