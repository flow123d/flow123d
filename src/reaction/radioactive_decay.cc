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
 * @file    radioactive_decay.cc
 * @brief   
 */

#include "reaction/radioactive_decay.hh"
#include "reaction/reaction_term.hh"
#include "reaction/first_order_reaction_base.hh"
#include "reaction/linear_ode_solver.hh"

#include "system/global_defs.h"
#include "mesh/mesh.h"
#include "input/factory.hh"
#include "input/accessors.hh"

#include "armadillo"

FLOW123D_FORCE_LINK_IN_CHILD(radioactiveDecay)


using namespace Input::Type;

const Record & RadioactiveDecay::get_input_type_product() {
	return Record("RadioactiveDecayProduct", "A record describing a product of a radioactive decay.")
		.allow_auto_conversion("name")
		.declare_key("name", String(), Default::obligatory(),
					 "The name of the product.")
		.declare_key("energy", Double(0.0), Default("0.0"),
					 "Not used at the moment! The released energy in MeV from the decay of the radionuclide into the product.")
		.declare_key("branching_ratio", Double(0.0), Default("1.0"),
					 "The branching ratio of the product when there is more than one."
					 "Considering only one product, the default ratio 1.0 is used."
					 "Its value must be positive. Further, the branching ratios of all products are normalized"
					 "by their sum, so the sum then gives 1.0 (this also resolves possible rounding errors).")
		.close();
}

const Record & RadioactiveDecay::get_input_type_single_decay() {
    return Record("Decay", "A model of a radioactive decay.")
		.declare_key("radionuclide", String(), Default::obligatory(),
					"The name of the parent radionuclide.")
		.declare_key("half_life", Double(0.0), Default::obligatory(),
					 "The half life of the parent radionuclide in seconds.")
		.declare_key("products", Array(RadioactiveDecay::get_input_type_product(), 1), Default::obligatory(),
					"An array of the decay products (daughters).")
		.close();
}

const Record & RadioactiveDecay::get_input_type() {
    return Record("RadioactiveDecay", "A model of a radioactive decay and possibly of a decay chain.")
        .derive_from( ReactionTerm::it_abstract_term() )
        .derive_from( ReactionTerm::it_abstract_mobile_term() )
        .derive_from( ReactionTerm::it_abstract_immobile_term() )
        .derive_from( ReactionTerm::it_abstract_reaction() )
		.declare_key("decays", Array( RadioactiveDecay::get_input_type_single_decay(), 1), Default::obligatory(),
					"An array of radioactive decays.")
		.close();
}

const int RadioactiveDecay::registrar =
		Input::register_class< RadioactiveDecay, Mesh &, Input::Record >("RadioactiveDecay") +
		RadioactiveDecay::get_input_type().size();



RadioactiveDecay::RadioactiveDecay(Mesh &init_mesh, Input::Record in_rec)
      : FirstOrderReactionBase(init_mesh, in_rec)
{
}

RadioactiveDecay::~RadioactiveDecay()
{
}


void RadioactiveDecay::initialize_from_input()
{
    Input::Array decay_array = input_record_.val<Input::Array>("decays");

    substance_ids_.resize( decay_array.size() );
    half_lives_.resize( decay_array.size() );
    bifurcation_.resize( decay_array.size() );

    unsigned int idx;
    int i_decay = 0;
    for (Input::Iterator<Input::Record> dec_it = decay_array.begin<Input::Record>(); 
         dec_it != decay_array.end(); ++dec_it, ++i_decay)
    {
        //read half-life
        half_lives_[i_decay] = dec_it->val<double>("half_life");

        //read radionuclide name and products array
        string radionuclide = dec_it->val<string>("radionuclide");
        Input::Array product_array = dec_it->val<Input::Array>("products");

        //resizing according to product count
        substance_ids_[i_decay].resize(product_array.size()+1);
        bifurcation_[i_decay].resize(product_array.size());
        
        //Reading products - setting substance ids, branching ratios.
        //iterating over products
        unsigned int i_product = 0;
        for (Input::Iterator<Input::Record> prod_it = product_array.begin<Input::Record>(); 
             prod_it != product_array.end(); ++prod_it, ++i_product)
        {   
            string product_name = prod_it->val<string>("name");
            idx = find_subst_name(product_name);
            if (idx < eq_data_base_->substances_.size())
                substance_ids_[i_decay][i_product+1] = idx;
            else THROW(ReactionTerm::ExcUnknownSubstance() 
                       << ReactionTerm::EI_Substance(product_name) 
                       << (*prod_it).ei_address());
            
            bifurcation_[i_decay][i_product] = prod_it->val<double>("branching_ratio");
        }
        
        // set radionuclide substance index
        idx = find_subst_name(radionuclide);
        if (idx < eq_data_base_->substances_.size())
            substance_ids_[i_decay][0] = idx;
        else THROW(ReactionTerm::ExcUnknownSubstance() 
                    << ReactionTerm::EI_Substance(radionuclide) 
                    << (*dec_it).ei_address());
            
        //Normalization of branching ratios in bifurcation vector by its norm.
        //sum:
        double sum=0;
        for(auto &b : bifurcation_[i_decay])
            sum += b;
        //Normalization:
        for(auto &b : bifurcation_[i_decay])
            b = b / sum;
    }
}


void RadioactiveDecay::assemble_ode_matrix(void )
{
    // create decay matrix
    reaction_matrix_ = arma::zeros(n_substances_, n_substances_);
    unsigned int reactant_index, product_index; //global indices of the substances
    double exponent;    //temporary variable k = ln(2)/t_half
    for (unsigned int i_decay = 0; i_decay < half_lives_.size(); i_decay++) {
        reactant_index = substance_ids_[i_decay][0];
        exponent = log(2) / half_lives_[i_decay];
        reaction_matrix_(reactant_index, reactant_index) = -exponent;
        
        for (unsigned int i_product = 1; i_product < substance_ids_[i_decay].size(); ++i_product){
            product_index = substance_ids_[i_decay][i_product];
            reaction_matrix_(product_index, reactant_index) = exponent * bifurcation_[i_decay][i_product-1];
        }
    }
}
