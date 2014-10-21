#include "reaction/linear_reaction.hh"
#include "reaction/reaction.hh"
#include "system/global_defs.h"

#include "la/distribution.hh"
#include "mesh/mesh.h"

using namespace std;
using namespace Input::Type;

Record LinearReaction::input_type_single_reaction
	= Record("Reaction", "Describes a single first order chemical reaction.")
	.declare_key("reactant", String(), Default::obligatory(),
				"The name of the reactant.")
//     .declare_key("half_life", Double(), Default::optional(),
//                 "The half life of the reactant in seconds")
    .declare_key("reaction_rate", Double(), Default::obligatory(),
                "The reaction rate coefficient of the first order reaction.")
    .declare_key("products", Array(String()), Default::obligatory(),
				"An array of the names of products.")
	.declare_key("branching_ratios", Array(Double()), Default("1.0"),   // default is one product, with ratio == 1.0
				"This is an array of branching ratios when more than one product is present. "
                "Considering only one product, the default ratio 1.0 is used."
                "The values are given as fractions in interval [0.0,1.0] and"
                "their sum must be 1.0; it is checked during input reading.");


Record LinearReaction::input_type
	= Record("FirstOrderReaction", "A model of first order chemical reactions (decompositions of a reactant into products).")
	.derive_from( ReactionTerm::input_type )
    .declare_key("reactions", Array( LinearReaction::input_type_single_reaction), Default::obligatory(),
                "An array of first order chemical reactions.")
    .declare_key("ode_solver", NumericalMethod::input_type, Default::optional(),
                 "Numerical solver for the system of first order ordinary differential equations coming from the model.");


LinearReaction::LinearReaction(Mesh &init_mesh, Input::Record in_rec)
      : LinearReactionBase(init_mesh, in_rec)
{
}

LinearReaction::~LinearReaction()
{
}

void LinearReaction::prepare_reaction_matrix(void )
{
    // create decay matrix
    reaction_matrix_ = zeros(n_substances_, n_substances_);
    unsigned int reactant_index, product_index; //global indices of the substances
    double exponent;    //temporary variable
    for (unsigned int i_decay = 0; i_decay < half_lives_.size(); i_decay++) {
        reactant_index = substance_ids_[i_decay][0];
        exponent = log(2) * time_->dt() / half_lives_[i_decay];
        reaction_matrix_(reactant_index, reactant_index) = -exponent;
        
        for (unsigned int i_product = 1; i_product < substance_ids_[i_decay].size(); ++i_product){
            product_index = substance_ids_[i_decay][i_product];
            reaction_matrix_(product_index, reactant_index) = exponent * bifurcation_[i_decay][i_product-1];
        }
    }
    //DBGMSG("reactions matrix prepared\n");
    //reaction_matrix_.print();
}

void LinearReaction::prepare_reaction_matrix_analytic(void)
{
    reaction_matrix_ = eye(n_substances_, n_substances_);
    
    unsigned int parent_idx, product_idx,   // global indices of substances
                 i_decay, i_product;        // local indices of substances
    double relative_timestep,   // exponent of 0.5
           temp_power;          // temporary power of 0.5
    
    // cycle over reactions/over rows/over parents
    for (i_decay = 0; i_decay < half_lives_.size(); i_decay++) {
        // setting diagonal elements
        parent_idx = substance_ids_[i_decay][0];
        relative_timestep = time_->dt() / half_lives_[i_decay];
        temp_power = pow(0.5, relative_timestep);
        reaction_matrix_(parent_idx,parent_idx) = temp_power;

        // cycle over products of specific reaction/row/parent
        for (i_product = 1; i_product < substance_ids_[i_decay].size(); ++i_product) {
            product_idx = substance_ids_[i_decay][i_product];
            reaction_matrix_(product_idx,parent_idx)
                                       = (1 - temp_power)* bifurcation_[i_decay][i_product-1];
        }
    }
}



//       raise warning if sum of ratios is not one
void LinearReaction::initialize_from_input()
{
    unsigned int idx;

	Input::Array decay_array = input_record_.val<Input::Array>("reactions");

	substance_ids_.resize( decay_array.size() );
	half_lives_.resize( decay_array.size() );
	bifurcation_.resize( decay_array.size() );

	int i_decay=0;
	for (Input::Iterator<Input::Record> dec_it = decay_array.begin<Input::Record>(); dec_it != decay_array.end(); ++dec_it, ++i_decay)
	{
// 		//half-lives determining part
// 		Input::Iterator<double> it_hl = dec_it->find<double>("half_life");
// 		if (it_hl) {
// 		   half_lives_[i_decay] = *it_hl;
// 		} else {
// 		   it_hl = dec_it->find<double>("kinetic");
// 		   if (it_hl) {
// 			   half_lives_[i_decay] = log(2)/(*it_hl);
// 		   } else {
// 		    xprintf(UsrErr, "Missing half-life or kinetic in the %d-th reaction.\n", i_decay);
// 		  }
// 		}
        
        //TODO:have array of reaction rates, not half lives
        //half-lives determining part
        double reaction_rate = dec_it->val<double>("reaction_rate");
        half_lives_[i_decay] = log(2)/(reaction_rate);
        
		//indices determining part
		string parent_name = dec_it->val<string>("reactant");
		Input::Array product_array = dec_it->val<Input::Array>("products");
		Input::Array ratio_array = dec_it->val<Input::Array>("branching_ratios"); // has default value [ 1.0 ]

		// substance_ids contains also parent id
		if (product_array.size() > 0)   substance_ids_[i_decay].resize( product_array.size()+1 );
		else			xprintf(UsrErr,"Empty array of products in the %d-th reaction.\n", i_decay);


		// set parent index
		idx = find_subst_name(parent_name);
		if (idx < substances_.size())	substance_ids_[i_decay][0] = idx;
		else                		xprintf(UsrErr,"Wrong name of parent substance in the %d-th reaction.\n", i_decay);

		// set products
		unsigned int i_product = 1;
		for(Input::Iterator<string> product_it = product_array.begin<string>(); product_it != product_array.end(); ++product_it, i_product++)
		{
			idx = find_subst_name(*product_it);
			if (idx < substances_.size())   substance_ids_[i_decay][i_product] = idx;
			else                    	xprintf(Warn,"Wrong name of %d-th product in the %d-th reaction.\n", i_product-1 , i_decay);
		}

		//bifurcation determining part
        if (ratio_array.size() == product_array.size() )   ratio_array.copy_to( bifurcation_[i_decay] );
        else            xprintf(UsrErr,"Number of branches %d has to match the number of products %d in the %d-th reaction.\n",
                                       ratio_array.size(), product_array.size(), i_decay);

	}
}

void LinearReaction::print_half_lives() {
    unsigned int i;

    xprintf(Msg, "\nHalf-lives are defined as:\n");
    for (i = 0; i < half_lives_.size(); i++) {
            xprintf(Msg, "parent_id: %d half_life: %f\n",substance_ids_[i][0], half_lives_[i]);
    }
}
