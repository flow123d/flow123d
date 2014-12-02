#include "reaction/first_order_reaction.hh"
#include "reaction/reaction_term.hh"
#include "reaction/linear_ode_solver.hh"

#include "system/global_defs.h"
#include "mesh/mesh.h"

using namespace Input::Type;

Record FirstOrderReaction::input_type_single_reaction
	= Record("Reaction", "Describes a single first order chemical reaction.")
	.declare_key("reactant", String(), Default::obligatory(),
				"The name of the reactant.")
    .declare_key("reaction_rate", Double(0.0), Default::obligatory(),
                "The reaction rate coefficient of the first order reaction.")
    .declare_key("products", Array(String(),1), Default::obligatory(),
				"An array of the names of products.")
	.declare_key("branching_ratios", Array(Double(0.0)), Default("1.0"),   // default is one product, with ratio == 1.0
                 "This is an array of branching ratios of the products when there is more than one present."
                 "Considering only one product, the default ratio 1.0 is used."
                 "Its value must be positive. Further, the branching ratios of all products are normalized" 
                 "by their sum, so the sum then gives 1.0 (this also resolves possible rouding errors).");


Record FirstOrderReaction::input_type
	= Record("FirstOrderReaction", "A model of first order chemical reactions (decompositions of a reactant into products).")
	.derive_from( ReactionTerm::input_type )
    .declare_key("reactions", Array( FirstOrderReaction::input_type_single_reaction), Default::obligatory(),
                "An array of first order chemical reactions.")
    .declare_key("ode_solver", LinearODESolverBase::input_type, Default::optional(),
                 "Numerical solver for the system of first order ordinary differential equations coming from the model.");


FirstOrderReaction::FirstOrderReaction(Mesh &init_mesh, Input::Record in_rec)
      : FirstOrderReactionBase(init_mesh, in_rec)
{
}

FirstOrderReaction::~FirstOrderReaction()
{
}

void FirstOrderReaction::assemble_ode_matrix(void )
{
    // create decay matrix
    reaction_matrix_ = arma::zeros(n_substances_, n_substances_);
    unsigned int reactant_index, product_index; //global indices of the substances
    double exponent;    //temporary variable for k
    for (unsigned int i_reaction = 0; i_reaction < reaction_rates_.size(); i_reaction++) {
        reactant_index = substance_ids_[i_reaction][0];
        exponent = reaction_rates_[i_reaction];
        reaction_matrix_(reactant_index, reactant_index) = -exponent;
        
        for (unsigned int i_product = 1; i_product < substance_ids_[i_reaction].size(); ++i_product){
            product_index = substance_ids_[i_reaction][i_product];
            reaction_matrix_(product_index, reactant_index) = exponent * bifurcation_[i_reaction][i_product-1];
        }
    }
}


void FirstOrderReaction::initialize_from_input()
{
    unsigned int idx;   //temporary variable, indexing substances

	Input::Array reactions_array = input_record_.val<Input::Array>("reactions");

	substance_ids_.resize( reactions_array.size() );
    reaction_rates_.resize( reactions_array.size() );
	bifurcation_.resize( reactions_array.size() );

	int i_reaction=0;
	for (Input::Iterator<Input::Record> dec_it = reactions_array.begin<Input::Record>(); 
         dec_it != reactions_array.end(); ++dec_it, ++i_reaction)
	{ 
        //read reaction rate
        reaction_rates_[i_reaction] = dec_it->val<double>("reaction_rate");
        
		//read reactant name, product names and branching ratios
		string reactant_name = dec_it->val<string>("reactant");
		Input::Array product_array = dec_it->val<Input::Array>("products");
		Input::Array ratio_array = dec_it->val<Input::Array>("branching_ratios"); // has default value [ 1.0 ]

		//resize substance_ids array
		substance_ids_[i_reaction].resize( product_array.size()+1 );

		// set reactant index
		idx = find_subst_name(reactant_name);
		if (idx < substances_.size())	
            substance_ids_[i_reaction][0] = idx;
		else THROW(ReactionTerm::ExcUnknownSubstance() 
                    << ReactionTerm::EI_Substance(reactant_name) 
                    << (*dec_it).ei_address());

		// set products
		unsigned int i_product = 1;
		for(Input::Iterator<string> product_it = product_array.begin<string>(); 
            product_it != product_array.end(); ++product_it, i_product++)
		{
			idx = find_subst_name(*product_it);
			if (idx < substances_.size())
                substance_ids_[i_reaction][i_product] = idx;
			else THROW(ReactionTerm::ExcUnknownSubstance() 
                        << ReactionTerm::EI_Substance(*product_it) 
                        << product_array.ei_address());
		}

		//bifurcation determining part
        if (ratio_array.size() == product_array.size() )   ratio_array.copy_to( bifurcation_[i_reaction] );
        else    xprintf(UsrErr,"Number of branches %d has to match the number of products %d in the %d-th reaction.\n",
                        ratio_array.size(), product_array.size(), i_reaction);

        //Normalization of branching ratios in bifurcation vector by its norm.
        //sum:
        double sum=0;
        for(auto &b : bifurcation_[i_reaction])
            sum += b;
        //Normalization:
        for(auto &b : bifurcation_[i_reaction])
            b = b / sum;
	}
}
