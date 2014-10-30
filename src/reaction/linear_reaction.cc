#include "reaction/linear_reaction.hh"
#include "reaction/reaction.hh"
#include "system/global_defs.h"

#include "reaction/linear_ode_solver.hh"
#include "la/distribution.hh"
#include "mesh/mesh.h"

using namespace std;
using namespace Input::Type;

Record LinearReaction::input_type_single_reaction
	= Record("Reaction", "Describes a single first order chemical reaction.")
	.declare_key("reactant", String(), Default::obligatory(),
				"The name of the reactant.")
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
    .declare_key("ode_solver", LinearODESolverBase::input_type, Default::optional(),
                 "Numerical solver for the system of first order ordinary differential equations coming from the model.");


LinearReaction::LinearReaction(Mesh &init_mesh, Input::Record in_rec)
      : LinearReactionBase(init_mesh, in_rec)
{
}

LinearReaction::~LinearReaction()
{
}

void LinearReaction::assemble_ode_matrix(void )
{
    // create decay matrix
    reaction_matrix_ = zeros(n_substances_, n_substances_);
    unsigned int reactant_index, product_index; //global indices of the substances
    double exponent;    //temporary variable for k
    for (unsigned int i_reaction = 0; i_reaction < reaction_rates_.size(); i_reaction++) {
        reactant_index = substance_ids_[i_reaction][0];
//         exponent = reaction_rates_[i_reaction] * time_->dt();
        exponent = reaction_rates_[i_reaction];
        reaction_matrix_(reactant_index, reactant_index) = -exponent;
        
        for (unsigned int i_product = 1; i_product < substance_ids_[i_reaction].size(); ++i_product){
            product_index = substance_ids_[i_reaction][i_product];
            reaction_matrix_(product_index, reactant_index) = exponent * bifurcation_[i_reaction][i_product-1];
        }
    }
}

// void LinearReaction::prepare_reaction_matrix_analytic(void)
// {
//     reaction_matrix_ = eye(n_substances_, n_substances_);
//     
//     unsigned int reactant_idx, product_idx,   // global indices of substances
//                  i_reaction, i_product;       // local indices of substances
//     double exponential; //temporary value for the exponential exp(-kt)
//     
//     // cycle over reactions/over rows/over parents
//     for (i_reaction = 0; i_reaction < reaction_rates_.size(); i_reaction++) {
//         // setting diagonal elements
//         reactant_idx = substance_ids_[i_reaction][0];
//         exponential = std::exp(- reaction_rates_[i_reaction] * time_->dt());
//         reaction_matrix_(reactant_idx,reactant_idx) = exponential;
// 
//         // cycle over products of specific reaction/row/reactant
//         for (i_product = 1; i_product < substance_ids_[i_reaction].size(); ++i_product) {
//             product_idx = substance_ids_[i_reaction][i_product];
//             reaction_matrix_(product_idx,reactant_idx)
//                                        = (1 - exponential)* bifurcation_[i_reaction][i_product-1];
//         }
//     }
// }


void LinearReaction::initialize_from_input()
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
		string parent_name = dec_it->val<string>("reactant");
		Input::Array product_array = dec_it->val<Input::Array>("products");
		Input::Array ratio_array = dec_it->val<Input::Array>("branching_ratios"); // has default value [ 1.0 ]

		// substance_ids contains also parent id
		if (product_array.size() > 0)   substance_ids_[i_reaction].resize( product_array.size()+1 );
		else    xprintf(UsrErr,"Empty array of products in the %d-th reaction.\n", i_reaction);


		// set parent index
		idx = find_subst_name(parent_name);
		if (idx < substances_.size())	substance_ids_[i_reaction][0] = idx;
		else    xprintf(UsrErr,"Unknown name of the reactant in the %d-th reaction.\n", i_reaction);

		// set products
		unsigned int i_product = 1;
		for(Input::Iterator<string> product_it = product_array.begin<string>(); 
            product_it != product_array.end(); ++product_it, i_product++)
		{
			idx = find_subst_name(*product_it);
			if (idx < substances_.size())   substance_ids_[i_reaction][i_product] = idx;
			else    xprintf(Warn,"Unknown name of the %d-th product in the %d-th reaction.\n", i_product-1 , i_reaction);
		}

		//bifurcation determining part
        if (ratio_array.size() == product_array.size() )   ratio_array.copy_to( bifurcation_[i_reaction] );
        else    xprintf(UsrErr,"Number of branches %d has to match the number of products %d in the %d-th reaction.\n",
                        ratio_array.size(), product_array.size(), i_reaction);

        //test the sum of branching ratios = 1.0
        double test_sum=0;
        for(auto &b : bifurcation_[i_reaction])
        {
            test_sum += b;
        }
        if(test_sum != 1.0)
            xprintf(UsrErr,"The sum of branching ratios %f in the %d-th reaction is not 1.0.\n",
                        test_sum, i_reaction);
	}
}
