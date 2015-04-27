#include "reaction/first_order_reaction.hh"
#include "reaction/reaction_term.hh"
#include "reaction/linear_ode_solver.hh"

#include "system/global_defs.h"
#include "mesh/mesh.h"

using namespace Input::Type;

Record FirstOrderReaction::input_type_reactant
    = Record("FirstOrderReactionReactant", "A record describing a reactant of a reaction.")
    .allow_auto_conversion("name")
    .declare_key("name", String(), Default::obligatory(), 
                 "The name of the reactant.")
    //.declare_key("stoichiometric_coefficient", Integer(0.0), Default::optional(1.0))   //in future
    ;
    
Record FirstOrderReaction::input_type_product 
    = Record("FirstOrderReactionProduct", "A record describing a product of a reaction.")
    .allow_auto_conversion("name")
    .declare_key("name", String(), Default::obligatory(), 
                 "The name of the product.")
    //.declare_key("stoichiometric_coefficient", Integer(0.0), Default::optional(1.0))   //in future
    .declare_key("branching_ratio", Double(0.0), Default("1.0"),
                 "The branching ratio of the product when there are more products.\n"
                 "The value must be positive. Further, the branching ratios of all products are normalized "
                 "in order to sum to one.\n"
                 "The default value 1.0, should only be used in the case of single product.");
    
Record FirstOrderReaction::input_type_single_reaction
	= Record("Reaction", "Describes a single first order chemical reaction.")
	.declare_key("reactants", Array(FirstOrderReaction::input_type_reactant,1), Default::obligatory(),
				"An array of reactants. Do not use array, reactions with only one reactant (decays) are implemented at the moment!")
    .declare_key("reaction_rate", Double(0.0), Default::obligatory(),
                "The reaction rate coefficient of the first order reaction.")
    .declare_key("products", Array(FirstOrderReaction::input_type_product,1), Default::obligatory(),
				"An array of products.")
	;


Record FirstOrderReaction::input_type
	= Record("FirstOrderReaction", "A model of first order chemical reactions (decompositions of a reactant into products).")
	.derive_from( ReactionTerm::get_input_type() )
    .declare_key("reactions", Array( FirstOrderReaction::input_type_single_reaction), Default::obligatory(),
                "An array of first order chemical reactions.")
    .declare_key("ode_solver", LinearODESolverBase::get_input_type(), Default::optional(),
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
	for (Input::Iterator<Input::Record> react_it = reactions_array.begin<Input::Record>(); 
         react_it != reactions_array.end(); ++react_it, ++i_reaction)
	{ 
        //read reaction rate
        reaction_rates_[i_reaction] = react_it->val<double>("reaction_rate");
        
		//read reactant name, product names and branching ratios
        Input::Array reactant_array = react_it->val<Input::Array>("reactants");
		Input::Array product_array = react_it->val<Input::Array>("products");

		//resize substance_ids array
		substance_ids_[i_reaction].resize( product_array.size()+1 );
        bifurcation_[i_reaction].resize(product_array.size());

        if(reactant_array.size() != 1)
            xprintf(UsrErr, "More than one reactant is not available at the moment.");
        
        //take only one reactant
        Input::Iterator<Input::Record> reactant_it = reactant_array.begin<Input::Record>();
        {
            string reactant_name = reactant_it->val<string>("name");
            idx = find_subst_name(reactant_name);
            if (idx < substances_.size())   
                substance_ids_[i_reaction][0] = idx;
            else THROW(ReactionTerm::ExcUnknownSubstance() 
                    << ReactionTerm::EI_Substance(reactant_name) 
                    << (*reactant_it).ei_address());
        }
		
		// set products
		unsigned int i_product = 0;
		for(Input::Iterator<Input::Record> product_it = product_array.begin<Input::Record>(); 
            product_it != product_array.end(); ++product_it, i_product++)
		{
            string product_name = product_it->val<string>("name");
			idx = find_subst_name(product_name);
			if (idx < substances_.size())
                substance_ids_[i_reaction][i_product+1] = idx;
			else THROW(ReactionTerm::ExcUnknownSubstance() 
                        << ReactionTerm::EI_Substance(product_name) 
                        << product_array.ei_address());
            
            bifurcation_[i_reaction][i_product] = product_it->val<double>("branching_ratio");
		}


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
