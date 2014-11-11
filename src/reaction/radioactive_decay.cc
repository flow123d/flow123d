#include "reaction/radioactive_decay.hh"
#include "reaction/reaction_term.hh"
#include "reaction/first_order_reaction_base.hh"

#include "system/global_defs.h"
#include "mesh/mesh.h"

#include "armadillo"

using namespace arma;
using namespace Input::Type;

Record RadioactiveDecay::input_type_single_decay
    = Record("Decay", "A model of a radioactive decay.")
    .declare_key("radionuclide", String(), Default::obligatory(),
                "The name of the parent radionuclide.")
    .declare_key("half_life", Double(), Default::obligatory(),
                 "The half life of the parent radionuclide in seconds.")
    .declare_key("products", Array(String()), Default::obligatory(),
                "An array of the decay products (daughters).")
    .declare_key("energy", Array(Double()), Default::optional(),
                "Not used now but it can be used for source in the thermal convection-diffusion model in future. "
                "An array of energy in MeV which is released from the parent radionuclide to a specific product.")
    .declare_key("branching_ratios", Array(Double()), Default("1.0"),  //default is one product, with ratio = 1.0
                "This is an array of branching ratios when more than one product is present. "
                "Considering only one product, the default ratio 1.0 is used. "
                "The values are given as fractions in interval [0.0,1.0] and "
                "their sum must be 1.0; it is checked during input reading.");

Record RadioactiveDecay::input_type
    = Record("RadioactiveDecay", "A model of a radioactive decay and possibly of a decay chain.")
    .derive_from( ReactionTerm::input_type )
    .declare_key("decays", Array( RadioactiveDecay::input_type_single_decay), Default::obligatory(),
                "An array of radioactive decays. They can make a chain.")
    .declare_key("ode_solver", LinearODESolverBase::input_type, Default::optional(),
                 "Numerical solver for the system of first order ordinary differential equations coming from the model.");



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
    int i_decay=0;
    for (Input::Iterator<Input::Record> dec_it = decay_array.begin<Input::Record>(); dec_it != decay_array.end(); ++dec_it, ++i_decay)
    {
        //read half-life
        half_lives_[i_decay] = dec_it->val<double>("half_life");

        //radionuclide name, product name array and branching ratio array
        string radionuclide = dec_it->val<string>("radionuclide");
        Input::Array product_array = dec_it->val<Input::Array>("products");
        Input::Array ratio_array = dec_it->val<Input::Array>("branching_ratios"); // has default value [ 1.0 ]

        // substance_ids contains also parent id
        if (product_array.size() > 0)   substance_ids_[i_decay].resize( product_array.size()+1 );
        else    xprintf(UsrErr,"Empty array of products in the %d-th reaction.\n", i_decay);


        // set radionuclide substance index
        idx = find_subst_name(radionuclide);
        if (idx < substances_.size())    substance_ids_[i_decay][0] = idx;
        else    xprintf(UsrErr,"Unknown name of the radionuclide in the %d-th reaction.\n", i_decay);

        // set products
        unsigned int i_product = 1;
        for(Input::Iterator<string> product_it = product_array.begin<string>(); product_it != product_array.end(); ++product_it, i_product++)
        {
            idx = find_subst_name(*product_it);
            if (idx < substances_.size())   substance_ids_[i_decay][i_product] = idx;
            else                        xprintf(Warn,"Unknown name of the %d-th product in the %d-th reaction.\n", i_product-1 , i_decay);
        }

        //set branching ratio array
        if (ratio_array.size() == product_array.size() )   ratio_array.copy_to( bifurcation_[i_decay] );
        else            xprintf(UsrErr,"Number of branches %d has to match number of products %d in the %d-th reaction.\n",
                                       ratio_array.size(), product_array.size(), i_decay);
        
        //test the sum of branching ratios = 1.0
        double test_sum=0;
        for(auto &b : bifurcation_[i_decay])
        {
            test_sum += b;
        }
        if(test_sum != 1.0)
            xprintf(UsrErr,"The sum of branching ratios %f in the %d-th reaction is not 1.0.\n",
                        test_sum, i_decay);
    }
}


void RadioactiveDecay::assemble_ode_matrix(void )
{
    // create decay matrix
    reaction_matrix_ = zeros(n_substances_, n_substances_);
    unsigned int reactant_index, product_index; //global indices of the substances
    double exponent;    //temporary variable k = ln(2)/t_half
    for (unsigned int i_decay = 0; i_decay < half_lives_.size(); i_decay++) {
        reactant_index = substance_ids_[i_decay][0];
        //exponent = log(2) * time_->dt() / half_lives_[i_decay];
        exponent = log(2) / half_lives_[i_decay];
        reaction_matrix_(reactant_index, reactant_index) = -exponent;
        
        for (unsigned int i_product = 1; i_product < substance_ids_[i_decay].size(); ++i_product){
            product_index = substance_ids_[i_decay][i_product];
            reaction_matrix_(product_index, reactant_index) = exponent * bifurcation_[i_decay][i_product-1];
        }
    }
}