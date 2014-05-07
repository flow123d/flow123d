#include <math.h>

#include "reaction/reaction.hh"
#include "reaction/linear_reaction.hh"
#include "reaction/pade_approximant.hh"
#include "system/global_defs.h"
#include "system/sys_profiler.hh"

#include "la/distribution.hh"
#include "mesh/mesh.h"

using namespace std;
using namespace Input::Type;

Record LinearReaction::input_type_one_decay_substep
	= Record("Substep", "Equation for reading information about radioactive decays.")
	.declare_key("parent", String(), Default::obligatory(),
				"Identifier of an isotope.")
    .declare_key("half_life", Double(), Default::optional(),
                "Half life of the parent substance.")
    .declare_key("kinetic", Double(), Default::optional(),
                "Kinetic constants describing first order reactions.")
    .declare_key("products", Array(String()), Default::obligatory(),
				"Identifies isotopes which decays parental atom to.")
	.declare_key("branch_ratios", Array(Double()), Default("1.0"),   // default is one product, with ratio == 1.0
				"Decay chain branching percentage.");


Record LinearReaction::input_type
	= Record("LinearReactions", "Information for a decision about the way to simulate radioactive decay.")
	.derive_from( ReactionTerm::input_type )
    .declare_key("decays", Array( LinearReaction::input_type_one_decay_substep ), Default::obligatory(),
                "Description of particular decay chain substeps.");




LinearReaction::LinearReaction(Mesh &init_mesh, Input::Record in_rec)
      : ReactionTerm(init_mesh, in_rec)
{
}

LinearReaction::~LinearReaction()
{
}

void LinearReaction::initialize()
{
    ASSERT(distribution_ != nullptr, "Distribution has not been set yet.\n");
    ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
    ASSERT_LESS(0, names_.size());
    
    n_substances_ = names_.size();
    initialize_from_input();

    // allocation
    prev_conc_.resize(n_substances_);
    reaction_matrix_.resize(n_substances_);
    for(unsigned int i_subst = 0; i_subst < n_substances_; i_subst++)
        reaction_matrix_[i_subst].resize(n_substances_);
}


void LinearReaction::zero_time_step()
{
  ASSERT(distribution_ != nullptr, "Distribution has not been set yet.\n");
  ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
  ASSERT_LESS(0, names_.size());

  //nothing is to be computed at zero_time_step
}


void LinearReaction::reset_reaction_matrix()
{
    unsigned int rows, cols;
    for(rows = 0; rows < n_substances_;rows++){
        for(cols = 0; cols < n_substances_; cols++){
         if(rows == cols)   reaction_matrix_[rows][cols] = 1.0;
         else               reaction_matrix_[rows][cols] = 0.0;
     }
    }
}


void LinearReaction::modify_reaction_matrix(void) //All the parameters are supposed to be known
{
    ASSERT(reaction_matrix_.size() > 0, "Reaction matrix is not allocated.\n");
    
    unsigned int index_par, i_decay, i_product;
    double rel_step;    

    for (i_decay = 0; i_decay < half_lives_.size(); i_decay++) {
        index_par = substance_ids_[i_decay][0];
        rel_step = time_->dt() / half_lives_[i_decay];
        reaction_matrix_[index_par][index_par] = pow(0.5, rel_step);

        for (i_product = 1; i_product < substance_ids_[i_decay].size(); ++i_product)
            reaction_matrix_[ substance_ids_[i_decay][i_product] ][index_par]
                                       = (1 - pow(0.5, rel_step))* bifurcation_[i_decay][i_product-1];
    }
    
    //print_half_lives();
    //print_reaction_matrix(); //just for control print
}

double **LinearReaction::compute_reaction(double **concentrations, int loc_el) //multiplication of concentrations array by reaction matrix
{
    unsigned int cols, rows;

    for(rows = 0; rows < names_.size(); rows++){
        prev_conc_[rows] = concentrations[rows][loc_el];
        concentrations[rows][loc_el] = 0.0;
    }
    
    for(rows = 0; rows < n_substances_; rows++){
        for(cols = 0; cols < n_substances_; cols++){
            concentrations[rows][loc_el] += reaction_matrix_[rows][cols]*prev_conc_[cols];
        }
    }

    return concentrations;
}

//       raise warning if sum of ratios is not one
void LinearReaction::initialize_from_input()
{
    unsigned int idx;

	Input::Array decay_array = input_record_.val<Input::Array>("decays");

	substance_ids_.resize( decay_array.size() );
	half_lives_.resize( decay_array.size() );
	bifurcation_.resize( decay_array.size() );

	int i_decay=0;
	for (Input::Iterator<Input::Record> dec_it = decay_array.begin<Input::Record>(); dec_it != decay_array.end(); ++dec_it, ++i_decay)
	{
		//half-lives determining part
		Input::Iterator<double> it_hl = dec_it->find<double>("half_life");
		if (it_hl) {
		   half_lives_[i_decay] = *it_hl;
		} else {
		   it_hl = dec_it->find<double>("kinetic");
		   if (it_hl) {
			   half_lives_[i_decay] = log(2)/(*it_hl);
		   } else {
		    xprintf(UsrErr, "Missing half-life or kinetic in the %d-th reaction.\n", i_decay);
		  }
		}

		//indices determining part
		string parent_name = dec_it->val<string>("parent");
		Input::Array product_array = dec_it->val<Input::Array>("products");
		Input::Array ratio_array = dec_it->val<Input::Array>("branch_ratios"); // has default value [ 1.0 ]

		// substance_ids contains also parent id
		if (product_array.size() > 0)   substance_ids_[i_decay].resize( product_array.size()+1 );
		else			xprintf(UsrErr,"Empty array of products in the %d-th reaction.\n", i_decay);


		// set parent index
		idx = find_subst_name(parent_name);
		if (idx < names_.size())	substance_ids_[i_decay][0] = idx;
		else                		xprintf(UsrErr,"Wrong name of parent substance in the %d-th reaction.\n", i_decay);

		// set products
		unsigned int i_product = 1;
		for(Input::Iterator<string> product_it = product_array.begin<string>(); product_it != product_array.end(); ++product_it, i_product++)
		{
			idx = find_subst_name(*product_it);
			if (idx < names_.size())   substance_ids_[i_decay][i_product] = idx;
			else                    	xprintf(Warn,"Wrong name of %d-th product in the %d-th reaction.\n", i_product-1 , i_decay);
		}

		//bifurcation determining part
        if (ratio_array.size() == product_array.size() )   ratio_array.copy_to( bifurcation_[i_decay] );
        else            xprintf(UsrErr,"Number of branches %d has to match number of products %d in the %d-th reaction.\n",
                                       ratio_array.size(), product_array.size(), i_decay);

	}
}

void LinearReaction::update_solution(void)
{
    //DBGMSG("LinearReactions - update solution\n");
    if(time_->is_changed_dt())
    {
        reset_reaction_matrix();
        modify_reaction_matrix();
    }

    START_TIMER("linear reaction step");
    
	for (unsigned int loc_el = 0; loc_el < distribution_->lsize(); loc_el++)
        this->compute_reaction(concentration_matrix_, loc_el);
    
    END_TIMER("linear reaction step");
}


void LinearReaction::print_reaction_matrix(void)
{
	unsigned int cols,rows;

	//DBGMSG("r mat: %p\n", reaction_matrix);
	if(reaction_matrix_.size() == n_substances_){
                if(time_ != NULL)
                  xprintf(Msg,"\ntime_step %f,Reaction matrix looks as follows:\n",time_->dt());
		for(rows = 0; rows < n_substances_; rows++){
			for(cols = 0; cols < n_substances_; cols++){
					xprintf(Msg,"%f\t",reaction_matrix_[rows][cols]);
			}
			xprintf(Msg,"\n");
		}
	}else{
		xprintf(Msg,"\nReaction matrix needs to be allocated.\n");
	}
	return;
}

void LinearReaction::print_half_lives() {
    unsigned int i;

    xprintf(Msg, "\nHalf-lives are defined as:\n");
    for (i = 0; i < half_lives_.size(); i++) {
            xprintf(Msg, "parent_id: %d half_life: %f\n",substance_ids_[i][0], half_lives_[i]);
    }
}

unsigned int LinearReaction::find_subst_name(const string &name)
{

    unsigned int k=0;
        for(; k < names_.size(); k++)
                if (name == names_[k]) return k;

        return k;
}
