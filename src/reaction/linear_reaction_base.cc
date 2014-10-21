#include "reaction/linear_reaction_base.hh"
#include "reaction/reaction.hh"

#include "reaction/pade_approximant.hh"

#include "system/global_defs.h"
#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "la/distribution.hh"

using namespace Input::Type;
using namespace arma;


LinearReactionBase::LinearReactionBase(Mesh &init_mesh, Input::Record in_rec)
    : ReactionTerm(init_mesh, in_rec)
{
    Input::Iterator<Input::AbstractRecord> num_it = input_record_.find<Input::AbstractRecord>("ode_solver");
    if ( num_it )
    {
        if (num_it->type() == PadeApproximant::input_type) 
        {
            pade_approximant_ = new PadeApproximant(*num_it);
            numerical_method_ = NumericalMethod::pade_approximant;
        }
    }
    else
        numerical_method_ = NumericalMethod::analytic;  //no numerical method, use analytic solution
}

LinearReactionBase::~LinearReactionBase()
{
}

void LinearReactionBase::initialize_from_input()
{
    xprintf(Warn, "The method initialize_from_input() should be reimplemented in descendants.");
}

void LinearReactionBase::initialize()
{
    ASSERT(distribution_ != nullptr, "Distribution has not been set yet.\n");
    ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
    ASSERT_LESS(0, substances_.size());
    
    n_substances_ = substances_.size();
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
    	molar_matrix_(i,i) = substances_[i].molar_mass();
    	molar_mat_inverse_(i,i) = 1./substances_[i].molar_mass();
    }
}


void LinearReactionBase::zero_time_step()
{
    ASSERT(distribution_ != nullptr, "Distribution has not been set yet.\n");
    ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
    ASSERT_LESS(0, substances_.size());

    //nothing is to be computed at zero_time_step
}

void LinearReactionBase::compute_reaction_matrix(void )
{
    switch(numerical_method_)
    {
    case NumericalMethod::analytic:
    	prepare_reaction_matrix_analytic();
    	break;

    case NumericalMethod::pade_approximant:
    	prepare_reaction_matrix();
    	pade_approximant_->approximate_matrix(reaction_matrix_);
    	break;

    default:
    	prepare_reaction_matrix_analytic();
    }
    
    //TODO: do this on the original matrix of the system
    // make scaling that takes into account different molar masses of substances
    reaction_matrix_ = molar_matrix_ * reaction_matrix_ * molar_mat_inverse_;
}


double **LinearReactionBase::compute_reaction(double **concentrations, int loc_el) //multiplication of concentrations array by reaction matrix
{      
    unsigned int rows;  // row in the concentration matrix, regards the substance index
    
    // save previous concentrations to column vector
    for(rows = 0; rows < n_substances_; rows++)
        prev_conc_(rows) = concentrations[rows][loc_el];
    
    // compute new concetrations R*c
    vec new_conc = reaction_matrix_ * prev_conc_;
    
    // save new concentrations to the concentration matrix
    for(rows = 0; rows < n_substances_; rows++)
        concentrations[rows][loc_el] = new_conc(rows);
 
    return concentrations;
}

//       raise warning if sum of ratios is not one
// void LinearReactionBase::initialize_from_input()
// {
//     unsigned int idx;
// 
//     Input::Array decay_array = input_record_.val<Input::Array>("decays");
// 
//     substance_ids_.resize( decay_array.size() );
//     half_lives_.resize( decay_array.size() );
//     bifurcation_.resize( decay_array.size() );
// 
//     int i_decay=0;
//     for (Input::Iterator<Input::Record> dec_it = decay_array.begin<Input::Record>(); dec_it != decay_array.end(); ++dec_it, ++i_decay)
//     {
//         //half-lives determining part
//         Input::Iterator<double> it_hl = dec_it->find<double>("half_life");
//         if (it_hl) {
//            half_lives_[i_decay] = *it_hl;
//         } else {
//            it_hl = dec_it->find<double>("kinetic");
//            if (it_hl) {
//                half_lives_[i_decay] = log(2)/(*it_hl);
//            } else {
//             xprintf(UsrErr, "Missing half-life or kinetic in the %d-th reaction.\n", i_decay);
//           }
//         }
// 
//         //indices determining part
//         string parent_name = dec_it->val<string>("parent");
//         Input::Array product_array = dec_it->val<Input::Array>("products");
//         Input::Array ratio_array = dec_it->val<Input::Array>("branch_ratios"); // has default value [ 1.0 ]
// 
//         // substance_ids contains also parent id
//         if (product_array.size() > 0)   substance_ids_[i_decay].resize( product_array.size()+1 );
//         else            xprintf(UsrErr,"Empty array of products in the %d-th reaction.\n", i_decay);
// 
// 
//         // set parent index
//         idx = find_subst_name(parent_name);
//         if (idx < names_.size())    substance_ids_[i_decay][0] = idx;
//         else                        xprintf(UsrErr,"Wrong name of parent substance in the %d-th reaction.\n", i_decay);
// 
//         // set products
//         unsigned int i_product = 1;
//         for(Input::Iterator<string> product_it = product_array.begin<string>(); product_it != product_array.end(); ++product_it, i_product++)
//         {
//             idx = find_subst_name(*product_it);
//             if (idx < names_.size())   substance_ids_[i_decay][i_product] = idx;
//             else                        xprintf(Warn,"Wrong name of %d-th product in the %d-th reaction.\n", i_product-1 , i_decay);
//         }
// 
//         //bifurcation determining part
//         if (ratio_array.size() == product_array.size() )   ratio_array.copy_to( bifurcation_[i_decay] );
//         else            xprintf(UsrErr,"Number of branches %d has to match number of products %d in the %d-th reaction.\n",
//                                        ratio_array.size(), product_array.size(), i_decay);
// 
//     }
// }

void LinearReactionBase::update_solution(void)
{
    //DBGMSG("LinearReactionBases - update solution\n");
    if(time_->is_changed_dt())
    {
        compute_reaction_matrix();
        //reaction_matrix_.print();
    }

    START_TIMER("linear reaction step");
    
    for (unsigned int loc_el = 0; loc_el < distribution_->lsize(); loc_el++)
        this->compute_reaction(concentration_matrix_, loc_el);
    
    END_TIMER("linear reaction step");
}


unsigned int LinearReactionBase::find_subst_name(const string &name)
{
    unsigned int k=0;
        for(; k < n_substances_; k++)
                if (name == substances_[k].name()) return k;

        return k;
}
