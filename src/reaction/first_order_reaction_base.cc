#include "reaction/first_order_reaction_base.hh"
#include "reaction/reaction_term.hh"

#include "reaction/linear_ode_solver.hh"
#include "reaction/pade_approximant.hh"
#include "reaction/linear_ode_analytic.hh"


#include "system/global_defs.h"
#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "la/distribution.hh"

using namespace Input::Type;


FirstOrderReactionBase::FirstOrderReactionBase(Mesh &init_mesh, Input::Record in_rec)
    : ReactionTerm(init_mesh, in_rec)
{
	FLOW123D_FORCE_LINK_IN_PARENT(padeApproximant);
	FLOW123D_FORCE_LINK_IN_PARENT(linearODEAnalytic);

	Input::Iterator<Input::AbstractRecord> num_it = input_record_.find<Input::AbstractRecord>("ode_solver");
    if ( num_it )
    {
    	linear_ode_solver_ = (*num_it).factory< LinearODESolverBase, Input::Record >(*num_it);
    }
    else    //default linear ode solver
    {
        linear_ode_solver_ = Input::Factory< LinearODESolverBase, Input::Record >::instance()->create("LinearODEAnalytic", Input::Record());
    }
}

FirstOrderReactionBase::~FirstOrderReactionBase()
{
}

void FirstOrderReactionBase::initialize()
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


void FirstOrderReactionBase::zero_time_step()
{
    ASSERT(distribution_ != nullptr, "Distribution has not been set yet.\n");
    ASSERT(time_ != nullptr, "Time governor has not been set yet.\n");
    ASSERT_LESS(0, substances_.size());

    assemble_ode_matrix();
    // make scaling that takes into account different molar masses of substances
    reaction_matrix_ = molar_matrix_ * reaction_matrix_ * molar_mat_inverse_;
    linear_ode_solver_->set_system_matrix(reaction_matrix_);
}


double **FirstOrderReactionBase::compute_reaction(double **concentrations, int loc_el) //multiplication of concentrations array by reaction matrix
{      
    unsigned int rows;  // row in the concentration matrix, regards the substance index
    arma::vec new_conc;
    
    // save previous concentrations to column vector
    for(rows = 0; rows < n_substances_; rows++)
        prev_conc_(rows) = concentrations[rows][loc_el];
    
    // compute new concetrations R*c
    linear_ode_solver_->update_solution(prev_conc_, new_conc);
    
    // save new concentrations to the concentration matrix
    for(rows = 0; rows < n_substances_; rows++)
        concentrations[rows][loc_el] = new_conc(rows);
 
    return concentrations;
}

void FirstOrderReactionBase::update_solution(void)
{
    //DBGMSG("FirstOrderReactionBases - update solution\n");
    if(time_->is_changed_dt())
    {
        linear_ode_solver_->set_step(time_->dt());
    }

    START_TIMER("linear reaction step");

    for (unsigned int loc_el = 0; loc_el < distribution_->lsize(); loc_el++)
        this->compute_reaction(concentration_matrix_, loc_el);
    
    END_TIMER("linear reaction step");
}


unsigned int FirstOrderReactionBase::find_subst_name(const string &name)
{
    unsigned int k=0;
        for(; k < n_substances_; k++)
                if (name == substances_[k].name()) return k;

        return k;
}
