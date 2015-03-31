
#include "reaction/linear_ode_analytic.hh"
#include "input/accessors.hh"
#include "input/factory.hh"
#include "system/sys_profiler.hh"

using namespace Input::Type;

Record LinearODEAnalytic::input_type
    = Record("LinearODEAnalytic", "Evaluate analytic solution of the system of ODEs.")
    .derive_from(LinearODESolverBase::input_type);
    
const int LinearODEAnalytic::registrar =
		Input::register_class< LinearODEAnalytic, Input::Record >("LinearODEAnalytic");

LinearODEAnalytic::LinearODEAnalytic(Input::Record in_rec)
{
}

LinearODEAnalytic::~LinearODEAnalytic()
{
}

void LinearODEAnalytic::update_solution(arma::vec& init_vector, arma::vec& output_vec)
{
    if(step_changed_)
    {
        compute_matrix();
        step_changed_ = false;
    }
    
    output_vec = solution_matrix_ * init_vector;
}

void LinearODEAnalytic::compute_matrix()
{
    START_TIMER("ODEAnalytic::compute_matrix");

    ASSERT(system_matrix_.n_cols == system_matrix_.n_rows, "Matrix is not square.");
    solution_matrix_.copy_size(system_matrix_);
    
    double exponential, temp;
    for(unsigned int i = 0; i < solution_matrix_.n_rows; i++)
    {
        exponential = std::exp(system_matrix_(i,i) * step_);
        for(unsigned int j = 0; j < solution_matrix_.n_cols; j++)
        {
            temp = 1.0;
            if( (i != j) && (system_matrix_(i,i) != 0.0) )
                temp = -system_matrix_(j,i)/system_matrix_(i,i);
            
            solution_matrix_(j,i) = (1-exponential)*temp;
        }
        solution_matrix_(i,i) = exponential;
    }
}
