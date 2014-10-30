
#include "reaction/linear_ode_solver.hh"

#include "armadillo"
#include "input/accessors.hh"

using namespace Input::Type;

AbstractRecord LinearODESolverBase::input_type
    = AbstractRecord("NumericalMethod", "Numerical method used in reaction computation.");

Record LinearODEAnalytic::input_type
    = Record("LinearODEAnalytic", "Evaluate analytic solution of the system of ODEs.")
    .derive_from(LinearODESolverBase::input_type);
    
LinearODESolverBase::LinearODESolverBase()
:step_(0), step_changed_(true)
{
}

LinearODESolverBase::~LinearODESolverBase()
{
}

void LinearODESolverBase::set_system_matrix(const mat& matrix)
{
    system_matrix_ = matrix;
    step_changed_ = true;
}

void LinearODESolverBase::set_step(double step)
{
    step_ = step;
    step_changed_ = true;
}


/********************************* LinearODEAnalytic implementation ****************************************/

LinearODEAnalytic::LinearODEAnalytic(Input::Record in_rec)
{
}

LinearODEAnalytic::~LinearODEAnalytic()
{
}

void LinearODEAnalytic::set_bifurcation_matrix(const mat& bifurcation)
{
    bifurcation_matrix_ = bifurcation;
}

void LinearODEAnalytic::update_solution(vec& init_vector, vec& output_vec)
{
    if(step_changed_)
    {
        compute_matrix();
        step_changed_ = false;
    }
    
    output_vec = solution_matrix_ * init_vector;
}

// void LinearODEAnalytic::update_solution(mat& init_vectors, mat& output_vecs, 
//                                         const std::vector< unsigned int >& mask)
// {
//     ASSERT(0, "Method must be implemented in the template class.");
// }

void LinearODEAnalytic::compute_matrix()
{
    ASSERT(system_matrix_.n_cols == system_matrix_.n_rows, "Matrix is not square.");
    ASSERT(bifurcation_matrix_.n_rows == system_matrix_.n_rows, "Bifurcation matrix dimensions are wrong.");
    ASSERT(bifurcation_matrix_.n_cols == system_matrix_.n_cols, "Bifurcation matrix dimensions are wrong.");
    solution_matrix_.copy_size(system_matrix_);
    
    double exponential;
    for(unsigned int i = 0; i < solution_matrix_.n_rows; i++)
    {
        exponential = std::exp(system_matrix_(i,i) * step_);
        for(unsigned int j = 0; j < solution_matrix_.n_cols; j++)
        {
            solution_matrix_(j,i) = (1-exponential)*bifurcation_matrix_(i,j);
        }
        solution_matrix_(i,i) = exponential;
    }
    
//     bifurcation_matrix_.print();
//     system_matrix_.print();
//     solution_matrix_.print();
}
