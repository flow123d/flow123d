
#include "reaction/linear_ode_solver.hh"

#include "armadillo"
#include "input/accessors.hh"

using namespace Input::Type;

AbstractRecord & LinearODESolverBase::get_input_type() {
	static AbstractRecord type = AbstractRecord("LinearODESolver",
			"Solver of a linear system of ODEs.");
	type.close();
	return type;
}
    
LinearODESolverBase::LinearODESolverBase()
:step_(0), step_changed_(true)
{
}

LinearODESolverBase::~LinearODESolverBase()
{
}

void LinearODESolverBase::set_system_matrix(const arma::mat& matrix)
{
    system_matrix_ = matrix;
    step_changed_ = true;
}

void LinearODESolverBase::set_step(double step)
{
    step_ = step;
    step_changed_ = true;
}
