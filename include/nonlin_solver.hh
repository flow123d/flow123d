/*
 * nonlin_solver.hh
 *
 *  Created on: Dec 27, 2011
 *      Author: jb
 */

#ifndef NONLIN_SOLVER_HH_
#define NONLIN_SOLVER_HH_

/**
 * Experimental Homotopy solver.
 *
 * We consider nonlinear vector function F(x,s) with
 * x unknown vector and s homotopy parameter. Where we
 * know solution  or can solve F(x,0)=0 by Newton solver and
 * F(x,1)=0 is final system we want to solve.
 *
 * 1) PRoblem of homotopy method is determination of step of parameter s
 * 2) Consider system  (F(x,s) , g(1-s) ) = 0 and solve by Newton method with LS
 *    where function g is monotone differentiable with nonzero derivative at 0.
 *
 *    at one iteration we solve:
 *
 *    (dF/dx(x^n)) \delta = -F(x^n) + dF/ds(x^n) * g(1-s^n)/g'(1-s^n) = 0
 *
 *   we find line search direction and then perform usual linesearch to update x
 *
 *   The function g should be large enough according to the remaining residual F(x,s).
 *   ...
 * 3) Split every step into
 *     homotopy step: (solving dF/dx(x^n) \delta = dF/ds(x^n) and line search for longest step upto given maximum tolerance)
 *     and Newton step (solving  dF/dx(x^n( \delta = -F(x^n) and linesearch for best decrease
 */

#include "base/parameter_handler.h"
#include <lac/sparse_ilu.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/solver_bicgstab.h>
#include <lac/precondition.h>
#include <limits>

using namespace dealii;

class NonlinSystemBase {
public:
    typedef Vector<double> VecType;
    typedef SparseMatrix<double> MatType;

    virtual void compute_jacobian(const VecType &x, double s, MatType &jac, bool symmetric) =0;
    virtual void compute_function(const VecType &x, double s, VecType &func) =0;
    virtual void compute_parameter_derivative(const VecType &x, double s, VecType &diff)=0;
    virtual VecType &get_left_scaling()=0; // possibly scale system from left to keep it symmetric
    virtual VecType &get_solution_vector()=0;
    virtual MatType &get_matrix()=0;
    virtual VecType &get_function_vector()=0;
    virtual ~NonlinSystemBase() {};
private:
};

class HomotopyNewton {
public:
    enum ConvergenceState { iterating =0, converged=1,
                            diverged_max_it=-1};
    typedef Vector<double> VecType;
    typedef SparseMatrix<double> MatType;

    class Params {
    public:
        double tol_final,    ///< Final tolerance at sucessful convergence.
               tol_max,      ///< Maximal error allowed after homotopy step.
               tol_step;     ///< Tolerance necessary to perform next homotopy step
        unsigned int max_it;

        double alpha; // decrease factor for Newton method
        double max_lambda; // maximum step in homotopy parameter
        double min_lambda; // minimum lambda in Newton step relative to relative step length
        bool check_jacobian_flag;

        // linear solver parameters
        double lin_rtol;
        double lin_atol;
        unsigned int lin_max_it;

        bool use_homotopy;
    };

    HomotopyNewton( ParameterHandler &prm);
    void set_system(NonlinSystemBase &sys);

    ConvergenceState solve();
    unsigned int get_iter() {return iter;}
    unsigned int get_cum_iter() {return cum_iter;}

private:
    void newton_step();
    void homotopy_step();
    void line_search();
    void linear_solve(VecType &sol, VecType &rhs);
    void check_jacobian();


    Params params;

    // work variables
    double s_param, lambda;
    double f_norm_save, f_norm;

    VecType *solution, *func;
    MatType *jac;
    NonlinSystemBase *system;

    ConvergenceState conv_state;

    // aux vectors
    VecType s_diff, solution_save, sol_step;

    SparseILU<double> precondition;

    unsigned int cum_lin_iter, iter, cum_iter, sym_iter;

    bool converged_homotopy;
    bool reuse_precondition_sparsity;
    bool symmetric_jac;



};


HomotopyNewton::HomotopyNewton(ParameterHandler &prm) {

    params.use_homotopy = prm.get_bool("nlin_use_homotopy");
    params.tol_final = prm.get_double("nlin_tol");
    params.tol_step = params.tol_final * prm.get_double("nlin_tol_step_mult");
    params.tol_max= params.tol_step * prm.get_double("nlin_tol_max_mult");
    Assert(params.tol_final < params.tol_step, ExcMessage("Tol_final should be smaller then tol_step."));
    Assert(params.tol_step < params.tol_max, ExcMessage("Tol_step should be smaller then tol_max."));
    params.max_it = prm.get_double("nlin_max_it");
    params.alpha = prm.get_double("nlin_alpha");
    params.max_lambda = prm.get_double("nlin_max_lambda");
    params.min_lambda = prm.get_double("nlin_min_lambda");
    params.check_jacobian_flag = prm.get_bool("nlin_check_jacobian");

    params.lin_rtol=prm.get_double("lin_rtol");
    params.lin_atol=prm.get_double("lin_atol");
    params.lin_max_it=prm.get_integer("lin_max_it");


    cum_lin_iter=0;
    cum_iter=0;

    reuse_precondition_sparsity=false;

}

void HomotopyNewton::set_system(NonlinSystemBase &sys)
{
    system = &sys;

    solution=&system->get_solution_vector();
    func=&system->get_function_vector();
    Assert( solution->size() == func->size(), ExcMessage("Solution and func size do not match.") );

    jac=&system->get_matrix();

    reuse_precondition_sparsity=false;

    // reinit auxiliary vectors
    s_diff.reinit(*solution);
    solution_save.reinit(*solution);
    sol_step.reinit(*solution);

    cum_lin_iter=0;
    cum_iter=0;
    symmetric_jac = false;


}


HomotopyNewton::ConvergenceState HomotopyNewton::solve() {

    s_param=0.0;
    iter=0;
    converged_homotopy=!params.use_homotopy;
    sym_iter=0;

    system->compute_function(*solution, s_param, *func);
    f_norm=func->l2_norm();

    double f_norm_save;

    while ( !converged_homotopy  || f_norm > params.tol_final ) {

        f_norm_save = f_norm;
        if (!converged_homotopy && f_norm < params.tol_step) {
            homotopy_step();
            if ( fabs(1.0 - s_param) < params.tol_final ) {
                converged_homotopy=true;
                s_param=1.0;
            }
            sym_iter=0;
            cout << "  "<<setw(10) << iter<< setw(15) << f_norm << " H" << setw(8) << lambda << " d: " << f_norm/f_norm_save << endl;
        }
        else {
            newton_step();
            cout << "  "<<setw(10) << iter<< setw(15) << f_norm << " N" << setw(8) << lambda << " d: " << f_norm/f_norm_save << endl;
        }
        iter++;
        cum_iter++;
        if (iter > params.max_it) {
            conv_state = diverged_max_it;
            return (conv_state);
        }
    }

    conv_state = converged;
    return (conv_state);
}

void HomotopyNewton::homotopy_step() {

    symmetric_jac = true;
    system->compute_jacobian(*solution, s_param, *jac, symmetric_jac);
    system->compute_parameter_derivative(*solution, s_param, s_diff);
    solution_save=*solution;
    f_norm_save=f_norm;


    linear_solve(sol_step, s_diff);
    // TODO: check sucesssful linear solve

    // linesearch to get maximum step but keep f_norm < tol_max


    lambda = 1.0 - s_param;
    lambda = min(lambda, params.max_lambda);
    while (true) {
     // compute F(sol - lambda * sol_step)
     solution->equ(1.0, solution_save, -lambda, sol_step);
     system->compute_function(*solution,s_param+lambda,*func);
     // TODO: check for solution in function domain

     f_norm=func->l2_norm();

     if (f_norm < params.tol_max ) { /* Sufficient reduction */
         goto theend;
     }
     lambda*=0.5;
    }

theend:
    s_param+=lambda;
    return;
}

void HomotopyNewton::newton_step() {
    //x_norm=solution->l2_norm();
    // TODO: check for NaN in function value

    // use non-symmetric jacobian only for for converged homotopy

    if (sym_iter<2) {
        sym_iter++;
        symmetric_jac = true;
    } else {
        symmetric_jac=false;
    }
    system->compute_jacobian(*solution, s_param, *jac, symmetric_jac);
    if (converged_homotopy && params.check_jacobian_flag) check_jacobian();

    linear_solve(sol_step, *func);
    //cout << "step: " << sol_step << endl;
    // TODO: check sucesssful linear solve

    line_search();
}

void HomotopyNewton::line_search() {
    // linesearch taken from PETSC
    /*
        Note that for line search purposes we work with with the related
        minimization problem:
           min  z(x):  R^n -> R,
        where z(x) = .5 * fnorm*fnorm, and fnorm = || f ||_2.
      */
    solution_save = *solution;
    f_norm_save = f_norm;


    //sol_new_norm = sol_new.l2_norm();
    // TODO: check for zero solution
    // TODO: check for max step size, scale back to max size

    // compute max of relative step
    // get min_lambda as min_lambda over relative l-inf norm of the step
    // so the min_lambda parameter is minimum relative l-inf norm of the final step
    double ratio_max=0.0;
    for(unsigned int i=0; i< solution->size();i++) ratio_max=max(sol_step(i)/(*solution)(i),ratio_max);
    double min_lambda = params.min_lambda / ratio_max;

    double init_slope;
    if (symmetric_jac) {
        const double eps = numeric_limits<double> :: epsilon() * 1024;
        // for Picard iterations use numerical derivative
        solution->equ(1.0, solution_save, -eps, sol_step);
        system->compute_function(*solution,s_param,*func);
        f_norm=func->l2_norm();
        init_slope=0.5*(f_norm *f_norm - f_norm_save*f_norm_save)/eps;
    } else {
        // -f(x) * M * sol_step = d/d_lambda (0.5 * f(x-ls) * f(x-ls)) at lambda=0
        init_slope=-jac->matrix_scalar_product(*func,sol_step);
    }
    if (init_slope >= 0.0) {
        cout << "positive slope: " << -init_slope << endl;
        // TODO: vetsinou se jedna o derivaci okolo nuly co s tim ?
        // ? pocitat numericky ?
    }
    if (init_slope > 0.0)  init_slope = -init_slope; // ??
    if (init_slope == 0.0) init_slope = -1.0;

    lambda     = 1.0;
    // compute F(sol - lambda * sol_step)
    solution->equ(1.0, solution_save, -lambda, sol_step);
    system->compute_function(*solution,s_param,*func);
    f_norm=func->l2_norm();

    if (.5*(f_norm)*(f_norm) <= .5*f_norm_save*f_norm_save + params.alpha*init_slope) { /* Sufficient reduction */
        return;
    }


      /* Fit points with quadratic */
      // ax^2 +bx +c; c = 0.5*f_norm_save ^2; b = init_slope; c = 0.5 f_norm^2 - b -c
      // min at: -b / 2a

      double lambdatemp = -init_slope/((f_norm)*(f_norm) - f_norm_save*f_norm_save - 2.0*init_slope);
      double lambdaprev = lambda;
      double gnormprev  = f_norm;
      if (lambdatemp > .5*lambda)  lambdatemp = .5*lambda;
      if (lambdatemp <= .1*lambda) lambda = .1*lambda;
      else                         lambda = lambdatemp;

      solution->equ(1.0, solution_save, -lambda, sol_step);
      system->compute_function(*solution,s_param,*func);
      f_norm=func->l2_norm();
      if (.5*(f_norm)*(f_norm) < .5*f_norm_save*f_norm_save + lambda*params.alpha * init_slope) { /* sufficient reduction */
          return;
      }

      /* Fit points with cubic */
      int count = 1;
      while (true) {
        if (lambda <= min_lambda) {
            cout << "LS failed!" << endl;
            return;
        }

        // ax^3 + bx^2 +cx +d; d = 05.f_norm_save^2; c = init_slope;
        // min at:
        // a l_prev +b = (0.5 gnormprev^2 - d - c l_prev)/ l_prev^2
        // a l + b = (0.5 f_nrom^2 -d - c l)/ l^2
        double t1 = .5*((f_norm)*(f_norm) - f_norm_save*f_norm_save) - lambda*init_slope;
        double t2 = .5*(gnormprev*gnormprev  - f_norm_save*f_norm_save) - lambdaprev*init_slope;
        double a  = (t1/(lambda*lambda) - t2/(lambdaprev*lambdaprev))/(lambda-lambdaprev);
        double b  = (-lambdaprev*t1/(lambda*lambda) + lambda*t2/(lambdaprev*lambdaprev))/(lambda-lambdaprev);

        // 2*determinant of optimality condition eq.
        double d  = b*b - 3*a*init_slope;
        if (d < 0.0) d = 0.0; //
        if (a == 0.0) {
          // quadratic case
          lambdatemp = -init_slope/(2.0*b);
        } else {
          // take positive root
          lambdatemp = (-b + sqrt(d))/(3.0*a);
        }
        lambdaprev = lambda;
        gnormprev  = f_norm;
        if (lambdatemp > .5*lambda)  lambdatemp = .5*lambda;
        if (lambdatemp <= .1*lambda) lambda     = .1*lambda;
        else                         lambda     = lambdatemp;

        solution->equ(1.0, solution_save, -lambda, sol_step);
        system->compute_function(*solution,s_param,*func);
        f_norm=func->l2_norm();
        if (.5*(f_norm)*(f_norm) < .5*f_norm_save*f_norm_save + lambda*params.alpha * init_slope) { /* sufficient reduction */
            return;
        }

        count++;
      }
}

void HomotopyNewton::linear_solve(VecType &sol, VecType &rhs)
{
  ReductionControl        solver_control (params.lin_max_it, params.lin_atol, params.lin_rtol);
  solver_control.log_history(true);

  if (true || ! symmetric_jac) {
      // precise jacobian can be non symetric
      SolverBicgstab<>        bcgs_solver (solver_control);

      // initalize preconditioner
      precondition.initialize(*jac,
          SparseILU<double>::AdditionalData(
                  0, //strengthen_diagonal
                  0, //  extra_off_diagonals
                  reuse_precondition_sparsity, //use_previous_sparsity
                  0 //    const SparsityPattern *    use_this_sparsity = 0
          ));
      reuse_precondition_sparsity=true;

      cout << "   " << "bcgs" << endl;
      bcgs_solver.solve (*jac, sol, rhs,
            precondition);

  } else {
      // assume symmetric jacobian
      SolverCG<>              cg_solver (solver_control);

      // scale system from left and solve
      VecType &scale_vec = system->get_left_scaling();
      for(unsigned int row=0; row < rhs.size(); ++row) {
          rhs(row)/=scale_vec(row);
          for(MatType::iterator col_it = jac->begin(row);
                  col_it != jac->end(row); ++col_it)
              col_it->value()/= scale_vec(row);
      }
      // initalize preconditioner
      precondition.initialize(*jac,
          SparseILU<double>::AdditionalData(
                  0, //strengthen_diagonal
                  0, //  extra_off_diagonals
                  reuse_precondition_sparsity, //use_previous_sparsity
                  0 //    const SparsityPattern *    use_this_sparsity = 0
          ));
      reuse_precondition_sparsity=true;

      cout << "   " << "cg" << endl;
      cg_solver.solve (*jac, sol, rhs,
            precondition);

      // scale back from left
      for(unsigned int row=0; row < rhs.size(); ++row) {
          rhs(row)*=scale_vec(row);
          for(MatType::iterator col_it = jac->begin(row);
                  col_it != jac->end(row); ++col_it)
              col_it->value()*= scale_vec(row);
      }
  }

  cum_lin_iter += solver_control.last_step();
/*std::cout << "CG: " << solver_control.last_step() << ", res: "
        <<solver_control.last_value() << " cum lin it: " << cum_lin_iter
        << std::endl;*/

  //solution.print(cout);
}

void HomotopyNewton::check_jacobian() {

    //MatType num_jac(jac->get_sparsity_pattern());
    solution_save=*solution;
    VecType &func_save = sol_step;
    func_save = *func;
    const double dx = 1.0E-10; // * f_norm / sqrt(solution->size());
    double diff, err=0;

    for(unsigned int col = 0; col < solution->size(); col++ ) {
        (*solution)(col)+=dx;
        system->compute_function(*solution, s_param, *func);
        for(unsigned int row= 0; row < solution->size(); row++) {
            diff = ((*func)(row) - func_save(row) ) / dx;
            if (fabs(diff) > 1.0E-20) {
                //if ( fabs((*jac)(row,col) - diff) > 1.0E-5*fabs(diff) ) {
                    //cout << "(" << row << ", " << col << ") " << (*jac)(row,col) <<" "<< diff << " " << (*jac)(row,col) - diff << endl;
                //}

                err =  max(err, fabs((*jac)(row,col) - diff) / max( fabs(diff),  fabs( (*jac)(row,col) ) ) );
            }
        }
        (*solution)(col)-=dx;
    }
    cout << "jac inf elem rel error: " << sqrt(err)/jac->frobenius_norm() << endl;
    *solution = solution_save;
}

#endif /* NONLIN_SOLVER_HH_ */
