/*
 * richards_lmh.hh
 *
 *  Created on: Oct 9, 2011
 *      Author: jb
 */

#ifndef RICHARDS_LMH_HH_
#define RICHARDS_LMH_HH_






#include <fstream>
#include <iostream>


#include <base/logstream.h>
//#include <base/function.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>

#include <dofs/dof_renumbering.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_dg_vector.h>
#include <base/polynomials_raviart_thomas.h>
#include <fe/fe_face.h>
#include <numerics/data_out.h>

#include <base/tensor_function.h>
#include <lac/sparse_ilu.h>
#include <lac/vector.h>
//#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

 #include <numerics/error_estimator.h>
 #include <grid/grid_refinement.h>
 #include <numerics/solution_transfer.h>

//#include <petsc.h>

using namespace dealii;


//#include <boost/assign.hpp>

//#include <schur_solver.hh> // eliminace rychlosti, zavest jako metodu rozhrani k solveru
//#include "direct_solver.hh"



#include <UMFPACK_system.hh>
#include <time.hh>
#include <spatial_functions.hh>
#include <hydro_functions.hh>
#include <output.hh>
#include <local_assembly.hh>
#include <nonlin_solver.hh>



/**
 *  Richards' equation solver.
 */



/***
 * rovnice Richards
 * zamyslene vstupy:
 * - hruba sit, pri konstrukci
 * - pocatecni podminka, set .. (i behem vypoctu)
 * - okrajove podminky (std:map pro okrajove indexy hrube site)
 * - hydrologicke funkce  (std:map pro materialove indexy hrube site)
 * - pocatecni a cilovy cas
 *  nepovine:
 * - scaling faktory, anisotropie (Function ..)
 *
 * metody pro pouziti:
 *  1) run - simulace na danem casovem intervalu
 *  2) get_solution, set_solution (pouze tlaky)
 *  3) assembly, .. pouziti nad danou casti globalni matice, v ramci externiho newtonovskeho resice
 *
 */



template <int dim>
class Richards_LMH : public NonlinSystemBase {
public:
    Richards_LMH(Triangulation<dim> &coarse_tria, ParameterHandler &prm, unsigned int order);
    //! Add boundary condition bc for the given boundary_index.
    //! boundary without assigned BC are homogeneous Dirichlet, but rather leads to exception


    void reinit(struct RichardsData<dim> *data);

    virtual void compute_jacobian(const VecType &x, double s, MatType &jac, bool symmetric);
    virtual void compute_function(const VecType &x, double s, VecType &func);
    virtual void compute_parameter_derivative(const VecType &x, double s, VecType &diff);
    virtual VecType &get_solution_vector()
            { return solution->phead;}
    virtual MatType &get_matrix()
            {return matrix;}
    virtual VecType &get_function_vector()
            {return solution->residual;}


//    void solve();

    void run();

    virtual ~Richards_LMH();

private:
    //enum block_names {traces=0, saturation=1};

    void assemble_system();
    void grid_refine();
    //void compute_errors () const;

    ParameterHandler &prm;

    //! numerics
    SolverTime time;

    int order; // order of FE

    Triangulation<dim> triangulation;

    FE_FaceQ<dim> trace_fe;
    FESystem<dim> fe;


    DoFHandler<dim> dof_handler;

    /// linear algebra objects
//    Vector<double> residual_source;
//    Vector<double> rhs;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> matrix;

    ConstraintMatrix     constraints;
    std::map<unsigned int, double> boundary_values;

    //! field vectors
//    Vector<double> saturation;
//    Vector<double> old_saturation;
//    Vector<double> m_mat_crit;
//    Vector<double> theta_diff;

//    BCOutput<dim>    bc_out;

    HomotopyNewton  *nlin_solver;
    RichardsData<dim> *richards_data;
    Solution<dim>    *solution;

    LocalAssembly<dim> local_assembly;
    FieldOutput<dim> field_output;
};

/**
 *  Construct the solver - discretization space and DOF handler
 */

template <int dim>
Richards_LMH<dim>::Richards_LMH (Triangulation<dim> &coarse_tria,ParameterHandler &prm, unsigned int order)
        :
                //bc(MAX_NUM_OF_DEALII_BOUNDARIES),
                prm(prm),
                time(prm),                  // start, end, dt_init
                order(order),

                triangulation(),
                trace_fe(order),
                fe(trace_fe,2),

                dof_handler (triangulation),             // dof_handler is initialized by empty tria.

                local_assembly(order),
                field_output(triangulation,order, local_assembly)
{




    // setup triangulation and dof handler
    cout << "tria: " << endl;
    triangulation.copy_triangulation(coarse_tria);

    cout<<"dofs_per_cell: "<<fe.dofs_per_cell<<endl;
    cout<<"dofs_per_face: "<<fe.dofs_per_face<<endl;

    dof_handler.distribute_dofs(trace_fe);
    constraints.clear();
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    constraints.close();

    CompressedSparsityPattern csp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler,csp, constraints);
    constraints.condense (csp);
    csp.compress();

    sparsity_pattern.copy_from(csp);

    std::cout << "active cells "<< triangulation.n_active_cells()
        << " total cells: " << triangulation.n_cells()
            << " dofs: " << dof_handler.n_dofs() << std::endl;

    field_output.reinit(prm);
    nlin_solver = new HomotopyNewton(prm);


}


template <int dim>
Richards_LMH<dim>::~Richards_LMH() {

}

/**
 * Reinit by richards data.
 */
template <int dim>
void Richards_LMH<dim>::reinit(struct RichardsData<dim> *data) {
    // setup support classes
    richards_data = data;
    solution = new Solution<dim>(richards_data);
    solution->reinit(dof_handler, false);
    local_assembly.set_data(solution);

    matrix.reinit(sparsity_pattern);
//    residual_source.reinit(rhs);
}



/**
 * Refine grid between time steps, this is far to be optimal, but is simplest
 * since we need not any grid and solution saving.
 * refinement and coarsing shall be done after sucesfull solution of the timestep.
 */
template <int dim>
void Richards_LMH<dim>::grid_refine()
{
   Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
   std::vector<bool> component_mask(dim+2, true);
   component_mask[dim]=false;
   

   // compute refinement indicators
   //KellyErrorEstimator<dim>::estimate (dof_handler,
   //                                   QGauss<dim-1>(order+1),
   //                                   typename FunctionMap<dim>::type(),    // Neuman boundary function
   //                                   solution,
   //                                   estimated_error_per_cell,
    //                                  component_mask);
/*
   { // compute error estimate based on jumps in pressure
    // Main cycle over the cells.
    typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell, ++sat_cell) {
        fe_values.reinit(cell);
        sat_fe_values.reinit(sat_cell);
        local_matrix = 0;
        local_rhs = 0;

    }
    }*/
  /* for(unsigned int i=0; i< estimated_error_per_cell.size(); i++) cout << i << " : " << estimated_error_per_cell(i) <<endl;*/
   /*GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   0.1, 0.4 , 1000 );*/
   /*GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   0, 0.33, 1000 );
    */
   /*
    triangulation.prepare_coarsening_and_refinement();
    SolutionTransfer<dim, BlockVector<double> > solution_transfer(dof_handler);
    solution_transfer.prepare_for_coarsening_and_refinement (solution);

    SolutionTransfer<dim, Vector<double> > sat_transfer(sat_dh);
    sat_transfer.prepare_for_coarsening_and_refinement (saturation);

    triangulation.execute_coarsening_and_refinement ();
    dof_handler.distribute_dofs (fe);
    sat_dh.distribute_dofs (sat_fe);
        std::cout << "active cells "<< triangulation.n_active_cells()
        << " total cells: " << triangulation.n_cells()
            << " dofs: " << dof_handler.n_dofs() << std::endl;

    DoFRenumbering::component_wise (dof_handler);
    linear_system.reinit(dof_handler);
    old_solution.reinit (linear_system.get_rhs());
    old_saturation.reinit(old_solution.block(head_b));

    solution_transfer.interpolate (solution, old_solution);
    solution.reinit(linear_system.get_rhs());
    solution=old_solution;
    linear_system.get_constraints().distribute (old_solution);
    linear_system.get_constraints().distribute (solution);

    sat_transfer.interpolate(saturation, old_saturation);
    saturation.reinit(solution.block(head_b));
    // saturation is discontinuous -> no constrains,
    // also we don't copy saturation, since it will be imediatelly updated from solution
    residual.reinit(linear_system.get_rhs());
*/
}

template <int dim>
void Richards_LMH<dim>::compute_jacobian(const VecType &x, double s, MatType &jac, bool symmetric)
{
    std::vector<unsigned int> local_dof_indices(dof_handler.get_fe().dofs_per_cell);
    std::map<unsigned int, double> boundary_values;

    jac=0;
    solution->update();

    // Main cycle over the cells.
    typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
            endc = dof_handler.end();
    for (; cell != endc; ++cell) {
        local_assembly.reinit(cell);
        local_assembly.apply_bc(boundary_values);
/*
        cout << "------------\n";
        local_matrix.print_formatted(cout);
        */

        cell->get_dof_indices(local_dof_indices);
        matrix.add(local_dof_indices, local_assembly.get_matrix());
//        rhs.add(local_dof_indices, local_assembly.get_rhs());


    }

    // apply boundary values to solution vector
    for(map<unsigned int,double>::iterator iter = boundary_values.begin();
        iter != boundary_values.end();
        ++iter) {

        solution->phead(iter->first) = iter->second;
    }

//    MatrixTools::apply_boundary_values(boundary_values,
//            matrix,
//            solution->phead,
//            rhs);


}

template <int dim>
void Richards_LMH<dim>::compute_function(const VecType &x, double s, VecType &func) {
    std::vector<unsigned int> local_dof_indices(dof_handler.get_fe().dofs_per_cell);

    func=0;
    solution->update();

    // Main cycle over the cells.
    typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
            endc = dof_handler.end();
    for (; cell != endc; ++cell) {
        local_assembly.reinit(cell);
        local_assembly.apply_bc(boundary_values);

        cell->get_dof_indices(local_dof_indices);
        func.add(local_dof_indices, local_assembly.get_rhs());

    }
    // apply boundary values to solution vector
    for(map<unsigned int,double>::iterator iter = boundary_values.begin();
        iter != boundary_values.end();
        ++iter) {

        func(iter->first) = 0.0;
    }

    // debugging residual output
    //field_output.output_fields(dof_handler, solution->phead, time.t(), true);
}

template <int dim>
void Richards_LMH<dim>::compute_parameter_derivative(const VecType &x, double s, VecType &diff)
{

}


/**
 *   Computation method.
 *
 */


template <int dim>
void Richards_LMH<dim>::run ()
{

    Vector<double> sol_last;//, sat_last;
    Vector<double> sol_new;
    double res_norm;
    double last_res_norm;
    double decrease, last_decrease;

    richards_data->set_time(time.t());
    local_assembly.set_dt(time.dt(), time.t());

    // set initial condition
    std::vector<unsigned int> face_dofs(trace_fe.dofs_per_face);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();

    for (; cell != dof_handler.end(); ++cell) {

        for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
            cell->face(face_no)->get_dof_indices(face_dofs);

            solution->phead(face_dofs[0])=richards_data->initial->value(cell->face(face_no)->barycenter());

        }
    }

    solution->timestep_update(time.dt()); // update nonlinearities
    field_output.output_fields(dof_handler, solution->phead, time.t(),false);

    compute_jacobian(solution->phead, 0.0, matrix, true); // apply boundary conditions to initial condition in order to get true initial residum



    time.inc();
    nlin_solver->set_system(*this);


  do {  // ---------------------------- Time loop
      //grid_refine();
      solution->timestep_update(time.dt());
renew_timestep:
      std::cout << "Timestep " << time.n_step() << " t= "<<time.t()<<", dt= " << time.dt() << endl;

      local_assembly.set_dt(time.dt(), time.t());
      richards_data->set_time(time.t());
      // TODO: set nonlinear tolerance relative to time step
      //res_norm*= time.end_t()/time.dt();


      HomotopyNewton::ConvergenceState state = nlin_solver->solve();

      if (state < 0) {
          // divergence
          if (time.reinc_time(0.3)) {
             solution->revert();
             goto renew_timestep;
         } else {
             cout << "DIVERGENCE OF SOLVER" <<endl;
             return;
         }
      }


      /**** SOME EXPERIMENTS WITH DT ADAPTIVITY
                  {
                  // vypada to, ze posun cela by mohla lepe zachytit norma zmeny reseni - integral
                  Vector<double> diff(saturation);
                  Vector<double> diff_integral(triangulation.n_active_cells());

                  VectorTools::integrate_difference(sat_dh, residual.block(head_b), ZeroFunction<dim>(1),diff_integral ,QGauss<dim>(order+1),VectorTools::L2_norm);
                  double res_diff=diff_integral.l2_norm();

                  diff=saturation; diff-=sat_last;
                  VectorTools::integrate_difference(sat_dh, diff, ZeroFunction<dim>(1),diff_integral ,QGauss<dim>(order+1),VectorTools::L2_norm);
                  double sat_diff=diff_integral.l2_norm();

                  diff=solution.block(head_b); diff-=head_last;
                  VectorTools::integrate_difference(sat_dh, diff, ZeroFunction<dim>(1),diff_integral ,QGauss<dim>(order+1),VectorTools::L2_norm);
                  double head_diff=diff_integral.l2_norm();

                  cout<< "res norm "<<res_norm << " decrease: "<< decrease <<endl;
                  cout<< iter<<" (r,dic,dr,ds,dh): "
                      <<res_norm <<" "
                      <<decrease <<" "
                      <<res_diff <<" "
                      <<sat_diff <<" "
                      <<head_diff << " " << endl;

                  }
      */



      // time error estimator
/*
      if (theta_diff.size() == 0) {
        theta_diff.reinit(saturation.size(),false);
      }

      double diff;
      double err=0;
      for(int i=0;i<saturation.size();++i) {
          diff=(saturation(i)-old_saturation(i))/time.dt();
          //cout << diff <<endl;
          err=std::max(err, diff - theta_diff(i) );
          theta_diff(i) = diff;
      }
      err*=time.dt()/2.0;

      cout << "time err est: " << err << endl;
*/


      field_output.output_fields(dof_handler, solution->phead, time.t(),false);

      // adapt time step
      double factor=1.0;
      if (nlin_solver->get_iter()<5) factor=1.0 / 0.9;
      else if (nlin_solver->get_iter()>13) factor=0.9;
      time.scale_time_step(factor);

      std::cout << "nonlin iter: "<< nlin_solver->get_iter() << " cum: "<< nlin_solver->get_cum_iter() << endl;
    }
  while (time.inc());

}


#endif /* RICHARDS_LMH_HH_ */
