/*
 * output.hh
 *
 *  Created on: Oct 26, 2011
 *      Author: jb
 */

#ifndef OUTPUT_HH_
#define OUTPUT_HH_


#include <fstream>
#include <iostream>

//#include <base/quadrature_lib.h>
//#include <base/logstream.h>
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
#include <lac/block_vector.h>
#include "base/parameter_handler.h"

#include "spatial_functions.hh"

//#include <base/tensor_function.h>

// #include <numerics/error_estimator.h>
// #include <grid/grid_refinement.h>
// #include <numerics/solution_transfer.h>

//#include <petsc.h>

#include "local_assembly.hh"

using namespace dealii;




/**
 * Output class.
 *
 * We need to reconstruct some output data per every cell idealy this should be done
 * by introducing our own child of DataOut_DoFData, and implement our own build_patches function.
 * But in order not to stuck in technicalities we just use another DoFHandler with all necessary fields.
 *
 * We may also use merge_patches method of DataOut_DoFData.
 */


template <int dim>
class FieldOutput {
public:
    FieldOutput(Triangulation<dim> &tria, unsigned int order, LocalAssembly<dim> &la);
    void reinit(ParameterHandler &prm);
    void output_fields(DoFHandler<dim> &solution_dh, Vector<double> &solution_vector, double time, bool force);
    void update_fields(DoFHandler<dim> &solution_dh, double time);
    ~FieldOutput();

private:

    unsigned int order;

    //FE_DGRaviartThomas<dim, dim> velocity_fe;
    //FE_DGQ<dim> pressure_trace_fe;
    //FE_DGQ<dim> pressure_fe;
    //FESystem<dim> fe;
    DoFHandler<dim> dh;
    BlockVector<double> out_vec;
    DataOut<dim> data_out;

    std::string file_name;
    unsigned int print_level;
    double print_time;
    double print_time_step;
    std::vector<unsigned int> blocks;
    LocalAssembly<dim> & local_assembly;

    //static const unsigned int n_output_components = dim + 4;
    double x_size;

    double bc_flux_total, last_volume, volume, init_volume, cum_bc_flux, head_var, q_var,norm_head_second_diff;
    std::ofstream bc_output;


};

template <int dim>
FieldOutput<dim>::FieldOutput(Triangulation<dim> &tria, unsigned int p_order, LocalAssembly<dim> &la)
: order(p_order),
  //velocity_fe(order),
  //pressure_trace_fe(order+2),
  //pressure_fe(order),
  //fe (velocity_fe,1,pressure_fe,4,pressure_trace_fe,1),
  dh(tria),

  print_level(0),
  print_time(0.0),
  print_time_step(0.1),
  local_assembly(la),
  bc_output("bc_output.out")
{
    cout << "fo construct " <<endl;


}

template <int dim>
void FieldOutput<dim>::reinit(ParameterHandler &prm)
{

    print_time_step=prm.get_double("print_time_step");
    x_size=prm.get_double("x_size");


    FESystem<dim> &fe = local_assembly.output_fe;
    dh.distribute_dofs(fe);
    DoFRenumbering::component_wise (dh);
    cout << "b sizE:" << fe.n_blocks() << endl;
    blocks.resize(fe.n_blocks());

    DoFTools::count_dofs_per_block (dh, blocks);

    out_vec.reinit(blocks.size());
    for (unsigned int i=0; i < blocks.size(); i++) {
        std::cout << blocks[i] << std::endl;
          out_vec.block(i).reinit(blocks[i]); // fast reinit, since
    }
    out_vec.collect_sizes();

    data_out.attach_dof_handler (dh);

    std::vector<std::string> names;

    names.push_back("flux");
    names.push_back("flux");

    names.push_back("half_flux");
    names.push_back("half_flux");

    names.push_back("pressure_head");
    names.push_back("saturation");
    names.push_back("p_error");
    names.push_back("q_error");

    names.push_back("post_phead");
    names.push_back("post_old_phead");
    names.push_back("residual");
    names.push_back("post_lambda");
    names.push_back("post_aux");

    std::vector<DataComponentInterpretation::DataComponentInterpretation> component_interpretation(2*dim,
            DataComponentInterpretation::component_is_part_of_vector);
    component_interpretation .push_back(DataComponentInterpretation::component_is_scalar);
    component_interpretation .push_back(DataComponentInterpretation::component_is_scalar);
    component_interpretation .push_back(DataComponentInterpretation::component_is_scalar);
    component_interpretation .push_back(DataComponentInterpretation::component_is_scalar);
    component_interpretation .push_back(DataComponentInterpretation::component_is_scalar);
    component_interpretation .push_back(DataComponentInterpretation::component_is_scalar);
    component_interpretation .push_back(DataComponentInterpretation::component_is_scalar);
    component_interpretation .push_back(DataComponentInterpretation::component_is_scalar);
    component_interpretation .push_back(DataComponentInterpretation::component_is_scalar);

    data_out.add_data_vector(out_vec, names, DataOut<dim>::type_automatic, component_interpretation);

    init_volume=-1;
    cum_bc_flux=0;
}

template <int dim>
FieldOutput<dim>::~FieldOutput()
{

}




template <int dim>
void FieldOutput<dim>::output_fields(DoFHandler<dim> &solution_dh, Vector<double> &solution_vector, double time, bool force)
{

  // update vector
  update_fields(solution_dh, time);

  // bc output
  if (cum_bc_flux == 0) {
      bc_output << " time   err%    loc_err%    bc_flux     volume      vol_diff" <<endl;
  }
  double error_local = (volume-last_volume) + bc_flux_total;
  cum_bc_flux +=bc_flux_total;
  double error = (volume - init_volume) + cum_bc_flux;
  bc_output << setw(12) << time
          << setw(12) << error/volume * 100
          << setw(12) << error_local/(volume - last_volume) * 100
          << setw(12) << bc_flux_total
          << setw(12) << volume
          << setw(12) << volume - last_volume<<endl;


  if (!force && time * 1.0000001 < print_time) return;
  std::cout << "PRINT time (" << print_level << "): " << time << "    variations: "<< head_var << " " << q_var << " " << norm_head_second_diff << std::endl;
  print_time += print_time_step;

  // file name
  std::ostringstream fns;
  fns << "./output/solution-" << std::setfill('0')<<std::setw(3) << print_level << ".vtk";
  file_name = fns.str();
  print_level++;

  //AssertDimension(out_vec.block(p_traces_bl).size(),  solution_vector.size());
  //out_vec.block(p_traces_bl) = solution_vector;

  data_out.build_patches (order+3);
  std::ofstream output (file_name.c_str());
  data_out.write_vtk (output);

}

template <int dim>
void FieldOutput<dim>::update_fields(DoFHandler<dim> &solution_dh, double time)
{
    typename DoFHandler<dim>::active_cell_iterator
        cell = dh.begin_active(),
        endc = dh.end(),
        sol_cell = solution_dh.begin_active();
    Vector<double> local_output(local_assembly.output_fe.dofs_per_cell);
    std::vector<unsigned int> local_dof_indices(local_assembly.output_fe.dofs_per_cell);

    local_assembly.solution->compute_head_second_diff();
    norm_head_second_diff = local_assembly.solution->head_second_diff.l2_norm();

    double l2_error=0;
    double l2_flux_error = 0;
    out_vec=0;
    bc_flux_total =0;
    last_volume=volume;
    volume =0;
    q_var=head_var =0;
    for (; cell != endc; ++cell, ++sol_cell) {
        local_assembly.reinit(sol_cell);
        local_assembly.output_evaluate(bc_flux_total, volume);
        local_assembly.compute_add_variation(q_var, head_var);


        l2_error += local_assembly.get_p_error() /x_size;
        l2_flux_error += local_assembly.get_q_error() /x_size;


        cell->get_dof_indices(local_dof_indices);
        out_vec.add(local_dof_indices, local_assembly.get_output_vector());





    }
    if (init_volume < 0) init_volume = volume;
    if ( local_assembly.solution->richards_data->has_exact_solution() ) {
        cout << "Time: " << time << " L2 error: " << sqrt(l2_error) << " " << sqrt(l2_flux_error) << endl;
    }
}

#endif /* OUTPUT_HH_ */
