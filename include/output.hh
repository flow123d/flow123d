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


//#include <base/tensor_function.h>

// #include <numerics/error_estimator.h>
// #include <grid/grid_refinement.h>
// #include <numerics/solution_transfer.h>

//#include <petsc.h>

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
    FieldOutput(Triangulation<dim> &tria, unsigned int order = 0);
    void reinit(ParameterHandler &prm);
    void output_fields(DoFHandler<dim> &solution_dh, Vector<double> &solution_vector, double time);
    ~FieldOutput();

private:
    enum block_index_names {
        velocity_bl=0,
        pressure_bl=1,
        saturation_bl=2,
        estimator_bl=3,
        p_traces_bl=4
    };

    FE_DGVector<PolynomialsRaviartThomas<dim>, dim> velocity_fe;
    FE_FaceQ<dim> pressure_trace_fe;
    FE_DGQ<dim> pressure_fe;
    FESystem<dim> fe;
    DoFHandler<dim> dh;
    BlockVector<double> out_vec;
    DataOut<dim> data_out;

    unsigned int order;
    std::string file_name;
    double print_time;
    unsigned int print_level;
    double print_time_step;
    std::vector<unsigned int> blocks;

    static const unsigned int n_output_components = dim + 4;

    void set_parameters(ParameterHandler &prm);

};

template <int dim>
FieldOutput<dim>::FieldOutput(Triangulation<dim> &tria, unsigned int p_order)
: order(p_order),
  velocity_fe(order,mapping_piola),
  pressure_fe(order),
  pressure_trace_fe(order),
  fe (velocity_fe,1,pressure_fe,3,pressure_trace_fe,1),
  dh(tria),
  print_level(0),
  print_time(0.0),
  print_time_step(0.1)
{}

template <int dim>
void FieldOutput<dim>::reinit(ParameterHandler &prm)
{
    set_parameters(prm);

    dh.distribute_dofs(fe);
    DoFRenumbering::component_wise (dh);
    blocks.resize(5);

    DoFTools::count_dofs_per_block (dh, blocks);

    out_vec.reinit(blocks.size());
    for (unsigned int i=0; i < blocks.size(); i++) {
        std::cout << blocks[i] << std::endl;
          out_vec.block(i).reinit(blocks[i]); // fast reinit, since
    }
    out_vec.collect_sizes();

    data_out.attach_dof_handler (dh);

    std::vector<std::string> names(dim, "flux");
    names.push_back("pressure");
    names.push_back("saturation");
    names.push_back("estimator");
    names.push_back("pressure_traces");

    std::vector<DataComponentInterpretation::DataComponentInterpretation> component_interpretation(dim,
            DataComponentInterpretation::component_is_part_of_vector);
    component_interpretation .push_back(DataComponentInterpretation::component_is_scalar);
    component_interpretation .push_back(DataComponentInterpretation::component_is_scalar);
    component_interpretation .push_back(DataComponentInterpretation::component_is_scalar);
    component_interpretation .push_back(DataComponentInterpretation::component_is_scalar);

    data_out.add_data_vector(out_vec, names, DataOut<dim>::type_automatic, component_interpretation);
}

template <int dim>
FieldOutput<dim>::~FieldOutput()
{

}


template <int dim>
void FieldOutput<dim>::set_parameters(ParameterHandler &prm)
{
    prm.declare_entry ("print_time_step", "0.1",
                        Patterns::Double(),
                        "Time step for filed output.");
    print_time_step=prm.get_double("print_time_step");

}


template <int dim>
void FieldOutput<dim>::output_fields(DoFHandler<dim> &solution_dh, Vector<double> &solution_vector, double time)
{

//  bc_out.output(time, dof_handler, sat_dh, solution, saturation);

  if (time < print_time) return;
  std::cout << "PRINT time (" << print_level << "): " << time << std::endl;
  print_time += print_time_step;

  // file name
  std::ostringstream fns;
  fns << "./output/solution-" << std::setfill('0')<<std::setw(3) << print_level << ".vtk";
  file_name = fns.str();
  print_level++;

  // update vector
  AssertDimension(out_vec.block(p_traces_bl).size(),  solution_vector.size());
  out_vec.block(p_traces_bl) = solution_vector;

  data_out.build_patches (order+1);
  std::ofstream output (file_name.c_str());
  data_out.write_vtk (output);

}
#endif /* OUTPUT_HH_ */
