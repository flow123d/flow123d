/* 
 * File:   output.hh
 * Author: jb
 *
 * Created on June 12, 2010, 7:44 AM
 */

/**
 * you can plot 4th column of bc_output by gnuplot command:
 *
 * gnuplot>  plot "./bc_output" using ($1):($4)
 * 
 */

#ifndef _OUTPUT_HH
#define	_OUTPUT_HH

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <numerics/vectors.h>
#include <fe/fe_values.h>

#include "time.hh"

using namespace std;
using namespace dealii;

template <int dim>
class BCOutput {
public:
    BCOutput( char * filename, unsigned int order)
      :  bc_output(filename), order(order)
      {
      // TODO: this is indexed by boundary indicators
      // should by indexed by application index of the boundary
      // when we implement some boundary descriptor
      bc_flux.reinit(4);
      bc_head.reinit(4);
      bc_surface.reinit(4);
      init_volume=-1;
      cum_bc_flux.reinit(4,0.0);
          bc_output<< "* time flux[i] head[i] ... volume error" << endl;

      }

    template <class DH, class  sat_DH>
    void output(SolverTime &time, DH &dh, sat_DH &sat_dh, BlockVector<double> &solution, Vector<double> &saturation) {

    QGauss<dim> quadrature_formula(order + 1);
    QGauss < dim - 1 > face_quadrature_formula(order + 1);

    FEValues<dim> sat_fe_values(sat_dh.get_fe(), quadrature_formula,
            update_values | update_JxW_values);
    FEFaceValues<dim> fe_face_values(dh.get_fe(), face_quadrature_formula,
            update_values | update_normal_vectors |
            update_quadrature_points | update_JxW_values);


    //const unsigned int dofs_per_cell = dh.get_fe().dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

//    std::vector<unsigned int> local_dof_indices(dofs_per_cell);
//    std::vector<unsigned int> face_dofs(dh.get_fe().dofs_per_face);

    std::vector<double> saturation_values(n_q_points);
    std::vector<Tensor <1,dim> > face_u_values(n_face_q_points);
    std::vector<double> face_p_values(n_face_q_points);
    std::vector<unsigned int> face_dofs(dh.get_fe().dofs_per_face);

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);

    bc_flux=0;
    bc_head=0;
    bc_surface=0;
    volume=0;

    // Main cycle over the cells.
    typename DoFHandler<dim>::active_cell_iterator
        cell = dh.begin_active(),
            endc = dh.end();
    typename DoFHandler<dim>::active_cell_iterator
        sat_cell = sat_dh.begin_active();
    for (; cell != endc; ++cell, ++sat_cell) {
        // volume integrals
        sat_fe_values.reinit(sat_cell);
        sat_fe_values.get_function_values(saturation, saturation_values);
        // cycle over quadrature points
        for (unsigned int q = 0; q < n_q_points; ++q) 
            volume += saturation_values[q] * sat_fe_values.JxW(q);

        if (cell->at_boundary()) {
            for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
                unsigned int b_ind = cell->face(face_no)->boundary_indicator();
                if (b_ind == 255) continue;
                fe_face_values.reinit(cell, face_no);

                // compute boundary integrals
                fe_face_values[velocities].get_function_values(solution, face_u_values);
                fe_face_values[pressure].get_function_values(solution, face_p_values);
                for (unsigned int q = 0; q < n_face_q_points; ++q) {
                    bc_head(b_ind) += face_p_values[q] * fe_face_values.JxW(q);
                    bc_flux(b_ind) += (face_u_values[q] * fe_face_values.normal_vector(q)) * fe_face_values.JxW(q);
                    bc_surface(b_ind) += fe_face_values.JxW(q);
                }
            }
        }
    }

     for(unsigned int i=0;i<4;i++) {
        if (bc_surface(i) != 0) {
            //bc_flux(i)/=bc_surface(i);
            bc_head(i)/=bc_surface(i);
            cum_bc_flux(i)+=bc_flux(i) * time.dt();
        }
     }
    // check of stability, negative velocity, min and max of pressure
    
    double u_min=*(min_element(solution.block(0).begin(),solution.block(0).end()));
    double u_max=*(max_element(solution.block(0).begin(),solution.block(0).end()));
    double h_min=*(min_element(solution.block(1).begin(),solution.block(1).end()));
    double h_max=*(max_element(solution.block(1).begin(),solution.block(1).end()));
/*
    for(int i=0;i<solution.block(0).size();++i) {
        if ((solution.block(0))(i)<u_min) u_min=solution.block(0).(i);
        if (solution.block(0)(i)>u_max) u_max=solution.block(0).(i);
    }
    double h_min=1000;
    double h_max=-1000;
    for(int i=0;i<solution.block(0).size();++i) {
        if (solution.block(1).(i)<h_min) h_min=solution.block(1).(i);
        if (solution.block(1).(i)>h_max) h_max=solution.block(1).(i);
    }
*/
    if (init_volume < 0.0) init_volume = volume;

    double error = 0;
    for (unsigned int i=0; i < cum_bc_flux.size();i++) {
        error+=cum_bc_flux(i);
        cout << i << " :BC " << cum_bc_flux(i) <<endl;
    }
    error=( (volume - init_volume)  + error  ) / volume;

      bc_output << setw(10) << time.t();
      bc_output << setw(13) << bc_flux(2);
      bc_output << setw(13) << bc_head(2);
      bc_output << setw(13) << bc_flux(3);
      bc_output << setw(13) << bc_head(3)
                << setw(13) << u_min
              << setw(13)<< u_max
              << setw(13)<< h_min
              << setw(13)<< h_max
                << setw(13) << volume
                << setw(13) << error << endl;


    }
private:
    std::ofstream bc_output;
    int order;

    // bc values
    Vector<double> bc_flux;
    Vector<double> bc_head;
    Vector<double> bc_surface;
    Vector<double> cum_bc_flux;

    double init_volume;
    double volume;

};

#endif	/* _OUTPUT_HH */

