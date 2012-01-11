/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id: quadrature.hh 1352 2011-09-23 14:14:47Z jan.stebel $
 * $Revision: 1352 $
 * $LastChangedBy: jan.stebel $
 * $LastChangedDate: 2011-09-23 16:14:47 +0200 (Fri, 23 Sep 2011) $
 *
 * @file
 * @brief Discontinuous Galerkin method for equation of transport with dispersion.
 *  @author Jan Stebel
 */

#include "transport_dg.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/mapping_p1.hh"
#include "solve.h"
#include "petscmat.h"
#include <armadillo>


TransportDG::TransportDG(TimeMarks & marks, Mesh & init_mesh, MaterialDatabase & material_database)
        : TransportBase(marks, init_mesh, material_database)
{
    // set up time governor
    time_=new TimeGovernor(
            0.0,
            OptGetDbl("Global", "Stop_time", "1.0"),
            *time_marks
            );

    time_->set_permanent_constrain(
            OptGetDbl("Global", "Time_step", "1.0"),
            OptGetDbl("Global", "Time_step", "1.0")
            );
    time_->fix_dt_until_mark();
    
    // set up solver
    solver = new Solver;
    solver_init(solver);

    /*
     * Read initial condition.
     *
     * Read boundary conditions.
     *
     * Distribute DOFs.
     */

    dof_handler2d = new DOFHandler<2>(mesh());
    fe2d = new FE_P_disc<1,2>;

    dof_handler2d->distribute_dofs(*fe2d);

    ls = new LinSys_MPIAIJ(dof_handler2d->n_global_dofs());

    // assemble system matrix
    ls->start_allocation();
    assemble();
    ls->start_add_assembly();
    assemble();
    ls->finalize();

}



void TransportDG::update_solution()
{
    time_->next_time();
    time_->view();

    // solve
    solve_system(solver, ls);

	// save solution
}



void TransportDG::get_solution_vector(double *& vector, unsigned int & size)
{
    vector = ls->get_solution_array();
    size   = ls->vec_lsize();
}



void TransportDG::get_parallel_solution_vector(Vec & vector)
{
    vector = ls->get_solution();
}



void TransportDG::set_velocity_field(Vec & velocity_vector)
{
	flux_vector = velocity_vector;
}



void TransportDG::output_data()
{
    double *solution = ls->get_solution_array();
    unsigned int ndofs = fe2d->n_dofs();
    unsigned int dof_indices[ndofs];
    ofstream f("transport_dg.dat");
    vec3 p;

    for (DOFHandler<2>::CellIterator cell = dof_handler2d->begin_cell(); cell != dof_handler2d->end_cell(); ++cell)
    {
        if (cell->dim != 2) continue;

        dof_handler2d->get_dof_indices(cell, dof_indices);

        for (int i=0; i<ndofs; i++)
        {
            p = cell->node[i]->point();
            f << p(0) << " "
              << p(1) << " "
              << p(2) << " "
              << solution[dof_indices[i]] << endl;
        }
        // save the first node again for nice drawing in Gnuplot
        p = cell->node[0]->point();
        f << p(0) << " "
          << p(1) << " "
          << p(2) << " "
          << solution[dof_indices[0]] << endl;
        f << endl << endl;
    }
    f.close();
}



void TransportDG::assemble()
{
    MappingP1<2,3> map;
    QGauss<2> q(2);
    FEValues<2,3> fe_values(map, q, *fe2d, update_values | update_gradients | update_JxW_values | update_quadrature_points);
    unsigned int ndofs = fe2d->n_dofs();
    unsigned int dof_indices[ndofs], side_dof_indices[2][ndofs];
    PetscScalar local_matrix[ndofs*ndofs], local_rhs[ndofs], local_flux_matrix[ndofs*ndofs];
    double gamma, omega[2], Kgradn[2][ndofs];
    double alpha = 1e3;

    QGauss<1> q1d(2);
    FESideValues<2,3> *fe_values_side[2];
    DOFHandler<2>::CellIterator cell = dof_handler2d->begin_cell();

    for (int sid=0; sid<2; sid++)
        fe_values_side[sid] = new FESideValues<2,3>(map, q1d, *fe2d, update_values | update_gradients | update_side_JxW_values | update_normal_vectors);



	// assemble integral over elements
    for (cell = dof_handler2d->begin_cell(); cell != dof_handler2d->end_cell(); ++cell)
    {
        if (cell->dim != 2) continue;

        fe_values.reinit(cell);

        dof_handler2d->get_dof_indices(cell, dof_indices);

        for (int i=0; i<ndofs; i++)
        {
            for (int j=0; j<ndofs; j++)
            {
                local_matrix[i*ndofs+j] = 0;
                for (int k=0; k<q.size(); k++)
                    local_matrix[i*ndofs+j] += arma::dot(fe_values.shape_grad(j,k),fe_values.shape_grad(i,k))
                                              *fe_values.JxW(k);
            }

            local_rhs[i] = 0;
            for (int k=0; k<q.size(); k++)
            {
                vec3 p = fe_values.point(k);
                local_rhs[i] += (p[0]*(1-p[0])+p[1]*(1-p[1]))*fe_values.shape_value(i,k)*fe_values.JxW(k);
            }
        }

        ls->set_values(ndofs, (int *)dof_indices, ndofs, (int *)dof_indices, local_matrix, local_rhs);
    }

    // assemble integral over sides
    for(EdgeFullIter edge = mesh().edge.begin(); edge!=mesh().edge.end(); ++edge)
    {
        ASSERT(edge->n_sides <= 2, "Don't know how to treat an edge with more than 2 sides.");

        bool skip_edge = false;
        for (int sid=0; sid<edge->n_sides; sid++)
            if (edge->side[sid]->dim != 1 || mesh().element.full_iter(edge->side[sid]->element)->dim != 2)
            {
                skip_edge = true;
                break;
            }
        if (skip_edge) continue;
    
        for (int sid=0; sid<edge->n_sides; sid++)
        {
            cell = mesh().element.full_iter(edge->side[sid]->element);
            dof_handler2d->get_dof_indices(cell, side_dof_indices[sid]);
            fe_values_side[sid]->reinit(cell, edge->side[0]);
        }

        double det = fe_values_side[0]->JxW(0)/q1d.weight(0);
        if (edge->n_sides == 1)
        {
            omega[0] = 1;
            gamma = alpha/det; // this should be simply alpha/determinant, but we do not have direct access to the side determinant
        }
        else
        {
            omega[0] = 0.5;
            omega[1] = 0.5;
            gamma = alpha*0.5/det; // this should be simply alpha/determinant, but we do not have direct access to the side determinant
        }

        // boundary condition
        for (int i=0; i<ndofs; i++)
        {
            local_rhs[i] = 0;
            if (edge->n_sides==1)
                for (int k=0; k<q1d.size(); k++)
                    local_rhs[i] += gamma*0*fe_values_side[0]->shape_value(i,k)
                                    *fe_values_side[0]->JxW(k);
        }

        // fluxes and penalty
        for (int s1=0; s1<edge->n_sides; s1++)
        {
            for (int s2=0; s2<edge->n_sides; s2++)
            {
                for (int i=0; i<ndofs; i++)
                    for (int j=0; j<ndofs; j++)
                    {
                        local_flux_matrix[i*ndofs+j] = 0;
                        for (int k=0; k<q1d.size(); k++)
                        {
                            local_flux_matrix[i*ndofs+j] += (
                                    ((s1==s2)?1.:-1.)*gamma*fe_values_side[s1]->shape_value(j,k)*fe_values_side[s2]->shape_value(i,k)
                                    +((s2==0)?-1.:1.)*dot(omega[s1]*fe_values_side[s1]->shape_grad(j,k),fe_values_side[0]->normal_vector(k))*fe_values_side[s2]->shape_value(i,k)
                                                            )*fe_values_side[0]->JxW(k);
                        }
                    }
                ls->set_values(ndofs, (int *)side_dof_indices[s2], ndofs, (int *)side_dof_indices[s1], local_flux_matrix, local_rhs);
            }
        }
    }
}

