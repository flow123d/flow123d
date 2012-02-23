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
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
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
#include "fem/fe_rt.hh"
#include "io/output.h"
#include "xio.h"
#include <iostream>
#include <iomanip>
#include "mesh/boundaries.h"


TransportDG::TransportDG(TimeMarks & marks, Mesh & init_mesh, MaterialDatabase & material_database)
        : TransportBase(marks, init_mesh, material_database),
          tol_flux_bc(1e-6),
          advection(1e0)
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
     * Set up physical parameters.
     */
    alphaL = OptGetDbl("Transport", "Longitudal_dispersivity", "20e-30");
    alphaT = OptGetDbl("Transport", "Transversal_dispersivity", "5e-30");
    molecular_diffusivity = OptGetDbl("Transport", "Molecular_diffusivity", "1e-9");
    tortuosity = 1;
//    printf("Dispersivities:\nDm = %f\nalphaL = %f\nalphaT = %f\n", molecular_diffusivity, alphaL, alphaT);


    /*
     * Read names of transported substances.
     */
    n_substances = OptGetInt("Transport", "N_substances", NULL );
    INPUT_CHECK(n_substances >= 1 ,"Number of substances must be positive.\n");
    read_subst_names();


    /*
     * Read boundary conditions.
     */
    bcv = new Vec[n_substances];
    bc_values = new double*[n_substances];
    for (int sbi=0; sbi<n_substances; sbi++)
    {
        bc_values[sbi] = new double[mesh_->n_boundaries()];
        VecCreateSeqWithArray(PETSC_COMM_SELF, mesh_->n_boundaries(), bc_values[sbi], &bcv[sbi]);
    }

    // read times of time dependent boundary condition and check the input files

    // bc marks do not influent choice of time step (not fixed_time type)
    bc_mark_type_ = time_marks->type_bc_change() | equation_mark_type_;
    std::vector<double> bc_times;
    OptGetDblArray("Transport", "bc_times", "", bc_times);
    FILE * f;
    std::string fname;
    if (bc_times.size() == 0) {
        // only one boundary condition, check filename and read bc condition
        bc_time_level = -1;
        fname = make_bc_file_name(-1);
        if ( !(f = xfopen(fname.c_str(), "rt")) ) {
            xprintf(UsrErr,"Missing file: %s", fname.c_str());
        }
        xfclose(f);
        read_bc_vector();
    } else {
        bc_time_level = 0;
        // set initial bc to zero
        for(unsigned int sbi=0; sbi < n_substances; ++sbi) VecZeroEntries(bcv[sbi]);
        // check files for bc time levels and set TimeMarks
        for(unsigned int i=0;i<bc_times.size();++i) {
            fname = make_bc_file_name(i);
            if ( !(f = xfopen(fname.c_str(), "rt")) ) {
                xprintf(UsrErr,"Missing file: %s", fname.c_str());
            }
            xfclose(f);

            time_marks->add(TimeMark(bc_times[i],bc_mark_type_));
        }
    }


    // distribute DOFs
    dof_handler1d = new DOFHandler<1,3>(mesh());
    dof_handler2d = new DOFHandler<2,3>(mesh());
    dof_handler3d = new DOFHandler<3,3>(mesh());
    fe1d = new FE_P_disc<1,1,3>;
    fe2d = new FE_P_disc<1,2,3>;
    fe3d = new FE_P_disc<1,3,3>;
    dof_handler1d->distribute_dofs(*fe1d);
    dof_handler2d->distribute_dofs(*fe2d, dof_handler1d->n_global_dofs());
    dof_handler3d->distribute_dofs(*fe3d, dof_handler1d->n_global_dofs() + dof_handler2d->n_global_dofs());

    // set up output class
    string output_file = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Transport", "Transport_out", "\\"));
    transport_output = new OutputTime(mesh_, output_file);
    output_solution.resize(n_substances);
    for (int i=0; i<n_substances; i++)
    {
        output_solution[i] = new double[mesh().node_vector.size()];
        transport_output->register_node_data<double>(substance_names[0], "M/L^3", output_solution[i], mesh().node_vector.size());
    }
    

    ls    = new LinSys_MPIAIJ(dof_handler1d->n_global_dofs() + dof_handler2d->n_global_dofs() + dof_handler3d->n_global_dofs());
    ls_dt = new LinSys_MPIAIJ(dof_handler1d->n_global_dofs() + dof_handler2d->n_global_dofs() + dof_handler3d->n_global_dofs());


    // assemble mass matrix
    ls_dt->start_allocation();
    assemble_mass_matrix<1>(dof_handler1d, fe1d);
    assemble_mass_matrix<2>(dof_handler2d, fe2d);
    assemble_mass_matrix<3>(dof_handler3d, fe3d);
    ls_dt->start_add_assembly();
    assemble_mass_matrix<1>(dof_handler1d, fe1d);
    assemble_mass_matrix<2>(dof_handler2d, fe2d);
    assemble_mass_matrix<3>(dof_handler3d, fe3d);
    ls_dt->finalize();
    mass_matrix = ls_dt->get_matrix();

    // preallocate system matrix
    ls->start_allocation();
    assemble<1>(dof_handler1d, fe1d);
    assemble<2>(dof_handler2d, fe2d);
    assemble<3>(dof_handler3d, fe3d);
    set_boundary_conditions<1>(dof_handler1d, fe1d);
    set_boundary_conditions<2>(dof_handler2d, fe2d);
    set_boundary_conditions<3>(dof_handler3d, fe3d);
    stiffness_matrix = 0;


    /*
     * Read initial condition.
     */
    read_initial_condition();

    // save initial state
    output_data();
}


TransportDG::~TransportDG()
{
    delete transport_output;
    delete time_;
    delete solver;
    delete ls;
    delete ls_dt;
    for (int i=0; i<n_substances; i++) delete[] output_solution[i];

    gamma.clear();
    substance_names.clear();
}



void TransportDG::update_solution()
{
    time_->next_time();
    time_->view();
    
/*    // Crank-Nicholson
    Vec old_rhs;
    if (time_->tlevel() == 1)
    {
        VecCreateSeq(PETSC_COMM_SELF, ls->size(), &rhs);
        VecSet(rhs, 0);
        VecDuplicate(rhs, &old_rhs);
        VecCopy(rhs, old_rhs);

    }
    else
    {
        VecDuplicate(rhs, &old_rhs);
        VecCopy(rhs, old_rhs);
        VecScale(old_rhs, 0.5*time_->dt());
        MatMultAdd(ls->get_matrix(), ls->get_solution(), old_rhs, old_rhs);
    }*/

    // assemble system matrix
    if (flux_changed || (bc_time_level != -1 && time_->is_current(bc_mark_type_)))
    {
        flux_changed = false;

        // possibly read boundary conditions
        if (bc_time_level != -1 && time_->is_current(bc_mark_type_)) read_bc_vector();

        ls->start_add_assembly();
        MatZeroEntries(ls->get_matrix());
        VecSet(ls->get_rhs(), 0);
        assemble<1>(dof_handler1d, fe1d);
        assemble<2>(dof_handler2d, fe2d);
        assemble<3>(dof_handler3d, fe3d);
        set_boundary_conditions<1>(dof_handler1d, fe1d);
        set_boundary_conditions<2>(dof_handler2d, fe2d);
        set_boundary_conditions<3>(dof_handler3d, fe3d);
        ls->finalize();

        if (stiffness_matrix == 0)
            MatConvert(ls->get_matrix(), MATSAME, MAT_INITIAL_MATRIX, &stiffness_matrix);
        else
            MatCopy(ls->get_matrix(), stiffness_matrix, DIFFERENT_NONZERO_PATTERN);
        VecDuplicate(ls->get_rhs(), &rhs);
        VecCopy(ls->get_rhs(), rhs);
    }



    /* Apply backward Euler time integration.
     *
     * Denoting A the stiffness matrix and M the mass matrix, the algebraic system at the k-th time level reads
     *
     *   (1/dt M + A)u^k = f + 1/dt M.u^{k-1}
     *
     * Hence we modify at each time level the right hand side:
     *
     *   f^k = f + 1/dt M u^{k-1},
     *
     * where f stands for the term stemming from the force and boundary conditions.
     * Accordingly, we set
     *
     *   A^k = A + 1/dt M.
     *
     */

    MatCopy(stiffness_matrix, ls->get_matrix(), DIFFERENT_NONZERO_PATTERN);
    MatAXPY(ls->get_matrix(), 1./time_->dt(), mass_matrix, SUBSET_NONZERO_PATTERN);
    Vec y;
    VecDuplicate(rhs, &y);
    MatMult(mass_matrix, ls->get_solution(), y);
    VecWAXPY(ls->get_rhs(), 1./time_->dt(), y, rhs);


    /* Apply Crank-Nicholson time integration.
     *
     */
/*    VecWAXPY(ls->get_rhs(), 0.5*time_->dt(), rhs, old_rhs);
    MatCopy(stiffness_matrix, ls->get_matrix(), DIFFERENT_NONZERO_PATTERN);
    MatAYPX(ls->get_matrix(), -0.5*time_->dt(), mass_matrix, SUBSET_NONZERO_PATTERN);*/


    // solve
    solve_system(solver, ls);

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
    // So far the velocity_vector contains zeros, so we ignore it.
    // Instead we use the value Side.flux.
	flux_vector = velocity_vector;

	flux_changed = true;
}



void TransportDG::output_data()
{
    double *solution;
    unsigned int dof_indices[max(fe1d->n_dofs(), max(fe2d->n_dofs(), fe3d->n_dofs()))];
    int n_nodes = mesh_->node_vector.size();
    int count[n_nodes];

    solution = ls->get_solution_array();

    // interpolate solution to a continuous field and write to output file


    for (int i=0; i< n_nodes; i++)
    {
        count[i] = 0;
        output_solution[0][i] = 0;
    }
    for (DOFHandler<1,3>::CellIterator cell = dof_handler1d->begin_cell(); cell != dof_handler1d->end_cell(); ++cell)
    {
        if (cell->dim != 1) continue;

        dof_handler1d->get_dof_indices(cell, dof_indices);

        for (int i=0; i<fe1d->n_dofs(); i++)
        {
            int nid = mesh_->node_vector.index(cell->node[i]);
            count[nid]++;
            output_solution[0][nid] += solution[dof_indices[i]];
        }
    }
    for (DOFHandler<2,3>::CellIterator cell = dof_handler2d->begin_cell(); cell != dof_handler2d->end_cell(); ++cell)
    {
        if (cell->dim != 2) continue;

        dof_handler2d->get_dof_indices(cell, dof_indices);

        for (int i=0; i<fe2d->n_dofs(); i++)
        {
            int nid = mesh_->node_vector.index(cell->node[i]);
            count[nid]++;
            output_solution[0][nid] += solution[dof_indices[i]];
        }
    }
    for (DOFHandler<3,3>::CellIterator cell = dof_handler3d->begin_cell(); cell != dof_handler3d->end_cell(); ++cell)
    {
        if (cell->dim != 3) continue;

        dof_handler3d->get_dof_indices(cell, dof_indices);

        for (int i=0; i<fe3d->n_dofs(); i++)
        {
            int nid = mesh_->node_vector.index(cell->node[i]);
            count[nid]++;
            output_solution[0][nid] += solution[dof_indices[i]];
        }
    }

    for (int i=0; i<n_nodes; i++)
        if (count[i] > 1)
            output_solution[0][i] /= count[i];

    transport_output->write_data(time_->t());
}



template<unsigned int dim>
void TransportDG::assemble_mass_matrix(DOFHandler<dim,3> *dh, FiniteElement<dim,3> *fe)
{
    MappingP1<dim,3> map;
    QGauss<dim> q(2);
    FEValues<dim,3> fe_values(map, q, *fe, update_values | update_JxW_values);
    unsigned int ndofs = fe->n_dofs();
    unsigned int dof_indices[ndofs];
    PetscScalar local_mass_matrix[ndofs*ndofs], local_rhs[ndofs];

    typename DOFHandler<dim,3>::CellIterator cell = dh->begin_cell();

    // assemble integral over elements
    for (cell = dh->begin_cell(); cell != dh->end_cell(); ++cell)
    {
        if (cell->dim != dim) continue;

        fe_values.reinit(cell);

        dh->get_dof_indices(cell, dof_indices);

        // assemble the local stiffness and mass matrix
        for (int i=0; i<ndofs; i++)
        {
            for (int j=0; j<ndofs; j++)
            {
                local_mass_matrix[i*ndofs+j] = 0;
                for (int k=0; k<q.size(); k++)
                    local_mass_matrix[i*ndofs+j] += fe_values.shape_value(j,k)*fe_values.shape_value(i,k)*fe_values.JxW(k);
            }
            local_rhs[i] = 0;
        }

        ls_dt->set_values(ndofs, (int *)dof_indices, ndofs, (int *)dof_indices, local_mass_matrix, local_rhs);
    }
}







template<unsigned int dim>
void TransportDG::assemble(DOFHandler<dim,3> *dh, FiniteElement<dim,3> *fe)
{
    // coefficient of DG method influencing the continuity across edges
    double alpha;
    MappingP1<dim,3> map;
    QGauss<dim> q(2);
    QGauss<dim-1> side_q(2);
    FE_RT0<dim,3> fe_rt;
    FEValues<dim,3> fv_rt(map, q, fe_rt, update_values | update_gradients);
    FEValues<dim,3> fe_values(map, q, *fe, update_values | update_gradients | update_JxW_values | update_quadrature_points);
    vector<FESideValues<dim,3>*> fe_values_side;
    FESideValues<dim,3> fsv_rt(map, side_q, fe_rt, update_values | update_normal_vectors | update_side_JxW_values);
    typename DOFHandler<dim,3>::CellIterator cell = dh->begin_cell();

    const unsigned int ndofs = fe->n_dofs();
    unsigned int dof_indices[ndofs];
    vector<unsigned int*> side_dof_indices;
    PetscScalar local_matrix[ndofs*ndofs], local_rhs[ndofs];
    vector<mat33> K;
    vector< vector<mat33> > side_K;
    vector<vec3> velocity;
    vector< vector<vec3> > side_velocity;
    double gamma_l;
    double omega[2];
    double transport_flux;
    double *solution = ls->get_solution_array();

    gamma.resize(mesh_->boundary.size());

    if (molecular_diffusivity != 0)
        alpha = 1./molecular_diffusivity;

	// assemble integral over elements
    for (cell = dh->begin_cell(); cell != dh->end_cell(); ++cell)
    {
        if (cell->dim != dim) continue;

        fe_values.reinit(cell);
        fv_rt.reinit(cell);
        
        calculate_velocity(cell, velocity, fv_rt);
        calculate_dispersivity_tensor(K, velocity);

        dh->get_dof_indices(cell, dof_indices);

        // assemble the local stiffness matrix
        for (int i=0; i<ndofs; i++)
        {
            for (int j=0; j<ndofs; j++)
            {
                local_matrix[i*ndofs+j] = 0;
                for (int k=0; k<q.size(); k++)
                {
                    local_matrix[i*ndofs+j] += (dot(K[k]*fe_values.shape_grad(j,k),fe_values.shape_grad(i,k))
                                               +dot(fe_values.shape_grad(j,k),velocity[k])*fe_values.shape_value(i,k)*advection
                                               )*fe_values.JxW(k);
                }

            }
        }

        ls->mat_set_values(ndofs, (int *)dof_indices, ndofs, (int *)dof_indices, local_matrix/*, local_rhs*/);
    }

    // assemble integral over sides
    for(EdgeFullIter edge = mesh().edge.begin(); edge!=mesh().edge.end(); ++edge)
    {
        bool skip_edge = false;
        for (int sid=0; sid<edge->n_sides; sid++)
            if (edge->side[sid]->dim != dim-1 || mesh().element.full_iter(edge->side[sid]->element)->dim != dim)
            {
                skip_edge = true;
                break;
            }
        if (skip_edge) continue;

        side_K.resize(edge->n_sides);
        side_velocity.resize(edge->n_sides);

        if (side_dof_indices.size() < edge->n_sides)
            for (int i=side_dof_indices.size(); i<edge->n_sides; i++)
                side_dof_indices.push_back(new unsigned int[ndofs]);

        if (fe_values_side.size() < edge->n_sides)
            for (int sid=fe_values_side.size(); sid<edge->n_sides; sid++)
                fe_values_side.push_back(new FESideValues<dim,3>(map, side_q, *fe, update_values | update_gradients | update_side_JxW_values | update_normal_vectors));
        
        for (int sid=0; sid<edge->n_sides; sid++)
        {
            cell = mesh().element.full_iter(edge->side[sid]->element);
            dh->get_dof_indices(cell, side_dof_indices[sid]);
            fe_values_side[sid]->reinit(cell, edge->side[0]);
            fsv_rt.reinit(cell, edge->side[0]);
            
            calculate_velocity(cell, side_velocity[sid], fsv_rt);
            calculate_dispersivity_tensor(side_K[sid], side_velocity[sid]);
        }

        // fluxes and penalty
        if (edge->n_sides == 1)
        {
            // set up the parameters for DG method
            set_DG_parameters(edge, 0, 0, side_q.size(), side_K, fe_values_side[0]->normal_vector(0), alpha, advection, gamma_l, omega, transport_flux);

            if (edge->side[0]->cond != 0) gamma[mesh_->boundary.full_iter(edge->side[0]->cond).id()] = gamma_l;

            for (int i=0; i<ndofs; i++)
                for (int j=0; j<ndofs; j++)
                {
                    local_matrix[i*ndofs+j] = 0;
                    for (int k=0; k<side_q.size(); k++)
                    {
                        // penalty enforcing continuity across edges (applied on interior and Dirichlet edges)
                        if (edge->side[0]->cond!=0 && edge->side[0]->flux < -tol_flux_bc)
                        {
                            local_matrix[i*ndofs+j] += gamma_l*fe_values_side[0]->shape_value(j,k)*fe_values_side[0]->shape_value(i,k)*fe_values_side[0]->JxW(k);
                        }
                        // flux due to transport (applied on interior edges)
//                        if (edge->side[0]->cond==0 || edge->side[0]->flux >= -tol_flux_bc)
//                            local_matrix[i*ndofs+j] -= advection*0.5*transport_flux*fe_values_side[0]->shape_value(j,k)*fe_values_side[0]->shape_value(i,k)*fe_values_side[0]->JxW(k);
                        // terms due to diffusion
                        local_matrix[i*ndofs+j] -= dot(side_K[0][k]*fe_values_side[0]->shape_grad(j,k),fe_values_side[0]->normal_vector(k))*fe_values_side[0]->shape_value(i,k)*fe_values_side[0]->JxW(k);
                    }
                }
            ls->mat_set_values(ndofs, (int *)side_dof_indices[0], ndofs, (int *)side_dof_indices[0], local_matrix);
        }
        else // if edge->n_sides != 1
        {
            for (int s1=0; s1<edge->n_sides; s1++)
            {
                for (int s2=s1+1; s2<edge->n_sides; s2++)
                {
                    // set up the parameters for DG method
                    set_DG_parameters(edge, s1, s2, side_q.size(), side_K, fe_values_side[0]->normal_vector(0), alpha, advection, gamma_l, omega, transport_flux);

                    if (edge->side[s1]->cond != 0) gamma[mesh_->boundary.full_iter(edge->side[s1]->cond).id()] = gamma_l;

                    int sd[2];
                    sd[0] = s1;
                    sd[1] = s2;
                    for (int m=0; m<2; m++)
                    {
                        for (int n=0; n<2; n++)
                        {
                            for (int i=0; i<ndofs; i++)
                                for (int j=0; j<ndofs; j++)
                                {
                                    local_matrix[i*ndofs+j] = 0;
                                    for (int k=0; k<side_q.size(); k++)
                                    {
                                        // penalty enforcing continuity across edges (applied on interior and Dirichlet edges)
                                        local_matrix[i*ndofs+j] += (m==n?1:-1)*gamma_l*fe_values_side[sd[m]]->shape_value(j,k)*fe_values_side[sd[n]]->shape_value(i,k)*fe_values_side[0]->JxW(k);
                                        // flux due to transport (applied on interior edges)
                                        local_matrix[i*ndofs+j] -= (m==0?1:-1)*advection*0.5*transport_flux*fe_values_side[sd[m]]->shape_value(j,k)*fe_values_side[sd[n]]->shape_value(i,k)*fe_values_side[0]->JxW(k);
                                        // terms due to diffusion
                                        local_matrix[i*ndofs+j] += ((n==0)?-1.:1.)*dot(omega[0]*side_K[sd[m]][k]*fe_values_side[sd[m]]->shape_grad(j,k),fe_values_side[0]->normal_vector(k))*fe_values_side[sd[n]]->shape_value(i,k)*fe_values_side[0]->JxW(k);
                                        local_matrix[i*ndofs+j] += ((m==0)?-1.:1.)*dot(omega[1]*side_K[sd[n]][k]*fe_values_side[sd[n]]->shape_grad(i,k),fe_values_side[0]->normal_vector(k))*fe_values_side[sd[m]]->shape_value(j,k)*fe_values_side[0]->JxW(k);
                                    }
                                }
                            ls->mat_set_values(ndofs, (int *)side_dof_indices[sd[n]], ndofs, (int *)side_dof_indices[sd[m]], local_matrix);
                        }
                    }
                }
            }
        }
    }

    for (int i=0; i<fe_values_side.size(); i++)
        delete fe_values_side[i];

    for (int i=0; i<side_dof_indices.size(); i++)
        delete[] side_dof_indices[i];
}


template<unsigned int dim>
void TransportDG::set_boundary_conditions(DOFHandler<dim,3> *dh, FiniteElement<dim,3> *fe)
{
    typename DOFHandler<dim,3>::CellIterator cell = dh->begin_cell();
    MappingP1<dim,3> map;
    QGauss<dim-1> side_q(2);
    FESideValues<dim,3> fe_values_side(map, side_q, *fe, update_values | update_side_JxW_values);
    unsigned int side_dof_indices[fe->n_dofs()];
    double local_rhs[fe->n_dofs()*fe->n_dofs()];

    for (BoundaryFullIter b = mesh_->boundary.begin(); b != mesh_->boundary.end(); ++b)
    {
        cell = mesh_->element.full_iter(b->side->element);
        if (cell->dim != dim || b->side->flux >= -tol_flux_bc) continue;

        fe_values_side.reinit(cell, b->side);
        dh->get_dof_indices(cell, side_dof_indices);

        for (int i=0; i<fe->n_dofs(); i++)
        {
            local_rhs[i] = 0;
            for (int k=0; k<side_q.size(); k++)
            {
                local_rhs[i] += (
                        +gamma[b.id()]*bc_values[0][b.id()]*fe_values_side.shape_value(i,k)
//                       -advection*0.5*min(b->side->flux,-tol_flux_bc)*bc_values[0][b.id()]*fe_values_side.shape_value(i,k)
                                )*fe_values_side.JxW(k);
            }
        }
        ls->rhs_set_values(fe->n_dofs(), (int *)side_dof_indices, local_rhs);
    }
}



template<unsigned int dim>
void TransportDG::calculate_velocity(typename DOFHandler<dim,3>::CellIterator cell, vector<vec3> &velocity, FEValuesBase<dim,3> &fv)
{
    std::map<Node*, int> node_nums;
    for (int i=0; i<cell->n_nodes; i++)
        node_nums[cell->node[i]] = i;

    velocity.resize(fv.n_points());

    for (int k=0; k<fv.n_points(); k++)
    {
        velocity[k].zeros();
        for (int sid=0; sid<cell->n_sides; sid++)
        {
            if (cell->side[sid]->dim != dim-1) continue;
            int num = dim*(dim+1)/2;
            for (int i=0; i<cell->side[sid]->n_nodes; i++)
                num -= node_nums[cell->side[sid]->node[i]];
            velocity[k] += fv.shape_vector(num,k)*cell->side[sid]->flux;
        }
    }
}




void TransportDG::calculate_dispersivity_tensor(vector<mat33> &K, vector<vec3> &velocity)
{
    double vnorm;

    K.resize(velocity.size());

    for (int k=0; k<velocity.size(); k++)
    {
        vnorm = norm(velocity[k], 2);

        if (fabs(vnorm)>1e-6)
            for (int i=0; i<3; i++)
                for (int j=0; j<3; j++)
                    K[k](i,j) = velocity[k][i]*velocity[k][j]/vnorm*(alphaL-alphaT);
        else
            K[k].zeros();

        K[k] += eye(3,3)*(molecular_diffusivity*tortuosity + alphaT*vnorm);
    }
}


void TransportDG::set_DG_parameters(const Edge *edge,
            const int s1,
            const int s2,
            const unsigned int n_points,
            const vector< vector<mat33> > &K,
            const vec3 &normal_vector,
            const double alpha,
            const double advection,
            double &gamma,
            double *omega,
            double &transport_flux)
{
    double delta[2];
    double h = 0;

    if (edge->side[s1]->dim == 0)
    {
        h = 1;
    }
    else
    {
        for (int i=0; i<edge->side[s1]->n_nodes; i++)
            for (int j=i+1; j<edge->side[s1]->n_nodes; j++)
                h = max(h, edge->side[s1]->node[i]->distance(*edge->side[s1]->node[j]));
    }

    // calculate the total in- and out-flux through the edge
    double pflux = 0, nflux = 0;
    for (int i=0; i<edge->n_sides; i++)
    {
        if (edge->side[i]->flux > 0)
            pflux += edge->side[i]->flux;
        else
            nflux += edge->side[i]->flux;
    }

    if (edge->side[s2]->flux > 0 && edge->side[s1]->flux < 0 && s1 < s2)
    {
        transport_flux = edge->side[s1]->flux*fabs(edge->side[s2]->flux/pflux);
    }
    else if (edge->side[s2]->flux < 0 && edge->side[s1]->flux > 0 && s1 < s2)
    {
        transport_flux = edge->side[s1]->flux*fabs(edge->side[s2]->flux/nflux);
    }
    else if (s1==s2)
    {
        transport_flux = edge->side[s1]->flux;
    }
    else
    {
        transport_flux = 0;
    }

    gamma = 0.5*advection*fabs(transport_flux);

    if (s1 == s2)
    {
        omega[0] = 1;

        // delta is set to the average value of Kn.n on the side
        delta[0] = 0;
        for (int k=0; k<n_points; k++)
            delta[0] += dot(K[s1][k]*normal_vector,normal_vector);
        delta[0] /= n_points;

        gamma += alpha/h*delta[0];
    }
    else
    {
        delta[0] = 0;
        delta[1] = 0;
        for (int k=0; k<n_points; k++)
        {
            delta[0] += dot(K[s1][k]*normal_vector,normal_vector)/n_points;
            delta[1] += dot(K[s2][k]*normal_vector,normal_vector)/n_points;
        }

        double delta_sum = delta[0] + delta[1];

        // the following has to be modified for the case of more than 2 sides per edge
        if (delta_sum > 1e-15)
        {
            omega[0] = delta[1]/delta_sum;
            omega[1] = delta[0]/delta_sum;
            gamma += alpha*delta[0]*delta[1]/(delta_sum*h);
        }
        else
            for (int i=0; i<2; i++) omega[i] = 0;
    }

}



std::string TransportDG::make_bc_file_name(int level)
{

    std::string bc_fname = IONameHandler::get_instance()->get_input_file_name(OptGetFileName("Transport",
            "Transport_BCD", "\\"));
    if (level >= 0 ) {
        std::stringstream name_str;
        name_str << bc_fname << "_" << setfill('0') << setw(3) << level;
        bc_fname = name_str.str();
    }

    return bc_fname;
}


void TransportDG::read_bc_vector()
{
    int sbi;

    xprintf(Msg, "Transport_DG: Reading BC (level = %d)\n",bc_time_level);

    int bcd_id, boundary_id, boundary_index;
    double bcd_conc;
    char line[LINE_SIZE]; // line of data file

    // make bc filename

    FILE *in = xfopen(make_bc_file_name(bc_time_level).c_str(), "rt");
    skip_to(in, "$Transport_BCD");
    xfgets(line, LINE_SIZE - 2, in);
    int n_bcd = atoi(xstrtok(line));
    for (int i_bcd = 0; i_bcd < n_bcd; i_bcd++) {
        xfgets(line, LINE_SIZE - 2, in);
        bcd_id = atoi(xstrtok(line)); // scratch transport bcd id
        boundary_id = atoi(xstrtok(NULL));
        //        DBGMSG("transp b. id: %d\n",boundary_id);
        boundary_index = mesh_->boundary.find_id(boundary_id).index();
        INPUT_CHECK(boundary_index >= 0,"Wrong boundary index %d for bcd id %d in transport bcd file!", boundary_id, bcd_id);
        for (sbi = 0; sbi < n_substances; sbi++) {
            bcd_conc = atof(xstrtok(NULL));
            VecSetValue(bcv[sbi], boundary_index, bcd_conc, INSERT_VALUES);
        }
    }
    xfclose(in);
    for (sbi = 0; sbi < n_substances; sbi++)
        VecAssemblyBegin(bcv[sbi]);
    //for (sbi = 0; sbi < n_substances; sbi++)
    //    VecZeroEntries(bcvcorr[sbi]);
    for (sbi = 0; sbi < n_substances; sbi++)
        VecAssemblyEnd(bcv[sbi]);

    // update source vectors
    // TODO: rather use Lazy dependency
//    if (bc_time_level != -1)
//        for (unsigned int sbi = 0; sbi < n_substances; sbi++) {
//            MatMult(bcm, bcv[sbi], bcvcorr[sbi]);
//        }

    // VecView(bcvcorr[0],PETSC_VIEWER_STDOUT_SELF);
    // getchar();
    bc_time_level++;

}


void TransportDG::read_initial_condition()
{
    FILE    *in;          // input file
    char     line[ LINE_SIZE ]; // line of data file
    int sbi,index, id, eid, i,n_concentrations;
    double value;

    std::string concentration_fname = IONameHandler::get_instance()->get_input_file_name(OptGetFileName("Transport", "Concentration", "\\"));
    in = xfopen( concentration_fname, "rt" );

    skip_to( in, "$Concentrations" );
    xfgets( line, LINE_SIZE - 2, in );
    n_concentrations = atoi( xstrtok( line) );
    INPUT_CHECK(!( n_concentrations < 1 ),"Number of concentrations < 1 in function read_concentration_list()\n");
    INPUT_CHECK(!( mesh_->n_elements() != n_concentrations),"Different number of elements and concentrations\n");

    unsigned int maxndofs = max(fe1d->n_dofs(), max(fe2d->n_dofs(), fe3d->n_dofs()));
    double values[maxndofs];
    unsigned int ndofs, dof_indices[maxndofs];

    for (i = 0; i < n_concentrations; i++)
    {
        xfgets( line, LINE_SIZE - 2, in );
        ASSERT(!(line == NULL),"NULL as argument of function parse_concentration_line()\n");
        id    = atoi( xstrtok( line) ); // TODO: id musi byt >0 nebo >=0 ???
        INPUT_CHECK(!( id < 0 ),"Id number of concentration must be > 0\n");
        eid   = atoi( xstrtok( NULL) );
        ElementFullIter cell = mesh_->element.find_id(eid);

        switch (cell->dim)
        {
        case 1:
            dof_handler1d->get_dof_indices(cell, dof_indices);
            ndofs = fe1d->n_dofs();
            break;
        case 2:
            dof_handler2d->get_dof_indices(cell, dof_indices);
            ndofs = fe2d->n_dofs();
            break;
        case 3:
            dof_handler3d->get_dof_indices(cell, dof_indices);
            ndofs = fe3d->n_dofs();
            break;
        }

        for( sbi = 0; sbi < n_substances; sbi++ )
        {
            value = atof( xstrtok( NULL) );
            for (int j=0; j<ndofs; j++)
                values[j] = value;
            if (sbi==0) VecSetValues(ls->get_solution(), ndofs, (int *)dof_indices, values, INSERT_VALUES);
        }
    }

    xfclose( in );
    VecAssemblyBegin(ls->get_solution());
    VecAssemblyEnd(ls->get_solution());
}


void TransportDG::read_subst_names()
{
    int sbi;
    char *line = OptGetStr("Transport", "Substances", "none");

    ASSERT(!( (n_substances < 1) || (line == NULL) ),"Bad parameter of function subst_names()\n");
    for (sbi = 0; sbi < n_substances; sbi++)
        substance_names.push_back(strtok((sbi == 0 ? line : NULL), " \t,;"));

    xfree(line);
}


