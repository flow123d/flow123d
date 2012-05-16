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

#include "petscmat.h"
#include <armadillo>
#include "xio.h"
#include "transport/transport_dg.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/mapping_p1.hh"
#include "la/solve.h"
#include "fem/fe_rt.hh"
#include "io/output.h"
#include "mesh/boundaries.h"
#include "system/par_distribution.hh"
#include "transport/transport_bc.hh"


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
    Dm = OptGetDbl("Transport", "Molecular_diffusivity", "1e-6");
    tortuosity = 1;
//    printf("Dispersivities:\nDm = %f\nalphaL = %f\nalphaT = %f\n", molecular_diffusivity, alphaL, alphaT);

    if (Dm != 0)
    	alpha = max(1e0,1e0/Dm);


    /*
     * Read names of transported substances.
     */
    n_substances = OptGetInt("Transport", "N_substances", NULL );
    INPUT_CHECK(n_substances >= 1 ,"Number of substances must be positive.\n");
    read_subst_names();

    sorption = OptGetBool("Transport", "Sorption", "no");
	dual_porosity = OptGetBool("Transport", "Dual_porosity", "no");
	mat_base->read_transport_materials(dual_porosity, sorption,n_substances);

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

	// distribute solution vectors on processors
	distr = new Distribution(Distribution::Block, dof_handler1d->n_global_dofs() + dof_handler2d->n_global_dofs() + dof_handler3d->n_global_dofs());

    /*
     * Read boundary conditions.
     */
    bc = new TransportBC(mesh_, n_substances);
    bc_mark_type_ = time_marks->type_bc_change() | equation_mark_type_;
    if (bc->get_times().size() > 0)
    {
    	for (int i=0; i<bc->get_times().size(); i++)
    		time_marks->add(TimeMark(bc->get_times()[i], bc_mark_type_));
    }


    // set up output class
    string output_file = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Transport", "Transport_out", "\\"));
    transport_output = new OutputTime(mesh_, output_file);
    output_solution.resize(n_substances);
    for (int i=0; i<n_substances; i++)
    {
        output_solution[i] = new double[mesh().node_vector.size()];
        transport_output->register_node_data<double>(substance_names[i], "M/L^3", output_solution[i], mesh().node_vector.size());
    }
    output_cell_solution.resize(n_substances);
	for (int i=0; i<n_substances; i++)
	{
		output_cell_solution[i] = new double[mesh().element.size()];
		transport_output->register_elem_data<double>(substance_names[i]+"_cell", "M/L^3", output_cell_solution[i], mesh().element.size());
	}

	// set time marks for writing the output
	output_mark_type = this->mark_type() | time_marks->type_fixed_time() | time_marks->type_output();
    time_marks->add_time_marks(0.0, OptGetDbl("Global", "Save_step", "1.0"), time_->end_time(), output_mark_type);
    

    ls    = new LinSys_MPIAIJ(distr->lsize());
    ls_dt = new LinSys_MPIAIJ(distr->lsize());


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
    assemble_stiffness_matrix();
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
    for (int i=0; i<n_substances; i++)
    {
    	delete[] output_solution[i];
    	delete[] output_cell_solution[i];
    }

    gamma.clear();
    substance_names.clear();
}



void TransportDG::update_solution()
{
    time_->next_time();
    time_->view();
    
    // assemble system matrix if velocity flux or boundary conditions changed
    if (flux_changed || (bc->get_time_level() != -1 && time_->is_current(bc_mark_type_)))
    {
        flux_changed = false;

        // possibly read boundary conditions
        if (bc->get_time_level() != -1 && time_->is_current(bc_mark_type_)) bc->read();

        // new fluxes can change the location of Neumann boundary,
        // thus stiffness matrix must be reassembled
        ls->start_add_assembly();
        MatZeroEntries(ls->get_matrix());
        VecSet(ls->get_rhs(), 0);
        assemble_stiffness_matrix();
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
    // ls->get_matrix() = 1/dt*mass_matrix + ls->get_matrix()
    MatAXPY(ls->get_matrix(), 1./time_->dt(), mass_matrix, SUBSET_NONZERO_PATTERN);
    Vec y;
    VecDuplicate(rhs, &y);
    // y = mass_matrix*ls->get_solution()
    MatMult(mass_matrix, ls->get_solution(), y);
    // ls->get_rhs() = 1/dt*y + rhs
    VecWAXPY(ls->get_rhs(), 1./time_->dt(), y, rhs);

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

    if (!time_->is_current(output_mark_type)) return;


    // gather the solution from all processors
    IS is;
	VecScatter output_scatter;
	int row_ids[distr->size()];
	Vec solution_vec[n_substances];

	for (int i=0; i<n_substances; i++)
		VecCreateSeq(PETSC_COMM_SELF, ls->size(), &solution_vec[i]);

	for (int i=0; i<distr->size(); i++)
		row_ids[i] = i;

	ISCreateGeneral(PETSC_COMM_SELF, ls->size(), row_ids, PETSC_COPY_VALUES, &is); //WithArray
	VecScatterCreate(ls->get_solution(), is, solution_vec[0], PETSC_NULL, &output_scatter);
	for (int sbi = 0; sbi < n_substances; sbi++)
	{
		VecScatterBegin(output_scatter, ls->get_solution(), solution_vec[sbi], INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(output_scatter, ls->get_solution(), solution_vec[sbi], INSERT_VALUES, SCATTER_FORWARD);
	}
	VecScatterDestroy(&(output_scatter));
	ISDestroy(&(is));

    VecGetArray(solution_vec[0], &solution);

    // interpolate solution to a continuous field and write to output file

    if (distr->myp() == 0)
    {
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

		for (DOFHandler<1,3>::CellIterator cell = dof_handler1d->begin_cell(); cell != dof_handler1d->end_cell(); ++cell)
		{
			if (cell->dim != 1) continue;

			dof_handler1d->get_dof_indices(cell, dof_indices);
			output_cell_solution[0][cell.index()] = 0;

			for (int i=0; i<fe1d->n_dofs(); i++)
				output_cell_solution[0][cell.index()] += solution[dof_indices[i]]/fe1d->n_dofs();
		}
		for (DOFHandler<2,3>::CellIterator cell = dof_handler2d->begin_cell(); cell != dof_handler2d->end_cell(); ++cell)
		{
			if (cell->dim != 2) continue;

			dof_handler2d->get_dof_indices(cell, dof_indices);
			output_cell_solution[0][cell.index()] = 0;

			for (int i=0; i<fe2d->n_dofs(); i++)
				output_cell_solution[0][cell.index()] += solution[dof_indices[i]]/fe2d->n_dofs();
		}
		for (DOFHandler<3,3>::CellIterator cell = dof_handler3d->begin_cell(); cell != dof_handler3d->end_cell(); ++cell)
		{
			if (cell->dim != 3) continue;

			dof_handler3d->get_dof_indices(cell, dof_indices);
			output_cell_solution[0][cell.index()] = 0;

			for (int i=0; i<fe3d->n_dofs(); i++)
				output_cell_solution[0][cell.index()] += solution[dof_indices[i]]/fe3d->n_dofs();
		}

		transport_output->write_data(time_->t());
    }
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
        	local_rhs[i] = 0;
            for (int j=0; j<ndofs; j++)
            {
                local_mass_matrix[i*ndofs+j] = 0;
                if (dof_indices[i] < distr->begin() || dof_indices[i] > distr->end()) continue;
                for (int k=0; k<q.size(); k++)
                    local_mass_matrix[i*ndofs+j] += cell->material->size*fe_values.shape_value(j,k)*fe_values.shape_value(i,k)*fe_values.JxW(k);
            }
        }

        ls_dt->set_values(ndofs, (int *)dof_indices, ndofs, (int *)dof_indices, local_mass_matrix, local_rhs);
    }
}




void TransportDG::assemble_stiffness_matrix()
{
	assemble_volume_integrals<1>(dof_handler1d, fe1d);
	assemble_fluxes_boundary<1>(dof_handler1d, 0, fe1d, 0);
	assemble_fluxes_element_element<1>(dof_handler1d, 0, fe1d, 0);
	assemble_fluxes_element_side<1>(dof_handler1d, 0, fe1d, 0);

	assemble_volume_integrals<2>(dof_handler2d, fe2d);
	assemble_fluxes_boundary<2>(dof_handler2d, dof_handler1d, fe2d, fe1d);
	assemble_fluxes_element_element<2>(dof_handler2d, dof_handler1d, fe2d, fe1d);
	assemble_fluxes_element_side<2>(dof_handler2d, dof_handler1d, fe2d, fe1d);

	assemble_volume_integrals<3>(dof_handler3d, fe3d);
    assemble_fluxes_boundary<3>(dof_handler3d, dof_handler2d, fe3d, fe2d);
    assemble_fluxes_element_element<3>(dof_handler3d, dof_handler2d, fe3d, fe2d);
    assemble_fluxes_element_side<3>(dof_handler3d, dof_handler2d, fe3d, fe2d);
}



template<unsigned int dim>
void TransportDG::assemble_volume_integrals(DOFHandler<dim,3> *dh, FiniteElement<dim,3> *fe)
{
    MappingP1<dim,3> map;
    QGauss<dim> q(2);
    FE_RT0<dim,3> fe_rt;
    FEValues<dim,3> fv_rt(map, q, fe_rt, update_values | update_gradients);
    FEValues<dim,3> fe_values(map, q, *fe, update_values | update_gradients | update_JxW_values | update_quadrature_points);
    typename DOFHandler<dim,3>::CellIterator cell = dh->begin_cell();
    vector<mat33> K;
    vector<vec3> velocity;
    vector<double> divergence;
    const unsigned int ndofs = fe->n_dofs();
    unsigned int dof_indices[ndofs];
    PetscScalar local_matrix[ndofs*ndofs];

	// assemble integral over elements
    for (cell = dh->begin_cell(); cell != dh->end_cell(); ++cell)
    {
        if (cell->dim != dim) continue;

        fe_values.reinit(cell);
        fv_rt.reinit(cell);
        
        calculate_velocity(cell, velocity, fv_rt);
        calculate_dispersivity_tensor(K, velocity);
        calculate_velocity_divergence(cell, divergence, fv_rt);

        dh->get_dof_indices(cell, dof_indices);

        // assemble the local stiffness matrix
        for (int i=0; i<ndofs; i++)
        {
            for (int j=0; j<ndofs; j++)
            {
                local_matrix[i*ndofs+j] = 0;
                if (dof_indices[i] < distr->begin() || dof_indices[i] > distr->end()) continue;
                for (int k=0; k<q.size(); k++)
                {
                    local_matrix[i*ndofs+j] += (dot(K[k]*fe_values.shape_grad(j,k),fe_values.shape_grad(i,k))
                                               +dot(fe_values.shape_grad(j,k),velocity[k])*fe_values.shape_value(i,k)*advection/cell->material->por_m
//                                               +divergence[k]*fe_values.shape_value(j,k)*fe_values.shape_value(i,k)*advection
                                               )*fe_values.JxW(k);
                }

            }
        }

        ls->mat_set_values(ndofs, (int *)dof_indices, ndofs, (int *)dof_indices, local_matrix);
    }
}




template<unsigned int dim>
void TransportDG::assemble_fluxes_element_element(DOFHandler<dim,3> *dh, DOFHandler<dim-1,3> *dh_sub, FiniteElement<dim,3> *fe, FiniteElement<dim-1,3> *fe_sub)
{
    MappingP1<dim,3> map;
    QGauss<dim-1> side_q(2);
    FE_RT0<dim,3> fe_rt;
    vector<FESideValues<dim,3>*> fe_values_side;
    FESideValues<dim,3> fsv_rt(map, side_q, fe_rt, update_values | update_normal_vectors | update_side_JxW_values);
    vector<FEValuesSpaceBase<3>*> fv_sb;
    typename DOFHandler<dim,3>::CellIterator cell = dh->begin_cell();
    const unsigned int ndofs = fe->n_dofs();
    vector<unsigned int*> side_dof_indices;
    PetscScalar local_matrix[ndofs*ndofs], local_rhs[ndofs];
    vector< vector<mat33> > side_K;
    vector< vector<vec3> > side_velocity;
    double gamma_l, omega[2], transport_flux;

    gamma.resize(mesh_->boundary.size());

    // assemble integral over sides
    for(Neighbour *nb = mesh_->neighbour; nb != NULL; nb = nb->next)
    {
        if (nb->n_sides < 2 ||
        	(nb->type != BB_E && nb->type != BB_EL) ||
        	nb->element[0]->dim != dim) continue;

        fv_sb.resize(nb->n_sides);
        side_K.resize(nb->n_sides);
        side_velocity.resize(nb->n_sides);

        if (side_dof_indices.size() < nb->n_sides)
            for (int i=side_dof_indices.size(); i<nb->n_sides; i++)
                side_dof_indices.push_back(new unsigned int[ndofs]);

        if (fe_values_side.size() < nb->n_sides)
            for (int sid=fe_values_side.size(); sid<nb->n_sides; sid++)
                fe_values_side.push_back(new FESideValues<dim,3>(map, side_q, *fe, update_values | update_gradients | update_side_JxW_values | update_normal_vectors));

		for (int sid=0; sid<nb->n_sides; sid++)
		{
			cell = mesh_->element.full_iter(nb->element[sid]);
			dh->get_dof_indices(cell, side_dof_indices[sid]);
			fe_values_side[sid]->reinit(cell, nb->side[0]);
			fsv_rt.reinit(cell, nb->side[0]);
			calculate_velocity(cell, side_velocity[sid], fsv_rt);
			calculate_dispersivity_tensor(side_K[sid], side_velocity[sid]);
			fv_sb[sid] = fe_values_side[sid];
		}


        // fluxes and penalty
		for (int s1=0; s1<nb->n_sides; s1++)
		{
			for (int s2=s1+1; s2<nb->n_sides; s2++)
			{
				vec3 nv = (nb->side[s1]==NULL)?(-fv_sb[s2]->normal_vector(0)):fv_sb[s1]->normal_vector(0);

				// set up the parameters for DG method
				set_DG_parameters(nb, s1, s2, side_q.size(), side_K, -fv_sb[1]->normal_vector(0), alpha, advection, gamma_l, omega, transport_flux);

				if (nb->side[s1] != NULL && nb->side[s1]->cond != 0) gamma[mesh_->boundary.full_iter(nb->side[s1]->cond).id()] = gamma_l;

				int sd[2];
				sd[0] = s1;
				sd[1] = s2;

				for (int m=0; m<2; m++)
				{
					for (int n=0; n<2; n++)
					{
						for (int i=0; i<fv_sb[sd[n]]->n_dofs(); i++)
						{
							for (int j=0; j<fv_sb[sd[m]]->n_dofs(); j++)
							{
								int index = i*fv_sb[sd[m]]->n_dofs()+j;
								local_matrix[index] = 0;
								if (side_dof_indices[sd[n]][i] < distr->begin() || side_dof_indices[sd[n]][i] > distr->end()) continue;
								for (int k=0; k<side_q.size(); k++)
								{
									// flux due to transport (applied on interior edges)
									local_matrix[index] -= (m==0?1:-1)*advection/nb->element[sd[m]]->material->por_m*0.5*transport_flux*fv_sb[sd[m]]->shape_value(j,k)*fv_sb[sd[n]]->shape_value(i,k)*fv_sb[0]->JxW(k);

									// penalty enforcing continuity across edges (applied on interior and Dirichlet edges)
									local_matrix[index] += (m==n?1:-1)*gamma_l*fv_sb[sd[m]]->shape_value(j,k)*fv_sb[sd[n]]->shape_value(i,k)*fv_sb[0]->JxW(k);

									// terms due to diffusion
									local_matrix[index] += ((n==0)?-1.:1.)*dot(omega[0]*side_K[sd[m]][k]*fv_sb[sd[m]]->shape_grad(j,k),nv)*fv_sb[sd[n]]->shape_value(i,k)*fv_sb[0]->JxW(k);
									local_matrix[index] += ((m==0)?-1.:1.)*dot(omega[1]*side_K[sd[n]][k]*fv_sb[sd[n]]->shape_grad(i,k),nv)*fv_sb[sd[m]]->shape_value(j,k)*fv_sb[0]->JxW(k);
								}
							}
						}
						ls->mat_set_values(fv_sb[sd[n]]->n_dofs(), (int *)side_dof_indices[sd[n]], fv_sb[sd[m]]->n_dofs(), (int *)side_dof_indices[sd[m]], local_matrix);
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
void TransportDG::assemble_fluxes_boundary(DOFHandler<dim,3> *dh, DOFHandler<dim-1,3> *dh_sub, FiniteElement<dim,3> *fe, FiniteElement<dim-1,3> *fe_sub)
{
    MappingP1<dim,3> map;
    MappingP1<dim-1,3> map_vb;
    QGauss<dim-1> side_q(2);
    FE_RT0<dim,3> fe_rt;
    FESideValues<dim,3> fe_values_side(map, side_q, *fe, update_values | update_gradients | update_side_JxW_values | update_normal_vectors);
    FESideValues<dim,3> fsv_rt(map, side_q, fe_rt, update_values | update_normal_vectors | update_side_JxW_values);
    typename DOFHandler<dim,3>::CellIterator cell = dh->begin_cell();
    const unsigned int ndofs = fe->n_dofs();
    unsigned int side_dof_indices[ndofs];
    PetscScalar local_matrix[ndofs*ndofs], local_rhs[ndofs];
    vector< vector<mat33> > side_K;
    vector< vector<vec3> > side_velocity;
    double gamma_l, omega[2], transport_flux;

    gamma.resize(mesh_->boundary.size());

    // assemble boundary integral
    for (EdgeFullIter edge = mesh_->edge.begin(); edge != mesh_->edge.end(); ++edge)
    {
        if (edge->n_sides != 1 || edge->side[0]->dim != dim-1) continue;

        side_K.resize(edge->n_sides);
        side_velocity.resize(edge->n_sides);

        cell = mesh().element.full_iter(edge->side[0]->element);
        dh->get_dof_indices(cell, side_dof_indices);
        fe_values_side.reinit(cell, edge->side[0]);
        fsv_rt.reinit(cell, edge->side[0]);

        calculate_velocity(cell, side_velocity[0], fsv_rt);
        calculate_dispersivity_tensor(side_K[0], side_velocity[0]);

        // set up the parameters for DG method
        set_DG_parameters_edge(edge, side_q.size(), side_K, fe_values_side.normal_vector(0), alpha, advection, gamma_l, omega);
        if (edge->side[0]->cond != 0) gamma[mesh_->boundary.full_iter(edge->side[0]->cond).id()] = gamma_l;

        // fluxes and penalty
        for (int i=0; i<ndofs; i++)
        {
            for (int j=0; j<ndofs; j++)
            {
                local_matrix[i*ndofs+j] = 0;
                if (side_dof_indices[i] < distr->begin() || side_dof_indices[i] > distr->end()) continue;
                for (int k=0; k<side_q.size(); k++)
                {
                    // penalty enforcing continuity across edges (applied on interior and Dirichlet edges)
                    if (edge->side[0]->cond!=0 && edge->side[0]->flux < -tol_flux_bc)
                        local_matrix[i*ndofs+j] += gamma_l*fe_values_side.shape_value(j,k)*fe_values_side.shape_value(i,k)*fe_values_side.JxW(k);
                }
            }
        }
        ls->mat_set_values(ndofs, (int *)side_dof_indices, ndofs, (int *)side_dof_indices, local_matrix);
    }
}



template<unsigned int dim>
void TransportDG::assemble_fluxes_element_side(DOFHandler<dim,3> *dh, DOFHandler<dim-1,3> *dh_sub, FiniteElement<dim,3> *fe, FiniteElement<dim-1,3> *fe_sub)
{
    MappingP1<dim,3> map;
    MappingP1<dim-1,3> map_vb;
    QGauss<dim-1> side_q(2);
    FEValues<dim-1,3> *fe_values_vb;
    vector<FESideValues<dim,3>*> fe_values_side;
    vector<FEValuesSpaceBase<3>*> fv_sb;
    typename DOFHandler<dim,3>::CellIterator cell = dh->begin_cell();
    const unsigned int ndofs = fe->n_dofs();
    vector<unsigned int*> side_dof_indices;
    PetscScalar local_matrix[ndofs*ndofs], local_rhs[ndofs];
    double transport_flux;

    if (dim > 1)
        fe_values_vb = new FEValues<dim-1,3>(map_vb, side_q, *fe_sub, update_values | update_gradients | update_JxW_values | update_quadrature_points);

    // assemble integral over sides
    for(Neighbour *nb = mesh_->neighbour; nb != NULL; nb = nb->next)
    {
        if (nb->n_sides < 2 ||
        	nb->type != VB_ES ||
        	nb->element[0]->dim != dim-1) continue;

        fv_sb.resize(nb->n_sides);

        if (side_dof_indices.size() < nb->n_sides)
            for (int i=side_dof_indices.size(); i<nb->n_sides; i++)
                side_dof_indices.push_back(new unsigned int[ndofs]);

        if (fe_values_side.size() < nb->n_sides)
            for (int sid=fe_values_side.size(); sid<nb->n_sides; sid++)
                fe_values_side.push_back(new FESideValues<dim,3>(map, side_q, *fe, update_values | update_gradients | update_side_JxW_values | update_normal_vectors));

		typename DOFHandler<dim-1,3>::CellIterator cell_sub = mesh_->element.full_iter(nb->element[0]);
		dh_sub->get_dof_indices(cell_sub, side_dof_indices[0]);
		fe_values_vb->reinit(cell_sub);
		fv_sb[0] = fe_values_vb;

		cell = mesh_->element.full_iter(nb->element[1]);
		dh->get_dof_indices(cell, side_dof_indices[1]);
		fe_values_side[1]->reinit(cell, nb->side[1]);
		fv_sb[1] = fe_values_side[1];


		// flux from the higher dimension to the lower one
		transport_flux = nb->side[1]->flux;

		// set transmission condition for dim-1
		for (int j=0; j<fv_sb[0]->n_dofs(); j++)
		{
			for (int i=0; i<fv_sb[0]->n_dofs(); i++)
			{
				int index = i*fv_sb[0]->n_dofs() + j;
				local_matrix[index] = 0;
				if (side_dof_indices[0][i] < distr->begin() || side_dof_indices[0][i] > distr->end()) continue;
				for (int k=0; k<side_q.size(); k++)
					local_matrix[index] += advection/nb->element[0]->material->por_m*fabs(transport_flux)*fv_sb[0]->shape_value(j,k)*fv_sb[0]->shape_value(i,k)*fv_sb[0]->JxW(k);
			}
		}
		ls->mat_set_values(fv_sb[0]->n_dofs(), (int *)side_dof_indices[0], fv_sb[0]->n_dofs(), (int *)side_dof_indices[0], local_matrix);
		for (int j=0; j<fv_sb[1]->n_dofs(); j++)
		{
			for (int i=0; i<fv_sb[0]->n_dofs(); i++)
			{
				int index = i*fv_sb[1]->n_dofs() + j;
				local_matrix[index] = 0;
				if (side_dof_indices[0][i] < distr->begin() || side_dof_indices[0][i] > distr->end()) continue;
				for (int k=0; k<side_q.size(); k++)
					local_matrix[index] -= advection/nb->element[1]->material->por_m*fabs(transport_flux)*fv_sb[1]->shape_value(j,k)*fv_sb[0]->shape_value(i,k)*fv_sb[0]->JxW(k);
			}
		}
		ls->mat_set_values(fv_sb[0]->n_dofs(), (int *)side_dof_indices[0], fv_sb[1]->n_dofs(), (int *)side_dof_indices[1], local_matrix);

		// set transmission condition for dim
		for (int j=0; j<fv_sb[0]->n_dofs(); j++)
		{
			for (int i=0; i<fv_sb[1]->n_dofs(); i++)
			{
				int index = i*fv_sb[0]->n_dofs() + j;
				local_matrix[index] = 0;
				if (side_dof_indices[1][i] < distr->begin() || side_dof_indices[1][i] > distr->end()) continue;
				for (int k=0; k<side_q.size(); k++)
					local_matrix[index] += advection/nb->element[0]->material->por_m*min(0.,transport_flux)*fv_sb[0]->shape_value(j,k)*fv_sb[1]->shape_value(i,k)*fv_sb[0]->JxW(k);
			}
		}
		ls->mat_set_values(fv_sb[1]->n_dofs(), (int *)side_dof_indices[1], fv_sb[0]->n_dofs(), (int *)side_dof_indices[0], local_matrix);
		for (int j=0; j<fv_sb[1]->n_dofs(); j++)
		{
			for (int i=0; i<fv_sb[1]->n_dofs(); i++)
			{
				int index = i*fv_sb[1]->n_dofs() + j;
				local_matrix[index] = 0;
				if (side_dof_indices[1][i] < distr->begin() || side_dof_indices[1][i] > distr->end()) continue;
				for (int k=0; k<side_q.size(); k++)
					local_matrix[index] -= advection/nb->element[1]->material->por_m*min(0.,transport_flux)*fv_sb[1]->shape_value(j,k)*fv_sb[1]->shape_value(i,k)*fv_sb[0]->JxW(k);
			}
		}
		ls->mat_set_values(fv_sb[1]->n_dofs(), (int *)side_dof_indices[1], fv_sb[1]->n_dofs(), (int *)side_dof_indices[1], local_matrix);
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
        if (b.id() < bc->distribution()->begin() || b.id() > bc->distribution()->end()) continue;

        fe_values_side.reinit(cell, b->side);
        dh->get_dof_indices(cell, side_dof_indices);

        for (int i=0; i<fe->n_dofs(); i++)
        {
            local_rhs[i] = 0;
            for (int k=0; k<side_q.size(); k++)
            {
                local_rhs[i] += (
                        +gamma[b.id()]*bc->get_array(0)[b.id()-bc->distribution()->begin()]
//                       -advection*0.5*min(b->side->flux,-tol_flux_bc)*bc_values[0][b.id()]
                                )*fe_values_side.shape_value(i,k)*fe_values_side.JxW(k);
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


template<unsigned int dim>
void TransportDG::calculate_velocity_divergence(typename DOFHandler<dim,3>::CellIterator cell, vector<double> &divergence, FEValuesBase<dim,3> &fv)
{
    std::map<Node*, int> node_nums;
    for (int i=0; i<cell->n_nodes; i++)
        node_nums[cell->node[i]] = i;

    divergence.resize(fv.n_points());

    for (int k=0; k<fv.n_points(); k++)
    {
        divergence[k] = 0;
        for (int sid=0; sid<cell->n_sides; sid++)
        {
            if (cell->side[sid]->dim != dim-1) continue;
            int num = dim*(dim+1)/2;
            for (int i=0; i<cell->side[sid]->n_nodes; i++)
                num -= node_nums[cell->side[sid]->node[i]];
            divergence[k] += trace(fv.shape_grad_vector(num,k))*cell->side[sid]->flux;
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

        K[k] += eye(3,3)*(Dm*tortuosity + alphaT*vnorm);
    }
}


void TransportDG::set_DG_parameters(const Neighbour *n,
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
    double fluxes[n->n_sides];
    Side *s = (n->side[s1]==NULL)?n->side[s2]:n->side[s1];

    // calculate the side diameter
    if (s->dim == 0)
    {
        h = 1;
    }
    else
    {
        for (int i=0; i<s->n_nodes; i++)
            for (int j=i+1; j<s->n_nodes; j++)
                h = max(h, s->node[i]->distance(*s->node[j]));
    }

    // calculate the total in- and out-flux through the edge
    if (n->side[0] == NULL)
    {
        fluxes[0] = -n->side[1]->flux;
    }
    else
    {
        fluxes[0] = n->side[0]->flux;
    }
    for (int i=1; i<n->n_sides; i++) fluxes[i] = n->side[i]->flux;
    double pflux = 0, nflux = 0;
    for (int i=0; i<n->n_sides; i++)
    {
        if (fluxes[i] > 0)
            pflux += fluxes[i];
        else
            nflux += fluxes[i];
    }

    // calculate the flux from s1 to s2
    if (fluxes[s2] > 0 && fluxes[s1] < 0 && s1 < s2)
    {
        transport_flux = fluxes[s1]*fabs(fluxes[s2]/pflux);
    }
    else if (fluxes[s2] < 0 && fluxes[s1] > 0 && s1 < s2)
    {
        transport_flux = fluxes[s1]*fabs(fluxes[s2]/nflux);
    }
    else if (s1==s2)
    {
        transport_flux = fluxes[s1];
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







void TransportDG::set_DG_parameters_edge(const Edge *edge,
            const unsigned int n_points,
            const vector< vector<mat33> > &K,
            const vec3 &normal_vector,
            const double alpha,
            const double advection,
            double &gamma,
            double *omega)
{
    double delta[2];
    double h = 0;
    Side *side = edge->side[0];

    // calculate the side diameter
    if (side->dim == 0)
    {
        h = 1;
    }
    else
    {
        for (int i=0; i<side->n_nodes; i++)
            for (int j=i+1; j<side->n_nodes; j++)
                h = max(h, side->node[i]->distance(*side->node[j]));
    }

    gamma = 0.5*advection*fabs(side->flux);

	omega[0] = 1;

	// delta is set to the average value of Kn.n on the side
	delta[0] = 0;
	for (int k=0; k<n_points; k++)
		delta[0] += dot(K[0][k]*normal_vector,normal_vector);
	delta[0] /= n_points;

	gamma += alpha/h*delta[0];
}










void TransportDG::read_initial_condition()
{
    FILE    *in;          // input file
    char     line[ LINE_SIZE ]; // line of data file
    int sbi,index, eid, i,n_concentrations;
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
        eid   = atoi( xstrtok( line ) );
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


