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
#include "system/xio.h"
#include "system/sys_profiler.hh"
#include "transport/transport_dg.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/mapping_p1.hh"
#include "la/solve.h"
#include "fem/fe_rt.hh"
#include "io/output.h"
#include "mesh/boundaries.h"
#include "la/distribution.hh"
#include "input/accessors.hh"
#include "flow/old_bcd.hh"


using namespace Input::Type;

Selection TransportDG::dg_variant_selection_input_type
	= Selection("DG_variant")
	.add_value(non_symmetric, "non-symmetric", "non-symmetric weighted interior penalty DG method")
	.add_value(incomplete,    "incomplete",    "incomplete weighted interior penalty DG method")
	.add_value(symmetric,     "symmetric",     "symmetric weighted interior penalty DG method");

Record TransportDG::input_type
	= Record("AdvectionDiffusion_DG", "DG solver for transport with diffusion.")
	.derive_from(TransportBase::input_type)
    .declare_key("solver", Solver::input_type, Default::obligatory(),
            "Linear solver for MH problem.")
    .declare_key("bc_data", Array(TransportDG::EqData().boundary_input_type()
    		.declare_key("old_boundary_file", IT::FileName::input(),
    				"Input file with boundary conditions (obsolete).")
    		.declare_key("bc_times", Array(Double()), Default::optional(),
    				"Times for changing the boundary conditions (obsolete).")
    		), IT::Default::obligatory(), "")
    .declare_key("bulk_data", Array(TransportDG::EqData().bulk_input_type()),
    		IT::Default::obligatory(), "")
    .declare_key("dg_variant", TransportDG::dg_variant_selection_input_type, Default("non-symmetric"),
    		"Variant of interior penalty discontinuous Galerkin method.")
    .declare_key("dg_order", Integer(0,2), Default("1"),
    		"Polynomial order for finite element in DG method (order 0 is suitable if there is no diffusion/dispersion).");




TransportDG::FEObjects::FEObjects(Mesh *mesh_, unsigned int fe_order)
{
	unsigned int q_order;

	switch (fe_order)
	{
	case 0:
		q_order = 0;
		fe1_ = new FE_P_disc<0,1,3>;
		fe2_ = new FE_P_disc<0,2,3>;
		fe3_ = new FE_P_disc<0,3,3>;
		break;

	case 1:
		q_order = 2;
		fe1_ = new FE_P_disc<1,1,3>;
		fe2_ = new FE_P_disc<1,2,3>;
		fe3_ = new FE_P_disc<1,3,3>;
		break;

	case 2:
		q_order = 4;
		fe1_ = new FE_P_disc<2,1,3>;
		fe2_ = new FE_P_disc<2,2,3>;
		fe3_ = new FE_P_disc<2,3,3>;
		break;

	default:
		xprintf(PrgErr, "Unsupported polynomial order %d for finite elements in TransportDG ", fe_order);
		break;
	}

	fe_rt1_ = new FE_RT0<1,3>;
	fe_rt2_ = new FE_RT0<2,3>;
	fe_rt3_ = new FE_RT0<3,3>;

	q0_ = new QGauss<0>(q_order);
	q1_ = new QGauss<1>(q_order);
	q2_ = new QGauss<2>(q_order);
	q3_ = new QGauss<3>(q_order);

	map0_ = new MappingP1<0,3>;
	map1_ = new MappingP1<1,3>;
	map2_ = new MappingP1<2,3>;
	map3_ = new MappingP1<3,3>;

	dh1_ = new DOFHandler<1,3>(*mesh_);
	dh2_ = new DOFHandler<2,3>(*mesh_);
	dh3_ = new DOFHandler<3,3>(*mesh_);

	dh1_->distribute_dofs(*fe1_);
	dh2_->distribute_dofs(*fe2_, dh1_->n_global_dofs());
	dh3_->distribute_dofs(*fe3_, dh1_->n_global_dofs() + dh2_->n_global_dofs());
}

TransportDG::FEObjects::~FEObjects()
{
	delete fe1_;
	delete fe2_;
	delete fe3_;
	delete fe_rt1_;
	delete fe_rt2_;
	delete fe_rt3_;
	delete q0_;
	delete q1_;
	delete q2_;
	delete q3_;
	delete map0_;
	delete map1_;
	delete map2_;
	delete map3_;
	delete dh1_;
	delete dh2_;
	delete dh3_;
}

template<> FiniteElement<0,3> *TransportDG::FEObjects::fe<0>() { return 0; }
template<> FiniteElement<1,3> *TransportDG::FEObjects::fe<1>() { return fe1_; }
template<> FiniteElement<2,3> *TransportDG::FEObjects::fe<2>() { return fe2_; }
template<> FiniteElement<3,3> *TransportDG::FEObjects::fe<3>() { return fe3_; }

template<> FiniteElement<0,3> *TransportDG::FEObjects::fe_rt<0>() { return 0; }
template<> FiniteElement<1,3> *TransportDG::FEObjects::fe_rt<1>() { return fe_rt1_; }
template<> FiniteElement<2,3> *TransportDG::FEObjects::fe_rt<2>() { return fe_rt2_; }
template<> FiniteElement<3,3> *TransportDG::FEObjects::fe_rt<3>() { return fe_rt3_; }

template<> Quadrature<0> *TransportDG::FEObjects::q<0>() { return q0_; }
template<> Quadrature<1> *TransportDG::FEObjects::q<1>() { return q1_; }
template<> Quadrature<2> *TransportDG::FEObjects::q<2>() { return q2_; }
template<> Quadrature<3> *TransportDG::FEObjects::q<3>() { return q3_; }

template<> Mapping<0,3> *TransportDG::FEObjects::map<0>() { return map0_; }
template<> Mapping<1,3> *TransportDG::FEObjects::map<1>() { return map1_; }
template<> Mapping<2,3> *TransportDG::FEObjects::map<2>() { return map2_; }
template<> Mapping<3,3> *TransportDG::FEObjects::map<3>() { return map3_; }

template<> DOFHandler<0,3> *TransportDG::FEObjects::dh<0>() { return 0; }
template<> DOFHandler<1,3> *TransportDG::FEObjects::dh<1>() { return dh1_; }
template<> DOFHandler<2,3> *TransportDG::FEObjects::dh<2>() { return dh2_; }
template<> DOFHandler<3,3> *TransportDG::FEObjects::dh<3>() { return dh3_; }



TransportDG::EqData::EqData() : TransportEqData("TransportDG")
{
	ADD_FIELD(disp_l, "Longitudal dispersivity (for each substance).", Default("0"));
	ADD_FIELD(disp_t, "Transversal dispersivity (for each substance).", Default("0"));
	ADD_FIELD(diff_m, "Molecular diffusivity (for each substance).", Default("0"));
	ADD_FIELD(sigma_c, "Coefficient of diffusive transfer through fractures (for each substance).", Default("0"));
	ADD_FIELD(dg_penalty, "Penalty parameter influencing the discontinuity of the solution (for each substance). "
			"Its default value 1 is sufficient in most cases. Higher value diminishes the inter-element jumps.", Default("1.0"));
}

RegionSet TransportDG::EqData::read_boundary_list_item(Input::Record rec) {
	// Base method EqDataBase::read_boundary_list_item must be called first!
	RegionSet domain = EqDataBase::read_boundary_list_item(rec);
    FilePath bcd_file;

    // read transport boundary conditions using old file format .tbc
    if (rec.opt_val("old_boundary_file", bcd_file) )
        OldBcdInput::instance()->read_transport(bcd_file, bc_conc);

    return domain;
}

TransportDG::TransportDG(Mesh & init_mesh, const Input::Record &in_rec)
        : TransportBase(init_mesh, in_rec),
          mass_matrix(0),
          flux_changed(true),
          allocation_done(false)
{
    time_ = new TimeGovernor(in_rec.val<Input::Record>("time"), equation_mark_type_);
    time_->fix_dt_until_mark();
    
    // set up solver
    solver = new Solver;
    solver_init(solver, in_rec.val<Input::AbstractRecord>("solver"));

    // Read names of transported substances.
    in_rec.val<Input::Array>("substances").copy_to(subst_names);
    n_subst = subst_names.size();

    // Set up physical parameters.
    data.set_mesh(&init_mesh);
    data.init_conc.set_n_comp(n_subst);
    data.bc_conc.set_n_comp(n_subst);
    data.diff_m.set_n_comp(n_subst);
    data.disp_l.set_n_comp(n_subst);
    data.disp_t.set_n_comp(n_subst);
    data.sigma_c.set_n_comp(n_subst);
    data.dg_penalty.set_n_comp(n_subst);
    data.sources_density.set_n_comp(n_subst);
    data.sources_sigma.set_n_comp(n_subst);
    data.sources_conc.set_n_comp(n_subst);
    data.init_from_input( in_rec.val<Input::Array>("bulk_data"), in_rec.val<Input::Array>("bc_data") );
    data.set_time(*time_);

    // sorption and dual_porosity is currently not used in TransportDG
    sorption = in_rec.val<bool>("sorption_enable");
    dual_porosity = in_rec.val<bool>("dual_porosity");

    // DG variant
    dg_variant = in_rec.val<DGVariant>("dg_variant");

    // DG stabilization parameters on boundary edges
    gamma.resize(n_subst);
    for (unsigned int sbi=0; sbi<n_subst; sbi++)
    	gamma[sbi].resize(mesh_->boundary_.size());


    // create finite element structures and distribute DOFs
    dg_order = in_rec.val<unsigned int>("dg_order");
    feo = new FEObjects(mesh_, dg_order);

    // distribute solution vectors on processors
    distr = new Distribution(Distribution::Block, feo->dh<1>()->n_global_dofs() + feo->dh<2>()->n_global_dofs() + feo->dh<3>()->n_global_dofs());


    // set up output class
    Input::Record output_rec = in_rec.val<Input::Record>("output");
    transport_output = OutputTime::output_stream(output_rec.val<Input::Record>("output_stream"));

    // allocate output arrays
    if (distr->myp() == 0)
    {
    	int n_corners = 0;
    	FOR_ELEMENTS(mesh_, elem)
    		n_corners += elem->dim()+1;
    	output_solution.resize(n_subst);
    	for (int i=0; i<n_subst; i++)
    	{
			output_solution[i] = new double[n_corners];
			for(int j=0; j<n_corners; j++)
				output_solution[i][j] = 0.0;
			OutputTime::register_corner_data<double>(mesh_, subst_names[i], "M/L^3",
					output_rec.val<Input::Record>("output_stream"), output_solution[i], n_corners);
		}
    }

    // set time marks for writing the output
    output_mark_type = this->mark_type() | time_->marks().type_fixed_time() | time_->marks().type_output();
    time_->marks().add_time_marks(0.0, output_rec.val<double>("save_step"), time_->end_time(), output_mark_type);
    
    // allocate matrix and vector structures
    ls    = new LinSys*[n_subst];
    ls_dt = new LinSys_MPIAIJ(distr->lsize());
    for (unsigned int sbi = 0; sbi < n_subst; sbi++)
    	ls[sbi] = new LinSys_MPIAIJ(distr->lsize());
    stiffness_matrix = new Mat[n_subst];
    rhs = new Vec[n_subst];


    // set initial conditions
    set_initial_condition();

    // save initial state
    output_data();
}


TransportDG::~TransportDG()
{
    //delete transport_output;
    delete time_;
    delete solver;
    delete ls_dt;

    if (distr->myp() == 0)
    {
		for (int i=0; i<n_subst; i++)
			delete[] output_solution[i];
    }

    for (int i=0; i<n_subst; i++) delete ls[i];
    delete[] ls;
    delete[] stiffness_matrix;
    delete[] rhs;
    delete feo;

    gamma.clear();
    subst_names.clear();
}


void TransportDG::set_eq_data(Field< 3, FieldValue<3>::Scalar > *cross_section)
{
  data.cross_section = cross_section;
}



void TransportDG::update_solution()
{
	START_TIMER("DG-ONE STEP");

	// calculate mass balance at initial time
	if (!allocation_done)
	{
		mass_balance();
	}

    time_->next_time();
    time_->view("TDG");
    
    START_TIMER("data reinit");
    data.set_time(*time_);
    END_TIMER("data reinit");
    
    // check first time assembly - needs preallocation
    if (!allocation_done)
    {
    	// preallocate mass matrix
    	ls_dt->start_allocation();
    	assemble_mass_matrix();
    	mass_matrix = NULL;

		// preallocate system matrix
		for (unsigned int i=0; i<n_subst; i++)
		{
			ls[i]->start_allocation();
			stiffness_matrix[i] = NULL;
			rhs[i] = NULL;
		}
		assemble_stiffness_matrix();
		set_sources();
		set_boundary_conditions();

		allocation_done = true;
    }

	// assemble mass matrix
	if (mass_matrix == NULL ||
		data.cross_section->changed() ||
		data.por_m.changed())
	{
		ls_dt->start_add_assembly();
		assemble_mass_matrix();
		ls_dt->finalize();
		mass_matrix = ls_dt->get_matrix();
	}

	// assemble stiffness matrix
    if (stiffness_matrix[0] == NULL ||
    	flux_changed ||
    	data.disp_l.changed() ||
    	data.disp_t.changed() ||
    	data.diff_m.changed() ||
    	data.sigma_c.changed() ||
    	data.dg_penalty.changed() ||
    	data.por_m.changed() ||
    	data.cross_section->changed())
    {
        // new fluxes can change the location of Neumann boundary,
        // thus stiffness matrix must be reassembled
    	for (unsigned int i=0; i<n_subst; i++)
    	{
    		ls[i]->start_add_assembly();
    		MatZeroEntries(ls[i]->get_matrix());
    	}
        assemble_stiffness_matrix();
        for (unsigned int i=0; i<n_subst; i++)
        {
        	ls[i]->finalize();

        	if (stiffness_matrix[i] == NULL)
        		MatConvert(ls[i]->get_matrix(), MATSAME, MAT_INITIAL_MATRIX, &stiffness_matrix[i]);
        	else
        		MatCopy(ls[i]->get_matrix(), stiffness_matrix[i], DIFFERENT_NONZERO_PATTERN);
        }
    }

    // assemble right hand side (due to sources and boundary conditions)
//    if (rhs[0] == NULL ||
//    	flux_changed ||
//    	data.bc_conc.changed() ||
//    	data.dg_penalty.changed() ||
//    	data.sources_conc.changed() ||
//    	data.sources_density.changed() ||
//    	data.sources_sigma.changed())
//    {
    	for (unsigned int i=0; i<n_subst; i++)
    	{
    		ls[i]->start_add_assembly();
    		VecSet(ls[i]->get_rhs(), 0);
    	}
    	set_sources();
    	set_boundary_conditions();
    	for (unsigned int i=0; i<n_subst; i++)
    	{
    		ls[i]->finalize();

    		VecDuplicate(ls[i]->get_rhs(), &rhs[i]);
    		VecCopy(ls[i]->get_rhs(), rhs[i]);
    	}
//    }

    flux_changed = false;


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

    for (unsigned int i=0; i<n_subst; i++)
    {
		MatCopy(stiffness_matrix[i], ls[i]->get_matrix(), DIFFERENT_NONZERO_PATTERN);
		// ls->get_matrix() = 1/dt*mass_matrix + ls->get_matrix()
		MatAXPY(ls[i]->get_matrix(), 1./time_->dt(), mass_matrix, SUBSET_NONZERO_PATTERN);
		Vec y;
		VecDuplicate(rhs[i], &y);
		// y = mass_matrix*ls->get_solution()
		MatMult(mass_matrix, ls[i]->get_solution(), y);
		// ls->get_rhs() = 1/dt*y + rhs
		VecWAXPY(ls[i]->get_rhs(), 1./time_->dt(), y, rhs[i]);

		//MatView( ls->get_matrix(), PETSC_VIEWER_STDOUT_SELF );

		VecDestroy(&y);

		//VecView( ls->get_rhs(), PETSC_VIEWER_STDOUT_SELF );
		// solve
		solve_system(solver, ls[i]);
    }

    mass_balance();

    //VecView( ls->get_solution(), PETSC_VIEWER_STDOUT_SELF );
    
    END_TIMER("DG-ONE STEP");
}




void TransportDG::get_solution_vector(double *& vector, unsigned int & size)
{}



void TransportDG::get_parallel_solution_vector(Vec & vector)
{}



void TransportDG::set_velocity_field(const MH_DofHandler &dh)
{
    // So far the velocity_vector contains zeros, so we ignore it.
    // Instead we use the value Side.flux.

    mh_dh = &dh;
	flux_changed = true;

}



void TransportDG::output_data()
{
    double *solution;
    unsigned int dof_indices[max(feo->fe<1>()->n_dofs(), max(feo->fe<2>()->n_dofs(), feo->fe<3>()->n_dofs()))];
    int n_nodes = mesh_->node_vector.size();
    int count[n_nodes];

    if (!time_->is_current(output_mark_type)) return;

    START_TIMER("DG-OUTPUT");

    // gather the solution from all processors
    IS is;
    VecScatter output_scatter;
    int row_ids[distr->size()];
	Vec solution_vec;

	for (int i=0; i<distr->size(); i++)
		row_ids[i] = i;

	for (int sbi=0; sbi<n_subst; sbi++)
	{
		VecCreateSeq(PETSC_COMM_SELF, ls[sbi]->size(), &solution_vec);

		ISCreateGeneral(PETSC_COMM_SELF, ls[sbi]->size(), row_ids, PETSC_COPY_VALUES, &is); //WithArray
		VecScatterCreate(ls[sbi]->get_solution(), is, solution_vec, PETSC_NULL, &output_scatter);
		VecScatterBegin(output_scatter, ls[sbi]->get_solution(), solution_vec, INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(output_scatter, ls[sbi]->get_solution(), solution_vec, INSERT_VALUES, SCATTER_FORWARD);
		VecScatterDestroy(&(output_scatter));
		ISDestroy(&(is));

		// on the main processor fill the output array and save to file
		if (distr->myp() == 0)
		{
			int corner_id = 0, node_id;

			// The solution is evaluated at vertices of each element.
			// For this reason we construct quadratures whose quadrature
			// points are placed at the vertices of the reference element.
			Quadrature<1> q1(2);
			Quadrature<2> q2(3);
			Quadrature<3> q3(4);

			for (int i=0; i<2; i++)
				q1.set_point(i, RefElement<1>::node_coords(i).subvec(0,0));
			for (int i=0; i<3; i++)
				q2.set_point(i, RefElement<2>::node_coords(i).subvec(0,1));
			for (int i=0; i<4; i++)
				q3.set_point(i, RefElement<3>::node_coords(i).subvec(0,2));

			FEValues<1,3> fv1(*feo->map<1>(), q1, *feo->fe<1>(), update_values);
			FEValues<2,3> fv2(*feo->map<2>(), q2, *feo->fe<2>(), update_values);
			FEValues<3,3> fv3(*feo->map<3>(), q3, *feo->fe<3>(), update_values);

			VecGetArray(solution_vec, &solution);
			FOR_ELEMENTS(mesh_, elem)
			{
				switch (elem->dim())
				{
				case 1:
					feo->dh<1>()->get_dof_indices(elem, dof_indices);
					fv1.reinit(elem);
					for (int k=0; k<q1.size(); k++)
					{
						output_solution[sbi][corner_id] = 0;
						for (int i=0; i<fv1.n_dofs(); i++)
							output_solution[sbi][corner_id] += solution[dof_indices[i]]*fv1.shape_value(i,k);
						corner_id++;
					}
					break;
				case 2:
					feo->dh<2>()->get_dof_indices(elem, dof_indices);
					fv2.reinit(elem);
					for (int k=0; k<q2.size(); k++)
					{
						output_solution[sbi][corner_id] = 0;
						for (int i=0; i<fv2.n_dofs(); i++)
							output_solution[sbi][corner_id] += solution[dof_indices[i]]*fv2.shape_value(i,k);
						corner_id++;
					}
					break;
				case 3:
					feo->dh<3>()->get_dof_indices(elem, dof_indices);
					fv3.reinit(elem);
					for (int k=0; k<q3.size(); k++)
					{
						output_solution[sbi][corner_id] = 0;
						for (int i=0; i<fv3.n_dofs(); i++)
							output_solution[sbi][corner_id] += solution[dof_indices[i]]*fv3.shape_value(i,k);
						corner_id++;
					}
					break;
				default:
					break;
				}

			}
		}

		VecDestroy(&solution_vec);
	}

	if (distr->myp() == 0)
	{
		if(transport_output) {
			xprintf(MsgLog, "transport DG: write_data()\n");
			transport_output->write_data(time_->t());
		}
	}

    END_TIMER("DG-OUTPUT");
}


void TransportDG::assemble_mass_matrix()
{
  START_TIMER("assemble_mass");
	assemble_mass_matrix<1>();
	assemble_mass_matrix<2>();
	assemble_mass_matrix<3>();
  END_TIMER("assemble_mass");
}


template<unsigned int dim>
void TransportDG::assemble_mass_matrix()
{
    FEValues<dim,3> fe_values(*feo->map<dim>(), *feo->q<dim>(), *feo->fe<dim>(), update_values | update_JxW_values | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
    unsigned int dof_indices[ndofs];
    PetscScalar local_mass_matrix[ndofs*ndofs], local_rhs[ndofs];
    vector<double> elem_csec(qsize), por_m(qsize);
    typename DOFHandler<dim,3>::CellIterator cell = feo->dh<dim>()->begin_cell();

    // assemble integral over elements
    for (cell = feo->dh<dim>()->begin_cell(); cell != feo->dh<dim>()->end_cell(); ++cell)
    {
        if (cell->dim() != dim) continue;

        fe_values.reinit(cell);
        feo->dh<dim>()->get_dof_indices(cell, dof_indices);
        ElementAccessor<3> ele_acc = cell->element_accessor();

        for (int k=0; k<qsize; k++)
        {
        	elem_csec[k] = data.cross_section->value(fe_values.point(k), ele_acc);
        	por_m[k]     = data.por_m.value(fe_values.point(k), ele_acc);
        }

        // assemble the local stiffness and mass matrix
        for (int i=0; i<ndofs; i++)
        {
        	local_rhs[i] = 0;
            for (int j=0; j<ndofs; j++)
            {
                local_mass_matrix[i*ndofs+j] = 0;
                if (dof_indices[i] < distr->begin() || dof_indices[i] > distr->end()) continue;
                for (int k=0; k<qsize; k++)
                    local_mass_matrix[i*ndofs+j] += elem_csec[k]*por_m[k]*fe_values.shape_value(j,k)*fe_values.shape_value(i,k)*fe_values.JxW(k);
            }
        }

        ls_dt->set_values(ndofs, (int *)dof_indices, ndofs, (int *)dof_indices, local_mass_matrix, local_rhs);
    }
}




void TransportDG::assemble_stiffness_matrix()
{
  START_TIMER("assemble_stiffness");
   START_TIMER("assemble_volume_integrals");
	assemble_volume_integrals<1>();
	assemble_volume_integrals<2>();
	assemble_volume_integrals<3>();
   END_TIMER("assemble_volume_integrals");

   START_TIMER("assemble_fluxes_boundary");
	assemble_fluxes_boundary<1>();
	assemble_fluxes_boundary<2>();
	assemble_fluxes_boundary<3>();
   END_TIMER("assemble_fluxes_boundary");

   START_TIMER("assemble_fluxes_elem_elem");
	assemble_fluxes_element_element<1>();
	assemble_fluxes_element_element<2>();
	assemble_fluxes_element_element<3>();
   END_TIMER("assemble_fluxes_elem_elem");

   START_TIMER("assemble_fluxes_elem_side");
	assemble_fluxes_element_side<1>();
	assemble_fluxes_element_side<2>();
    assemble_fluxes_element_side<3>();
   END_TIMER("assemble_fluxes_elem_side");
  END_TIMER("assemble_stiffness");
}




template<unsigned int dim>
void TransportDG::assemble_volume_integrals()
{
    FEValues<dim,3> fv_rt(*feo->map<dim>(), *feo->q<dim>(), *feo->fe_rt<dim>(),
    		update_values | update_gradients);
    FEValues<dim,3> fe_values(*feo->map<dim>(), *feo->q<dim>(), *feo->fe<dim>(),
    		update_values | update_gradients | update_JxW_values | update_quadrature_points);
    typename DOFHandler<dim,3>::CellIterator cell = feo->dh<dim>()->begin_cell();
    arma::mat33 K;
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
    unsigned int dof_indices[ndofs];
    vector<arma::vec3> velocity(qsize);
    vector<arma::vec> Dm(qsize), alphaL(qsize), alphaT(qsize);
    vector<double> divergence, por_m(qsize), csection(qsize);

    PetscScalar local_matrix[ndofs*ndofs];

	// assemble integral over elements
    for (cell = feo->dh<dim>()->begin_cell(); cell != feo->dh<dim>()->end_cell(); ++cell)
    {
        if (cell->dim() != dim) continue;

        fe_values.reinit(cell);
        fv_rt.reinit(cell);
        ElementAccessor<3> ele_acc = cell->element_accessor();
        feo->dh<dim>()->get_dof_indices(cell, dof_indices);
        
        for (int k=0; k<qsize; k++)
        {
        	Dm[k]       = data.diff_m.value(fe_values.point(k), ele_acc);
        	alphaL[k]   = data.disp_l.value(fe_values.point(k), ele_acc);
        	alphaT[k]   = data.disp_t.value(fe_values.point(k), ele_acc);
        	por_m[k]    = data.por_m.value(fe_values.point(k), ele_acc);
        	csection[k] = data.cross_section->value(fe_values.point(k), ele_acc);
        }

        calculate_velocity(cell, velocity, fv_rt);
        calculate_velocity_divergence(cell, divergence, fv_rt);

        // assemble the local stiffness matrix
        for (unsigned int sbi=0; sbi<n_subst; sbi++)
        {
        	for (int i=0; i<ndofs; i++)
        		for (int j=0; j<ndofs; j++)
        			local_matrix[i*ndofs+j] = 0;

        	for (int k=0; k<qsize; k++)
        	{
        		calculate_dispersivity_tensor(K, velocity[k], Dm[k][sbi], alphaL[k][sbi], alphaT[k][sbi], por_m[k], csection[k]);
        		double por_times_csection_times_JxW = por_m[k]*csection[k]*fe_values.JxW(k);

        		for (int i=0; i<ndofs; i++)
        		{
        			if (dof_indices[i] < distr->begin() || dof_indices[i] > distr->end()) continue;

        			arma::vec3 Kt_grad_i = K.t()*fe_values.shape_grad(i,k);
        			arma::vec3 velocity_times_value_i_times_JxW = velocity[k]*(fe_values.shape_value(i,k)*fe_values.JxW(k));
        			double div_times_value_i_times_JxW = divergence[k]*fe_values.shape_value(i,k)*fe_values.JxW(k);

        			for (int j=0; j<ndofs; j++)
        			{
						local_matrix[i*ndofs+j] += arma::dot(Kt_grad_i, fe_values.shape_grad(j,k))*por_times_csection_times_JxW
												   +dot(fe_values.shape_grad(j,k),velocity_times_value_i_times_JxW)
												   +div_times_value_i_times_JxW*fe_values.shape_value(j,k);
					}

				}
			}
			ls[sbi]->mat_set_values(ndofs, (int *)dof_indices, ndofs, (int *)dof_indices, local_matrix);
        }
    }
}


void TransportDG::set_sources()
{
  START_TIMER("assemble_sources");
	set_sources<1>();
	set_sources<2>();
	set_sources<3>();
  END_TIMER("assemble_sources");
}

template<unsigned int dim>
void TransportDG::set_sources()
{
    FEValues<dim,3> fe_values(*feo->map<dim>(), *feo->q<dim>(), *feo->fe<dim>(),
    		update_values | update_JxW_values | update_quadrature_points);
    typename DOFHandler<dim,3>::CellIterator cell = feo->dh<dim>()->begin_cell();
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
    vector<arma::vec> conc(qsize), sources_conc(qsize), sources_density(qsize), sources_sigma(qsize);
    unsigned int dof_indices[ndofs];
    PetscScalar local_rhs[ndofs];
    double conc_diff, source;

    for (unsigned int k=0; k<qsize; k++)
    	conc[k].resize(n_subst);

	// assemble integral over elements
    for (cell = feo->dh<dim>()->begin_cell(); cell != feo->dh<dim>()->end_cell(); ++cell)
    {
        if (cell->dim() != dim) continue;

        fe_values.reinit(cell);
        feo->dh<dim>()->get_dof_indices(cell, dof_indices);

        for (int k=0; k<qsize; k++)
        {
        	sources_conc[k]  = data.sources_conc.value(fe_values.point(k), cell->element_accessor());
        	sources_density[k]  = data.sources_density.value(fe_values.point(k), cell->element_accessor());
        	sources_sigma[k] = data.sources_sigma.value(fe_values.point(k), cell->element_accessor());
        	for (unsigned int sbi=0; sbi<n_subst; sbi++)
        	{
        		conc[k][sbi] = 0;
        		// TODO: Calculation of concentration has to take into account all values on the element.
        		// The following code will be wrong in case that dofs on some element are distributed
        		// to more than one process.
        		for (unsigned int i=0; i<ndofs; i++)
        			if (dof_indices[i] >= distr->begin() && dof_indices[i] <= distr->end())
        				conc[k][sbi] += ls[sbi]->get_solution_array()[dof_indices[i] - distr->begin()]*fe_values.shape_value(i,k);
        	}
        }

        // assemble the local stiffness matrix
        for (unsigned int sbi=0; sbi<n_subst; sbi++)
        {
        	for (int i=0; i<ndofs; i++)
        		local_rhs[i] = 0;

        	// compute sources
        	for (unsigned int k=0; k<qsize; k++)
        	{
        		conc_diff = sources_conc[k][sbi] - conc[k][sbi];
        		if (conc_diff > 0.0)
        			source = (sources_density[k][sbi] + conc_diff*sources_sigma[k][sbi])*fe_values.JxW(k);
        		else
        			source = sources_density[k][sbi]*fe_values.JxW(k);

        		for (int i=0; i<ndofs; i++)
        		{
        			if (dof_indices[i] < distr->begin() || dof_indices[i] > distr->end()) continue;

        			local_rhs[i] += source*fe_values.shape_value(i,k);
            	}
            }
        	ls[sbi]->rhs_set_values(ndofs, (int *)dof_indices, local_rhs);
        }
    }
}




template<unsigned int dim>
void TransportDG::assemble_fluxes_element_element()
{
    vector<FESideValues<dim,3>*> fe_values;
    FESideValues<dim,3> fsv_rt(*feo->map<dim>(), *feo->q<dim-1>(), *feo->fe_rt<dim>(),
    		update_values | update_normal_vectors | update_side_JxW_values);
    typename DOFHandler<dim,3>::CellIterator cell = feo->dh<dim>()->begin_cell();
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim-1>()->size();
    vector<unsigned int*> side_dof_indices;
    PetscScalar local_matrix[ndofs*ndofs], local_rhs[ndofs];
    vector<vector<arma::mat33> > side_K(qsize);
    vector<vector<arma::vec3> > side_velocity;
    vector<vector<double> > por_m(qsize), csection(qsize);
    vector<vector<arma::vec> > Dm(qsize), alphaL(qsize), alphaT(qsize);
    vector<arma::vec> dg_penalty;
    double gamma_l, omega[2], transport_flux;

    // assemble integral over sides
    FOR_EDGES( mesh_, edg )
    {
        if (edg->n_sides < 2 || edg->side(0)->element()->dim() != dim) continue;

        side_velocity.resize(edg->n_sides);
        for (unsigned int k=0; k<qsize; k++)
        {
        	por_m[k].resize(edg->n_sides);
        	csection[k].resize(edg->n_sides);
        	Dm[k].resize(edg->n_sides);
        	alphaL[k].resize(edg->n_sides);
        	alphaT[k].resize(edg->n_sides);
        	side_K[k].resize(edg->n_sides);
        }
        dg_penalty.resize(edg->n_sides);

        for (int i=side_dof_indices.size(); i<edg->n_sides; i++)
        	side_dof_indices.push_back(new unsigned int[ndofs]);

        for (int sid=fe_values.size(); sid<edg->n_sides; sid++)
        	fe_values.push_back(new FESideValues<dim,3>(*feo->map<dim>(), *feo->q<dim-1>(), *feo->fe<dim>(),
        			update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points));

		for (int sid=0; sid<edg->n_sides; sid++)
		{
			cell = mesh_->element.full_iter(edg->side(sid)->element());
			ElementAccessor<3> ele_acc = cell->element_accessor();
			feo->dh<dim>()->get_dof_indices(cell, side_dof_indices[sid]);
			fe_values[sid]->reinit(cell, edg->side(sid)->el_idx());
			fsv_rt.reinit(cell, edg->side(sid)->el_idx());
			calculate_velocity(cell, side_velocity[sid], fsv_rt);

			for (int k=0; k<qsize; k++)
			{
				por_m[k][sid] = data.por_m.value(fe_values[sid]->point(k), ele_acc);
				csection[k][sid] = data.cross_section->value(fe_values[sid]->point(k), ele_acc);
				Dm[k][sid] = data.diff_m.value(fe_values[sid]->point(k), ele_acc);
				alphaL[k][sid] = data.disp_l.value(fe_values[sid]->point(k), ele_acc);
				alphaT[k][sid] = data.disp_t.value(fe_values[sid]->point(k), ele_acc);
			}
			dg_penalty[sid] = data.dg_penalty.value(cell->centre(), ele_acc);
		}


        // fluxes and penalty
		for (unsigned int sbi=0; sbi<n_subst; sbi++)
		{
			for (unsigned int k=0; k<qsize; k++)
				for (unsigned int sid=0; sid<edg->n_sides; sid++)
					calculate_dispersivity_tensor(side_K[k][sid], side_velocity[sid][k], Dm[k][sid][sbi], alphaL[k][sid][sbi], alphaT[k][sid][sbi], por_m[k][sid], csection[k][sid]);

			for (int s1=0; s1<edg->n_sides; s1++)
			{
				for (int s2=s1+1; s2<edg->n_sides; s2++)
				{
					ASSERT(edg->side(s1)->valid(), "Invalid side of edge.");
					arma::vec3 nv = fe_values[s1]->normal_vector(0);

					// set up the parameters for DG method
					set_DG_parameters_edge(*edg, s1, s2, side_K, fe_values[0]->normal_vector(0), dg_penalty[s1][sbi], dg_penalty[s2][sbi], gamma_l, omega, transport_flux);

					int sd[2];
					sd[0] = s1;
					sd[1] = s2;

#define AVERAGE(i,k,side_id)  (fe_values[sd[side_id]]->shape_value(i,k)*0.5)
#define WAVERAGE(i,k,side_id) (por_m[k][sd[side_id]]*csection[k][sd[side_id]]*arma::dot(side_K[k][sd[side_id]]*fe_values[sd[side_id]]->shape_grad(i,k),nv)*omega[side_id])
#define JUMP(i,k,side_id)     ((side_id==0?1:-1)*fe_values[sd[side_id]]->shape_value(i,k))

					// For selected pair of elements:
					for (int m=0; m<2; m++)
					{
						for (int n=0; n<2; n++)
						{
							for (int i=0; i<fe_values[sd[n]]->n_dofs(); i++)
								for (int j=0; j<fe_values[sd[m]]->n_dofs(); j++)
									local_matrix[i*fe_values[sd[m]]->n_dofs()+j] = 0;

							for (int k=0; k<qsize; k++)
							{
								double flux_times_JxW = transport_flux*fe_values[0]->JxW(k);
								double gamma_times_JxW = gamma_l*fe_values[0]->JxW(k);

								for (int i=0; i<fe_values[sd[n]]->n_dofs(); i++)
								{
									if (side_dof_indices[sd[n]][i] < distr->begin() || side_dof_indices[sd[n]][i] > distr->end()) continue;

									double flux_JxW_avg_i = flux_times_JxW*AVERAGE(i,k,n);
									double gamma_JxW_jump_i = gamma_times_JxW*JUMP(i,k,n);
									double JxW_jump_i = fe_values[0]->JxW(k)*JUMP(i,k,n);
									double JxW_var_wavg_i = fe_values[0]->JxW(k)*WAVERAGE(i,k,n)*dg_variant;

									for (int j=0; j<fe_values[sd[m]]->n_dofs(); j++)
									{
										int index = i*fe_values[sd[m]]->n_dofs()+j;

										// flux due to transport (applied on interior edges) (average times jump)
										local_matrix[index] -= flux_JxW_avg_i*JUMP(j,k,m);

										// penalty enforcing continuity across edges (applied on interior and Dirichlet edges) (jump times jump)
										local_matrix[index] += gamma_JxW_jump_i*JUMP(j,k,m);

										// terms due to diffusion
										local_matrix[index] -= WAVERAGE(j,k,m)*JxW_jump_i;
										local_matrix[index] -= JUMP(j,k,m)*JxW_var_wavg_i;
									}
								}
							}
							ls[sbi]->mat_set_values(fe_values[sd[n]]->n_dofs(), (int *)side_dof_indices[sd[n]], fe_values[sd[m]]->n_dofs(), (int *)side_dof_indices[sd[m]], local_matrix);
						}
					}
#undef AVERAGE
#undef WAVERAGE
#undef JUMP
				}
			}
		}
    }

    for (int i=0; i<fe_values.size(); i++)
        delete fe_values[i];

    for (int i=0; i<side_dof_indices.size(); i++)
        delete[] side_dof_indices[i];
}


template<unsigned int dim>
void TransportDG::assemble_fluxes_boundary()
{
    FESideValues<dim,3> fe_values_side(*feo->map<dim>(), *feo->q<dim-1>(), *feo->fe<dim>(),
    		update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
    FESideValues<dim,3> fsv_rt(*feo->map<dim>(), *feo->q<dim-1>(), *feo->fe_rt<dim>(),
    		update_values | update_normal_vectors | update_side_JxW_values);
    typename DOFHandler<dim,3>::CellIterator cell = feo->dh<dim>()->begin_cell();
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim-1>()->size();
    unsigned int side_dof_indices[ndofs];
    PetscScalar local_matrix[ndofs*ndofs], local_rhs[ndofs];
    vector<arma::mat33> side_K(qsize);
    vector<arma::vec3> side_velocity;
    vector<double> por_m(qsize), csection(qsize);
    vector<arma::vec> Dm(qsize), alphaL(qsize), alphaT(qsize);
    arma::vec dg_penalty;
    double gamma_l;

    // assemble boundary integral
    FOR_BOUNDARIES(mesh_, b)
    {
        if (b->side()->dim() != dim-1) continue;

        if (b->side()->cond() == NULL || mh_dh->side_flux( *(b->side()) ) >= -mh_dh->precision()) continue;

        cell = mesh().element.full_iter(b->side()->element());
        ElementAccessor<3> ele_acc = cell->element_accessor();
        feo->dh<dim>()->get_dof_indices(cell, side_dof_indices);
        fe_values_side.reinit(cell, b->side()->el_idx());
        fsv_rt.reinit(cell, b->side()->el_idx());

        calculate_velocity(cell, side_velocity, fsv_rt);
        for (int k=0; k<qsize; k++)
        {
        	Dm[k] = data.diff_m.value(fe_values_side.point(k), ele_acc);
        	alphaL[k] = data.disp_l.value(fe_values_side.point(k), ele_acc);
        	alphaT[k] = data.disp_t.value(fe_values_side.point(k), ele_acc);
        	por_m[k]  = data.por_m.value(fe_values_side.point(k), ele_acc);
        	csection[k] = data.cross_section->value(fe_values_side.point(k), ele_acc);
        }
        dg_penalty = data.dg_penalty.value(cell->centre(), ele_acc);
        for (unsigned int sbi=0; sbi<n_subst; sbi++)
        {
        	for (unsigned int k=0; k<qsize; k++)
        		calculate_dispersivity_tensor(side_K[k], side_velocity[k], Dm[k][sbi], alphaL[k][sbi], alphaT[k][sbi], por_m[k], csection[k]);

			// set up the parameters for DG method
			set_DG_parameters_boundary(b->side(), side_K, fe_values_side.normal_vector(0), dg_penalty[sbi], gamma_l);
			if (b->side()->cond() != 0)
				gamma[sbi][b->side()->cond_idx()] = gamma_l;

			// fluxes and penalty
			for (int i=0; i<ndofs; i++)
			{
				for (int j=0; j<ndofs; j++)
				{
					local_matrix[i*ndofs+j] = 0;
					if (side_dof_indices[i] < distr->begin() || side_dof_indices[i] > distr->end()) continue;
					for (int k=0; k<qsize; k++)
					{
						// penalty enforcing continuity across edges (applied on interior and Dirichlet edges)
						local_matrix[i*ndofs+j] += gamma_l*fe_values_side.shape_value(j,k)*fe_values_side.shape_value(i,k)*fe_values_side.JxW(k);
					}
				}
			}
			ls[sbi]->mat_set_values(ndofs, (int *)side_dof_indices, ndofs, (int *)side_dof_indices, local_matrix);
        }
    }
}


template<unsigned int dim>
void TransportDG::assemble_fluxes_element_side()
{
	if (dim == 1) return;

    FEValues<dim-1,3> fe_values_vb(*feo->map<dim-1>(), *feo->q<dim-1>(), *feo->fe<dim-1>(),
    		update_values | update_gradients | update_JxW_values | update_quadrature_points);
    FESideValues<dim,3> fe_values_side(*feo->map<dim>(), *feo->q<dim-1>(), *feo->fe<dim>(),
    		update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
    vector<FEValuesSpaceBase<3>*> fv_sb(2);
    typename DOFHandler<dim,3>::CellIterator cell = feo->dh<dim>()->begin_cell();
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim-1>()->size();
    unsigned int side_dof_indices[2*ndofs], n_dofs[2];
    PetscScalar local_matrix[4*ndofs*ndofs];
    double transport_flux, comm_flux[2][2], por_m[2][qsize], csection[qsize];
    arma::vec sigma[qsize];

    // index 0 = element with lower dimension,
    // index 1 = side of element with higher dimension
    fv_sb[0] = &fe_values_vb;
    fv_sb[1] = &fe_values_side;

    // assemble integral over sides
    FOR_NEIGHBOURS(mesh_ , nb)
    {
        // skip neighbours of different dimension
        if (nb->element()->dim() != dim-1) continue;

		typename DOFHandler<dim-1,3>::CellIterator cell_sub = mesh_->element.full_iter(nb->element());
		feo->dh<dim-1>()->get_dof_indices(cell_sub, side_dof_indices);
		fe_values_vb.reinit(cell_sub);
		n_dofs[0] = fv_sb[0]->n_dofs();

		cell = nb->side()->element();
		feo->dh<dim>()->get_dof_indices(cell, side_dof_indices+n_dofs[0]);
		fe_values_side.reinit(cell, nb->side()->el_idx());
		n_dofs[1] = fv_sb[1]->n_dofs();

		// flux from the higher dimension to the lower one
		transport_flux = mh_dh->side_flux( *(nb->side()) )/nb->side()->measure();
		for (unsigned int k=0; k<qsize; k++)
		{
			por_m[0][k] = data.por_m.value(fe_values_vb.point(k), nb->element()->element_accessor());
			por_m[1][k] = data.por_m.value(fe_values_side.point(k), cell->element_accessor());
			csection[k] = data.cross_section->value(fe_values_side.point(k), cell->element_accessor());
			sigma[k]    = data.sigma_c.value(fe_values_vb.point(k), nb->element()->element_accessor());
		}

		for (unsigned int sbi=0; sbi<n_subst; sbi++)
		{
			for (int i=0; i<n_dofs[0]+n_dofs[1]; i++)
				for (int j=0; j<n_dofs[0]+n_dofs[1]; j++)
					local_matrix[i*(n_dofs[0]+n_dofs[1])+j] = 0;

			// set transmission conditions
			for (int k=0; k<qsize; k++)
			{
				/* The communication flux has two parts:
				 * - "diffusive" term containing sigma
				 * - "advective" term representing usual upwind
				 *
				 * The last term comm_flux[1][1] is the total mass flux,
				 * hence there is an extra contribution "-transport_flux".
				 */
				comm_flux[0][0] =  (csection[k]*por_m[0][k]*sigma[k][sbi]-min(0.,transport_flux))*fv_sb[0]->JxW(k);
				comm_flux[0][1] = -(csection[k]*por_m[0][k]*sigma[k][sbi]-min(0.,transport_flux))*fv_sb[0]->JxW(k);
				comm_flux[1][0] = -(csection[k]*por_m[1][k]*sigma[k][sbi]+max(0.,transport_flux))*fv_sb[0]->JxW(k);
				comm_flux[1][1] =  (csection[k]*por_m[1][k]*sigma[k][sbi]+max(0.,transport_flux)-transport_flux)*fv_sb[0]->JxW(k);

				for (int n=0; n<2; n++)
				{
					for (int i=0; i<n_dofs[n]; i++)
					{
						if (side_dof_indices[i+n*n_dofs[0]] < distr->begin() || side_dof_indices[i+n*n_dofs[0]] > distr->end()) continue;
						for (int m=0; m<2; m++)
						{
							for (int j=0; j<n_dofs[m]; j++)
							{
								int index = (i+n*n_dofs[0])*(n_dofs[0]+n_dofs[1]) + m*n_dofs[0] + j;
								local_matrix[index] += comm_flux[m][n]*fv_sb[m]->shape_value(j,k)*fv_sb[n]->shape_value(i,k);
							}
						}
					}
				}
			}
			ls[sbi]->mat_set_values(n_dofs[0]+n_dofs[1], (int *)side_dof_indices, n_dofs[0]+n_dofs[1], (int *)side_dof_indices, local_matrix);
		}
    }
}




void TransportDG::set_boundary_conditions()
{
  START_TIMER("assemble_bc");
	set_boundary_conditions<1>();
	set_boundary_conditions<2>();
	set_boundary_conditions<3>();
  END_TIMER("assemble_bc");
}


template<unsigned int dim>
void TransportDG::set_boundary_conditions()
{
    typename DOFHandler<dim,3>::CellIterator cell = feo->dh<dim>()->begin_cell();
    FESideValues<dim,3> fe_values_side(*feo->map<dim>(), *feo->q<dim-1>(), *feo->fe<dim>(),
    		update_values | update_side_JxW_values | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim-1>()->size();
    unsigned int side_dof_indices[ndofs];
    double local_rhs[ndofs];
    vector<arma::vec> bc_values(qsize);
    int rank;

    //TODO: distribution of BC elements
    // BC are now set only by process 0
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank != 0) return;

    FOR_BOUNDARIES(mesh_, b)
    {
    	cell = b->side()->element();
    	if (cell->dim() != dim) continue;

    	// skip outflow boundaries
    	if (mh_dh->side_flux( *(b->side()) ) >= -mh_dh->precision()) continue;

        ElementAccessor<3> ele_acc = b->element_accessor();
        fe_values_side.reinit(cell, b->side()->el_idx());
        feo->dh<dim>()->get_dof_indices(cell, side_dof_indices);

        for (unsigned int k=0; k<qsize; k++)
        	bc_values[k] = data.bc_conc.value(fe_values_side.point(k), ele_acc);

        for (unsigned int sbi=0; sbi<n_subst; sbi++)
        {
        	for (int i=0; i<ndofs; i++) local_rhs[i] = 0;

        	for (int k=0; k<qsize; k++)
        	{
        		double bc_term = gamma[sbi][b->side()->cond_idx()]*bc_values[k][sbi]*fe_values_side.JxW(k);
        		for (int i=0; i<ndofs; i++)
					local_rhs[i] += bc_term*fe_values_side.shape_value(i,k);
			}
			ls[sbi]->rhs_set_values(ndofs, (int *)side_dof_indices, local_rhs);
        }
    }
}



// TODO: The detection of side number from SideIter
// in TransportDG::calculate_velocity() and TransportDG::calculate_divergence()
// should be done with help of class RefElement. This however requires
// that the MH_DofHandler uses the node/side ordering defined in
// the respective RefElement.
template<unsigned int dim>
void TransportDG::calculate_velocity(const typename DOFHandler<dim,3>::CellIterator &cell, vector<arma::vec3> &velocity, FEValuesBase<dim,3> &fv)
{
    std::map<const Node*, int> node_nums;
    for (int i=0; i<cell->n_nodes(); i++)
        node_nums[cell->node[i]] = i;

    velocity.resize(fv.n_points());

    for (int k=0; k<fv.n_points(); k++)
    {
        velocity[k].zeros();
        for (int sid=0; sid<cell->n_sides(); sid++)
        {
            if (cell->side(sid)->dim() != dim-1) continue;
            int num = dim*(dim+1)/2;
            for (int i=0; i<cell->side(sid)->n_nodes(); i++)
                num -= node_nums[cell->side(sid)->node(i)];
            velocity[k] += fv.shape_vector(num,k) * mh_dh->side_flux( *(cell->side(sid)) );
        }
    }
}


template<unsigned int dim>
void TransportDG::calculate_velocity_divergence(const typename DOFHandler<dim,3>::CellIterator &cell, vector<double> &divergence, FEValuesBase<dim,3> &fv)
{
    std::map<const Node*, int> node_nums;
    for (int i=0; i<cell->n_nodes(); i++)
        node_nums[cell->node[i]] = i;

    divergence.resize(fv.n_points());

    for (int k=0; k<fv.n_points(); k++)
    {
        divergence[k] = 0;
        for (int sid=0; sid<cell->n_sides(); sid++)
        {
            if (cell->side(sid)->dim() != dim-1) continue;
            int num = dim*(dim+1)/2;
            for (int i=0; i<cell->side(sid)->n_nodes(); i++)
                num -= node_nums[cell->side(sid)->node(i)];
            divergence[k] += arma::trace(fv.shape_grad_vector(num,k)) * mh_dh->side_flux( *(cell->side(sid)) );
        }
    }
}




void TransportDG::calculate_dispersivity_tensor(arma::mat33 &K, const arma::vec3 &velocity,
		double Dm, double alphaL, double alphaT, double porosity, double cross_cut)
{
    double vnorm = arma::norm(velocity, 2);

	if (fabs(vnorm) > sqrt(numeric_limits<double>::epsilon()))
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				K(i,j) = velocity[i]*velocity[j]/(vnorm*vnorm)*(alphaL-alphaT) + alphaT*(i==j?1:0);
	else
		K.zeros();

	K = (vnorm*porosity*cross_cut)*K + (Dm*pow(porosity, 1./3))*arma::eye(3,3);
}


void TransportDG::set_DG_parameters_edge(const Edge &edg,
            const int s1,
            const int s2,
            const vector< vector<arma::mat33> > &K,
            const arma::vec3 &normal_vector,
            const double alpha1,
            const double alpha2,
            double &gamma,
            double *omega,
            double &transport_flux)
{
    double delta[2];
    double h = 0;
    double fluxes[edg.n_sides];
    double local_alpha = 0;

    ASSERT(edg.side(s1)->valid(), "Invalid side of an edge.");
    SideIter s = edg.side(s1);

    // calculate the side diameter
    if (s->dim() == 0)
    {
        h = 1;
    }
    else
    {
        for (int i=0; i<s->n_nodes(); i++)
            for (int j=i+1; j<s->n_nodes(); j++)
                h = max(h, s->node(i)->distance(*s->node(j)));
    }

    // calculate the total in- and out-flux through the edge
    for (int i=0; i<edg.n_sides; i++) fluxes[i] = mh_dh->side_flux( *(edg.side(i)) )/edg.side(i)->measure();
    double pflux = 0, nflux = 0;
    for (int i=0; i<edg.n_sides; i++)
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

    gamma = 0.5*fabs(transport_flux);


    // determine local DG penalty

    local_alpha = max(alpha1, alpha2);

    if (s1 == s2)
    {
        omega[0] = 1;

        // delta is set to the average value of Kn.n on the side
        delta[0] = 0;
        for (int k=0; k<K.size(); k++)
            delta[0] += dot(K[k][s1]*normal_vector,normal_vector);
        delta[0] /= K.size();

        gamma += local_alpha/h*(delta[0]);
    }
    else
    {
        delta[0] = 0;
        delta[1] = 0;
        for (int k=0; k<K.size(); k++)
        {
            delta[0] += dot(K[k][s1]*normal_vector,normal_vector);
            delta[1] += dot(K[k][s2]*normal_vector,normal_vector);
        }
        delta[0] /= K.size();
        delta[1] /= K.size();

        double delta_sum = delta[0] + delta[1];

        if (delta_sum > numeric_limits<double>::epsilon())
        {
            omega[0] = delta[1]/delta_sum;
            omega[1] = delta[0]/delta_sum;
            gamma += local_alpha/h*(delta[0]*delta[1]/delta_sum);
        }
        else
            for (int i=0; i<2; i++) omega[i] = 0;
    }
}







void TransportDG::set_DG_parameters_boundary(const SideIter side,
            const vector<arma::mat33> &K,
            const arma::vec3 &normal_vector,
            const double alpha,
            double &gamma)
{
    double delta[2];
    double h = 0;

    // calculate the side diameter
    if (side->dim() == 0)
    {
        h = 1;
    }
    else
    {
        for (int i=0; i<side->n_nodes(); i++)
            for (int j=i+1; j<side->n_nodes(); j++)
                h = max(h, side->node(i)->distance( *side->node(j) ));
    }

    gamma = 0.5*fabs( mh_dh->side_flux( *(side) )/side->measure() );

	// delta is set to the average value of Kn.n on the side
	delta[0] = 0;
	for (int k=0; k<K.size(); k++)
		delta[0] += dot(K[k]*normal_vector,normal_vector);
	delta[0] /= K.size();

	gamma += alpha/h*delta[0];
}






void TransportDG::set_initial_condition()
{
	for (unsigned int sbi=0; sbi<n_subst; sbi++)
		ls[sbi]->start_allocation();
	prepare_initial_condition<1>();
	prepare_initial_condition<2>();
	prepare_initial_condition<3>();

	for (unsigned int sbi=0; sbi<n_subst; sbi++)
		ls[sbi]->start_add_assembly();
	prepare_initial_condition<1>();
	prepare_initial_condition<2>();
	prepare_initial_condition<3>();

	for (unsigned int sbi=0; sbi<n_subst; sbi++)
	{
		ls[sbi]->finalize();
		solve_system(solver, ls[sbi]);
	}
}

template<unsigned int dim>
void TransportDG::prepare_initial_condition()
{
	FEValues<dim,3> fe_values(*feo->map<dim>(), *feo->q<dim>(), *feo->fe<dim>(),
			update_values | update_JxW_values | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
    unsigned int dof_indices[ndofs], nid;
    double matrix[ndofs*ndofs], rhs[ndofs];
    arma::vec init_values[qsize];

    FOR_ELEMENTS(mesh_, elem)
    {
    	if (elem->dim() != dim) continue;

    	ElementAccessor<3> ele_acc = elem->element_accessor();
    	feo->dh<dim>()->get_dof_indices(elem, dof_indices);
    	fe_values.reinit(elem);

    	for (unsigned int k=0; k<qsize; k++)
    		init_values[k] = data.init_conc.value(fe_values.point(k), ele_acc);

    	for (unsigned int sbi=0; sbi<n_subst; sbi++)
    	{
    		for (unsigned int i=0; i<ndofs; i++)
    		{
    			rhs[i] = 0;
    			for (unsigned int j=0; j<ndofs; j++)
    				matrix[i*ndofs+j] = 0;
    		}

    		for (unsigned int k=0; k<qsize; k++)
    		{
    			double rhs_term = init_values[k](sbi)*fe_values.JxW(k);

    			for (unsigned int i=0; i<ndofs; i++)
    			{
    				if (dof_indices[i] < distr->begin() || dof_indices[i] > distr->end()) continue;

    				for (unsigned int j=0; j<ndofs; j++)
    					matrix[i*ndofs+j] += fe_values.shape_value(i,k)*fe_values.shape_value(j,k)*fe_values.JxW(k);

    				rhs[i] += fe_values.shape_value(i,k)*rhs_term;
    			}
    		}
    		ls[sbi]->set_values(ndofs, (int *)dof_indices, ndofs, (int *)dof_indices, matrix, rhs);
    	}
    }
}



void TransportDG::calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance)
{
	calc_fluxes<1>(bcd_balance, bcd_plus_balance, bcd_minus_balance);
	calc_fluxes<2>(bcd_balance, bcd_plus_balance, bcd_minus_balance);
	calc_fluxes<3>(bcd_balance, bcd_plus_balance, bcd_minus_balance);
}

template<unsigned int dim>
void TransportDG::calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance)
{
	FESideValues<dim,3> fe_values(*feo->map<dim>(), *feo->q<dim-1>(), *feo->fe<dim>(),
			update_values | update_gradients | update_side_JxW_values | update_quadrature_points);
	FESideValues<dim,3> fsv_rt(*feo->map<dim>(), *feo->q<dim-1>(), *feo->fe_rt<dim>(),
			update_values | update_normal_vectors | update_side_JxW_values);
	const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim-1>()->size();
	unsigned int dof_indices[ndofs];
	DOFHandler<1,3>::CellIterator cell = feo->dh<dim>()->begin_cell();
	vector<double> por_m(qsize), csection(qsize);
	vector<arma::vec3> side_velocity(qsize);
	double conc, mass_flux, water_flux;
	arma::vec3 c_grad;
	arma::vec Dm[qsize], alphaL[qsize], alphaT[qsize];
	arma::mat33 D;

    FOR_BOUNDARIES(mesh_, bcd) {

    	if (bcd->side()->dim() != dim-1) continue;

        cell = bcd->side()->element();
        ElementAccessor<3> ele_acc = cell->element_accessor();

		water_flux = mh_dh->side_flux(*(bcd->side()))/bcd->side()->measure();

		fe_values.reinit(cell, bcd->side()->el_idx());
		fsv_rt.reinit(cell, bcd->side()->el_idx());
		feo->dh<dim>()->get_dof_indices(cell, dof_indices);

		calculate_velocity(cell, side_velocity, fsv_rt);
		for (unsigned int k=0; k<qsize; k++)
		{
			por_m[k] = data.por_m.value(fe_values.point(k), ele_acc);
			csection[k] = data.cross_section->value(fe_values.point(k), ele_acc);
			Dm[k]     = data.diff_m.value(fe_values.point(k), ele_acc);
			alphaL[k] = data.disp_l.value(fe_values.point(k), ele_acc);
			alphaT[k] = data.disp_t.value(fe_values.point(k), ele_acc);
		}
		for (unsigned int sbi=0; sbi<n_subst; sbi++)
		{
			mass_flux = 0;

			for (unsigned int k=0; k<qsize; k++)
			{
				calculate_dispersivity_tensor(D, side_velocity[k], Dm[k][sbi], alphaL[k][sbi], alphaT[k][sbi], por_m[k], csection[k]);
				conc = 0;
				c_grad.zeros();
				for (unsigned int i=0; i<ndofs; i++)
				{
					if (dof_indices[i] < distr->begin() || dof_indices[i] > distr->end()) continue;
					conc += fe_values.shape_value(i,k)*ls[sbi]->get_solution_array()[dof_indices[i]-distr->begin()];
					c_grad += fe_values.shape_grad(i,k)*ls[sbi]->get_solution_array()[dof_indices[i]-distr->begin()];
				}

				mass_flux += (-csection[k]*por_m[k]*dot(D*c_grad,fsv_rt.normal_vector(k)) + water_flux*conc)*fe_values.JxW(k);
			}

			Region r = bcd->region();
			if (! r.is_valid()) xprintf(Msg, "Invalid region, ele % d, edg: % d\n", bcd->bc_ele_idx_, bcd->edge_idx_);
			unsigned int bc_region_idx = r.boundary_idx();
			bcd_balance[sbi][bc_region_idx] += mass_flux;

			if (mass_flux > 0) bcd_plus_balance[sbi][bc_region_idx] += mass_flux;
			else bcd_minus_balance[sbi][bc_region_idx] += mass_flux;
		}
    }

}

void TransportDG::calc_elem_sources(vector<vector<double> > &mass, vector<vector<double> > &src_balance)
{
	calc_elem_sources<1>(mass, src_balance);
	calc_elem_sources<2>(mass, src_balance);
	calc_elem_sources<3>(mass, src_balance);
}

template<unsigned int dim>
void TransportDG::calc_elem_sources(vector<vector<double> > &mass, vector<vector<double> > &src_balance)
{
	FEValues<dim,3> fe_values(*feo->map<dim>(), *feo->q<dim>(), *feo->fe<dim>(),
			update_values | update_JxW_values | update_quadrature_points);
	const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
	unsigned int dof_indices[ndofs];
	vector<arma::vec> sources_conc(qsize), sources_density(qsize), sources_sigma(qsize);
	vector<double> por_m(qsize), csection(qsize);
	double mass_sum, sources_sum, conc, conc_diff;

	FOR_ELEMENTS(mesh_, elem)
	{
		if (elem->dim() != dim) continue;

		ElementAccessor<3> ele_acc = elem->element_accessor();
		fe_values.reinit(elem);
		feo->dh<dim>()->get_dof_indices(elem, dof_indices);

		for (unsigned int k=0; k<qsize; k++)
		{
			por_m[k] = data.por_m.value(fe_values.point(k), ele_acc);
			csection[k] = data.cross_section->value(fe_values.point(k), ele_acc);
			sources_conc[k]    = data.sources_conc.value(fe_values.point(k), ele_acc);
			sources_density[k] = data.sources_density.value(fe_values.point(k), ele_acc);
			sources_sigma[k]   = data.sources_sigma.value(fe_values.point(k), ele_acc);
		}

		for (unsigned int sbi=0; sbi<n_subst; sbi++)
		{
			mass_sum = 0;
			sources_sum = 0;

			for (unsigned int k=0; k<qsize; k++)
			{
				conc = 0;
				for (unsigned int i=0; i<ndofs; i++)
				{
					// TODO: Correct calculation of actual concentration (see note in TransportDG::set_sources()).
					if (dof_indices[i] < distr->begin() || dof_indices[i] > distr->end()) continue;
					conc += fe_values.shape_value(i,k)*ls[sbi]->get_solution_array()[dof_indices[i]-distr->begin()];
				}

				mass_sum += por_m[k]*csection[k]*conc*fe_values.JxW(k);

				conc_diff = sources_conc[k][sbi] - conc;
				if ( conc_diff > 0.0)
					sources_sum += (sources_density[k][sbi] + conc_diff*sources_sigma[k][sbi])*fe_values.JxW(k);
				else
					sources_sum += sources_density[k][sbi]*fe_values.JxW(k);
			}

			Region r = elem->element_accessor().region();
			if (! r.is_valid()) xprintf(Msg, "Invalid region, ele % d\n", elem.index());
			unsigned int region_idx = r.bulk_idx();
			mass[sbi][region_idx] += mass_sum;
			src_balance[sbi][region_idx] += sources_sum;
		}
	}
}











