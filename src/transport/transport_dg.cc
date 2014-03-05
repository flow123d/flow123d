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
#include <boost/foreach.hpp>
#include "system/xio.h"
#include "system/sys_profiler.hh"
#include "transport/transport_dg.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_values.hh"
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "io/output.h"
#include "mesh/boundaries.h"
#include "la/distribution.hh"
#include "input/accessors.hh"
#include "flow/darcy_flow_mh.hh"
#include "la/linsys_PETSC.hh"
#include "mesh/partitioning.hh"
#include "transport/advection_diffusion_model.hh"
#include "transport/concentration_model.hh"
#include "transport/heat_model.hh"

using namespace Input::Type;

template<class Model>
Selection TransportDG<Model>::dg_variant_selection_input_type
	= Selection("DG_variant", "Type of penalty term.")
	.add_value(non_symmetric, "non-symmetric", "non-symmetric weighted interior penalty DG method")
	.add_value(incomplete,    "incomplete",    "incomplete weighted interior penalty DG method")
	.add_value(symmetric,     "symmetric",     "symmetric weighted interior penalty DG method");

template<class Model>
Selection TransportDG<Model>::EqData::bc_type_selection =
              Selection("TransportDG_EqData_bc_Type")
               .add_value(inflow, "inflow", "Dirichlet BC on inflow and homogeneous Neumann BC on outflow.")
               .add_value(dirichlet, "dirichlet")
               .add_value(neumann, "neumann")
               .add_value(robin, "robin");

template<class Model>
Record TransportDG<Model>::input_type
	= Model::get_input_type("DG", "DG solver")
    .declare_key("solver", LinSys_PETSC::input_type, Default::obligatory(),
            "Linear solver for MH problem.")
    .declare_key("bc_data", Array(TransportDG<Model>::EqData().boundary_input_type()), IT::Default::obligatory(), "")
    .declare_key("bulk_data", Array(TransportDG<Model>::EqData().bulk_input_type()),
    		IT::Default::obligatory(), "")
    .declare_key("dg_variant", TransportDG<Model>::dg_variant_selection_input_type, Default("non-symmetric"),
    		"Variant of interior penalty discontinuous Galerkin method.")
    .declare_key("dg_order", Integer(0,3), Default("1"),
    		"Polynomial order for finite element in DG method (order 0 is suitable if there is no diffusion/dispersion).");



FEObjects::FEObjects(Mesh *mesh_, unsigned int fe_order)
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

	case 3:
		q_order = 6;
		fe1_ = new FE_P_disc<3,1,3>;
		fe2_ = new FE_P_disc<3,2,3>;
		fe3_ = new FE_P_disc<3,3,3>;
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

	dh_ = new DOFHandlerMultiDim(*mesh_);

	dh_->distribute_dofs(*fe1_, *fe2_, *fe3_);
}


FEObjects::~FEObjects()
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
	delete dh_;
}

template<> FiniteElement<0,3> *FEObjects::fe<0>() { return 0; }
template<> FiniteElement<1,3> *FEObjects::fe<1>() { return fe1_; }
template<> FiniteElement<2,3> *FEObjects::fe<2>() { return fe2_; }
template<> FiniteElement<3,3> *FEObjects::fe<3>() { return fe3_; }

template<> FiniteElement<0,3> *FEObjects::fe_rt<0>() { return 0; }
template<> FiniteElement<1,3> *FEObjects::fe_rt<1>() { return fe_rt1_; }
template<> FiniteElement<2,3> *FEObjects::fe_rt<2>() { return fe_rt2_; }
template<> FiniteElement<3,3> *FEObjects::fe_rt<3>() { return fe_rt3_; }

template<> Quadrature<0> *FEObjects::q<0>() { return q0_; }
template<> Quadrature<1> *FEObjects::q<1>() { return q1_; }
template<> Quadrature<2> *FEObjects::q<2>() { return q2_; }
template<> Quadrature<3> *FEObjects::q<3>() { return q3_; }

template<> Mapping<0,3> *FEObjects::mapping<0>() { return map0_; }
template<> Mapping<1,3> *FEObjects::mapping<1>() { return map1_; }
template<> Mapping<2,3> *FEObjects::mapping<2>() { return map2_; }
template<> Mapping<3,3> *FEObjects::mapping<3>() { return map3_; }

DOFHandlerMultiDim *FEObjects::dh() { return dh_; }


template<class Model>
TransportDG<Model>::EqData::EqData() : Model::ModelEqData()
{
	ADD_FIELD(sigma_c, "Coefficient of diffusive transfer through fractures (for each substance).", Default("1"));
	ADD_FIELD(dg_penalty, "Penalty parameter influencing the discontinuity of the solution (for each substance). "
			"Its default value 1 is sufficient in most cases. Higher value diminishes the inter-element jumps.", Default("1.0"));

    ADD_FIELD(bc_type,"Boundary condition type, possible values: inflow, dirichlet, neumann, robin.", Default("inflow") );
    bc_type.set_selection(&bc_type_selection);

//    std::vector<FieldEnum> list; list.push_back(neumann);
    ADD_FIELD(bc_flux,"Flux in Neumann boundary condition.", Default("0.0"));
    ADD_FIELD(bc_robin_sigma,"Conductivity coefficient in Robin boundary condition.", Default("0.0"));
}


template<class Model>
TransportDG<Model>::TransportDG(Mesh & init_mesh, const Input::Record &in_rec)
        : TransportBase(init_mesh, in_rec),
          mass_matrix(0),
          allocation_done(false)
{
	// Check that Model is derived from AdvectionDiffusionModel.
	static_assert(std::is_base_of<AdvectionDiffusionModel, Model>::value, "");

    time_ = new TimeGovernor(in_rec.val<Input::Record>("time"), equation_mark_type_);
    time_->fix_dt_until_mark();

    // Read names of transported substances.
    Model::set_component_names(subst_names_, in_rec);
    n_subst_ = subst_names_.size();

    Input::Iterator<Input::Record> it = in_rec.find<Input::Record>("mass_balance");
    if (it)
    	mass_balance_ = new MassBalance(this, *it);

    // Set up physical parameters.
    data_.set_mesh(&init_mesh);
    data_.bc_type.set_n_comp(n_subst_);
    data_.bc_flux.set_n_comp(n_subst_);
    data_.bc_robin_sigma.set_n_comp(n_subst_);
    data_.sigma_c.set_n_comp(n_subst_);
    data_.dg_penalty.set_n_comp(n_subst_);
    Model::init_data(n_subst_);
    data_.init_from_input( in_rec.val<Input::Array>("bulk_data"), in_rec.val<Input::Array>("bc_data") );
    data_.set_time(*time_);

    // DG variant
    dg_variant = in_rec.val<DGVariant>("dg_variant");

    // DG stabilization parameters on boundary edges
    gamma.resize(n_subst_);
    for (int sbi=0; sbi<n_subst_; sbi++)
    	gamma[sbi].resize(mesh_->boundary_.size());


    // create finite element structures and distribute DOFs
    dg_order = in_rec.val<unsigned int>("dg_order");
    feo = new FEObjects(mesh_, dg_order);
    DBGMSG("TDG: solution size %d\n", feo->dh()->n_global_dofs());

    // Resize coefficient arrays
    int qsize = max(feo->q<0>()->size(), max(feo->q<1>()->size(), max(feo->q<2>()->size(), feo->q<3>()->size())));
    int max_edg_sides = max(mesh_->max_edge_sides(1), max(mesh_->max_edge_sides(2), mesh_->max_edge_sides(3)));
    mm_coef.resize(qsize);
    ad_coef.resize(n_subst_);
    dif_coef.resize(n_subst_);
    for (int sbi=0; sbi<n_subst_; sbi++)
    {
      ad_coef[sbi].resize(qsize);
      dif_coef[sbi].resize(qsize);
    }
    ad_coef_edg.resize(max_edg_sides);
    dif_coef_edg.resize(max_edg_sides);
    for (int sd=0; sd<max_edg_sides; sd++)
    {
    	ad_coef_edg[sd].resize(n_subst_);
    	dif_coef_edg[sd].resize(n_subst_);
    	for (int sbi=0; sbi<n_subst_; sbi++)
    	{
    		ad_coef_edg[sd][sbi].resize(qsize);
    		dif_coef_edg[sd][sbi].resize(qsize);
    	}
    }

    // set up output class
    Input::Record output_rec = in_rec.val<Input::Record>("output");
    transport_output = OutputTime::output_stream(output_rec.val<Input::Record>("output_stream"));

    // allocate output arrays
    if (feo->dh()->el_ds()->myp() == 0)
    {
    	int n_corners = 0;
    	FOR_ELEMENTS(mesh_, elem)
    		n_corners += elem->dim()+1;
    	output_solution.resize(n_subst_);
    	for (int i=0; i<n_subst_; i++)
    	{
			output_solution[i] = new double[n_corners];
			for(int j=0; j<n_corners; j++)
				output_solution[i][j] = 0.0;
			OutputTime::register_corner_data<double>(mesh_, subst_names_[i], "M/L^3",
					output_rec.val<Input::Record>("output_stream"), output_solution[i], n_corners);
		}
    }

    // set time marks for writing the output
    output_mark_type = this->mark_type() | time_->marks().type_fixed_time() | time_->marks().type_output();
    time_->marks().add_time_marks(0.0, output_rec.val<double>("save_step"), time_->end_time(), output_mark_type);

    // allocate matrix and vector structures
    
    ls    = new LinSys*[n_subst_];
    ls_dt = new LinSys_PETSC(feo->dh()->distr());
    ( (LinSys_PETSC *)ls_dt )->set_from_input( in_rec.val<Input::Record>("solver") );
//    ( (LinSys_PETSC *)ls_dt )->set_initial_guess_nonzero();
    for (int sbi = 0; sbi < n_subst_; sbi++) {
    	ls[sbi] = new LinSys_PETSC(feo->dh()->distr());
    	( (LinSys_PETSC *)ls[sbi] )->set_from_input( in_rec.val<Input::Record>("solver") );
    	ls[sbi]->set_solution(NULL);
    }
    stiffness_matrix = new Mat[n_subst_];
    rhs = new Vec[n_subst_];


    // set initial conditions
    set_initial_condition();
    for (int sbi = 0; sbi < n_subst_; sbi++)
    	( (LinSys_PETSC *)ls[sbi] )->set_initial_guess_nonzero();

}

template<class Model>
TransportDG<Model>::~TransportDG()
{
    delete time_;
    delete ls_dt;

    if (feo->dh()->el_ds()->myp() == 0)
    {
		for (int i=0; i<n_subst_; i++)
			delete[] output_solution[i];
    }

    for (int i=0; i<n_subst_; i++) delete ls[i];
    delete[] ls;
    delete[] stiffness_matrix;
    delete[] rhs;
    delete feo;
    if (mass_balance_ != NULL) delete mass_balance_;

    gamma.clear();
    subst_names_.clear();
}



template<class Model>
void TransportDG<Model>::update_solution()
{
	START_TIMER("DG-ONE STEP");

	// save solution and calculate mass balance at initial time
	if (!allocation_done)
	{
		output_data();
	}

    time_->next_time();
    time_->view("TDG");
    
    START_TIMER("data reinit");
    data_.set_time(*time_);
    END_TIMER("data reinit");
    
    // check first time assembly - needs preallocation
    if (!allocation_done)
    {
    	// preallocate mass matrix
    	ls_dt->start_allocation();
    	assemble_mass_matrix();
    	mass_matrix = NULL;

		// preallocate system matrix
		for (int i=0; i<n_subst_; i++)
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
		mass_matrix_changed())
	{
		ls_dt->start_add_assembly();
		assemble_mass_matrix();
		ls_dt->finish_assembly();
		mass_matrix = ls_dt->get_matrix();
	}

	// assemble stiffness matrix
    if (stiffness_matrix[0] == NULL ||
    	stiffness_matrix_changed())
    {
        // new fluxes can change the location of Neumann boundary,
        // thus stiffness matrix must be reassembled
    	for (int i=0; i<n_subst_; i++)
    	{
    		ls[i]->start_add_assembly();
    		ls[i]->mat_zero_entries();
    	}
        assemble_stiffness_matrix();
        for (int i=0; i<n_subst_; i++)
        {
        	ls[i]->finish_assembly();

        	if (stiffness_matrix[i] == NULL)
        		MatConvert(ls[i]->get_matrix(), MATSAME, MAT_INITIAL_MATRIX, &stiffness_matrix[i]);
        	else
        		MatCopy(ls[i]->get_matrix(), stiffness_matrix[i], DIFFERENT_NONZERO_PATTERN);
        }
    }

    // assemble right hand side (due to sources and boundary conditions)
    if (rhs[0] == NULL ||
    	rhs_changed())
    {
    	for (int i=0; i<n_subst_; i++)
    	{
    		ls[i]->start_add_assembly();
    		ls[i]->rhs_zero_entries();
    	}
    	set_sources();
    	set_boundary_conditions();
    	for (int i=0; i<n_subst_; i++)
    	{
    		ls[i]->finish_assembly();

    		VecDuplicate(ls[i]->get_rhs(), &rhs[i]);
    		VecCopy(ls[i]->get_rhs(), rhs[i]);
    	}
    }

    Model::flux_changed = false;


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
    Mat m;
    START_TIMER("solve");
    for (int i=0; i<n_subst_; i++)
    {
    	MatConvert(stiffness_matrix[i], MATSAME, MAT_INITIAL_MATRIX, &m);
		MatAXPY(m, 1./time_->dt(), mass_matrix, SUBSET_NONZERO_PATTERN);
		ls[i]->set_matrix(m, DIFFERENT_NONZERO_PATTERN);
		Vec y,w;
		VecDuplicate(rhs[i], &y);
		VecDuplicate(rhs[i], &w);
		MatMult(mass_matrix, ls[i]->get_solution(), y);
		VecWAXPY(w, 1./time_->dt(), y, rhs[i]);
		ls[i]->set_rhs(w);

		VecDestroy(&y);
		VecDestroy(&w);
		MatDestroy(&m);

		ls[i]->solve();
    }
    END_TIMER("solve");

    if (mass_balance() != NULL)
    	mass_balance()->calculate(time_->t());

    END_TIMER("DG-ONE STEP");
}



template<class Model>
void TransportDG<Model>::get_solution_vector(double *& vector, unsigned int & size)
{}


template<class Model>
void TransportDG<Model>::get_parallel_solution_vector(Vec & vector)
{}


template<class Model>
void TransportDG<Model>::set_velocity_field(const MH_DofHandler &dh)
{
    // So far the velocity_vector contains zeros, so we ignore it.
    // Instead we use the value Side.flux.

    mh_dh = &dh;
	Model::flux_changed = true;

}


template<class Model>
void TransportDG<Model>::output_data()
{
    double *solution;
    unsigned int dof_indices[max(feo->fe<1>()->n_dofs(), max(feo->fe<2>()->n_dofs(), feo->fe<3>()->n_dofs()))];

    if (!time_->is_current(output_mark_type)) return;

    START_TIMER("DG-OUTPUT");

    // gather the solution from all processors
    IS is;
    VecScatter output_scatter;
    int row_ids[feo->dh()->n_global_dofs()];
	Vec solution_vec;

	for (int i=0; i<feo->dh()->n_global_dofs(); i++)
		row_ids[i] = i;

	for (int sbi=0; sbi<n_subst_; sbi++)
	{
		VecCreateSeq(PETSC_COMM_SELF, ls[sbi]->size(), &solution_vec);

		ISCreateGeneral(PETSC_COMM_SELF, ls[sbi]->size(), row_ids, PETSC_COPY_VALUES, &is); //WithArray
		VecScatterCreate(ls[sbi]->get_solution(), is, solution_vec, PETSC_NULL, &output_scatter);
		VecScatterBegin(output_scatter, ls[sbi]->get_solution(), solution_vec, INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(output_scatter, ls[sbi]->get_solution(), solution_vec, INSERT_VALUES, SCATTER_FORWARD);
		VecScatterDestroy(&(output_scatter));
		ISDestroy(&(is));

		// on the main processor fill the output array and save to file
		if (feo->dh()->el_ds()->myp() == 0)
		{
			int corner_id = 0; //, node_id;

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

			FEValues<1,3> fv1(*feo->mapping<1>(), q1, *feo->fe<1>(), update_values);
			FEValues<2,3> fv2(*feo->mapping<2>(), q2, *feo->fe<2>(), update_values);
			FEValues<3,3> fv3(*feo->mapping<3>(), q3, *feo->fe<3>(), update_values);

			VecGetArray(solution_vec, &solution);
			FOR_ELEMENTS(mesh_, elem)
			{
				feo->dh()->get_dof_indices(elem, dof_indices);
				switch (elem->dim())
				{
				case 1:
					fv1.reinit(elem);
					for (unsigned int k=0; k<q1.size(); k++)
					{
						output_solution[sbi][corner_id] = 0;
						for (unsigned int i=0; i<fv1.n_dofs(); i++)
							output_solution[sbi][corner_id] += solution[dof_indices[i]]*fv1.shape_value(i,k);
						corner_id++;
					}
					break;
				case 2:
					fv2.reinit(elem);
					for (unsigned int k=0; k<q2.size(); k++)
					{
						output_solution[sbi][corner_id] = 0;
						for (unsigned int i=0; i<fv2.n_dofs(); i++)
							output_solution[sbi][corner_id] += solution[dof_indices[i]]*fv2.shape_value(i,k);
						corner_id++;
					}
					break;
				case 3:
					fv3.reinit(elem);
					for (unsigned int k=0; k<q3.size(); k++)
					{
						output_solution[sbi][corner_id] = 0;
						for (unsigned int i=0; i<fv3.n_dofs(); i++)
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

	if(transport_output) {
		xprintf(MsgLog, "transport DG: write_data()\n");
		transport_output->write_data(time_->t());
	}
	if (mass_balance() != NULL)
		mass_balance()->output(time_->t());

    END_TIMER("DG-OUTPUT");
}


template<class Model>
void TransportDG<Model>::assemble_mass_matrix()
{
  START_TIMER("assemble_mass");
	assemble_mass_matrix<1>();
	assemble_mass_matrix<2>();
	assemble_mass_matrix<3>();
  END_TIMER("assemble_mass");
}


template<class Model> template<unsigned int dim>
void TransportDG<Model>::assemble_mass_matrix()
{
    FEValues<dim,3> fe_values(*feo->mapping<dim>(), *feo->q<dim>(), *feo->fe<dim>(), update_values | update_JxW_values | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
    unsigned int dof_indices[ndofs];
    PetscScalar local_mass_matrix[ndofs*ndofs];

    // assemble integral over elements
    for (int i_cell=0; i_cell<feo->dh()->el_ds()->lsize(); i_cell++)
    {
    	typename DOFHandlerBase::CellIterator cell = mesh_->element(feo->dh()->el_index(i_cell));
        if (cell->dim() != dim) continue;

        fe_values.reinit(cell);
        feo->dh()->get_dof_indices(cell, dof_indices);
        ElementAccessor<3> ele_acc = cell->element_accessor();

        Model::compute_mass_matrix_coefficient(fe_values.point_list(), ele_acc, mm_coef);

        // assemble the local stiffness and mass matrix
        for (unsigned int i=0; i<ndofs; i++)
        {
            for (unsigned int j=0; j<ndofs; j++)
            {
                local_mass_matrix[i*ndofs+j] = 0;
                for (unsigned int k=0; k<qsize; k++)
                    local_mass_matrix[i*ndofs+j] += mm_coef[k]*fe_values.shape_value(j,k)*fe_values.shape_value(i,k)*fe_values.JxW(k);
            }
        }

        ls_dt->mat_set_values(ndofs, (int *)dof_indices, ndofs, (int *)dof_indices, local_mass_matrix);
    }
}




template<class Model>
void TransportDG<Model>::assemble_stiffness_matrix()
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



template<class Model>
template<unsigned int dim>
void TransportDG<Model>::assemble_volume_integrals()
{
    FEValues<dim,3> fv_rt(*feo->mapping<dim>(), *feo->q<dim>(), *feo->fe_rt<dim>(),
    		update_values | update_gradients);
    FEValues<dim,3> fe_values(*feo->mapping<dim>(), *feo->q<dim>(), *feo->fe<dim>(),
    		update_values | update_gradients | update_JxW_values | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
    unsigned int dof_indices[ndofs];
    vector<arma::vec3> velocity(qsize);
    vector<arma::vec> sources_sigma(qsize);
    PetscScalar local_matrix[ndofs*ndofs];

	// assemble integral over elements
    for (int i_cell=0; i_cell<feo->dh()->el_ds()->lsize(); i_cell++)
    {
    	typename DOFHandlerBase::CellIterator cell = mesh_->element(feo->dh()->el_index(i_cell));
        if (cell->dim() != dim) continue;

        fe_values.reinit(cell);
        fv_rt.reinit(cell);
        ElementAccessor<3> ele_acc = cell->element_accessor();
        feo->dh()->get_dof_indices(cell, dof_indices);

        calculate_velocity(cell, velocity, fv_rt);
        Model::compute_advection_diffusion_coefficients(fe_values.point_list(), velocity, ele_acc, ad_coef, dif_coef);
        Model::compute_sources_sigma(fe_values.point_list(), ele_acc, sources_sigma);

        // assemble the local stiffness matrix
        for (int sbi=0; sbi<n_subst_; sbi++)
        {
        	for (unsigned int i=0; i<ndofs; i++)
        		for (unsigned int j=0; j<ndofs; j++)
        			local_matrix[i*ndofs+j] = 0;

        	for (unsigned int k=0; k<qsize; k++)
        	{
        		for (unsigned int i=0; i<ndofs; i++)
        		{
        			arma::vec3 Kt_grad_i = dif_coef[sbi][k].t()*fe_values.shape_grad(i,k);
        			double ad_dot_grad_i = arma::dot(ad_coef[sbi][k], fe_values.shape_grad(i,k));

        			for (unsigned int j=0; j<ndofs; j++)
						local_matrix[i*ndofs+j] += (arma::dot(Kt_grad_i, fe_values.shape_grad(j,k))
												   -fe_values.shape_value(j,k)*ad_dot_grad_i
												   +sources_sigma[k][sbi]*fe_values.shape_value(j,k)*fe_values.shape_value(i,k))*fe_values.JxW(k);
				}
			}
			ls[sbi]->mat_set_values(ndofs, (int *)dof_indices, ndofs, (int *)dof_indices, local_matrix);
        }
    }
}


template<class Model>
void TransportDG<Model>::set_sources()
{
  START_TIMER("assemble_sources");
	set_sources<1>();
	set_sources<2>();
	set_sources<3>();
  END_TIMER("assemble_sources");
}

template<class Model>
template<unsigned int dim>
void TransportDG<Model>::set_sources()
{
    FEValues<dim,3> fe_values(*feo->mapping<dim>(), *feo->q<dim>(), *feo->fe<dim>(),
    		update_values | update_JxW_values | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
    vector<arma::vec> sources_conc(qsize), sources_density(qsize), sources_sigma(qsize);
    unsigned int dof_indices[ndofs];
    PetscScalar local_rhs[ndofs];
    double source;

	// assemble integral over elements
    for (int i_cell=0; i_cell<feo->dh()->el_ds()->lsize(); i_cell++)
    {
    	typename DOFHandlerBase::CellIterator cell = mesh_->element(feo->dh()->el_index(i_cell));
        if (cell->dim() != dim) continue;

        fe_values.reinit(cell);
        feo->dh()->get_dof_indices(cell, dof_indices);

        Model::compute_source_coefficients(fe_values.point_list(), cell->element_accessor(), sources_conc, sources_density, sources_sigma);

        // assemble the local stiffness matrix
        for (int sbi=0; sbi<n_subst_; sbi++)
        {
        	for (unsigned int i=0; i<ndofs; i++)
        		local_rhs[i] = 0;

        	// compute sources
        	for (unsigned int k=0; k<qsize; k++)
        	{
        		source = (sources_density[k][sbi] + sources_conc[k][sbi]*sources_sigma[k][sbi])*fe_values.JxW(k);

        		for (unsigned int i=0; i<ndofs; i++)
        			local_rhs[i] += source*fe_values.shape_value(i,k);
            }
        	ls[sbi]->rhs_set_values(ndofs, (int *)dof_indices, local_rhs);
        }
    }
}



template<class Model>
template<unsigned int dim>
void TransportDG<Model>::assemble_fluxes_element_element()
{
    vector<FESideValues<dim,3>*> fe_values;
    FESideValues<dim,3> fsv_rt(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe_rt<dim>(),
    		update_values);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim-1>()->size(),
    		n_max_sides = ad_coef_edg.size();
    vector<unsigned int*> side_dof_indices;
    PetscScalar local_matrix[ndofs*ndofs];
    vector<vector<arma::vec3> > side_velocity(n_max_sides);
    vector<arma::vec> dg_penalty(n_max_sides);
    double gamma_l, omega[2], transport_flux;

    for (int sid=0; sid<n_max_sides; sid++)
    {
    	side_dof_indices.push_back(new unsigned int[ndofs]);
    	fe_values.push_back(new FESideValues<dim,3>(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe<dim>(),
    			update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points));
    }

    // assemble integral over sides
    for (int iedg=0; iedg<feo->dh()->n_loc_edges(); iedg++)
    {
    	Edge *edg = &mesh_->edges[feo->dh()->edge_index(iedg)];
        if (edg->n_sides < 2 || edg->side(0)->element()->dim() != dim) continue;

		for (int sid=0; sid<edg->n_sides; sid++)
		{
			typename DOFHandlerBase::CellIterator cell = edg->side(sid)->element();
			ElementAccessor<3> ele_acc = cell->element_accessor();
			feo->dh()->get_dof_indices(cell, side_dof_indices[sid]);
			fe_values[sid]->reinit(cell, edg->side(sid)->el_idx());
			fsv_rt.reinit(cell, edg->side(sid)->el_idx());
			calculate_velocity(cell, side_velocity[sid], fsv_rt);
			Model::compute_advection_diffusion_coefficients(fe_values[sid]->point_list(), side_velocity[sid], ele_acc, ad_coef_edg[sid], dif_coef_edg[sid]);
			dg_penalty[sid] = data_.dg_penalty.value(cell->centre(), ele_acc);
		}

        // fluxes and penalty
		for (int sbi=0; sbi<n_subst_; sbi++)
		{
			vector<double> fluxes(edg->n_sides);
			for (unsigned int sid=0; sid<edg->n_sides; sid++)
			{
				fluxes[sid] = 0;
				for (unsigned int k=0; k<qsize; k++)
					fluxes[sid] += arma::dot(ad_coef_edg[sid][sbi][k], fe_values[sid]->normal_vector(k))*fe_values[sid]->JxW(k);
				fluxes[sid] /= edg->side(sid)->measure();
			}

			for (int s1=0; s1<edg->n_sides; s1++)
			{
				for (int s2=s1+1; s2<edg->n_sides; s2++)
				{
					ASSERT(edg->side(s1)->valid(), "Invalid side of edge.");
					if (!feo->dh()->el_is_local(edg->side(s1)->element().index())
							&& !feo->dh()->el_is_local(edg->side(s2)->element().index())) continue;

					arma::vec3 nv = fe_values[s1]->normal_vector(0);

					// set up the parameters for DG method
					set_DG_parameters_edge(*edg, s1, s2, qsize, dif_coef_edg[s1][sbi], dif_coef_edg[s2][sbi], fluxes, fe_values[0]->normal_vector(0), dg_penalty[s1][sbi], dg_penalty[s2][sbi], gamma_l, omega, transport_flux);

					int sd[2];
					sd[0] = s1;
					sd[1] = s2;

#define AVERAGE(i,k,side_id)  (fe_values[sd[side_id]]->shape_value(i,k)*0.5)
#define WAVERAGE(i,k,side_id) (arma::dot(dif_coef_edg[sd[side_id]][sbi][k]*fe_values[sd[side_id]]->shape_grad(i,k),nv)*omega[side_id])
#define JUMP(i,k,side_id)     ((side_id==0?1:-1)*fe_values[sd[side_id]]->shape_value(i,k))

					// For selected pair of elements:
					for (int n=0; n<2; n++)
					{
						if (!feo->dh()->el_is_local(edg->side(sd[n])->element().index())) continue;

						for (int m=0; m<2; m++)
						{
							for (unsigned int i=0; i<fe_values[sd[n]]->n_dofs(); i++)
								for (unsigned int j=0; j<fe_values[sd[m]]->n_dofs(); j++)
									local_matrix[i*fe_values[sd[m]]->n_dofs()+j] = 0;

							for (unsigned int k=0; k<qsize; k++)
							{
								double flux_times_JxW = transport_flux*fe_values[0]->JxW(k);
								double gamma_times_JxW = gamma_l*fe_values[0]->JxW(k);

								for (unsigned int i=0; i<fe_values[sd[n]]->n_dofs(); i++)
								{
									double flux_JxW_jump_i = flux_times_JxW*JUMP(i,k,n);
									double gamma_JxW_jump_i = gamma_times_JxW*JUMP(i,k,n);
									double JxW_jump_i = fe_values[0]->JxW(k)*JUMP(i,k,n);
									double JxW_var_wavg_i = fe_values[0]->JxW(k)*WAVERAGE(i,k,n)*dg_variant;

									for (unsigned int j=0; j<fe_values[sd[m]]->n_dofs(); j++)
									{
										int index = i*fe_values[sd[m]]->n_dofs()+j;

										// flux due to transport (applied on interior edges) (average times jump)
										local_matrix[index] += flux_JxW_jump_i*AVERAGE(j,k,m);

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

    for (unsigned int i=0; i<n_max_sides; i++)
    {
        delete fe_values[i];
        delete[] side_dof_indices[i];
    }
}


template<class Model>
template<unsigned int dim>
void TransportDG<Model>::assemble_fluxes_boundary()
{
    FESideValues<dim,3> fe_values_side(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe<dim>(),
    		update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
    FESideValues<dim,3> fsv_rt(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe_rt<dim>(),
    		update_values);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim-1>()->size();
    unsigned int side_dof_indices[ndofs];
    PetscScalar local_matrix[ndofs*ndofs];
    vector<arma::vec3> side_velocity;
    vector<arma::vec> robin_sigma(qsize);
    arma::vec dg_penalty;
    double gamma_l;

    // assemble boundary integral
    for (int iedg=0; iedg<feo->dh()->n_loc_edges(); iedg++)
    {
    	Edge *edg = &mesh_->edges[feo->dh()->edge_index(iedg)];
    	if (edg->n_sides > 1) continue;
    	// check spatial dimension
    	if (edg->side(0)->dim() != dim-1) continue;
    	// skip edges lying not on the boundary
    	if (edg->side(0)->cond() == NULL) continue;

    	SideIter side = edg->side(0);
        typename DOFHandlerBase::CellIterator cell = side->element();
        ElementAccessor<3> ele_acc = cell->element_accessor();
        feo->dh()->get_dof_indices(cell, side_dof_indices);
        fe_values_side.reinit(cell, side->el_idx());
        fsv_rt.reinit(cell, side->el_idx());

        calculate_velocity(cell, side_velocity, fsv_rt);
        Model::compute_advection_diffusion_coefficients(fe_values_side.point_list(), side_velocity, ele_acc, ad_coef, dif_coef);
        dg_penalty = data_.dg_penalty.value(cell->centre(), ele_acc);
        arma::uvec bc_type = data_.bc_type.value(side->cond()->element()->centre(), side->cond()->element_accessor());
        data_.bc_robin_sigma.value_list(fe_values_side.point_list(), side->cond()->element_accessor(), robin_sigma);

        for (int sbi=0; sbi<n_subst_; sbi++)
        {
        	for (unsigned int i=0; i<ndofs; i++)
        		for (unsigned int j=0; j<ndofs; j++)
        			local_matrix[i*ndofs+j] = 0;

			// On Neumann boundaries we have only term from integrating by parts the advective term,
			// on Dirichlet boundaries we additionally apply the penalty which enforces the prescribed value.
			double side_flux = 0;
			for (unsigned int k=0; k<qsize; k++)
				side_flux += arma::dot(ad_coef[sbi][k], fe_values_side.normal_vector(k))*fe_values_side.JxW(k);
			double transport_flux = side_flux/side->measure();

			if ((bc_type[sbi] == EqData::dirichlet) || (bc_type[sbi] == EqData::inflow ))
			{
				// set up the parameters for DG method
				set_DG_parameters_boundary(side, qsize, dif_coef[sbi], transport_flux, fe_values_side.normal_vector(0), dg_penalty[sbi], gamma_l);
				gamma[sbi][side->cond_idx()] = gamma_l;
				if (bc_type[sbi] == EqData::dirichlet || side_flux < -mh_dh->precision())
					transport_flux += gamma_l;
			}

			// fluxes and penalty
			for (unsigned int k=0; k<qsize; k++)
			{
				double flux_times_JxW;
				if (bc_type[sbi] == EqData::robin)
					flux_times_JxW = (transport_flux + robin_sigma[k][sbi])*fe_values_side.JxW(k);
				else
					flux_times_JxW = transport_flux*fe_values_side.JxW(k);

				for (unsigned int i=0; i<ndofs; i++)
					for (unsigned int j=0; j<ndofs; j++)
					{
					    // flux due to advection and penalty
						local_matrix[i*ndofs+j] += flux_times_JxW*fe_values_side.shape_value(i,k)*fe_values_side.shape_value(j,k);

						// flux due to diffusion (only on dirichlet and inflow boundary)
						if (bc_type[sbi] == EqData::dirichlet || (bc_type[sbi] == EqData::inflow && side_flux < -mh_dh->precision()))
							local_matrix[i*ndofs+j] -= (arma::dot(dif_coef[sbi][k]*fe_values_side.shape_grad(j,k),fe_values_side.normal_vector(k))*fe_values_side.shape_value(i,k)
									+ arma::dot(dif_coef[sbi][k]*fe_values_side.shape_grad(i,k),fe_values_side.normal_vector(k))*fe_values_side.shape_value(j,k)*dg_variant
									)*fe_values_side.JxW(k);
					}
			}
			ls[sbi]->mat_set_values(ndofs, (int *)side_dof_indices, ndofs, (int *)side_dof_indices, local_matrix);
        }
    }
}


template<class Model>
template<unsigned int dim>
void TransportDG<Model>::assemble_fluxes_element_side()
{

	if (dim == 1) return;
    FEValues<dim-1,3> fe_values_vb(*feo->mapping<dim-1>(), *feo->q<dim-1>(), *feo->fe<dim-1>(),
    		update_values | update_gradients | update_JxW_values | update_quadrature_points);
    FESideValues<dim,3> fe_values_side(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe<dim>(),
    		update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
    FESideValues<dim,3> fsv_rt(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe_rt<dim>(),
       		update_values);
    FEValues<dim-1,3> fv_rt(*feo->mapping<dim-1>(), *feo->q<dim-1>(), *feo->fe_rt<dim-1>(),
       		update_values);

    vector<FEValuesSpaceBase<3>*> fv_sb(2);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs();    // number of local dofs
    const unsigned int qsize = feo->q<dim-1>()->size();     // number of quadrature points
    unsigned int side_dof_indices[2*ndofs], n_dofs[2];
	vector<arma::vec3> velocity_higher, velocity_lower;
	vector<arma::vec> sigma_c(qsize);
	vector<double> csection_lower(qsize), csection_higher(qsize), mm_coef_lower(qsize), mm_coef_higher(qsize);
    PetscScalar local_matrix[4*ndofs*ndofs];
    double comm_flux[2][2];

    // index 0 = element with lower dimension,
    // index 1 = side of element with higher dimension
    fv_sb[0] = &fe_values_vb;
    fv_sb[1] = &fe_values_side;

    // assemble integral over sides
    for (int inb=0; inb<feo->dh()->n_loc_nb(); inb++)
    {
    	Neighbour *nb = &mesh_->vb_neighbours_[feo->dh()->nb_index(inb)];
        // skip neighbours of different dimension
        if (nb->element()->dim() != dim-1) continue;

		typename DOFHandlerBase::CellIterator cell_sub = mesh_->element.full_iter(nb->element());
		feo->dh()->get_dof_indices(cell_sub, side_dof_indices);
		fe_values_vb.reinit(cell_sub);
		n_dofs[0] = fv_sb[0]->n_dofs();

		typename DOFHandlerBase::CellIterator cell = nb->side()->element();
		feo->dh()->get_dof_indices(cell, side_dof_indices+n_dofs[0]);
		fe_values_side.reinit(cell, nb->side()->el_idx());
		n_dofs[1] = fv_sb[1]->n_dofs();

		// Element id's for testing if they belong to local partition.
		int element_id[2];
		element_id[0] = cell_sub.index();
		element_id[1] = cell.index();

		fsv_rt.reinit(cell, nb->side()->el_idx());
		fv_rt.reinit(cell_sub);
		calculate_velocity(cell, velocity_higher, fsv_rt);
		calculate_velocity(cell_sub, velocity_lower, fv_rt);
		Model::compute_advection_diffusion_coefficients(fe_values_vb.point_list(), velocity_lower, cell_sub->element_accessor(), ad_coef_edg[0], dif_coef_edg[0]);
		Model::compute_advection_diffusion_coefficients(fe_values_vb.point_list(), velocity_higher, cell->element_accessor(), ad_coef_edg[1], dif_coef_edg[1]);
		Model::compute_mass_matrix_coefficient(fe_values_vb.point_list(), cell_sub->element_accessor(), mm_coef_lower);
		Model::compute_mass_matrix_coefficient(fe_values_vb.point_list(), cell->element_accessor(), mm_coef_higher);
		data_.cross_section->value_list(fe_values_vb.point_list(), cell_sub->element_accessor(), csection_lower);
		data_.cross_section->value_list(fe_values_vb.point_list(), cell->element_accessor(), csection_higher);
		data_.sigma_c.value_list(fe_values_vb.point_list(), cell_sub->element_accessor(), sigma_c);

		for (int sbi=0; sbi<n_subst_; sbi++)
		{
			for (unsigned int i=0; i<n_dofs[0]+n_dofs[1]; i++)
				for (unsigned int j=0; j<n_dofs[0]+n_dofs[1]; j++)
					local_matrix[i*(n_dofs[0]+n_dofs[1])+j] = 0;

			// set transmission conditions
			for (unsigned int k=0; k<qsize; k++)
			{
				/* The communication flux has two parts:
				 * - "diffusive" term containing sigma
				 * - "advective" term representing usual upwind
				 *
				 * The calculation differs from the reference manual, since ad_coef and dif_coef have different meaning
				 * than b and A in the manual.
				 * In calculation of sigma there appears one more csection_lower in the denominator.
				 */
				double sigma = sigma_c[k][sbi]*arma::dot(dif_coef_edg[0][sbi][k]*fe_values_side.normal_vector(k),fe_values_side.normal_vector(k))*
				        2*csection_higher[k]*csection_higher[k]/(csection_lower[k]*csection_lower[k]);
				double transport_flux = arma::dot(ad_coef_edg[1][sbi][k], fe_values_side.normal_vector(k));
				double por_lower_over_higher = mm_coef_lower[k]*csection_higher[k]/(mm_coef_higher[k]*csection_lower[k]);

				comm_flux[0][0] =  (sigma-min(0.,transport_flux*por_lower_over_higher))*fv_sb[0]->JxW(k);
				comm_flux[0][1] = -(sigma-min(0.,transport_flux*por_lower_over_higher))*fv_sb[0]->JxW(k);
				comm_flux[1][0] = -(sigma+max(0.,transport_flux))*fv_sb[0]->JxW(k);
				comm_flux[1][1] =  (sigma+max(0.,transport_flux))*fv_sb[0]->JxW(k);

				for (int n=0; n<2; n++)
				{
					if (!feo->dh()->el_is_local(element_id[n])) continue;

					for (unsigned int i=0; i<n_dofs[n]; i++)
						for (int m=0; m<2; m++)
							for (unsigned int j=0; j<n_dofs[m]; j++)
								local_matrix[(i+n*n_dofs[0])*(n_dofs[0]+n_dofs[1]) + m*n_dofs[0] + j] +=
								        comm_flux[m][n]*fv_sb[m]->shape_value(j,k)*fv_sb[n]->shape_value(i,k);
				}
			}
			ls[sbi]->mat_set_values(n_dofs[0]+n_dofs[1], (int *)side_dof_indices, n_dofs[0]+n_dofs[1], (int *)side_dof_indices, local_matrix);
		}
    }

}






template<class Model>
void TransportDG<Model>::set_boundary_conditions()
{
  START_TIMER("assemble_bc");
	set_boundary_conditions<1>();
	set_boundary_conditions<2>();
	set_boundary_conditions<3>();
  END_TIMER("assemble_bc");
}


template<class Model>
template<unsigned int dim>
void TransportDG<Model>::set_boundary_conditions()
{
    FESideValues<dim,3> fe_values_side(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe<dim>(),
    		update_values | update_gradients | update_normal_vectors | update_side_JxW_values | update_quadrature_points);
    FESideValues<dim,3> fsv_rt(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe_rt<dim>(),
           		update_values);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim-1>()->size();
    unsigned int side_dof_indices[ndofs];
    double local_rhs[ndofs];
    vector<arma::vec> bc_values(qsize), bc_fluxes(qsize), bc_sigma(qsize);
	vector<arma::vec3> velocity;

    for (int i=0; i<qsize; i++)
    {
    	bc_values[i].resize(n_subst_);
    	bc_fluxes[i].resize(n_subst_);
    	bc_sigma[i].resize(n_subst_);
    }

    for (int iedg=0; iedg<feo->dh()->n_loc_edges(); iedg++)
    {
    	Edge *edg = &mesh_->edges[feo->dh()->edge_index(iedg)];
    	if (edg->n_sides > 1) continue;
    	if (edg->side(0)->dim() != dim-1) continue;
    	// skip edges lying not on the boundary
    	if (edg->side(0)->cond() == NULL) continue;

    	SideIter side = edg->side(0);
    	typename DOFHandlerBase::CellIterator cell = mesh().element.full_iter(side->element());
        ElementAccessor<3> ele_acc = side->cond()->element_accessor();

        arma::uvec bc_type = data_.bc_type.value(side->cond()->element()->centre(), ele_acc);

        fe_values_side.reinit(cell, side->el_idx());
		fsv_rt.reinit(cell, side->el_idx());
		calculate_velocity(cell, velocity, fsv_rt);

        Model::compute_advection_diffusion_coefficients(fe_values_side.point_list(), velocity, side->element()->element_accessor(), ad_coef, dif_coef);
        Model::compute_dirichlet_bc(fe_values_side.point_list(), ele_acc, bc_values);
        data_.bc_flux.value_list(fe_values_side.point_list(), ele_acc, bc_fluxes);
        data_.bc_robin_sigma.value_list(fe_values_side.point_list(), ele_acc, bc_sigma);

        feo->dh()->get_dof_indices(cell, side_dof_indices);

        for (int sbi=0; sbi<n_subst_; sbi++)
        {
        	for (unsigned int i=0; i<ndofs; i++) local_rhs[i] = 0;

        	for (unsigned int k=0; k<qsize; k++)
        	{
        		double bc_term = 0;
        		arma::vec3 bc_grad;
        		bc_grad.zeros();
                if ((bc_type[sbi] == EqData::inflow && mh_dh->side_flux( *edg->side(0) ) < -mh_dh->precision())
                		|| (bc_type[sbi] == EqData::dirichlet))
                {
                	bc_term = gamma[sbi][side->cond_idx()]*bc_values[k][sbi]*fe_values_side.JxW(k);
                	bc_grad = -bc_values[k][sbi]*fe_values_side.JxW(k)*dg_variant*(arma::trans(dif_coef[sbi][k])*fe_values_side.normal_vector(k));
                }
                else if (bc_type[sbi] == EqData::neumann)
                {
                	bc_term = bc_fluxes[k][sbi]*fe_values_side.JxW(k);
                }
                else if (bc_type[sbi] == EqData::robin)
                {
                	bc_term = bc_sigma[k][sbi]*bc_values[k][sbi]*fe_values_side.JxW(k);
                }

        		for (unsigned int i=0; i<ndofs; i++)
					local_rhs[i] += bc_term*fe_values_side.shape_value(i,k)
							+ arma::dot(bc_grad,fe_values_side.shape_grad(i,k));
			}
			ls[sbi]->rhs_set_values(ndofs, (int *)side_dof_indices, local_rhs);
        }
    }
}



// TODO: The detection of side number from SideIter
// in TransportDG::calculate_velocity()
// should be done with help of class RefElement. This however requires
// that the MH_DofHandler uses the node/side ordering defined in
// the respective RefElement.
template<class Model>
template<unsigned int dim>
void TransportDG<Model>::calculate_velocity(const typename DOFHandlerBase::CellIterator &cell, vector<arma::vec3> &velocity, FEValuesBase<dim,3> &fv)
{
    std::map<const Node*, int> node_nums;
    for (unsigned int i=0; i<cell->n_nodes(); i++)
        node_nums[cell->node[i]] = i;

    velocity.resize(fv.n_points());

    for (unsigned int k=0; k<fv.n_points(); k++)
    {
        velocity[k].zeros();
        for (unsigned int sid=0; sid<cell->n_sides(); sid++)
        {
            if (cell->side(sid)->dim() != dim-1) continue;
            int num = dim*(dim+1)/2;
            for (unsigned int i=0; i<cell->side(sid)->n_nodes(); i++)
                num -= node_nums[cell->side(sid)->node(i)];
            velocity[k] += fv.shape_vector(num,k) * mh_dh->side_flux( *(cell->side(sid)) );
        }
    }
}





template<class Model>
void TransportDG<Model>::set_DG_parameters_edge(const Edge &edg,
            const int s1,
            const int s2,
            const int K_size,
            const vector<arma::mat33> &K1,
            const vector<arma::mat33> &K2,
            const vector<double> &fluxes,
            const arma::vec3 &normal_vector,
            const double alpha1,
            const double alpha2,
            double &gamma,
            double *omega,
            double &transport_flux)
{
    double delta[2];
    double h = 0;
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
        for (unsigned int i=0; i<s->n_nodes(); i++)
            for (unsigned int j=i+1; j<s->n_nodes(); j++)
                h = max(h, s->node(i)->distance(*s->node(j)));
    }

    // calculate the total in- and out-flux through the edge
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
        transport_flux = fluxes[s1]*fabs(fluxes[s2]/pflux);
    else if (fluxes[s2] < 0 && fluxes[s1] > 0 && s1 < s2)
        transport_flux = fluxes[s1]*fabs(fluxes[s2]/nflux);
    else if (s1==s2)
        transport_flux = fluxes[s1];
    else
        transport_flux = 0;

    gamma = 0.5*fabs(transport_flux);


    // determine local DG penalty
    local_alpha = max(alpha1, alpha2);

    if (s1 == s2)
    {
        omega[0] = 1;

        // delta is set to the average value of Kn.n on the side
        delta[0] = 0;
        for (unsigned int k=0; k<K_size; k++)
            delta[0] += dot(K1[k]*normal_vector,normal_vector);
        delta[0] /= K_size;

        gamma += local_alpha/h*(delta[0]);
    }
    else
    {
        delta[0] = 0;
        delta[1] = 0;
        for (unsigned int k=0; k<K_size; k++)
        {
            delta[0] += dot(K1[k]*normal_vector,normal_vector);
            delta[1] += dot(K2[k]*normal_vector,normal_vector);
        }
        delta[0] /= K_size;
        delta[1] /= K_size;

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






template<class Model>
void TransportDG<Model>::set_DG_parameters_boundary(const SideIter side,
		    const int K_size,
            const vector<arma::mat33> &K,
            const double flux,
            const arma::vec3 &normal_vector,
            const double alpha,
            double &gamma)
{
    double delta = 0, h = 0;

    // calculate the side diameter
    if (side->dim() == 0)
    {
        h = 1;
    }
    else
    {
        for (unsigned int i=0; i<side->n_nodes(); i++)
            for (unsigned int j=i+1; j<side->n_nodes(); j++)
                h = max(h, side->node(i)->distance( *side->node(j) ));
    }

	// delta is set to the average value of Kn.n on the side
	for (unsigned int k=0; k<K_size; k++)
		delta += dot(K[k]*normal_vector,normal_vector);
	delta /= K_size;

	gamma = 0.5*fabs(flux) + alpha/h*delta;
}





template<class Model>
void TransportDG<Model>::set_initial_condition()
{
	START_TIMER("set_init_cond");
	for (int sbi=0; sbi<n_subst_; sbi++)
		ls[sbi]->start_allocation();
	prepare_initial_condition<1>();
	prepare_initial_condition<2>();
	prepare_initial_condition<3>();

	for (int sbi=0; sbi<n_subst_; sbi++)
		ls[sbi]->start_add_assembly();
	prepare_initial_condition<1>();
	prepare_initial_condition<2>();
	prepare_initial_condition<3>();

	for (int sbi=0; sbi<n_subst_; sbi++)
	{
		ls[sbi]->finish_assembly();
		ls[sbi]->solve();
	}
	END_TIMER("set_init_cond");
}

template<class Model>
template<unsigned int dim>
void TransportDG<Model>::prepare_initial_condition()
{
	FEValues<dim,3> fe_values(*feo->mapping<dim>(), *feo->q<dim>(), *feo->fe<dim>(),
			update_values | update_JxW_values | update_quadrature_points);
    const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
    unsigned int dof_indices[ndofs];
    double matrix[ndofs*ndofs], rhs[ndofs];
    std::vector<arma::vec> init_values(qsize);

    for (int i=0; i<qsize; i++)
    	init_values[i].resize(n_subst_);

    for (int i_cell=0; i_cell<feo->dh()->el_ds()->lsize(); i_cell++)
    {
    	typename DOFHandlerBase::CellIterator elem = mesh_->element(feo->dh()->el_index(i_cell));
    	if (elem->dim() != dim) continue;

    	ElementAccessor<3> ele_acc = elem->element_accessor();
    	feo->dh()->get_dof_indices(elem, dof_indices);
    	fe_values.reinit(elem);

   		Model::compute_init_cond(fe_values.point_list(), ele_acc, init_values);

    	for (int sbi=0; sbi<n_subst_; sbi++)
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
    				for (unsigned int j=0; j<ndofs; j++)
    					matrix[i*ndofs+j] += fe_values.shape_value(i,k)*fe_values.shape_value(j,k)*fe_values.JxW(k);

    				rhs[i] += fe_values.shape_value(i,k)*rhs_term;
    			}
    		}
    		ls[sbi]->set_values(ndofs, (int *)dof_indices, ndofs, (int *)dof_indices, matrix, rhs);
    	}
    }
}


template<class Model>
void TransportDG<Model>::calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance)
{
	calc_fluxes<1>(bcd_balance, bcd_plus_balance, bcd_minus_balance);
	calc_fluxes<2>(bcd_balance, bcd_plus_balance, bcd_minus_balance);
	calc_fluxes<3>(bcd_balance, bcd_plus_balance, bcd_minus_balance);
}

template<class Model>
template<unsigned int dim>
void TransportDG<Model>::calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance)
{
	FESideValues<dim,3> fe_values(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe<dim>(),
			update_values | update_gradients | update_normal_vectors | update_side_JxW_values | update_quadrature_points);
	FESideValues<dim,3> fsv_rt(*feo->mapping<dim>(), *feo->q<dim-1>(), *feo->fe_rt<dim>(),
			update_values);
	const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim-1>()->size();
	unsigned int dof_indices[ndofs];
	vector<arma::vec3> side_velocity(qsize);
    vector<arma::vec> bc_values(qsize);
	arma::vec3 conc_grad;

    for (int i=0; i<qsize; i++)
    	bc_values[i].resize(n_subst_);

    for (int iedg=0; iedg<feo->dh()->n_loc_edges(); iedg++)
    {
    	Edge *edg = &mesh_->edges[feo->dh()->edge_index(iedg)];
    	if (edg->n_sides > 1) continue;
    	if (edg->side(0)->dim() != dim-1) continue;
    	// skip edges lying not on the boundary
    	if (edg->side(0)->cond() == NULL) continue;

    	SideIter side = edg->side(0);
    	DOFHandlerMultiDim::CellIterator cell = side->element();
        ElementAccessor<3> ele_acc = cell->element_accessor();
        Region r = side->cond()->element_accessor().region();
        if (! r.is_valid()) xprintf(Msg, "Invalid region, ele % d, edg: % d\n", side->cond()->bc_ele_idx_, side->cond()->edge_idx_);
        unsigned int bc_region_idx = r.boundary_idx();

		fe_values.reinit(cell, side->el_idx());
		fsv_rt.reinit(cell, side->el_idx());
		feo->dh()->get_dof_indices(cell, dof_indices);
		calculate_velocity(cell, side_velocity, fsv_rt);
		Model::compute_advection_diffusion_coefficients(fe_values.point_list(), side_velocity, ele_acc, ad_coef, dif_coef);
		Model::compute_dirichlet_bc(fe_values.point_list(), side->cond()->element_accessor(), bc_values);
		arma::uvec bc_type = data_.bc_type.value(side->cond()->element()->centre(), side->cond()->element_accessor());

		for (int sbi=0; sbi<n_subst_; sbi++)
		{
			double mass_flux = 0;
			double water_flux = 0;
			for (unsigned int k=0; k<qsize; k++)
				water_flux += arma::dot(ad_coef[sbi][k], fe_values.normal_vector(k))*fe_values.JxW(k);
			water_flux /= side->measure();

			for (unsigned int k=0; k<qsize; k++)
			{
				double conc = 0;
				conc_grad.zeros();
				for (unsigned int i=0; i<ndofs; i++)
				{
					conc += fe_values.shape_value(i,k)*ls[sbi]->get_solution_array()[dof_indices[i]-feo->dh()->loffset()];
					conc_grad += fe_values.shape_grad(i,k)*ls[sbi]->get_solution_array()[dof_indices[i]-feo->dh()->loffset()];
				}

				// flux due to advection
				mass_flux += water_flux*conc*fe_values.JxW(k);

				if (bc_type[sbi] == EqData::dirichlet || (bc_type[sbi] == EqData::inflow && water_flux*side->measure() < -mh_dh->precision()))
				{
					// flux due to diffusion
					mass_flux -= arma::dot(dif_coef[sbi][k]*conc_grad,fe_values.normal_vector(k))*fe_values.JxW(k);

					// the penalty term has to be added otherwise the mass balance will not hold
					mass_flux -= gamma[sbi][side->cond_idx()]*(bc_values[k][sbi] - conc)*fe_values.JxW(k);
				}

			}
			bcd_balance[sbi][bc_region_idx] += mass_flux;
			if (mass_flux >= 0) bcd_plus_balance[sbi][bc_region_idx] += mass_flux;
			else bcd_minus_balance[sbi][bc_region_idx] += mass_flux;
		}
    }

}

template<class Model>
void TransportDG<Model>::calc_elem_sources(vector<vector<double> > &mass, vector<vector<double> > &src_balance)
{
	calc_elem_sources<1>(mass, src_balance);
	calc_elem_sources<2>(mass, src_balance);
	calc_elem_sources<3>(mass, src_balance);
}

template<class Model>
template<unsigned int dim>
void TransportDG<Model>::calc_elem_sources(vector<vector<double> > &mass, vector<vector<double> > &src_balance)
{
	FEValues<dim,3> fe_values(*feo->mapping<dim>(), *feo->q<dim>(), *feo->fe<dim>(),
			update_values | update_JxW_values | update_quadrature_points);
	const unsigned int ndofs = feo->fe<dim>()->n_dofs(), qsize = feo->q<dim>()->size();
	unsigned int dof_indices[ndofs];
	vector<arma::vec> sources_conc(qsize), sources_density(qsize), sources_sigma(qsize);
	double mass_sum, sources_sum, conc, conc_diff;

    for (int i_cell=0; i_cell<feo->dh()->el_ds()->lsize(); i_cell++)
    {
    	typename DOFHandlerBase::CellIterator elem = mesh_->element(feo->dh()->el_index(i_cell));
		if (elem->dim() != dim) continue;

		Region r = elem->element_accessor().region();
		if (! r.is_valid()) xprintf(Msg, "Invalid region, ele % d\n", elem.index());
		unsigned int region_idx = r.bulk_idx();

		ElementAccessor<3> ele_acc = elem->element_accessor();
		fe_values.reinit(elem);
		feo->dh()->get_dof_indices(elem, dof_indices);

		Model::compute_mass_matrix_coefficient(fe_values.point_list(), ele_acc, mm_coef);
		Model::compute_source_coefficients(fe_values.point_list(), ele_acc, sources_conc, sources_density, sources_sigma);

		for (int sbi=0; sbi<n_subst_; sbi++)
		{
			mass_sum = 0;
			sources_sum = 0;

			for (unsigned int k=0; k<qsize; k++)
			{
				conc = 0;
				for (unsigned int i=0; i<ndofs; i++)
					conc += fe_values.shape_value(i,k)*ls[sbi]->get_solution_array()[dof_indices[i]-feo->dh()->loffset()];

				mass_sum += mm_coef[k]*conc*fe_values.JxW(k);

				conc_diff = sources_conc[k][sbi] - conc;
				sources_sum += (sources_density[k][sbi] + conc_diff*sources_sigma[k][sbi])*fe_values.JxW(k);
			}

			mass[sbi][region_idx] += mass_sum;
			src_balance[sbi][region_idx] += sources_sum;
		}
	}
}





template class TransportDG<ConcentrationTransportModel>;
template class TransportDG<HeatTransferModel>;




