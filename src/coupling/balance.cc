/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    balance.cc
 * @ingroup transport
 * @brief   Mass balance
 */

#include <iostream>
#include <iomanip>

#include "system/system.hh"
#include "system/sys_profiler.hh"

#include <petscmat.h>
#include "mesh/mesh.h"
#include "io/output_time_set.hh"
#include "coupling/balance.hh"
#include "fields/unit_si.hh"
#include "tools/time_governor.hh"

using namespace Input::Type;

const Selection & Balance::get_format_selection_input_type() {
	return Selection("Balance_output_format", "Format of output file for balance.")
		.add_value(Balance::legacy, "legacy", "Legacy format used by previous program versions.")
		.add_value(Balance::txt, "txt", "Excel format with tab delimiter.")
		.add_value(Balance::gnuplot, "gnuplot", "Format compatible with GnuPlot datafile with fixed column width.")
		.close();
}

const Record & Balance::get_input_type() {
	return Record("Balance", "Balance of a conservative quantity, boundary fluxes and sources.")
	    .declare_key("times", OutputTimeSet::get_input_type(), Default("[]"), "" )
	    .declare_key("add_output_times", Bool(), Default("true"), "Add all output times of the balanced equation to the balance output times set. "
	            "Note that this is not the time set of the output stream.")
		//.declare_key("balance_on", Bool(), Default("true"), "Balance is computed if the value is true.")
		.declare_key("format", Balance::get_format_selection_input_type(), Default("\"txt\""), "Format of output file.")
		.declare_key("cumulative", Bool(), Default("false"), "Compute cumulative balance over time. "
				"If true, then balance is calculated at each computational time step, which can slow down the program.")
		.declare_key("file", FileName::output(), Default::read_time("File name generated from the balanced quantity: <quantity_name>_balance.*"), "File name for output of balance.")
		.close();
}

/*
std::shared_ptr<Balance> Balance::make_balance(
        const std::string &file_prefix,
        const Mesh *mesh,
        const Input::Record &in_rec,
        TimeGovernor &tg)
{
    auto ptr = make_shared<Balance>(file_prefix, mesh, in_rec, tg);
    auto &marks = TimeGovernor::marks();
    auto balance_output_type = tg.equation_mark_type() | TimeGovernor::marks().type_balance_output();
    if (marks.begin(balance_output_type) == marks.end(balance_output_type) ) {
        // no balance output time => force balance off
        ptr.reset();
    }
    return ptr;
}*/



Balance::Balance(const std::string &file_prefix, const Mesh *mesh)
	: 	  file_prefix_(file_prefix),
	      mesh_(mesh),
	  	  initial_time_(),
	  	  last_time_(),
	  	  initial_(true),
	  	  allocation_done_(false),
	  	  output_line_counter_(0),
          balance_on_(true),
          output_yaml_header_(false)

{
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
    ASSERT_PTR(mesh_);
}



Balance::~Balance()
{
	if (rank_ == 0) {
		output_.close();
		output_yaml_.close();
	}
	if (! allocation_done_) return;

	for (unsigned int c=0; c<quantities_.size(); ++c)
	{
		chkerr(MatDestroy(&(region_mass_matrix_[c])));
		chkerr(MatDestroy(&(be_flux_matrix_[c])));
		chkerr(MatDestroy(&(region_source_matrix_[c])));
		chkerr(MatDestroy(&(region_source_rhs_[c])));
		chkerr(VecDestroy(&(be_flux_vec_[c])));
		chkerr(VecDestroy(&(region_source_vec_[c])));
	}
	delete[] region_mass_matrix_;
	delete[] be_flux_matrix_;
	delete[] be_flux_vec_;
	delete[] region_source_matrix_;
	delete[] region_source_rhs_;
	delete[] region_source_vec_;

	chkerr(MatDestroy(&region_be_matrix_));
	chkerr(VecDestroy(&ones_));
	chkerr(VecDestroy(&ones_be_));
}


void Balance::init_from_input(
        const Input::Record &in_rec,
        TimeGovernor &tg)
{
    ASSERT(! allocation_done_);

    auto &marks = TimeGovernor::marks();

    output_mark_type_ = tg.equation_mark_type() | marks.type_output(),
    balance_output_type_ = tg.equation_fixed_mark_type() | marks.type_balance_output();

    cumulative_ = in_rec.val<bool>("cumulative");
    output_format_ = in_rec.val<OutputFormat>("format");

    OutputTimeSet time_set;
    time_set.read_from_input( in_rec.val<Input::Array>("times"), tg, balance_output_type_);
    add_output_times_ = in_rec.val<bool>("add_output_times");

    if (rank_ == 0) {
        // set default value by output_format_
        std::string default_file_name;
        switch (output_format_)
        {
        case txt:
            default_file_name = file_prefix_ + "_balance.txt";
            break;
        case gnuplot:
            default_file_name = file_prefix_ + "_balance.dat";
            break;
        case legacy:
            default_file_name = file_prefix_ + "_balance.txt";
            break;
        }



        balance_output_file_ = in_rec.val<FilePath>("file", FilePath(default_file_name, FilePath::output_file));
    }

}


void Balance::units(const UnitSI &unit)
{
    ASSERT(! allocation_done_);

    units_ = unit;
    units_.undef(false);
}

unsigned int Balance::add_quantity(const string &name)
{
    ASSERT(! allocation_done_);

	Quantity q(quantities_.size(), name);
	quantities_.push_back(q);

	return q.index_;
}


std::vector<unsigned int> Balance::add_quantities(const std::vector<string> &names)
{
    vector<unsigned int> indices;
	for (auto name : names)
		indices.push_back(add_quantity(name));
	return indices;
}


void Balance::allocate(unsigned int n_loc_dofs,
		unsigned int max_dofs_per_boundary)
{
    ASSERT(! allocation_done_);
    n_loc_dofs_ = n_loc_dofs;
    max_dofs_per_boundary_ = max_dofs_per_boundary;
}

void Balance::lazy_initialize()
{
    if (allocation_done_) return;

    auto &marks = TimeGovernor::marks();
    if (add_output_times_)
        marks.add_to_type_all(output_mark_type_, balance_output_type_);
    // if there are no balance marks turn balance off
    if (marks.begin(balance_output_type_) == marks.end(balance_output_type_) )
    {
        balance_on_ = false;
        cumulative_ = false;
        return;
    }

	// Max. number of regions to which a single dof can contribute.
	// TODO: estimate or compute this number directly (from mesh or dof handler).
	const int n_bulk_regs_per_dof = min(10, (int)mesh_->region_db().bulk_size());
	const unsigned int n_quant = quantities_.size();
	const unsigned int n_bdr_reg = mesh_->region_db().boundary_size();
	const unsigned int n_blk_reg = mesh_->region_db().bulk_size();


	// construct vector of regions of boundary edges
    for (unsigned int loc_el = 0; loc_el < mesh_->get_el_ds()->lsize(); loc_el++)
    {
        Element *elm = mesh_->element(mesh_->get_el_4_loc()[loc_el]);
        if (elm->boundary_idx_ != nullptr)
        {
            FOR_ELEMENT_SIDES(elm,si)
            {
                Boundary *b = elm->side(si)->cond();
                if (b != nullptr)
                	be_regions_.push_back(b->region().boundary_idx());
            }
        }
    }



	fluxes_    .resize(n_quant, vector<double>(n_bdr_reg, 0));
	fluxes_in_ .resize(n_quant, vector<double>(n_bdr_reg, 0));
	fluxes_out_.resize(n_quant, vector<double>(n_bdr_reg, 0));

	masses_     .resize(n_quant, vector<double>(n_blk_reg, 0));
	sources_    .resize(n_quant, vector<double>(n_blk_reg, 0));
	sources_in_ .resize(n_quant, vector<double>(n_blk_reg, 0));
	sources_out_.resize(n_quant, vector<double>(n_blk_reg, 0));

	sum_fluxes_     .resize(n_quant, 0);
	sum_fluxes_in_  .resize(n_quant, 0);
	sum_fluxes_out_ .resize(n_quant, 0);
	sum_masses_     .resize(n_quant, 0);
	sum_sources_    .resize(n_quant, 0);
	sum_sources_in_ .resize(n_quant, 0);
	sum_sources_out_.resize(n_quant, 0);

	if (cumulative_)
	{
		initial_mass_      .resize(n_quant, 0);
		integrated_fluxes_ .resize(n_quant, 0);
		integrated_sources_.resize(n_quant, 0);
		increment_sources_.resize(n_quant, 0);
		increment_fluxes_.resize(n_quant, 0);
	}



	region_mass_matrix_ = new Mat[n_quant];
	be_flux_matrix_ = new Mat[n_quant];
	region_source_matrix_ = new Mat[n_quant];
	region_source_rhs_ = new Mat[n_quant];
    region_mass_vec_ = new Vec[n_quant];
	be_flux_vec_ = new Vec[n_quant];
	region_source_vec_ = new Vec[n_quant];

	for (unsigned int c=0; c<n_quant; ++c)
	{
		chkerr(MatCreateAIJ(PETSC_COMM_WORLD,
				n_loc_dofs_,
				(rank_==0)?mesh_->region_db().bulk_size():0,
				PETSC_DECIDE,
				PETSC_DECIDE,
				(rank_==0)?n_bulk_regs_per_dof:0,
				0,
				(rank_==0)?0:n_bulk_regs_per_dof,
				0,
				&(region_mass_matrix_[c])));

        chkerr(MatCreateAIJ(PETSC_COMM_WORLD,
               n_loc_dofs_,
               (rank_==0)?mesh_->region_db().bulk_size():0,
               PETSC_DECIDE,
               PETSC_DECIDE,
               (rank_==0)?n_bulk_regs_per_dof:0,
               0,
               (rank_==0)?0:n_bulk_regs_per_dof,
               0,
               &(region_source_matrix_[c])));

       chkerr(MatCreateAIJ(PETSC_COMM_WORLD,
               n_loc_dofs_,
               (rank_==0)?mesh_->region_db().bulk_size():0,
               PETSC_DECIDE,
               PETSC_DECIDE,
               (rank_==0)?n_bulk_regs_per_dof:0,
               0,
               (rank_==0)?0:n_bulk_regs_per_dof,
               0,
               &(region_source_rhs_[c])));

       chkerr(MatCreateAIJ(PETSC_COMM_WORLD,
               be_regions_.size(),
               n_loc_dofs_,
               PETSC_DECIDE,
               PETSC_DECIDE,
               max_dofs_per_boundary_,
               0,
               0,
               0,
               &(be_flux_matrix_[c])));
        
        chkerr(VecCreateMPI(PETSC_COMM_WORLD,
                (rank_==0)?mesh_->region_db().bulk_size():0,
                PETSC_DECIDE,
                &(region_mass_vec_[c])));

		chkerr(VecCreateMPI(PETSC_COMM_WORLD,
				be_regions_.size(),
				PETSC_DECIDE,
				&(be_flux_vec_[c])));

		chkerr(VecCreateMPI(PETSC_COMM_WORLD,
				(rank_==0)?mesh_->region_db().bulk_size():0,
				PETSC_DECIDE,
				&(region_source_vec_[c])));
	}

	chkerr(MatCreateAIJ(PETSC_COMM_WORLD,
			be_regions_.size(),
			(rank_==0)?mesh_->region_db().boundary_size():0,
			PETSC_DECIDE,
			PETSC_DECIDE,
			(rank_==0)?1:0,
			0,
			(rank_==0)?0:1,
			0,
			&region_be_matrix_
			));
	chkerr(VecGetOwnershipRange(be_flux_vec_[0], &be_offset_, NULL));
	for (unsigned int loc_el=0; loc_el<be_regions_.size(); ++loc_el)
	{
	    chkerr(MatSetValue(region_be_matrix_,
				be_offset_+loc_el,
				be_regions_[loc_el],
				1,
				INSERT_VALUES));
	}
	chkerr(MatAssemblyBegin(region_be_matrix_, MAT_FINAL_ASSEMBLY));
	chkerr(MatAssemblyEnd(region_be_matrix_, MAT_FINAL_ASSEMBLY));

	double *ones_array;
	chkerr(VecCreateMPI(PETSC_COMM_WORLD,
			n_loc_dofs_,
			PETSC_DECIDE,
			&ones_));
	chkerr(VecGetArray(ones_, &ones_array));
	fill_n(ones_array, n_loc_dofs_, 1);
	chkerr(VecRestoreArray(ones_, &ones_array));

	chkerr(VecCreateMPI(PETSC_COMM_WORLD,
			be_regions_.size(),
			PETSC_DECIDE,
			&ones_be_));
	chkerr(VecGetArray(ones_be_, &ones_array));
	fill_n(ones_array, be_regions_.size(), 1);
	chkerr(VecRestoreArray(ones_be_, &ones_array));


    if (rank_ == 0) {
        balance_output_file_.open_stream(output_);
        // set file name of YAML output
        string yaml_file_name = file_prefix_ + "_balance.yaml";
        FilePath(yaml_file_name, FilePath::output_file).open_stream(output_yaml_);
    }

    allocation_done_ = true;
}



bool Balance::is_current(const TimeStep &step)
{
    if (! balance_on_) return false;

    auto &marks = TimeGovernor::marks();
    bool res =  marks.current(step, balance_output_type_) != marks.end(balance_output_type_);

    //cout << "flag: " << res << " time: " << step.end() << " type:" << hex << balance_output_type_ << endl;
    return res;
}

void Balance::start_mass_assembly(unsigned int quantity_idx)
{
    lazy_initialize();
    if (! balance_on_) return;
	chkerr(MatZeroEntries(region_mass_matrix_[quantity_idx]));
    chkerr(VecZeroEntries(region_mass_vec_[quantity_idx]));
}


void Balance::start_flux_assembly(unsigned int quantity_idx)
{
    lazy_initialize();
    if (! balance_on_) return;
	chkerr(MatZeroEntries(be_flux_matrix_[quantity_idx]));
	chkerr(VecZeroEntries(be_flux_vec_[quantity_idx]));
}


void Balance::start_source_assembly(unsigned int quantity_idx)
{
    lazy_initialize();
    if (! balance_on_) return;
	chkerr(MatZeroEntries(region_source_matrix_[quantity_idx]));
	chkerr(MatZeroEntries(region_source_rhs_[quantity_idx]));
	chkerr(VecZeroEntries(region_source_vec_[quantity_idx]));
}


void Balance::finish_mass_assembly(unsigned int quantity_idx)
{
	ASSERT(allocation_done_);
    if (! balance_on_) return;

	chkerr(MatAssemblyBegin(region_mass_matrix_[quantity_idx], MAT_FINAL_ASSEMBLY));
	chkerr(MatAssemblyEnd(region_mass_matrix_[quantity_idx], MAT_FINAL_ASSEMBLY));
    chkerr(VecAssemblyBegin(region_mass_vec_[quantity_idx]));
    chkerr(VecAssemblyEnd(region_mass_vec_[quantity_idx]));
}

void Balance::finish_flux_assembly(unsigned int quantity_idx)
{
    ASSERT(allocation_done_);
    if (! balance_on_) return;

    chkerr(MatAssemblyBegin(be_flux_matrix_[quantity_idx], MAT_FINAL_ASSEMBLY));
	chkerr(MatAssemblyEnd(be_flux_matrix_[quantity_idx], MAT_FINAL_ASSEMBLY));
	chkerr(VecAssemblyBegin(be_flux_vec_[quantity_idx]));
	chkerr(VecAssemblyEnd(be_flux_vec_[quantity_idx]));
}

void Balance::finish_source_assembly(unsigned int quantity_idx)
{
    ASSERT(allocation_done_);
    if (! balance_on_) return;

    chkerr(MatAssemblyBegin(region_source_matrix_[quantity_idx], MAT_FINAL_ASSEMBLY));
	chkerr(MatAssemblyEnd(region_source_matrix_[quantity_idx], MAT_FINAL_ASSEMBLY));
	chkerr(MatAssemblyBegin(region_source_rhs_[quantity_idx], MAT_FINAL_ASSEMBLY));
	chkerr(MatAssemblyEnd(region_source_rhs_[quantity_idx], MAT_FINAL_ASSEMBLY));
	chkerr(MatMultTranspose(region_source_rhs_[quantity_idx], ones_, region_source_vec_[quantity_idx]));
}




void Balance::add_mass_matrix_values(unsigned int quantity_idx,
		unsigned int region_idx,
		const vector<int> &dof_indices,
		const vector<double> &values)
{
    ASSERT_DBG(allocation_done_);
    if (! balance_on_) return;

	PetscInt reg_array[1] = { (int)region_idx };

	chkerr_assert(MatSetValues(region_mass_matrix_[quantity_idx],
			dof_indices.size(),
			&(dof_indices[0]),
			1,
			reg_array,
			&(values[0]),
			ADD_VALUES));
}


void Balance::add_flux_matrix_values(unsigned int quantity_idx,
		unsigned int elem_idx,
		const vector<int> &dof_indices,
		const vector<double> &values)
{
    ASSERT_DBG(allocation_done_);
    if (! balance_on_) return;

	PetscInt elem_array[1] = { int(be_offset_+elem_idx) };
	chkerr_assert(MatSetValues(be_flux_matrix_[quantity_idx],
			1,
			elem_array,
			dof_indices.size(),
			&(dof_indices[0]),
			&(values[0]),
			ADD_VALUES));
}


void Balance::add_source_matrix_values(unsigned int quantity_idx,
		unsigned int region_idx,
		const vector<int> &dof_indices,
		const vector<double> &values)
{
    ASSERT_DBG(allocation_done_);
    if (! balance_on_) return;

	PetscInt reg_array[1] = { (int)region_idx };

	chkerr_assert(MatSetValues(region_source_matrix_[quantity_idx],
			dof_indices.size(),
			&(dof_indices[0]),
			1,
			reg_array,
			&(values[0]),
			ADD_VALUES));
}


void Balance::add_mass_vec_value(unsigned int quantity_idx,
        unsigned int region_idx,
        double value)
{
  chkerr_assert(VecSetValue(region_mass_vec_[quantity_idx],
            region_idx,
            value,
            ADD_VALUES));
}


void Balance::add_flux_vec_value(unsigned int quantity_idx,
		unsigned int elem_idx,
		double value)
{
    ASSERT_DBG(allocation_done_);
    if (! balance_on_) return;

    chkerr_assert(VecSetValue(be_flux_vec_[quantity_idx],
			be_offset_+elem_idx,
			value,
			ADD_VALUES));
}


void Balance::add_source_vec_values(unsigned int quantity_idx,
		unsigned int region_idx,
		const vector<int> &dof_indices,
		const vector<double> &values)
{
    ASSERT_DBG(allocation_done_);
    if (! balance_on_) return;

	PetscInt reg_array[1] = { (int)region_idx };

	chkerr_assert(MatSetValues(region_source_rhs_[quantity_idx],
			dof_indices.size(),
			&(dof_indices[0]),
			1,
			reg_array,
			&(values[0]),
			ADD_VALUES));
}


void Balance::add_cumulative_source(unsigned int quantity_idx, double source)
{
    ASSERT_DBG(allocation_done_);
    if (!cumulative_) return;

    if (rank_ == 0)
        increment_sources_[quantity_idx] += source;
}


void Balance::calculate_cumulative_sources(unsigned int quantity_idx,
		const Vec &solution,
		double dt)
{
    ASSERT_DBG(allocation_done_);
	if (!cumulative_) return;

	Vec bulk_vec;

	chkerr(VecCreateMPIWithArray(PETSC_COMM_WORLD,
			1,
			(rank_==0)?mesh_->region_db().bulk_size():0,
			PETSC_DECIDE,
			&(sources_[quantity_idx][0]),
			&bulk_vec));

	// compute sources on bulk regions: S'.u + s
	chkerr(VecZeroEntries(bulk_vec));
	chkerr(MatMultTransposeAdd(region_source_matrix_[quantity_idx],
	        solution, region_source_vec_[quantity_idx], bulk_vec));

	double sum_sources;
	chkerr(VecSum(bulk_vec, &sum_sources));
	VecDestroy(&bulk_vec);

	if (rank_ == 0)
		// sum sources in one step
		increment_sources_[quantity_idx] += sum_sources*dt;
}


void Balance::calculate_cumulative_fluxes(unsigned int quantity_idx,
		const Vec &solution,
		double dt)
{
    ASSERT_DBG(allocation_done_);
    if (!cumulative_) return;

	Vec boundary_vec;

	chkerr(VecCreateMPIWithArray(PETSC_COMM_WORLD,
			1,
			(rank_==0)?mesh_->region_db().boundary_size():0,
			PETSC_DECIDE,
			&(fluxes_[quantity_idx][0]),
			&boundary_vec));

	// compute fluxes on boundary regions: R'.(F.u + f)
	chkerr(VecZeroEntries(boundary_vec));
	Vec temp;
	chkerr(VecDuplicate(ones_be_, &temp));
	chkerr(MatMultAdd(be_flux_matrix_[quantity_idx], solution, be_flux_vec_[quantity_idx], temp));
	// Since internally we keep outgoing fluxes, we change sign
	// to write to output _incoming_ fluxes.
	chkerr(VecScale(temp, -1));
	chkerr(MatMultTranspose(region_be_matrix_, temp, boundary_vec));
	chkerr(VecDestroy(&temp));

	double sum_fluxes;
	chkerr(VecSum(boundary_vec, &sum_fluxes));
	chkerr(VecDestroy(&boundary_vec));

	if (rank_ == 0)
		// sum fluxes in one step
		increment_fluxes_[quantity_idx] += sum_fluxes*dt;
}


void Balance::calculate_mass(unsigned int quantity_idx,
		const Vec &solution,
		vector<double> &output_array)
{
    ASSERT_DBG(allocation_done_);
    if (! balance_on_) return;

    Vec bulk_vec;

	chkerr(VecCreateMPIWithArray(PETSC_COMM_WORLD,
			1,
			(rank_==0)?mesh_->region_db().bulk_size():0,
			PETSC_DECIDE,
			&(output_array[0]),
			&bulk_vec));

	// compute mass on regions: M'.u
	chkerr(VecZeroEntries(bulk_vec));
	chkerr(MatMultTransposeAdd(region_mass_matrix_[quantity_idx], 
                               solution, 
                               region_mass_vec_[quantity_idx], 
                               bulk_vec));
	VecDestroy(&bulk_vec);
}


void Balance::calculate_source(unsigned int quantity_idx,
		const Vec &solution)
{
    ASSERT_DBG(allocation_done_);
    if (! balance_on_) return;

    Vec bulk_vec;

	chkerr(VecCreateMPIWithArray(PETSC_COMM_WORLD,
			1,
			(rank_==0)?mesh_->region_db().bulk_size():0,
			PETSC_DECIDE,
			&(sources_[quantity_idx][0]),
			&bulk_vec));

	// compute sources on bulk regions: S'.u + s
	chkerr(VecZeroEntries(bulk_vec));
	chkerr(MatMultTransposeAdd(region_source_matrix_[quantity_idx],
			solution,
			region_source_vec_[quantity_idx],
			bulk_vec));

	// compute positive/negative sources
	int lsize;
	Vec mat_r, rhs_r;
	const double *sol_array, *mat_array, *rhs_array;
	chkerr(VecGetLocalSize(solution, &lsize));
	chkerr(VecDuplicate(solution, &mat_r));
	chkerr(VecDuplicate(solution, &rhs_r));
	chkerr(VecGetArrayRead(solution, &sol_array));
	for (unsigned int r=0; r<mesh_->region_db().bulk_size(); ++r)
	{
		MatGetColumnVector(region_source_matrix_[quantity_idx], mat_r, r);
		MatGetColumnVector(region_source_rhs_[quantity_idx], rhs_r, r);

		VecGetArrayRead(mat_r, &mat_array);
		VecGetArrayRead(rhs_r, &rhs_array);

		sources_in_[quantity_idx][r] = 0;
		sources_out_[quantity_idx][r] = 0;
		for (int i=0; i<lsize; ++i)
		{
			double f = mat_array[i]*sol_array[i] + rhs_array[i];
			if (f > 0) sources_in_[quantity_idx][r] += f;
			else sources_out_[quantity_idx][r] += f;
		}

		VecRestoreArrayRead(mat_r, &mat_array);
		VecRestoreArrayRead(rhs_r, &rhs_array);
	}
	chkerr(VecRestoreArrayRead(solution, &sol_array));
	VecDestroy(&rhs_r);
	VecDestroy(&mat_r);
	VecDestroy(&bulk_vec);
}


void Balance::calculate_flux(unsigned int quantity_idx,
		const Vec &solution)
{
    ASSERT_DBG(allocation_done_);
    if (! balance_on_) return;

    Vec boundary_vec;

	chkerr(VecCreateMPIWithArray(PETSC_COMM_WORLD, 1,
	        (rank_==0)?mesh_->region_db().boundary_size():0,
	        PETSC_DECIDE, &(fluxes_[quantity_idx][0]), &boundary_vec));

	// compute fluxes on boundary regions: R'.(F.u + f)
	chkerr(VecZeroEntries(boundary_vec));
	Vec temp;
	chkerr(VecDuplicate(ones_be_, &temp));
	chkerr(MatMultAdd(be_flux_matrix_[quantity_idx], solution, be_flux_vec_[quantity_idx], temp));
	// Since internally we keep outgoing fluxes, we change sign
	// to write to output _incoming_ fluxes.
	chkerr(VecScale(temp, -1));
	chkerr(MatMultTranspose(region_be_matrix_, temp, boundary_vec));

	// compute positive/negative fluxes
	fluxes_in_[quantity_idx].assign(mesh_->region_db().boundary_size(), 0);
	fluxes_out_[quantity_idx].assign(mesh_->region_db().boundary_size(), 0);
	const double *flux_array;
	int lsize;
	chkerr(VecGetArrayRead(temp, &flux_array));
	chkerr(VecGetLocalSize(temp, &lsize));
	for (int e=0; e<lsize; ++e)
	{
		if (flux_array[e] < 0)
			fluxes_out_[quantity_idx][be_regions_[e]] += flux_array[e];
		else
			fluxes_in_[quantity_idx][be_regions_[e]] += flux_array[e];
	}
	chkerr(VecRestoreArrayRead(temp, &flux_array));
	VecDestroy(&temp);
	VecDestroy(&boundary_vec);
}





void Balance::output(double time)
{
    ASSERT_DBG(allocation_done_);
    if (! balance_on_) return;

	// gather results from processes and sum them up
	const unsigned int n_quant = quantities_.size();
	const unsigned int n_blk_reg = mesh_->region_db().bulk_size();
	const unsigned int n_bdr_reg = mesh_->region_db().boundary_size();
	const int buf_size = n_quant*2*n_blk_reg + n_quant*2*n_bdr_reg;
	double sendbuffer[buf_size], recvbuffer[buf_size];
	for (unsigned int qi=0; qi<n_quant; qi++)
	{
		for (unsigned int ri=0; ri<n_blk_reg; ri++)
		{
			sendbuffer[qi*2*n_blk_reg +           + ri] = sources_in_[qi][ri];
			sendbuffer[qi*2*n_blk_reg + n_blk_reg + ri] = sources_out_[qi][ri];
		}
		for (unsigned int ri=0; ri<n_bdr_reg; ri++)
		{
			sendbuffer[n_quant*2*n_blk_reg + qi*2*n_bdr_reg +           + ri] = fluxes_in_[qi][ri];
			sendbuffer[n_quant*2*n_blk_reg + qi*2*n_bdr_reg + n_bdr_reg + ri] = fluxes_out_[qi][ri];
		}
	}
	MPI_Reduce(&sendbuffer,recvbuffer,buf_size,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);
	// for other than 0th process update last_time and finish,
	// on process #0 sum balances over all regions and calculate
	// cumulative balance over time.
	if (rank_ == 0)
	{
		// update balance vectors
		for (unsigned int qi=0; qi<n_quant; qi++)
		{
			for (unsigned int ri=0; ri<n_blk_reg; ri++)
			{
				sources_in_[qi][ri]  = recvbuffer[qi*2*n_blk_reg +           + ri];
				sources_out_[qi][ri] = recvbuffer[qi*2*n_blk_reg + n_blk_reg + ri];
			}
			for (unsigned int ri=0; ri<n_bdr_reg; ri++)
			{
				fluxes_in_[qi][ri]  = recvbuffer[n_quant*2*n_blk_reg + qi*2*n_bdr_reg +           + ri];
				fluxes_out_[qi][ri] = recvbuffer[n_quant*2*n_blk_reg + qi*2*n_bdr_reg + n_bdr_reg + ri];
			}
		}
	}


	// The convention for input/output of fluxes is that positive means inward.
	// Therefore in the following code we switch sign of fluxes.
	if (rank_ == 0)
	{
		sum_fluxes_.assign(n_quant, 0);
		sum_fluxes_in_.assign(n_quant, 0);
		sum_fluxes_out_.assign(n_quant, 0);
		sum_masses_.assign(n_quant, 0);
		sum_sources_.assign(n_quant, 0);
		sum_sources_in_.assign(n_quant, 0);
		sum_sources_out_.assign(n_quant, 0);

		// sum all boundary fluxes
		const RegionSet & b_set = mesh_->region_db().get_region_set(".BOUNDARY");
		for( RegionSet::const_iterator reg = b_set.begin(); reg != b_set.end(); ++reg)
		{
			for (unsigned int qi=0; qi<n_quant; qi++)
			{
				sum_fluxes_[qi]     += fluxes_    [qi][reg->boundary_idx()];
				sum_fluxes_in_[qi]  += fluxes_in_ [qi][reg->boundary_idx()];
				sum_fluxes_out_[qi] += fluxes_out_[qi][reg->boundary_idx()];
			}
		}

		// sum all volume sources
		const RegionSet & bulk_set = mesh_->region_db().get_region_set("BULK");
		for( RegionSet::const_iterator reg = bulk_set.begin(); reg != bulk_set.end(); ++reg)
		{
			for (unsigned int qi=0; qi<n_quant; qi++)
			{
				sum_masses_[qi] += masses_[qi][reg->bulk_idx()];
				sum_sources_[qi] += sources_[qi][reg->bulk_idx()];
				sum_sources_in_[qi] += sources_in_[qi][reg->bulk_idx()];
				sum_sources_out_[qi] += sources_out_[qi][reg->bulk_idx()];
			}
		}

		// cumulative balance over time
		if (cumulative_)
		{
			// save initial time and mass
			if (initial_)
			{
				initial_time_ = time;
				last_time_ = initial_time_;
				for (unsigned int qi=0; qi<n_quant; qi++)
					initial_mass_[qi] = sum_masses_[qi];
				initial_ = false;
			}

			for (unsigned int qi=0; qi<n_quant; qi++)
			{
				integrated_fluxes_[qi] += increment_fluxes_[qi];
				integrated_sources_[qi] += increment_sources_[qi];
			}
		}
	}

	last_time_ = time;


	// perform actual output
	switch (output_format_)
	{
	case txt:
		output_csv(time, '\t', "");
		break;
	case gnuplot:
		output_csv(time, ' ', "#", 30);
		break;
	case legacy:
		output_legacy(time);
		break;
	}
	// output in YAML format
	output_yaml(time);

	if (rank_ == 0)
	{
		sum_fluxes_.assign(n_quant, 0);
		sum_sources_.assign(n_quant, 0);
		increment_fluxes_.assign(n_quant, 0);
		increment_sources_.assign(n_quant, 0);
	}
}


void Balance::output_legacy(double time)
{
	// write output only on process #0
	if (rank_ != 0) return;

	const unsigned int n_quant = quantities_.size();

	// print the head of mass balance file
	unsigned int c = 6; //column number without label
	unsigned int w = 14;  //column width
	unsigned int wl = 2*(w-5)+7;  //label column width
	string bc_head_format = "# %-*s%-*s%-*s%-*s%-*s%-*s\n",
		   bc_format = "%*s%-*d%-*s%-*s%-*g%-*g%-*g\n",
		   bc_total_format = "# %-*s%-*s%-*g%-*g%-*g\n";

	output_ << "# " << setw((w*c+wl-14)/2) << setfill('-') << "--"
			<< " MASS BALANCE "
	     	<< setw((w*c+wl-14)/2) << setfill('-') << "" << endl
			<< "# Time: " << time << "\n\n\n";

	// header for table of boundary fluxes
	output_ << "# Mass flux through boundary [M/T]:\n# "
			<< setiosflags(ios::left) << setfill(' ')
	        << setw(w)  << "[boundary_id]"
	        << setw(wl) << "[label]"
	        << setw(w)  << "[substance]"
	        << setw(w)  << "[total flux]"
	        << setw(w)  << "[outward flux]"
	        << setw(w)  << "[inward flux]"
	        << endl;

	// draw long line
	output_ << "# " << setw(w*c+wl) << setfill('-') << "" << setfill(' ') << endl;

	// print mass fluxes over boundaries
	const RegionSet & b_set = mesh_->region_db().get_region_set(".BOUNDARY");
	for( RegionSet::const_iterator reg = b_set.begin(); reg != b_set.end(); ++reg) {
		for (unsigned int qi=0; qi<n_quant; qi++) {
			output_ << setw(2)  << ""
					<< setw(w)  << (int)reg->id()
					<< setw(wl) << reg->label().c_str()
					<< setw(w)  << quantities_[qi].name_.c_str()
					<< setw(w)  << fluxes_[qi][reg->boundary_idx()]
					<< setw(w)  << fluxes_out_[qi][reg->boundary_idx()]
					<< setw(w)  << fluxes_in_[qi][reg->boundary_idx()]
					<< endl;
		}
	}

	// draw long line
	output_ << "# " << setw(w*c+wl) << setfill('-') << "" << setfill(' ') << endl;

	// total boundary balance
	for (unsigned int qi=0; qi<n_quant; qi++)
		output_ << "# " << setiosflags(ios::left)
				<< setw(w+wl) << "Total mass flux of substance [M/T]"
				<< setw(w)    << quantities_[qi].name_.c_str()
				<< setw(w)    << sum_fluxes_[qi]
				<< setw(w)    << sum_fluxes_out_[qi]
				<< setw(w)    << sum_fluxes_in_[qi]
				<< endl;
	output_ << "\n\n";


	// header for table of volume sources and masses
	string src_head_format = "# %-*s%-*s%-*s%-*s%-*s\n",
		   src_format = "%*s%-*d%-*s%-*s%-*g%-*g\n",
		   src_total_format = "# %-*s%-*s%-*g%-*g\n";
	output_ << "# Mass [M] and sources [M/T] on regions:\n"
			<< "# " << setiosflags(ios::left)
			<< setw(w)  << "[region_id]"
			<< setw(wl) << "[label]"
			<< setw(w)  << "[substance]"
			<< setw(w)  << "[total_mass]"
			<< setw(w)  << "[total_source]"
			<< endl;

	// draw long line
	output_ << "# " << setw(w*c+wl) << setfill('-') << "" << setfill(' ') << endl;

	// print  balance of volume sources and masses
	const RegionSet & bulk_set = mesh_->region_db().get_region_set("BULK");
	for( RegionSet::const_iterator reg = bulk_set.begin(); reg != bulk_set.end(); ++reg)
	{
		for (unsigned int qi=0; qi<n_quant; qi++)
		{
			output_ << setw(2)  << ""
					<< setw(w)  << (int)reg->id()
					<< setw(wl) << reg->label().c_str()
					<< setw(w)  << quantities_[qi].name_.c_str()
					<< setw(w)  << masses_[qi][reg->bulk_idx()]
					<< setw(w)  << sources_[qi][reg->bulk_idx()]
					<< endl;
		}
	}

	// draw long line
	output_ << "# " << setw(w*c+wl) << setfill('-') << "" << setfill(' ') << endl;

	// total sources balance
	for (unsigned int qi=0; qi<n_quant; qi++)
		output_ << "# " << setiosflags(ios::left) << setw(w+wl) << "Total mass [M] and sources [M/T]"
				<< setw(w) << quantities_[qi].name_.c_str()
				<< setw(w) << sum_masses_[qi]
				<< setw(w) << sum_sources_[qi]
				<< endl;

	if (cumulative_)
	{
		// Print cumulative sources
		output_ << "\n\n# Cumulative mass balance on time interval ["
				<< setiosflags(ios::left) << initial_time_ << ","
				<< setiosflags(ios::left) << time << "]\n"
				<< "# Initial mass [M] + sources integrated over time [M] - flux integrated over time [M] = current mass [M]\n"
				<< "# " << setiosflags(ios::left)
				<< setw(w) << "[substance]"
				<< setw(w) << "[A=init. mass]"
				<< setw(w) << "[B=source]"
				<< setw(w) << "[C=flux]"
				<< setw(w) << "[A+B-C]"
				<< setw(w) << "[D=curr. mass]"
				<< setw(w) << "[A+B-C-D=err.]"
				<< setw(w) << "[rel. error]"
				<< endl;

		for (unsigned int qi=0; qi<n_quant; qi++)
		{
			double denominator = max(fabs(initial_mass_[qi]+integrated_sources_[qi]-integrated_fluxes_[qi]),fabs(sum_masses_[qi]));
			output_ << "  " << setiosflags(ios::left)
					<< setw(w) << quantities_[qi].name_.c_str()
					<< setw(w) << initial_mass_[qi]
					<< setw(w) << integrated_sources_[qi]
					<< setw(w) << integrated_fluxes_[qi]
					<< setw(w) << initial_mass_[qi]+integrated_sources_[qi]-integrated_fluxes_[qi]
					<< setw(w) << sum_masses_[qi]
					<< setw(w) << initial_mass_[qi]+integrated_sources_[qi]-integrated_fluxes_[qi]-sum_masses_[qi]
					<< setw(w) << fabs(initial_mass_[qi]+integrated_sources_[qi]+integrated_fluxes_[qi]-sum_masses_[qi])/(denominator==0?1:denominator)
					<< endl;
		}
	}

	output_ << endl << endl;
}


std::string Balance::csv_zero_vals(unsigned int cnt, char delimiter)
{
	std::stringstream ss;
	for (unsigned int i=0; i<cnt; i++) ss << format_csv_val(0, delimiter);
	return ss.str();
}


void Balance::output_csv(double time, char delimiter, const std::string& comment_string, unsigned int repeat)
{
	// write output only on process #0
	if (rank_ != 0) return;

	const unsigned int n_quant = quantities_.size();

	// print data header only on first line
	if (repeat==0 && output_line_counter_==0) format_csv_output_header(delimiter, comment_string);

	// print sources and masses over bulk regions
	const RegionSet & bulk_set = mesh_->region_db().get_region_set("BULK");
	for( RegionSet::const_iterator reg = bulk_set.begin(); reg != bulk_set.end(); ++reg)
	{
		for (unsigned int qi=0; qi<n_quant; qi++)
		{
			// print data header (repeat header after every "repeat" lines)
			if (repeat && (output_line_counter_%repeat == 0)) format_csv_output_header(delimiter, comment_string);

			output_ << format_csv_val(time, delimiter, true)
					<< format_csv_val(reg->label(), delimiter)
					<< format_csv_val(quantities_[qi].name_, delimiter)
					<< csv_zero_vals(3, delimiter)
					<< format_csv_val(masses_[qi][reg->bulk_idx()], delimiter)
					<< format_csv_val(sources_[qi][reg->bulk_idx()], delimiter)
					<< format_csv_val(sources_in_[qi][reg->bulk_idx()], delimiter)
					<< format_csv_val(sources_out_[qi][reg->bulk_idx()], delimiter)
					<< csv_zero_vals(5, delimiter) << endl;
			++output_line_counter_;
		}
	}

	// print mass fluxes over boundaries
	const RegionSet & b_set = mesh_->region_db().get_region_set(".BOUNDARY");
	for( RegionSet::const_iterator reg = b_set.begin(); reg != b_set.end(); ++reg)
	{
		for (unsigned int qi=0; qi<n_quant; qi++) {
			// print data header (repeat header after every "repeat" lines)
			if (repeat && (output_line_counter_%repeat == 0)) format_csv_output_header(delimiter, comment_string);

			output_ << format_csv_val(time, delimiter, true)
					<< format_csv_val(reg->label(), delimiter)
					<< format_csv_val(quantities_[qi].name_, delimiter)
					<< format_csv_val(fluxes_[qi][reg->boundary_idx()], delimiter)
					<< format_csv_val(fluxes_in_[qi][reg->boundary_idx()], delimiter)
					<< format_csv_val(fluxes_out_[qi][reg->boundary_idx()], delimiter)
					<< csv_zero_vals(9, delimiter) << endl;
			++output_line_counter_;
		}
	}

	if (cumulative_)
	{
		for (unsigned int qi=0; qi<n_quant; qi++)
		{
			// print data header (repeat header after every "repeat" lines)
			if (repeat && (output_line_counter_%repeat == 0)) format_csv_output_header(delimiter, comment_string);

			double error = sum_masses_[qi] - (initial_mass_[qi] + integrated_sources_[qi] + integrated_fluxes_[qi]);
			output_ << format_csv_val(time, delimiter, true)
					<< format_csv_val("ALL", delimiter)
					<< format_csv_val(quantities_[qi].name_, delimiter)
					<< format_csv_val(sum_fluxes_[qi], delimiter)
					<< format_csv_val(sum_fluxes_in_[qi], delimiter)
					<< format_csv_val(sum_fluxes_out_[qi], delimiter)
					<< format_csv_val(sum_masses_[qi], delimiter)
					<< format_csv_val(sum_sources_[qi], delimiter)
					<< format_csv_val(sum_sources_in_[qi], delimiter)
					<< format_csv_val(sum_sources_out_[qi], delimiter)
					<< format_csv_val(increment_fluxes_[qi], delimiter)
					<< format_csv_val(increment_sources_[qi], delimiter)
					<< format_csv_val(integrated_fluxes_[qi], delimiter)
					<< format_csv_val(integrated_sources_[qi], delimiter)
					<< format_csv_val(error, delimiter) << endl;
			++output_line_counter_;
		}
	}

}


void Balance::format_csv_output_header(char delimiter, const std::string& comment_string)
{
	std::stringstream ss;
	if (delimiter == ' ') {
		ss << setw(output_column_width-comment_string.size()) << "\"time\"";
	} else {
		ss << "\"time\"";
	}

	output_ << comment_string << ss.str()
			<< format_csv_val("region", delimiter)
			<< format_csv_val("quantity [" + units_.format_text() + "]", delimiter)
			<< format_csv_val("flux", delimiter)
			<< format_csv_val("flux_in", delimiter)
			<< format_csv_val("flux_out", delimiter)
			<< format_csv_val("mass", delimiter)
			<< format_csv_val("source", delimiter)
			<< format_csv_val("source_in", delimiter)
			<< format_csv_val("source_out", delimiter)
			<< format_csv_val("flux_increment", delimiter)
			<< format_csv_val("source_increment", delimiter)
			<< format_csv_val("flux_cumulative", delimiter)
			<< format_csv_val("source_cumulative", delimiter)
			<< format_csv_val("error", delimiter)
			<< endl;
}

std::string Balance::format_csv_val(std::string val, char delimiter, bool initial)
{
	std::stringstream ss;
	std::replace( val.begin(), val.end(), '\"', '\'');

	if (!initial) ss << delimiter;
	if (delimiter == ' ') {
		std::stringstream sval;
		sval << "\"" << val << "\"";
		ss << " " << setw(output_column_width-1) << sval.str();
	} else {
		ss << "\"" << val << "\"";
	}

	return ss.str();
}

std::string Balance::format_csv_val(double val, char delimiter, bool initial)
{
	std::stringstream ss;

	if (!initial) ss << delimiter;
	if (delimiter == ' ') {
		ss << " " << setw(output_column_width-1) << val;
	} else {
		ss << val;
	}
	return ss.str();
}

void Balance::output_yaml(double time)
{
	// write output only on process #0
	if (rank_ != 0) return;

	const unsigned int n_quant = quantities_.size();

	// print data header only once
	if (!output_yaml_header_) {
		output_yaml_ << "column_names: [ flux,  flux_in,  flux_out,  mass,  source,  source_in,  source_out,  flux_increment,  "
			<< "source_increment,  flux_cumulative,  source_cumulative,  error ]"  << endl;
		output_yaml_ << "data:" << endl;
		output_yaml_header_ = true;
	}

	output_yaml_ << setfill(' ');

	// print sources and masses over bulk regions
	const RegionSet & bulk_set = mesh_->region_db().get_region_set("BULK");
	for( RegionSet::const_iterator reg = bulk_set.begin(); reg != bulk_set.end(); ++reg)
	{
		for (unsigned int qi=0; qi<n_quant; qi++)
		{
			output_yaml_ << "  - time: " << time << endl;
			output_yaml_ << setw(4) << "" << "region: " << reg->label() << endl;
			output_yaml_ << setw(4) << "" << "quantity: " << quantities_[qi].name_ << endl;
			output_yaml_ << setw(4) << "" << "data: " << "[ 0, 0, 0, " << masses_[qi][reg->bulk_idx()] << ", "
						 << sources_[qi][reg->bulk_idx()] << ", " << sources_in_[qi][reg->bulk_idx()] << ", "
						 << sources_out_[qi][reg->bulk_idx()] << ", 0, 0, 0, 0, 0 ]" << endl;
		}
	}

	// print mass fluxes over boundaries
	const RegionSet & b_set = mesh_->region_db().get_region_set(".BOUNDARY");
	for( RegionSet::const_iterator reg = b_set.begin(); reg != b_set.end(); ++reg)
	{
		for (unsigned int qi=0; qi<n_quant; qi++) {
			output_yaml_ << "  - time: " << time << endl;
			output_yaml_ << setw(4) << "" << "region: " << reg->label() << endl;
			output_yaml_ << setw(4) << "" << "quantity: " << quantities_[qi].name_ << endl;
			output_yaml_ << setw(4) << "" << "data: " << "[ " << fluxes_[qi][reg->boundary_idx()] << ", "
					     << fluxes_in_[qi][reg->boundary_idx()] << ", " << fluxes_out_[qi][reg->boundary_idx()]
					     << ", 0, 0, 0, 0, 0, 0, 0, 0, 0 ]" << endl;
		}
	}

	if (cumulative_)
	{
		for (unsigned int qi=0; qi<n_quant; qi++)
		{
			double error = sum_masses_[qi] - (initial_mass_[qi] + integrated_sources_[qi] + integrated_fluxes_[qi]);
			output_yaml_ << "  - time: " << time << endl;
			output_yaml_ << setw(4) << "" << "region: ALL" << endl;
			output_yaml_ << setw(4) << "" << "quantity: " << quantities_[qi].name_ << endl;
			output_yaml_ << setw(4) << "" << "data: " << "[ " << sum_fluxes_[qi] << ", "
					     << sum_fluxes_in_[qi] << ", " << sum_fluxes_out_[qi] << ", "
						 << sum_masses_[qi] << ", " << sum_sources_[qi] << ", "
						 << sum_sources_in_[qi] << ", " << sum_sources_out_[qi] << ", "
						 << increment_fluxes_[qi] << ", " << increment_sources_[qi] << ", "
						 << integrated_fluxes_[qi] << ", " << integrated_sources_[qi] << ", "
						 << error << " ]" << endl;
		}
	}
}






