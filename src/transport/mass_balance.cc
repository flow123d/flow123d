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
 * @ingroup transport
 * @brief  Mass balance
 *
 *
 */

#include <iostream>
#include <iomanip>

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/xio.h"

#include <petscmat.h>
#include "mesh/mesh.h"

#include "transport/mass_balance.hh"

using namespace Input::Type;

Selection Balance::format_selection_input_type
	= Selection("Balance_output_format", "Format of output file for balance.")
	.add_value(Balance::legacy, "legacy", "Legacy format used by previous program versions.")
	.add_value(Balance::csv, "csv", "Excel format with semicolon delimiter.")
	.add_value(Balance::gnuplot, "gnuplot", "Format compatible with GnuPlot datafile with fix column width.");

Record MassBalance::input_type
	= Record("MassBalance", "Balance of mass, boundary fluxes and sources for transport of substances.")
	.declare_key("mass_balance_on", Bool(), Default("true"), "Balance is computed if the value is true.")
	.declare_key("format", Balance::format_selection_input_type, Default("csv"), "Format of output file.")
	.declare_key("cumulative", Bool(), Default("false"), "Compute cumulative balance over time. "
			"If true, then balance is calculated at each computational time step, which can slow down the program.")
	.declare_key("file", FileName::output(), Default::read_time("FileName mass_balance.*"), "File name for output of mass balance.")
	.allow_auto_conversion("mass_balance_on")
;




MassBalance::MassBalance(EquationForMassBalance *eq, const Input::Record &in_rec)
	: equation_(eq), initial_time(0), last_time(-1), initial(true)
{
	balance_output_file = xfopen( in_rec.val<FilePath>("file"), "wt");

	cumulative = in_rec.val<bool>("cumulative");

    const int n_bcd_reg_ = equation_->region_db()->boundary_size();
    const int n_blk_reg_ = equation_->region_db()->bulk_size();

	bcd_balance.resize(eq->n_substances(), vector<double>(n_bcd_reg_, 0));
	bcd_plus_balance.resize(eq->n_substances(), vector<double>(n_bcd_reg_, 0));
	bcd_minus_balance.resize(eq->n_substances(), vector<double>(n_bcd_reg_, 0));
	mass.resize(eq->n_substances(), vector<double>(n_blk_reg_, 0));
	src_balance.resize(eq->n_substances(), vector<double>(n_blk_reg_, 0));

	bcd_total_balance.resize(eq->n_substances(), 0.);
	bcd_total_inflow.resize(eq->n_substances(), 0.);
	bcd_total_outflow.resize(eq->n_substances(), 0.);
	mass_total.resize(eq->n_substances(), 0.);
	src_total_balance.resize(eq->n_substances(), 0.);

	initial_mass.resize(eq->n_substances(), 0.);
	integrated_sources.resize(eq->n_substances(), 0.);
	integrated_fluxes.resize(eq->n_substances(), 0.);
}

MassBalance::~MassBalance()
{
	if (balance_output_file != NULL) xfclose(balance_output_file);
}


void MassBalance::calculate(double time) {
    // return if we already calculated at the given time
    if (last_time == time) return;

    const int n_bcd_reg_ = equation_->region_db()->boundary_size();
    const int n_blk_reg_ = equation_->region_db()->bulk_size();
    const int n_subst = equation_->n_substances();

    for (int i=0; i<n_subst; i++)
    {
    	for (int j=0; j<n_bcd_reg_; j++)
    	{
    		bcd_balance[i][j] = 0;
    		bcd_plus_balance[i][j] = 0;
    		bcd_minus_balance[i][j] = 0;
    	}
    	for (int j=0; j<n_blk_reg_; j++)
    	{
    		mass[i][j] = 0;
    		src_balance[i][j] = 0;
    	}
		bcd_total_balance[i] = 0;
		bcd_total_outflow[i] = 0;
		bcd_total_inflow[i]  = 0;
		mass_total[i] = 0;
		src_total_balance[i] = 0;
    }

    // compute fluxes, mass and volume sources
    equation_->calc_fluxes(bcd_balance, bcd_plus_balance, bcd_minus_balance);
	equation_->calc_elem_sources(mass, src_balance);


	// gather results from processes and sum them up
	int buf_size = n_subst*(3*n_bcd_reg_ + 2*n_blk_reg_);
	double sendbuffer[buf_size], recvbuffer[buf_size];
	for (int i=0; i<n_subst; i++)
	{
		for (int j=0; j<n_bcd_reg_; j++)
		{
			sendbuffer[i*3*n_bcd_reg_+            j] = bcd_balance[i][j];
			sendbuffer[i*3*n_bcd_reg_+  n_bcd_reg_+j] = bcd_plus_balance[i][j];
			sendbuffer[i*3*n_bcd_reg_+2*n_bcd_reg_+j] = bcd_minus_balance[i][j];
		}
		for (int j=0; j<n_blk_reg_; j++)
		{
			sendbuffer[n_subst*3*n_bcd_reg_+i*2*n_blk_reg_           +j] = mass[i][j];
			sendbuffer[n_subst*3*n_bcd_reg_+i*2*n_blk_reg_+n_blk_reg_+j] = src_balance[i][j];
		}
	}
	MPI_Reduce(&sendbuffer,recvbuffer,buf_size,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD);

	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	// for other than 0th process update last_time and finish,
	// on process #0 sum balances over all regions and calculate
	// cumulative balance over time.
	if (rank != 0)
	{
		last_time = time;
		return;
	}


	// update balance vectors
	for (int i=0; i<n_subst; i++)
	{
		for (int j=0; j<n_bcd_reg_; j++)
		{
			bcd_balance[i][j]       = recvbuffer[i*3*n_bcd_reg_+           j];
			bcd_plus_balance[i][j]  = recvbuffer[i*3*n_bcd_reg_+  n_bcd_reg_+j];
			bcd_minus_balance[i][j] = recvbuffer[i*3*n_bcd_reg_+2*n_bcd_reg_+j];
		}
		for (int j=0; j<n_blk_reg_; j++)
		{
			mass[i][j]        = recvbuffer[n_subst*3*n_bcd_reg_+i*2*n_blk_reg_         +j];
			src_balance[i][j] = recvbuffer[n_subst*3*n_bcd_reg_+i*2*n_blk_reg_+n_blk_reg_+j];
		}
	}


	// sum all boundary fluxes
	const RegionSet & b_set = equation_->region_db()->get_region_set("BOUNDARY");
	for( RegionSet::const_iterator reg = b_set.begin(); reg != b_set.end(); ++reg) {
		for (int sbi=0; sbi<n_subst; sbi++) {
			bcd_total_balance[sbi] += bcd_balance[sbi][reg->boundary_idx()];
			bcd_total_outflow[sbi] += bcd_plus_balance[sbi][reg->boundary_idx()];
			bcd_total_inflow[sbi] += bcd_minus_balance[sbi][reg->boundary_idx()];
		}
	}

	// sum all volume sources
	const RegionSet & bulk_set = equation_->region_db()->get_region_set("BULK");
	for( RegionSet::const_iterator reg = bulk_set.begin(); reg != bulk_set.end(); ++reg)
	{
		for (int sbi=0; sbi<n_subst; sbi++)
		{
			mass_total[sbi] += mass[sbi][reg->bulk_idx()];
			src_total_balance[sbi] += src_balance[sbi][reg->bulk_idx()];
		}
	}

	if (!cumulative) return;

    // cumulative balance over time

	// quantities need for the balance over time interval
	static vector<double> last_sources(n_subst, 0.),
			last_fluxes(n_subst, 0.);

	// save initial time and mass
	if (initial)
	{
		initial_time = time;
		last_time = initial_time;
		for (int i=0; i<n_subst; i++)
			initial_mass[i] = mass_total[i];
		initial = false;
	}


	// sum sources and fluxes according to the time integration scheme
	// TODO: Think if we really need TimeIntegrationScheme or the cummulative
	// quantities can be calculated in the same way for both explicit and implicit
	// methods.
	switch (equation_->time_scheme())
	{
	case EquationForMassBalance::explicit_euler:
//		for (int i=0; i<n_subst; i++)
//		{
//			integrated_sources[i] += last_sources[i]*(time-last_time);
//			integrated_fluxes[i] += last_fluxes[i]*(time-last_time);
//		}
//		break;
	case EquationForMassBalance::implicit_euler:
		for (int i=0; i<n_subst; i++)
		{
			integrated_sources[i] += src_total_balance[i]*(time-last_time);
			integrated_fluxes[i] += bcd_total_balance[i]*(time-last_time);
		}
		break;
	case EquationForMassBalance::crank_nicholson:
		for (int i=0; i<n_subst; i++)
		{
			integrated_sources[i] += (last_sources[i]+src_total_balance[i])*0.5*(time-last_time);
			integrated_fluxes[i] += (last_fluxes[i]+bcd_total_balance[i])*0.5*(time-last_time);
		}
		break;
	default:
		break;
	}

	last_time = time;
	for (int i=0; i<n_subst; i++)
	{
		last_sources[i] = src_total_balance[i];
		last_fluxes[i] = bcd_total_balance[i];
	}
}


void MassBalance::output(double time)
{
	// calculate balances for the given time
	if (last_time != time) calculate(time);
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	// write output only on process #0
	if (rank != 0) return;

	const int n_subst = equation_->n_substances();

	// print the head of mass balance file
	unsigned int c = 6; //column number without label
	unsigned int w = 14;  //column width
	unsigned int wl = 2*(w-5)+7;  //label column width
	stringstream s; //helpful stringstream
	string bc_head_format = "# %-*s%-*s%-*s%-*s%-*s%-*s\n",
		   bc_format = "%*s%-*d%-*s%-*s%-*g%-*g%-*g\n",
		   bc_total_format = "# %-*s%-*s%-*g%-*g%-*g\n";
	s << setw((w*c+wl-14)/2) << setfill('-') << "--"; //drawing half line
	fprintf(balance_output_file,"# %s MASS BALANCE %s\n",s.str().c_str(), s.str().c_str());
	fprintf(balance_output_file,"# Time: %f\n\n\n",time);

	// header for table of boundary fluxes
	fprintf(balance_output_file,"# Mass flux through boundary [M/T]:\n");
	fprintf(balance_output_file,bc_head_format.c_str(),w,"[boundary_id]",wl,"[label]",
							w,"[substance]",w,"[total flux]",w,"[outward flux]",w,"[inward flux]");
	s.clear();
	s.str(std::string());
	s << setw(w*c+wl) << setfill('-') << "-";
	fprintf(balance_output_file,"# %s\n",s.str().c_str());  //drawing long line

	// print mass fluxes over boundaries
	const RegionSet & b_set = equation_->region_db()->get_region_set("BOUNDARY");
	for( RegionSet::const_iterator reg = b_set.begin(); reg != b_set.end(); ++reg) {
		for (int sbi=0; sbi<n_subst; sbi++) {
			fprintf(balance_output_file, bc_format.c_str(),2,"",w,reg->id(),wl,reg->label().c_str(),
					w, equation_->substances().names()[sbi].c_str(),w, bcd_balance[sbi][reg->boundary_idx()],
					w, bcd_plus_balance[sbi][reg->boundary_idx()],
					w, bcd_minus_balance[sbi][reg->boundary_idx()]);
		}
	}
	// total boundary balance
	fprintf(balance_output_file,"# %s\n",s.str().c_str());  // drawing long line
	for (int sbi=0; sbi<n_subst; sbi++)
		fprintf(balance_output_file, bc_total_format.c_str(),w+wl,"Total mass flux of substance [M/T]",
				w,equation_->substances().names()[sbi].c_str(),w,bcd_total_balance[sbi], w, bcd_total_outflow[sbi], w, bcd_total_inflow[sbi]);
	fprintf(balance_output_file, "\n\n");


	// header for table of volume sources and masses
	string src_head_format = "# %-*s%-*s%-*s%-*s%-*s\n",
		   src_format = "%*s%-*d%-*s%-*s%-*g%-*g\n",
		   src_total_format = "# %-*s%-*s%-*g%-*g\n";
	fprintf(balance_output_file,"# Mass [M] and sources [M/T] on regions:\n");   //head
	fprintf(balance_output_file,src_head_format.c_str(),w,"[region_id]",wl,"[label]",
							w,"[substance]",w,"[total_mass]",w,"[total_source]");
	fprintf(balance_output_file,"# %s\n",s.str().c_str());  //long line

	// print  balance of volume sources and masses
	const RegionSet & bulk_set = equation_->region_db()->get_region_set("BULK");
	for( RegionSet::const_iterator reg = bulk_set.begin(); reg != bulk_set.end(); ++reg)
	{
		for (int sbi=0; sbi<n_subst; sbi++)
		{
			fprintf(balance_output_file, src_format.c_str(), 2,"", w, reg->id(), wl,
					reg->label().c_str(), w, equation_->substances().names()[sbi].c_str(),
					w,mass[sbi][reg->bulk_idx()],
					w,src_balance[sbi][reg->bulk_idx()]);
		}
	}
	// total sources balance
	fprintf(balance_output_file,"# %s\n",s.str().c_str());  //drawing long line
	for (int sbi=0; sbi<n_subst; sbi++)
		fprintf(balance_output_file, src_total_format.c_str(),w+wl,"Total mass [M] and sources [M/T]",
				w,equation_->substances().names()[sbi].c_str(),
				w,mass_total[sbi],
				w,src_total_balance[sbi]);

	if (cumulative)
	{
		// Print cumulative sources
		fprintf(balance_output_file, "\n\n# Cumulative mass balance on time interval [%-g,%-g]\n"
				"# Initial mass [M] + sources integrated over time [M] - flux integrated over time [M] = current mass [M]\n"
				"# %-*s%-*s%-*s%-*s%-*s%-*s%-*s%-*s\n",
				initial_time, time,
				w,"[substance]",
				w,"[A=init. mass]",
				w,"[B=source]",
				w,"[C=flux]",
				w,"[A+B-C]",
				w,"[D=curr. mass]",
				w,"[A+B-C-D=err.]",
				w,"[rel. error]");

		for (int i=0; i<n_subst; i++)
		{
			double denominator = max(fabs(initial_mass[i]+integrated_sources[i]-integrated_fluxes[i]),fabs(mass_total[i]));
			fprintf(balance_output_file, "  %-*s%-*g%-*g%-*g%-*g%-*g%-*g%-*g\n",
					w,equation_->substances().names()[i].c_str(),
					w,initial_mass[i],
					w,integrated_sources[i],
					w,integrated_fluxes[i],
					w,initial_mass[i]+integrated_sources[i]-integrated_fluxes[i],
					w,mass_total[i],
					w,initial_mass[i]+integrated_sources[i]-integrated_fluxes[i]-mass_total[i],
					w,fabs(initial_mass[i]+integrated_sources[i]-integrated_fluxes[i]-mass_total[i])/(denominator==0?1:denominator));
		}
	}

	fprintf(balance_output_file, "\n\n");
}















Balance::Balance(const std::vector<unsigned int> &elem_regions,
		const RegionDB *region_db,
		const Input::Record &in_rec)
	: 	  output_format_(legacy),
	  	  regions_(*region_db),
	  	  be_regions_(elem_regions),
	  	  initial_time_(),
	  	  last_time_(),
	  	  initial_(true),
	  	  allocation_done_(false),
	  	  output_line_counter_(0)

{
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);

	cumulative_ = in_rec.val<bool>("cumulative");
	output_format_ = in_rec.val<OutputFormat>("format");

	if (rank_ == 0) {
		// set default value by output_format_
		std::string default_file_name;
		switch (output_format_)
		{
		case csv:
			default_file_name = "mass_balance.csv";
			break;
		case gnuplot:
			default_file_name = "mass_balance.dat";
			break;
		case legacy:
			default_file_name = "mass_balance.txt";
			break;
		}

		output_.open(string(in_rec.val<FilePath>("file", FilePath(default_file_name, FilePath::output_file))).c_str());
	}
}


Balance::~Balance()
{
	if (rank_ == 0) output_.close();
	for (unsigned int c=0; c<quantities_.size(); ++c)
	{
		MatDestroy(&(region_mass_matrix_[c]));
		MatDestroy(&(be_flux_matrix_[c]));
		MatDestroy(&(region_source_matrix_[c]));
		MatDestroy(&(region_source_rhs_[c]));
		VecDestroy(&(be_flux_vec_[c]));
		VecDestroy(&(region_source_vec_[c]));
	}
	delete[] region_mass_matrix_;
	delete[] be_flux_matrix_;
	delete[] be_flux_vec_;
	delete[] region_source_matrix_;
	delete[] region_source_rhs_;
	delete[] region_source_vec_;

	MatDestroy(&region_be_matrix_);
	VecDestroy(&ones_);
	VecDestroy(&ones_be_);
}


unsigned int Balance::add_quantity(const string &name)
{
	ASSERT(!allocation_done_, "Attempt to add quantity after allocation.");

	Quantity q(quantities_.size(), name);
	quantities_.push_back(q);

	return q.index_;
}


std::vector<unsigned int> Balance::add_quantities(const std::vector<string> &names)
{
	ASSERT(!allocation_done_, "Attempt to add quantity after allocation.");

	vector<unsigned int> indices;

	for (auto name : names)
		indices.push_back(add_quantity(name));

	return indices;
}


void Balance::allocate(unsigned int n_loc_dofs,
		unsigned int max_dofs_per_boundary)
{
	ASSERT(!allocation_done_, "Attempt to allocate Balance object multiple times.");
	// Max. number of regions to which a single dof can contribute.
	// TODO: estimate or compute this number directly (from mesh or dof handler).
	const int n_bulk_regs_per_dof = min(10, (int)regions_.bulk_size());
	const unsigned int n_quant = quantities_.size();
	const unsigned int n_bdr_reg = regions_.boundary_size();
	const unsigned int n_blk_reg = regions_.bulk_size();


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
	}



	region_mass_matrix_ = new Mat[n_quant];
	be_flux_matrix_ = new Mat[n_quant];
	region_source_matrix_ = new Mat[n_quant];
	region_source_rhs_ = new Mat[n_quant];
	be_flux_vec_ = new Vec[n_quant];
	region_source_vec_ = new Vec[n_quant];

	for (unsigned int c=0; c<n_quant; ++c)
	{
		MatCreateAIJ(PETSC_COMM_WORLD,
				n_loc_dofs,
				(rank_==0)?regions_.bulk_size():0,
				PETSC_DECIDE,
				PETSC_DECIDE,
				(rank_==0)?n_bulk_regs_per_dof:0,
				0,
				(rank_==0)?0:n_bulk_regs_per_dof,
				0,
				&(region_mass_matrix_[c]));

		MatCreateAIJ(PETSC_COMM_WORLD,
				be_regions_.size(),
				n_loc_dofs,
				PETSC_DECIDE,
				PETSC_DECIDE,
				max_dofs_per_boundary,
				0,
				0,
				0,
				&(be_flux_matrix_[c]));

		MatCreateAIJ(PETSC_COMM_WORLD,
				n_loc_dofs,
				(rank_==0)?regions_.bulk_size():0,
				PETSC_DECIDE,
				PETSC_DECIDE,
				(rank_==0)?n_bulk_regs_per_dof:0,
				0,
				(rank_==0)?0:n_bulk_regs_per_dof,
				0,
				&(region_source_matrix_[c]));

		MatCreateAIJ(PETSC_COMM_WORLD,
				n_loc_dofs,
				(rank_==0)?regions_.bulk_size():0,
				PETSC_DECIDE,
				PETSC_DECIDE,
				(rank_==0)?n_bulk_regs_per_dof:0,
				0,
				(rank_==0)?0:n_bulk_regs_per_dof,
				0,
				&(region_source_rhs_[c]));

		VecCreateMPI(PETSC_COMM_WORLD,
				be_regions_.size(),
				PETSC_DECIDE,
				&(be_flux_vec_[c]));

		VecCreateMPI(PETSC_COMM_WORLD,
				(rank_==0)?regions_.bulk_size():0,
				PETSC_DECIDE,
				&(region_source_vec_[c]));
	}

	MatCreateAIJ(PETSC_COMM_WORLD,
			be_regions_.size(),
			(rank_==0)?regions_.boundary_size():0,
			PETSC_DECIDE,
			PETSC_DECIDE,
			(rank_==0)?1:0,
			0,
			(rank_==0)?0:1,
			0,
			&region_be_matrix_
			);
	VecGetOwnershipRange(be_flux_vec_[0], &be_offset_, NULL);
	for (unsigned int loc_el=0; loc_el<be_regions_.size(); ++loc_el)
	{
		MatSetValue(region_be_matrix_,
				be_offset_+loc_el,
				be_regions_[loc_el],
				1,
				INSERT_VALUES);
	}
	MatAssemblyBegin(region_be_matrix_, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(region_be_matrix_, MAT_FINAL_ASSEMBLY);

	double *ones_array;
	VecCreateMPI(PETSC_COMM_WORLD,
			n_loc_dofs,
			PETSC_DECIDE,
			&ones_);
	VecGetArray(ones_, &ones_array);
	fill_n(ones_array, n_loc_dofs, 1);
	VecRestoreArray(ones_, &ones_array);

	VecCreateMPI(PETSC_COMM_WORLD,
			be_regions_.size(),
			PETSC_DECIDE,
			&ones_be_);
	VecGetArray(ones_be_, &ones_array);
	fill_n(ones_array, be_regions_.size(), 1);
	VecRestoreArray(ones_be_, &ones_array);

	allocation_done_ = true;
}


void Balance::start_mass_assembly(unsigned int quantity_idx)
{
	ASSERT(allocation_done_, "Balance structures are not allocated!");
	MatZeroEntries(region_mass_matrix_[quantity_idx]);
}


void Balance::start_flux_assembly(unsigned int quantity_idx)
{
	ASSERT(allocation_done_, "Balance structures are not allocated!");
	MatZeroEntries(be_flux_matrix_[quantity_idx]);
	VecZeroEntries(be_flux_vec_[quantity_idx]);
}


void Balance::start_source_assembly(unsigned int quantity_idx)
{
	ASSERT(allocation_done_, "Balance structures are not allocated!");
	MatZeroEntries(region_source_matrix_[quantity_idx]);
	MatZeroEntries(region_source_rhs_[quantity_idx]);
	VecZeroEntries(region_source_vec_[quantity_idx]);
}


void Balance::finish_mass_assembly(unsigned int quantity_idx)
{
	ASSERT(allocation_done_, "Balance structures are not allocated!");
	MatAssemblyBegin(region_mass_matrix_[quantity_idx], MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(region_mass_matrix_[quantity_idx], MAT_FINAL_ASSEMBLY);
}

void Balance::finish_flux_assembly(unsigned int quantity_idx)
{
	ASSERT(allocation_done_, "Balance structures are not allocated!");
	MatAssemblyBegin(be_flux_matrix_[quantity_idx], MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(be_flux_matrix_[quantity_idx], MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(be_flux_vec_[quantity_idx]);
	VecAssemblyEnd(be_flux_vec_[quantity_idx]);
}

void Balance::finish_source_assembly(unsigned int quantity_idx)
{
	ASSERT(allocation_done_, "Balance structures are not allocated!");
	MatAssemblyBegin(region_source_matrix_[quantity_idx], MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(region_source_matrix_[quantity_idx], MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(region_source_rhs_[quantity_idx], MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(region_source_rhs_[quantity_idx], MAT_FINAL_ASSEMBLY);
	MatMultTranspose(region_source_rhs_[quantity_idx], ones_, region_source_vec_[quantity_idx]);
}




void Balance::add_mass_matrix_values(unsigned int quantity_idx,
		unsigned int region_idx,
		const vector<int> &dof_indices,
		const vector<double> &values)
{
	PetscInt reg_array[1] = { (int)region_idx };

	MatSetValues(region_mass_matrix_[quantity_idx],
			dof_indices.size(),
			&(dof_indices[0]),
			1,
			reg_array,
			&(values[0]),
			ADD_VALUES);
}


void Balance::add_flux_matrix_values(unsigned int quantity_idx,
		unsigned int elem_idx,
		const vector<int> &dof_indices,
		const vector<double> &values)
{
	PetscInt elem_array[1] = { int(be_offset_+elem_idx) };

	MatSetValues(be_flux_matrix_[quantity_idx],
			1,
			elem_array,
			dof_indices.size(),
			&(dof_indices[0]),
			&(values[0]),
			ADD_VALUES);
}


void Balance::add_source_matrix_values(unsigned int quantity_idx,
		unsigned int region_idx,
		const vector<int> &dof_indices,
		const vector<double> &values)
{
	PetscInt reg_array[1] = { (int)region_idx };

	MatSetValues(region_source_matrix_[quantity_idx],
			dof_indices.size(),
			&(dof_indices[0]),
			1,
			reg_array,
			&(values[0]),
			ADD_VALUES);
}


void Balance::add_flux_vec_value(unsigned int quantity_idx,
		unsigned int elem_idx,
		double value)
{
	VecSetValue(be_flux_vec_[quantity_idx],
			be_offset_+elem_idx,
			value,
			ADD_VALUES);
}


void Balance::add_source_rhs_values(unsigned int quantity_idx,
		unsigned int region_idx,
		const vector<int> &dof_indices,
		const vector<double> &values)
{
	PetscInt reg_array[1] = { (int)region_idx };

	MatSetValues(region_source_rhs_[quantity_idx],
			dof_indices.size(),
			&(dof_indices[0]),
			1,
			reg_array,
			&(values[0]),
			ADD_VALUES);
}


void Balance::calculate_cumulative_sources(unsigned int quantity_idx,
		const Vec &solution,
		double dt)
{
	if (!cumulative_) return;

	ASSERT(allocation_done_, "Balance structures are not allocated!");

	Vec bulk_vec;

	VecCreateMPIWithArray(PETSC_COMM_WORLD,
			1,
			(rank_==0)?regions_.bulk_size():0,
			PETSC_DECIDE,
			&(sources_[quantity_idx][0]),
			&bulk_vec);

	// compute sources on bulk regions: S'.u + s
	VecZeroEntries(bulk_vec);
	MatMultTransposeAdd(region_source_matrix_[quantity_idx], solution, region_source_vec_[quantity_idx], bulk_vec);

	double sum_sources;
	VecSum(bulk_vec, &sum_sources);
	VecDestroy(&bulk_vec);

	if (rank_ == 0)
		// sum sources in time interval
		integrated_sources_[quantity_idx] += sum_sources*dt;
}


void Balance::calculate_cumulative_fluxes(unsigned int quantity_idx,
		const Vec &solution,
		double dt)
{
	if (!cumulative_) return;

	ASSERT(allocation_done_, "Balance structures are not allocated!");

	Vec boundary_vec;

	VecCreateMPIWithArray(PETSC_COMM_WORLD,
			1,
			(rank_==0)?regions_.boundary_size():0,
			PETSC_DECIDE,
			&(fluxes_[quantity_idx][0]),
			&boundary_vec);

	// compute fluxes on boundary regions: R'.(F.u + f)
	VecZeroEntries(boundary_vec);
	Vec temp;
	VecDuplicate(ones_be_, &temp);
	MatMultAdd(be_flux_matrix_[quantity_idx], solution, be_flux_vec_[quantity_idx], temp);
	MatMultTranspose(region_be_matrix_, temp, boundary_vec);
	VecDestroy(&temp);

	double sum_fluxes;
	VecSum(boundary_vec, &sum_fluxes);
	VecDestroy(&boundary_vec);

	if (rank_ == 0)
		// sum fluxes in time interval
		integrated_fluxes_[quantity_idx] += sum_fluxes*dt;
}


void Balance::calculate_mass(unsigned int quantity_idx,
		const Vec &solution,
		vector<double> &output_array)
{
	ASSERT(allocation_done_, "Balance structures are not allocated!");
	Vec bulk_vec;

	VecCreateMPIWithArray(PETSC_COMM_WORLD,
			1,
			(rank_==0)?regions_.bulk_size():0,
			PETSC_DECIDE,
			&(output_array[0]),
			&bulk_vec);

	// compute mass on regions: M'.u
	VecZeroEntries(bulk_vec);
	MatMultTranspose(region_mass_matrix_[quantity_idx], solution, bulk_vec);
	VecDestroy(&bulk_vec);
}


void Balance::calculate_source(unsigned int quantity_idx,
		const Vec &solution)
{
	ASSERT(allocation_done_, "Balance structures are not allocated!");
	Vec bulk_vec;

	VecCreateMPIWithArray(PETSC_COMM_WORLD,
			1,
			(rank_==0)?regions_.bulk_size():0,
			PETSC_DECIDE,
			&(sources_[quantity_idx][0]),
			&bulk_vec);

	// compute sources on bulk regions: S'.u + s
	VecZeroEntries(bulk_vec);
	MatMultTransposeAdd(region_source_matrix_[quantity_idx],
			solution,
			region_source_vec_[quantity_idx],
			bulk_vec);

	// compute positive/negative sources
	int lsize;
	Vec mat_r, rhs_r;
	const double *sol_array, *mat_array, *rhs_array;
	VecGetLocalSize(solution, &lsize);
	VecDuplicate(solution, &mat_r);
	VecDuplicate(solution, &rhs_r);
	VecGetArrayRead(solution, &sol_array);
	for (unsigned int r=0; r<regions_.bulk_size(); ++r)
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
			if (f > 0) sources_out_[quantity_idx][r] += f;
			else sources_in_[quantity_idx][r] += f;
		}

		VecRestoreArrayRead(mat_r, &mat_array);
		VecRestoreArrayRead(rhs_r, &rhs_array);
	}
	VecRestoreArrayRead(solution, &sol_array);
	VecDestroy(&rhs_r);
	VecDestroy(&mat_r);
	VecDestroy(&bulk_vec);
}


void Balance::calculate_flux(unsigned int quantity_idx,
		const Vec &solution)
{
	ASSERT(allocation_done_, "Balance structures are not allocated!");
	Vec boundary_vec;

	VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, (rank_==0)?regions_.boundary_size():0, PETSC_DECIDE, &(fluxes_[quantity_idx][0]), &boundary_vec);

	// compute fluxes on boundary regions: R'.(F.u + f)
	VecZeroEntries(boundary_vec);
	Vec temp;
	VecDuplicate(ones_be_, &temp);
	MatMultAdd(be_flux_matrix_[quantity_idx], solution, be_flux_vec_[quantity_idx], temp);
	MatMultTranspose(region_be_matrix_, temp, boundary_vec);

	// compute positive/negative fluxes
	fluxes_in_[quantity_idx].assign(regions_.boundary_size(), 0);
	fluxes_out_[quantity_idx].assign(regions_.boundary_size(), 0);
	const double *flux_array;
	int lsize;
	VecGetArrayRead(temp, &flux_array);
	VecGetLocalSize(temp, &lsize);
	for (int e=0; e<lsize; ++e)
	{
		if (flux_array[e] > 0)
			fluxes_out_[quantity_idx][be_regions_[e]] += flux_array[e];
		else
			fluxes_in_[quantity_idx][be_regions_[e]] += flux_array[e];
	}
	VecRestoreArrayRead(temp, &flux_array);
	VecDestroy(&temp);
	VecDestroy(&boundary_vec);
}

void Balance::add_cumulative_source(unsigned int quantity_idx, double source)
{
	if (rank_ == 0)
		integrated_sources_[quantity_idx] += source;
}





void Balance::output(double time)
{
	ASSERT(allocation_done_, "Balance structures are not allocated!");

	// gather results from processes and sum them up
	const unsigned int n_quant = quantities_.size();
	const unsigned int n_blk_reg = regions_.bulk_size();
	const unsigned int n_bdr_reg = regions_.boundary_size();
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
		const RegionSet & b_set = regions_.get_region_set("BOUNDARY");
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
		const RegionSet & bulk_set = regions_.get_region_set("BULK");
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
		}
	}

	last_time_ = time;


	// perform actual output
	switch (output_format_)
	{
	case csv:
		output_csv(time, ';', "");
		break;
	case gnuplot:
		output_csv(time, ' ', "#", 30);
		break;
	case legacy:
		output_legacy(time);
		break;
	}

	if (rank_ == 0)
	{
		sum_fluxes_.assign(n_quant, 0);
		sum_sources_.assign(n_quant, 0);
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
	const RegionSet & b_set = regions_.get_region_set("BOUNDARY");
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
	const RegionSet & bulk_set = regions_.get_region_set("BULK");
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
					<< setw(w) << fabs(initial_mass_[qi]+integrated_sources_[qi]-integrated_fluxes_[qi]-sum_masses_[qi])/(denominator==0?1:denominator)
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
	const RegionSet & bulk_set = regions_.get_region_set("BULK");
	for( RegionSet::const_iterator reg = bulk_set.begin(); reg != bulk_set.end(); ++reg)
	{
		for (unsigned int qi=0; qi<n_quant; qi++)
		{
			// print data header (repeat header after every "repeat" lines)
			if (repeat && (output_line_counter_%repeat == 0)) format_csv_output_header(delimiter, comment_string);

			output_ << format_csv_val(time, delimiter)
					<< format_csv_val(reg->label(), delimiter)
					<< format_csv_val(quantities_[qi].name_, delimiter)
					<< csv_zero_vals(3, delimiter)
					<< format_csv_val(masses_[qi][reg->bulk_idx()], delimiter)
					<< format_csv_val(sources_[qi][reg->bulk_idx()], delimiter)
					<< format_csv_val(sources_in_[qi][reg->bulk_idx()], delimiter)
					<< format_csv_val(sources_out_[qi][reg->bulk_idx()], delimiter)
					<< csv_zero_vals(3, delimiter) << endl;
			++output_line_counter_;
		}
	}

	// print mass fluxes over boundaries
	const RegionSet & b_set = regions_.get_region_set("BOUNDARY");
	for( RegionSet::const_iterator reg = b_set.begin(); reg != b_set.end(); ++reg)
	{
		for (unsigned int qi=0; qi<n_quant; qi++) {
			// print data header (repeat header after every "repeat" lines)
			if (repeat && (output_line_counter_%repeat == 0)) format_csv_output_header(delimiter, comment_string);

			output_ << format_csv_val(time, delimiter)
					<< format_csv_val(reg->label(), delimiter)
					<< format_csv_val(quantities_[qi].name_, delimiter)
					<< format_csv_val(fluxes_[qi][reg->boundary_idx()], delimiter)
					<< format_csv_val(fluxes_in_[qi][reg->boundary_idx()], delimiter)
					<< format_csv_val(fluxes_in_[qi][reg->boundary_idx()], delimiter)
					<< csv_zero_vals(7, delimiter) << endl;
			++output_line_counter_;
		}
	}

	if (cumulative_)
	{
		for (unsigned int qi=0; qi<n_quant; qi++)
		{
			// print data header (repeat header after every "repeat" lines)
			if (repeat && (output_line_counter_%repeat == 0)) format_csv_output_header(delimiter, comment_string);

			double error = sum_masses_[qi] - (initial_mass_[qi] + integrated_sources_[qi] - integrated_fluxes_[qi]);
			output_ << format_csv_val(time, delimiter)
					<< format_csv_val("ALL", delimiter)
					<< format_csv_val(quantities_[qi].name_, delimiter)
					<< format_csv_val(sum_fluxes_[qi], delimiter)
					<< format_csv_val(sum_fluxes_in_[qi], delimiter)
					<< format_csv_val(sum_fluxes_out_[qi], delimiter)
					<< format_csv_val(sum_masses_[qi], delimiter)
					<< format_csv_val(sum_sources_[qi], delimiter)
					<< format_csv_val(sum_sources_in_[qi], delimiter)
					<< format_csv_val(sum_sources_out_[qi], delimiter)
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
		ss << delimiter << "\"time\"";
	}

	output_ << comment_string << ss.str()
			<< format_csv_val("region", delimiter)
			<< format_csv_val("quantity", delimiter)
			<< format_csv_val("flux", delimiter)
			<< format_csv_val("flux_in", delimiter)
			<< format_csv_val("flux_out", delimiter)
			<< format_csv_val("mass", delimiter)
			<< format_csv_val("source", delimiter)
			<< format_csv_val("source_in", delimiter)
			<< format_csv_val("source_out", delimiter)
			<< format_csv_val("integrated_flux", delimiter)
			<< format_csv_val("integrated_source", delimiter)
			<< format_csv_val("error", delimiter)
			<< endl;
}

std::string Balance::format_csv_val(std::string val, char delimiter)
{
	std::stringstream ss;
	std::replace( val.begin(), val.end(), '\"', '\'');
	if (delimiter == ' ') {
		std::stringstream sval;
		sval << "\"" << val << "\"";
		ss << " " << setw(output_column_width-1) << sval.str();
	} else {
		ss << delimiter << "\"" << val << "\"";
	}

	return ss.str();
}

std::string Balance::format_csv_val(double val, char delimiter)
{
	std::stringstream ss;
	if (delimiter == ' ') {
		ss << " " << setw(output_column_width-1) << val;
	} else {
		ss << delimiter << val;
	}
	return ss.str();
}





