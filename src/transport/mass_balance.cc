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
	.add_value(Balance::legacy, "legacy", "Legacy format used by previous program versions.");

Record MassBalance::input_type
	= Record("MassBalance", "Balance of mass, boundary fluxes and sources for transport of substances.")
	.declare_key("format", Balance::format_selection_input_type, Default("legacy"), "Format of output file.")
	.declare_key("cumulative", Bool(), Default("false"), "Compute cumulative balance over time. "
			"If true, then balance is calculated at each computational time step, which can slow down the program.")
	.declare_key("file", FileName::output(), Default("mass_balance.txt"), "File name for output of mass balance.")
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















Balance::Balance(const std::vector<std::string> &quantities,
		const std::vector<std::string> &components,
		const RegionDB *region_db,
		const bool cumulative,
		const std::string file)

	: 	  output_format_(legacy),
	  	  quantities_(quantities),
	  	  components_(components),
	  	  regions_(*region_db),
	  	  initial_time_(),
	  	  last_time_(),
	  	  initial_(true),
	  	  cumulative_(cumulative)

{
	const unsigned int n_bdr_reg = regions_.boundary_size();
	const unsigned int n_blk_reg = regions_.bulk_size();
	const unsigned int n_comp = components_.size();
	const unsigned int n_quant = quantities_.size();

	fluxes_    .resize(n_quant, vector<vector<double> >(n_comp, vector<double>(n_bdr_reg, 0)));
	fluxes_in_ .resize(n_quant, vector<vector<double> >(n_comp, vector<double>(n_bdr_reg, 0)));
	fluxes_out_.resize(n_quant, vector<vector<double> >(n_comp, vector<double>(n_bdr_reg, 0)));

	masses_     .resize(n_quant, vector<vector<double> >(n_comp, vector<double>(n_blk_reg, 0)));
	sources_    .resize(n_quant, vector<vector<double> >(n_comp, vector<double>(n_blk_reg, 0)));
	sources_in_ .resize(n_quant, vector<vector<double> >(n_comp, vector<double>(n_blk_reg, 0)));
	sources_out_.resize(n_quant, vector<vector<double> >(n_comp, vector<double>(n_blk_reg, 0)));

	sum_fluxes_     .resize(n_quant, 0);
	sum_fluxes_in_  .resize(n_quant, 0);
	sum_fluxes_out_ .resize(n_quant, 0);
	sum_masses_     .resize(n_quant, 0);
	sum_sources_    .resize(n_quant, 0);
	sum_sources_in_ .resize(n_quant, 0);
	sum_sources_out_.resize(n_quant, 0);

	if (cumulative_)
	{
		initial_mass_.resize(n_quant, 0);

		integrated_fluxes_. resize(n_quant, 0);
		integrated_sources_.resize(n_quant, 0);
	}

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);

	if (rank_ == 0)
		output_.open(file.c_str());
}


Balance::~Balance()
{
	if (rank_ == 0) output_.close();
	delete[] region_mass_matrix_;
	delete[] region_flux_matrix_;
	delete[] region_flux_rhs_;
	delete[] region_flux_vec_;
	delete[] region_source_matrix_;
	delete[] region_source_rhs_;
	delete[] region_source_vec_;
}


void Balance::allocate_matrices(unsigned int lsize)
{
	const int n_bulk_regs_per_dof = min(10, (int)regions_.bulk_size());
	const int n_boundary_regs_per_dof = min(10, (int)regions_.boundary_size());
	const unsigned int n_comp = quantities_.size()*components_.size();

	region_mass_matrix_ = new Mat[n_comp];
	region_flux_matrix_ = new Mat[n_comp];
	region_source_matrix_ = new Mat[n_comp];
	region_flux_rhs_ = new Mat[n_comp];
	region_source_rhs_ = new Mat[n_comp];
	region_flux_vec_ = new Vec[n_comp];
	region_source_vec_ = new Vec[n_comp];

	for (unsigned int c=0; c<n_comp; ++c)
	{
		MatCreateAIJ(PETSC_COMM_WORLD,
				lsize,
				(rank_==0)?regions_.bulk_size():0,
				PETSC_DECIDE,
				PETSC_DECIDE,
				(rank_==0)?n_bulk_regs_per_dof:0,
				0,
				(rank_==0)?0:n_bulk_regs_per_dof,
				0,
				&(region_mass_matrix_[c]));

		MatCreateAIJ(PETSC_COMM_WORLD,
				lsize,
				(rank_==0)?regions_.boundary_size():0,
				PETSC_DECIDE,
				PETSC_DECIDE,
				(rank_==0)?n_boundary_regs_per_dof:0,
				0,
				(rank_==0)?0:n_boundary_regs_per_dof,
				0,
				&(region_flux_matrix_[c]));

		MatCreateAIJ(PETSC_COMM_WORLD,
				lsize,
				(rank_==0)?regions_.bulk_size():0,
				PETSC_DECIDE,
				PETSC_DECIDE,
				(rank_==0)?n_bulk_regs_per_dof:0,
				0,
				(rank_==0)?0:n_bulk_regs_per_dof,
				0,
				&(region_source_matrix_[c]));

		MatCreateAIJ(PETSC_COMM_WORLD,
				lsize,
				(rank_==0)?regions_.boundary_size():0,
				PETSC_DECIDE,
				PETSC_DECIDE,
				(rank_==0)?n_boundary_regs_per_dof:0,
				0,
				(rank_==0)?0:n_boundary_regs_per_dof,
				0,
				&(region_flux_rhs_[c]));

		MatCreateAIJ(PETSC_COMM_WORLD,
				lsize,
				(rank_==0)?regions_.bulk_size():0,
				PETSC_DECIDE,
				PETSC_DECIDE,
				(rank_==0)?n_bulk_regs_per_dof:0,
				0,
				(rank_==0)?0:n_bulk_regs_per_dof,
				0,
				&(region_source_rhs_[c]));

		VecCreateMPI(PETSC_COMM_WORLD,
				(rank_==0)?regions_.boundary_size():0,
				PETSC_DECIDE,
				&(region_flux_vec_[c]));

		VecCreateMPI(PETSC_COMM_WORLD,
				(rank_==0)?regions_.bulk_size():0,
				PETSC_DECIDE,
				&(region_source_vec_[c]));
	}

	double *ones_array = new double[lsize];
	fill_n(ones_array, lsize, 1);
	VecCreateMPIWithArray(PETSC_COMM_WORLD,
			1,
			lsize,
			PETSC_DECIDE,
			ones_array,
			&ones_);
	VecAssemblyBegin(ones_);
	VecAssemblyEnd(ones_);
}


void Balance::set_mass_matrix_values(unsigned int quantity_idx,
		unsigned int component_idx,
		unsigned int region_idx,
		int n_dofs,
		int *dof_indices,
		double *values)
{
	PetscInt reg_array[1] = { (int)region_idx };

	MatSetValues(region_mass_matrix_[quantity_idx*components_.size()+component_idx],
			n_dofs,
			dof_indices,
			1,
			reg_array,
			values,
			ADD_VALUES);
}


void Balance::set_flux_matrix_values(unsigned int quantity_idx,
		unsigned int component_idx,
		unsigned int region_idx,
		int n_dofs,
		int *dof_indices,
		double *values)
{
	PetscInt reg_array[1] = { (int)region_idx };

	MatSetValues(region_flux_matrix_[quantity_idx*components_.size()+component_idx],
			n_dofs,
			dof_indices,
			1,
			reg_array,
			values,
			ADD_VALUES);

}


void Balance::set_source_matrix_values(unsigned int quantity_idx,
		unsigned int component_idx,
		unsigned int region_idx,
		int n_dofs,
		int *dof_indices,
		double *values)
{
	PetscInt reg_array[1] = { (int)region_idx };

	MatSetValues(region_source_matrix_[quantity_idx*components_.size()+component_idx],
			n_dofs,
			dof_indices,
			1,
			reg_array,
			values,
			ADD_VALUES);
}


void Balance::set_flux_rhs_values(unsigned int quantity_idx,
		unsigned int component_idx,
		unsigned int region_idx,
		int n_dofs,
		int *dof_indices,
		double *values)
{
	PetscInt reg_array[1] = { (int)region_idx };

	MatSetValues(region_flux_rhs_[quantity_idx*components_.size()+component_idx],
			n_dofs,
			dof_indices,
			1,
			reg_array,
			values,
			ADD_VALUES);
}


void Balance::set_source_rhs_values(unsigned int quantity_idx,
		unsigned int component_idx,
		unsigned int region_idx,
		int n_dofs,
		int *dof_indices,
		double *values)
{
	PetscInt reg_array[1] = { (int)region_idx };

	MatSetValues(region_source_rhs_[quantity_idx*components_.size()+component_idx],
			n_dofs,
			dof_indices,
			1,
			reg_array,
			values,
			ADD_VALUES);
}


void Balance::calculate_cumulative(unsigned int quantity_idx,
		unsigned int component_idx,
		const Vec &solution,
		double dt)
{
	if (!cumulative_) return;

	const unsigned int n_comp = components_.size();
	Vec bulk_vec, boundary_vec;

	VecCreateMPIWithArray(PETSC_COMM_WORLD,
			1,
			(rank_==0)?regions_.bulk_size():0,
			PETSC_DECIDE,
			&(sources_[quantity_idx][component_idx][0]),
			&bulk_vec);

	VecCreateMPIWithArray(PETSC_COMM_WORLD,
			1,
			(rank_==0)?regions_.boundary_size():0,
			PETSC_DECIDE,
			&(fluxes_[quantity_idx][component_idx][0]),
			&boundary_vec);

	// compute sources on bulk regions: S'.u + s
	VecZeroEntries(bulk_vec);
	MatMultTransposeAdd(region_source_matrix_[quantity_idx*n_comp+component_idx], solution, region_source_vec_[quantity_idx*n_comp+component_idx], bulk_vec);

	// compute fluxes on boundary regions: F'.u + f
	VecZeroEntries(boundary_vec);
	MatMultTransposeAdd(region_flux_matrix_[quantity_idx*n_comp+component_idx], solution, region_flux_vec_[quantity_idx*n_comp+component_idx], boundary_vec);

	if (rank_ == 0)
	{
		const unsigned int n_quant = quantities_.size();

		for (unsigned int qi=0; qi<n_quant; qi++)
		{
			sum_fluxes_[qi] = 0;
			sum_sources_[qi] = 0;
		}

		// sum all boundary fluxes
		const RegionSet & b_set = regions_.get_region_set("BOUNDARY");
		for( RegionSet::const_iterator reg = b_set.begin(); reg != b_set.end(); ++reg)
			for (unsigned int qi=0; qi<n_quant; qi++)
				for (unsigned int ci=0; ci<n_comp; ci++)
					sum_fluxes_[qi] += fluxes_[qi][ci][reg->boundary_idx()];

		// sum all volume sources
		const RegionSet & bulk_set = regions_.get_region_set("BULK");
		for( RegionSet::const_iterator reg = bulk_set.begin(); reg != bulk_set.end(); ++reg)
			for (unsigned int qi=0; qi<n_quant; qi++)
				for (unsigned int ci=0; ci<n_comp; ci++)
					sum_sources_[qi] += sources_[qi][ci][reg->bulk_idx()];

		// sum sources and fluxes in time interval
		integrated_sources_[quantity_idx] += sum_sources_[quantity_idx]*dt;
		integrated_fluxes_[quantity_idx] += sum_fluxes_[quantity_idx]*dt;
	}

}


void Balance::calculate_mass(unsigned int quantity_idx,
		unsigned int component_idx,
		const Vec &solution)
{
	const unsigned int n_comp = components_.size();
	Vec bulk_vec;

	VecCreateMPIWithArray(PETSC_COMM_WORLD,
			1,
			(rank_==0)?regions_.bulk_size():0,
			PETSC_DECIDE,
			&(masses_[quantity_idx][component_idx][0]),
			&bulk_vec);

	// compute mass on regions: M'.u
	VecZeroEntries(bulk_vec);
	MatMultTranspose(region_mass_matrix_[quantity_idx*n_comp+component_idx], solution, bulk_vec);
}


void Balance::calculate_source(unsigned int quantity_idx,
		unsigned int component_idx,
		const Vec &solution)
{
	const unsigned int n_comp = components_.size();
	Vec bulk_vec;

	VecCreateMPIWithArray(PETSC_COMM_WORLD,
			1,
			(rank_==0)?regions_.bulk_size():0,
			PETSC_DECIDE,
			&(sources_[quantity_idx][component_idx][0]),
			&bulk_vec);

	// compute sources on bulk regions: S'.u + s
	VecZeroEntries(bulk_vec);
	MatMultTransposeAdd(region_source_matrix_[quantity_idx*n_comp+component_idx],
			solution,
			region_source_vec_[quantity_idx*n_comp+component_idx],
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
		MatGetColumnVector(region_source_matrix_[quantity_idx*n_comp+component_idx], mat_r, r);
		MatGetColumnVector(region_source_rhs_[quantity_idx*n_comp+component_idx], rhs_r, r);

		VecGetArrayRead(mat_r, &mat_array);
		VecGetArrayRead(rhs_r, &rhs_array);

		sources_in_[quantity_idx][component_idx][r] = 0;
		sources_out_[quantity_idx][component_idx][r] = 0;
		for (int i=0; i<lsize; ++i)
		{
			double f = mat_array[i]*sol_array[i] + rhs_array[i];
			if (f > 0) sources_out_[quantity_idx][component_idx][r] += f;
			else sources_in_[quantity_idx][component_idx][r] += f;
		}
	}

	// gather results from processes and sum them up
	const unsigned int n_quant = quantities_.size();
	const unsigned int n_blk_reg = regions_.bulk_size();
	const int buf_size = n_quant*n_comp*2*n_blk_reg;
	double sendbuffer[buf_size], recvbuffer[buf_size];
	for (unsigned int qi=0; qi<n_quant; qi++)
	{
		for (unsigned int ci=0; ci<n_comp; ci++)
		{
			for (unsigned int ri=0; ri<n_blk_reg; ri++)
			{
				sendbuffer[(qi*n_comp+ci)*2*n_blk_reg +           + ri] = sources_in_[qi][ci][ri];
				sendbuffer[(qi*n_comp+ci)*2*n_blk_reg + n_blk_reg + ri] = sources_out_[qi][ci][ri];
			}
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
			for (unsigned int ci=0; ci<n_comp; ci++)
			{
				for (unsigned int ri=0; ri<n_blk_reg; ri++)
				{
					sources_in_[qi][ci][ri]  = recvbuffer[(qi*n_comp+ci)*2*n_blk_reg +           + ri];
					sources_out_[qi][ci][ri] = recvbuffer[(qi*n_comp+ci)*2*n_blk_reg + n_blk_reg + ri];
				}
			}
		}

	}
}


void Balance::calculate_flux(unsigned int quantity_idx,
		unsigned int component_idx,
		const Vec &solution)
{
	const unsigned int n_comp = components_.size();
	Vec boundary_vec;

	VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, (rank_==0)?regions_.boundary_size():0, PETSC_DECIDE, &(fluxes_[quantity_idx][component_idx][0]), &boundary_vec);

	// compute fluxes on boundary regions: F'.u + f
	VecZeroEntries(boundary_vec);
	MatMultTransposeAdd(region_flux_matrix_[quantity_idx*n_comp+component_idx], solution, region_flux_vec_[quantity_idx*n_comp+component_idx], boundary_vec);

	// compute positive/negative fluxes
	int lsize;
	Vec mat_r, rhs_r;
	const double *sol_array, *mat_array, *rhs_array;
	VecGetLocalSize(solution, &lsize);
	VecDuplicate(solution, &mat_r);
	VecDuplicate(solution, &rhs_r);
	VecGetArrayRead(solution, &sol_array);
	for (unsigned int r=0; r<regions_.boundary_size(); ++r)
	{
		MatGetColumnVector(region_flux_matrix_[quantity_idx*n_comp+component_idx], mat_r, r);
		MatGetColumnVector(region_flux_rhs_[quantity_idx*n_comp+component_idx], rhs_r, r);

		VecGetArrayRead(mat_r, &mat_array);
		VecGetArrayRead(rhs_r, &rhs_array);

		fluxes_in_[quantity_idx][component_idx][r] = 0;
		fluxes_out_[quantity_idx][component_idx][r] = 0;
		for (int i=0; i<lsize; ++i)
		{
			double f = mat_array[i]*sol_array[i] + rhs_array[i];
			if (f > 0) fluxes_out_[quantity_idx][component_idx][r] += f;
			else fluxes_in_[quantity_idx][component_idx][r] += f;
		}
	}

	// gather results from processes and sum them up
	const unsigned int n_quant = quantities_.size();
	const unsigned int n_bdr_reg = regions_.boundary_size();
	const int buf_size = n_quant*n_comp*2*n_bdr_reg;
	double sendbuffer[buf_size], recvbuffer[buf_size];
	for (unsigned int qi=0; qi<n_quant; qi++)
	{
		for (unsigned int ci=0; ci<n_comp; ci++)
		{
			for (unsigned int ri=0; ri<n_bdr_reg; ri++)
			{
				sendbuffer[(qi*n_comp+ci)*2*n_bdr_reg +           + ri] = fluxes_in_[qi][ci][ri];
				sendbuffer[(qi*n_comp+ci)*2*n_bdr_reg + n_bdr_reg + ri] = fluxes_out_[qi][ci][ri];
			}
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
			for (unsigned int ci=0; ci<n_comp; ci++)
			{
				for (unsigned int ri=0; ri<n_bdr_reg; ri++)
				{
					fluxes_in_[qi][ci][ri]  = recvbuffer[(qi*n_comp+ci)*2*n_bdr_reg +           + ri];
					fluxes_out_[qi][ci][ri] = recvbuffer[(qi*n_comp+ci)*2*n_bdr_reg + n_bdr_reg + ri];
				}
			}
		}

	}
}





void Balance::output(double time)
{
	if (rank_ == 0)
	{
		const unsigned int n_comp = components_.size();
		const unsigned int n_quant = quantities_.size();

		for (unsigned int qi=0; qi<n_quant; qi++)
		{
			sum_fluxes_[qi] = 0;
			sum_fluxes_in_[qi] = 0;
			sum_fluxes_out_[qi] = 0;
			sum_masses_[qi] = 0;
			sum_sources_[qi] = 0;
			sum_sources_in_[qi] = 0;
			sum_sources_out_[qi] = 0;
		}

		// sum all boundary fluxes
		const RegionSet & b_set = regions_.get_region_set("BOUNDARY");
		for( RegionSet::const_iterator reg = b_set.begin(); reg != b_set.end(); ++reg)
		{
			for (unsigned int qi=0; qi<n_quant; qi++)
			{
				for (unsigned int ci=0; ci<n_comp; ci++)
				{
					sum_fluxes_[qi]     += fluxes_    [qi][ci][reg->boundary_idx()];
					sum_fluxes_in_[qi]  += fluxes_in_ [qi][ci][reg->boundary_idx()];
					sum_fluxes_out_[qi] += fluxes_out_[qi][ci][reg->boundary_idx()];
				}
			}
		}

		// sum all volume sources
		const RegionSet & bulk_set = regions_.get_region_set("BULK");
		for( RegionSet::const_iterator reg = bulk_set.begin(); reg != bulk_set.end(); ++reg)
		{
			for (unsigned int qi=0; qi<n_quant; qi++)
			{
				for (unsigned int ci=0; ci<n_comp; ci++)
				{
					sum_masses_[qi] += masses_[qi][ci][reg->bulk_idx()];
					sum_sources_[qi] += sources_[qi][ci][reg->bulk_idx()];
					sum_sources_in_[qi] += sources_in_[qi][ci][reg->bulk_idx()];
					sum_sources_out_[qi] += sources_out_[qi][ci][reg->bulk_idx()];
				}
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


//			// sum sources and fluxes in time interval
//			for (unsigned int qi=0; qi<n_quant; qi++)
//			{
//				integrated_sources_[qi] += sum_sources_[qi]*(time-last_time_);
//				integrated_fluxes_[qi] += sum_fluxes_[qi]*(time-last_time_);
//			}
		}
	}

	last_time_ = time;



	// perform actual output
	switch (output_format_)
	{
	case csv:
	case gnuplot:
	case legacy:
		output_legacy(time);
		break;
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
					<< setw(w)  << quantities_[qi].c_str()
					<< setw(w)  << fluxes_[qi][0][reg->boundary_idx()]
					<< setw(w)  << fluxes_out_[qi][0][reg->boundary_idx()]
					<< setw(w)  << fluxes_in_[qi][0][reg->boundary_idx()]
					<< endl;
		}
	}

	// draw long line
	output_ << "# " << setw(w*c+wl) << setfill('-') << "" << setfill(' ') << endl;

	// total boundary balance
	for (unsigned int qi=0; qi<n_quant; qi++)
		output_ << "# " << setiosflags(ios::left)
				<< setw(w+wl) << "Total mass flux of substance [M/T]"
				<< setw(w)    << quantities_[qi].c_str()
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
					<< setw(w)  << quantities_[qi].c_str()
					<< setw(w)  << masses_[qi][0][reg->bulk_idx()]
					<< setw(w)  << sources_[qi][0][reg->bulk_idx()]
					<< endl;
		}
	}

	// draw long line
	output_ << "# " << setw(w*c+wl) << setfill('-') << "" << setfill(' ') << endl;

	// total sources balance
	for (unsigned int qi=0; qi<n_quant; qi++)
		output_ << "# " << setiosflags(ios::left) << setw(w+wl) << "Total mass [M] and sources [M/T]"
				<< setw(w) << quantities_[qi].c_str()
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
					<< setw(w) << quantities_[qi].c_str()
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






