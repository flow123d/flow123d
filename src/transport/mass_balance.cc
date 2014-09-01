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

Record MassBalance::input_type
	= Record("MassBalance", "Balance of mass, boundary fluxes and sources for transport of substances.")
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
		//DBGMSG("writing reg->idx() and id() and boundary_idx(): %d\t%d\t%d\n", reg->idx(), reg->id(), reg->boundary_idx());
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
					w, equation_->substance_names()[sbi].c_str(),w, bcd_balance[sbi][reg->boundary_idx()],
					w, bcd_plus_balance[sbi][reg->boundary_idx()],
					w, bcd_minus_balance[sbi][reg->boundary_idx()]);
		}
	}
	// total boundary balance
	fprintf(balance_output_file,"# %s\n",s.str().c_str());  // drawing long line
	for (int sbi=0; sbi<n_subst; sbi++)
		fprintf(balance_output_file, bc_total_format.c_str(),w+wl,"Total mass flux of substance [M/T]",
				w,equation_->substance_names()[sbi].c_str(),w,bcd_total_balance[sbi], w, bcd_total_outflow[sbi], w, bcd_total_inflow[sbi]);
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
					reg->label().c_str(), w, equation_->substance_names()[sbi].c_str(),
					w,mass[sbi][reg->bulk_idx()],
					w,src_balance[sbi][reg->bulk_idx()]);
		}
	}
	// total sources balance
	fprintf(balance_output_file,"# %s\n",s.str().c_str());  //drawing long line
	for (int sbi=0; sbi<n_subst; sbi++)
		fprintf(balance_output_file, src_total_format.c_str(),w+wl,"Total mass [M] and sources [M/T]",
				w,equation_->substance_names()[sbi].c_str(),
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
					w,equation_->substance_names()[i].c_str(),
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


