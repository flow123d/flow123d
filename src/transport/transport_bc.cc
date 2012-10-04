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
 * @brief Reader and storage for transport boundary conditions.
 *  @author Jan Stebel
 */


#include "io_namehandler.hh"

#include "transport/transport_bc.hh"
#include <iostream>
#include <iomanip>



TransportBC::TransportBC(Mesh *mesh, int n_subst, const Input::Record &in_rec)
	: mesh_(mesh),
	  n_substances(n_subst)
{
	int ierr, np, rank;

	MPI_Barrier(PETSC_COMM_WORLD);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &np);

	distr = new Distribution(Distribution::Localized, mesh_->n_boundaries());

	bcv = new Vec[n_subst];
	bc = new double*[n_subst];

	for (int sbi=0; sbi<n_subst; sbi++)
	{
		bc[sbi] = new double[distr->lsize(rank)];
		ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, distr->lsize(rank), mesh->n_boundaries(), bc[sbi], &bcv[sbi]);
	}

	fname = in_rec.val<FilePath>("boundary_file");

	Input::Iterator<Input::Array> bc_it = in_rec.find<Input::Array>("bc_times");
	if (bc_it) bc_it->copy_to(bc_times);

	FILE * f;
	string fn;

	if (bc_times.size() == 0) {
	    DBGMSG("one bc file \n");
		// only one boundary condition, check filename and read bc condition
		bc_time_level = -1;
		fn = make_bc_file_name(-1);
		if ( !(f = xfopen(fn.c_str(), "rt")) ) {
			xprintf(UsrErr,"Missing file: %s", fn.c_str());
		}
		xfclose(f);
		read();
	} else {
		bc_time_level = 0;
		// set initial bc to zero
		for(unsigned int sbi=0; sbi < n_substances; ++sbi) VecZeroEntries(bcv[sbi]);
		// check files for bc time levels
		for(unsigned int i=0;i<bc_times.size();++i) {
			fn = make_bc_file_name(i);
			if ( !(f = xfopen(fn.c_str(), "rt")) ) {
				xprintf(UsrErr,"Missing file: %s", fn.c_str());
			}
			xfclose(f);
		}
	}
}

TransportBC::~TransportBC()
{
	for (int i=0; i<n_substances; i++)
	{
		VecDestroy(&bcv[i]);
		delete[] bc[i];
	}
	delete[] bcv;
	delete[] bc;
}


void TransportBC::read()
{
	int rank, sbi;

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if (rank == 0)
	{
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
	}
	for (sbi = 0; sbi < n_substances; sbi++)
		VecAssemblyBegin(bcv[sbi]);
	//for (sbi = 0; sbi < n_substances; sbi++)
	//    VecZeroEntries(bcvcorr[sbi]);
	for (sbi = 0; sbi < n_substances; sbi++)
		VecAssemblyEnd(bcv[sbi]);

	if (bc_time_level != -1) bc_time_level++;
}

string TransportBC::make_bc_file_name(int level)
{
    string bc_fname = fname;

	if (level >= 0 )
	{
        stringstream name_str;
        name_str << fname << "_" << setfill('0') << setw(3) << level;
        bc_fname = name_str.str();
    }

    return bc_fname;
}


Vec &TransportBC::get_vector(int sbi)
{
	return bcv[sbi];
}

double *TransportBC::get_array(int sbi)
{
	return bc[sbi];
}

const vector<double> &TransportBC::get_times()
{
	return bc_times;
}

int TransportBC::get_time_level()
{
	return bc_time_level;
}

const Distribution *TransportBC::distribution()
{
	return distr;
}





