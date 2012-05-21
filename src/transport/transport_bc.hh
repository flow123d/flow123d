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

#ifndef TRANSPORT_BC_HH_
#define TRANSPORT_BC_HH_

#include <petscmat.h>
#include "mesh/mesh.h"
#include "system/par_distribution.hh"


/**
 * @brief Reader class for transport boundary conditions.
 *
 * The class handles reading from files and storing the transport boundary conditions.
 *
 */
class TransportBC
{
public:

	/**
	 * @brief Constructor.
	 * @param mesh The mesh.
	 * @param n_subst Number of substances.
	 */
	TransportBC(Mesh *mesh, int n_subst);

	/// Destructor.
	~TransportBC();

	/// Reads the boundary conditions for the next time level.
	void read();

	/**
	 * @brief Returns the PETSc vector with b.c. for the given substance.
	 * @param sbi Substance number.
	 */
	Vec &get_vector(int sbi);

	/**
	 * @brief Returns the array of b.c.
	 * @param sbi Substance number.
	 */
	double *get_array(int sbi);

	/// Return the vector of time instants when boundary conditions change.
	const std::vector<double> &get_times();

	/// Returns the current time level.
	int get_time_level();

	/// Returns the distribution of @p bcv onto processors.
	const Distribution *distribution();


private:

	/**
	 * @brief Creates the input file name for the given time level.
	 * @param level Time level.
	 */
	std::string make_bc_file_name(int level);

	/// Boundary condition vector (each for one substance).
	Vec *bcv;

	/// Boundary condition array (each for one substance).
	double **bc;

	/// The mesh.
	Mesh *mesh_;

	/// Number of substances.
	const int n_substances;

	/// Distribution of the boundary elements.
	Distribution *distr;

	/// Times of change of boundary conditions.
	std::vector<double> bc_times;

	/// Current time level.
	int bc_time_level;

	/// Input file name (or its base) for reading b.c.
	std::string fname;
};


#endif
