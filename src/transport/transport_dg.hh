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
 * $Id: quadrature.hh 1352 2011-09-23 14:14:47Z jan.stebel $
 * $Revision: 1352 $
 * $LastChangedBy: jan.stebel $
 * $LastChangedDate: 2011-09-23 16:14:47 +0200 (Fri, 23 Sep 2011) $
 *
 * @file
 * @brief Discontinuous Galerkin method for equation of transport with dispersion.
 *  @author Jan Stebel
 */

#ifndef TRANSPORT_DG_HH_
#define TRANSPORT_DG_HH_

#include "transport_operator_splitting.hh"
#include "la_linsys.hh"
#include "fem/dofhandler.hh"
#include "fem/finite_element.hh"



//using namespace arma;

class TransportDG : public TransportBase
{
public:

    TransportDG(TimeMarks &marks,  Mesh &init_mesh, MaterialDatabase &material_database);

	void update_solution();

	void get_solution_vector(double * &vector, unsigned int &size);

	void get_parallel_solution_vector(Vec &vector);

	void set_velocity_field(Vec &velocity_vector);

	void output_data();

private:

	/*
	 * Structures for:
	 * system matrix (sparse), system r.h.s., solution vector
	 * dofhandler, etc.
	 *
	 */
	void assemble();

	/**
	 * Vector of fluxes across element edges.
     */
	Vec flux_vector;

	Vec solution;

	LinSys *ls;

	struct Solver *solver;

	DOFHandler<2> *dof_handler2d;

	FiniteElement<2> *fe2d;


};


#endif /* TRANSPORT_DG_HH_ */
