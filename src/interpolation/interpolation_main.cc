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
 * $Id: interpolation_main.cc 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

//#include "system/sys_profiler.hh"
#include "system/system.hh"
#include "mesh/msh_gmshreader.h"
//#include "new_mesh/mesh.hh"
#include "mesh/mesh.h"
#include "new_mesh/bounding_interval_hierarchy.hh"
#include <armadillo>


int main(int argc, char **argv) {

	const std::string& file_name = "../tests/01_steady_flow_123d/input/test1.msh";

	MeshReader* meshReader = new GmshMeshReader();

	const std::string& ngh_fname = "../tests/01_steady_flow_123d/input/test1.ngh";
	const std::string& bcd_fname = "../tests/01_steady_flow_123d/input/test1.fbc";
	Mesh* mesh = new Mesh(ngh_fname, bcd_fname);
	meshReader->read(file_name, mesh);
	BoundingIntevalHierachy* bihTree = new BoundingIntevalHierachy(mesh);

	arma::vec3 point = arma::vec3("0.45 0.46 0.47");
	std::vector<BoundingBox *> searchedElements;
	bihTree->get_element(point, searchedElements);

	xprintf(Msg, " - searched elements ids: ");

	for (std::vector<BoundingBox *>::iterator tmp = searchedElements.begin(); tmp!=searchedElements.end(); tmp++)
	{
		BoundingBox* b = *tmp;

		if (b->contains_point(point)) xprintf(Msg, "%d ", b->getId());
	}
	xprintf(Msg, "\n");

	xprintf(Msg, " - interpolation_main executed\n");

	getchar();
}
