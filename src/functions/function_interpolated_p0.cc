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
 * $Id: function_interpolated_p0.cc 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

#include "function_interpolated_p0.hh"
#include "system/system.hh"
#include "mesh/msh_gmshreader.h"

template <int dim>
FunctionInterpolatedP0<dim>::FunctionInterpolatedP0() {

}


template <int dim>
void FunctionInterpolatedP0<dim>::set_element(ElementFullIter &element){
	element_ = element;
}

template <int dim>
void FunctionInterpolatedP0<dim>::set_source_of_interpolation(const std::string & mesh_file,
		const std::string & raw_output) {

	const std::string& ngh_fname = "";
	const std::string& bcd_fname = "";
	MeshReader* meshReader = new GmshMeshReader();

	mesh_ = new Mesh(ngh_fname, bcd_fname);
	meshReader->read(mesh_file, mesh_);

	FILE* file = xfopen(raw_output.c_str(), "rt");
	char line[ LINE_SIZE ];

	skip_to(file, "$FlowField");
	xfgets(line, LINE_SIZE - 2, file); //time
	xfgets(line, LINE_SIZE - 2, file); //count of elements
	int numElements = atoi(xstrtok(line));

	for (int i = 0; i < numElements; ++i) {
	    int id;
	    double pressure;

		xfgets(line, LINE_SIZE - 2, file);

	    //get element ID, presure
	    id = atoi(xstrtok(line));
	    pressure = atoi(xstrtok(NULL));
	    //ElementFullIter ele = mesh_->element.find_id(id);
	    //set pressure of element


	}

	printf("Count of elements: %s\n", line);

	xfclose(file);

}
