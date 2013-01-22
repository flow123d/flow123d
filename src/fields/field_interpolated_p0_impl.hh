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


#ifndef FIELD_INTERPOLATED_P0_IMPL_HH_
#define FIELD_INTERPOLATED_P0_IMPL_HH_

#include "field_interpolated_p0.hh"
#include "system/system.hh"
#include "mesh/msh_gmshreader.h"
#include "new_mesh/bih_tree.hh"
#include "new_mesh/ngh/include/intersection.h"
#include "new_mesh/ngh/include/point.h"
#include "system/sys_profiler.hh"
#include "boost/lexical_cast.hpp"
#include "system/tokenizer.hh"
#include "system/xio.h"


namespace it = Input::Type;
template <int spacedim, class Value>
it::Record FieldInterpolatedP0<spacedim, Value>::input_type
    = it::Record("FieldInterpolatedP0", "Field given by P0 data on another mesh. Currently defined only on boundary.")
	.derive_from(FieldBase<spacedim, Value>::input_type)
	// TODO: use mesh record here, but make it possibly of the same simplicity (the file name only)
	.declare_key("mesh", it::FileName::input(),it::Default::obligatory(),
			"File with the mesh from which we interpolate. (currently only GMSH supported)")
	// TODO: allow interpolation from VTK files (contains also mesh), and our own format of raw data, that includes:
	// mesh, dof_handler, and dof values
	.declare_key("raw_data", it::FileName::input(), it::Default::obligatory(),
			"File with raw output from flow calculation. Currently we can interpolate only pressure.");



template <int spacedim, class Value>
FieldInterpolatedP0<spacedim, Value>::FieldInterpolatedP0(const unsigned int n_comp)
: FieldBase<spacedim, Value>(n_comp)
{}




template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::init_from_input(const Input::Record &rec) {
    set_source_of_interpolation(
            rec.val<FilePath>("mesh"),
            rec.val<FilePath>("raw_data")
            );
}




template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::set_source_of_interpolation(const FilePath & mesh_file,
		const FilePath & raw_output) {

	// read mesh, create tree
    {
       mesh_ = new Mesh();
	   GmshMeshReader reader(mesh_file);
	   reader.read_mesh( mesh_);
	   // no call to mesh->setup_topology, we need only elements, no connectivity
    }
	bihTree_ = new BIHTree(mesh_);

	// read pressures
	Tokenizer tok(mesh_file);

	FILE* raw_output_file = xfopen( string(raw_output).c_str(), "rt");
	read_pressures(raw_output_file);
	xfclose(raw_output_file);
	//read_element_data_from_gmsh(tok, "xx");
}


template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::read_pressures(FILE* raw_output) {
	xprintf(Msg, " - FieldInterpolatedP0->read_pressures(FILE* raw_output)\n");

	int numElements;
	char line[ LINE_SIZE ];

	skip_to(raw_output, "$FlowField");
	xfgets(line, LINE_SIZE - 2, raw_output); //time
	xfgets(line, LINE_SIZE - 2, raw_output); //count of elements
	numElements = atoi(xstrtok(line));
	pressures_.reserve(numElements);
	for (int  i = 0; i < numElements; ++i) pressures_.push_back(0.0);
	xprintf(Msg, " - Reading pressures...");

	for (int i = 0; i < numElements; ++i) {
		int id;
		double pressure;

		xfgets(line, LINE_SIZE - 2, raw_output);

		//get element ID, presure
		id = atoi(xstrtok(line));
		pressure = atof(xstrtok(NULL));
		ElementFullIter ele = mesh_->element.find_id(id);
		pressures_[ ele.index() ] = pressure;
	}

	xprintf(Msg, " %d values of pressure read. O.K.\n", pressures_.size());
}


template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::calculate_triangle_pressure(TTriangle &element) {
	double total_measure, measure;
	BoundingBox elementBoundingBox = element.get_bounding_box();
	
	TIntersectionType iType;
	TTetrahedron tetrahedron;

	START_TIMER("find_elements_2D");
	((BIHTree *)bihTree_)->find_bounding_box(elementBoundingBox, searchedElements_);
	END_TIMER("find_elements_2D");

	total_measure = 0.0;
	pressure_ = 0.0;

        START_TIMER("compute_pressure_2D");
    ADD_CALLS( searchedElements_.size());
	for (std::vector<unsigned int>::iterator it = searchedElements_.begin(); it!=searchedElements_.end(); it++)
	{
                int idx = *it;
                ElementFullIter ele = mesh_->element( idx );
                if (ele->dim() == 3) {
                        createTetrahedron(ele, tetrahedron);
                        GetIntersection(element, tetrahedron, iType, measure);
                        if (iType == area) {
                                pressure_ += pressures_[ idx ] * measure;
                                total_measure += measure;
                        }
                } else {
                        xprintf(Err, "Dimension of source element must be 3!\n");
                }
        }
        pressure_ /= total_measure;
        END_TIMER("compute_pressure_2D");
}



template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::calculate_abscissa_pressure(TAbscissa &element) {
        double total_measure, measure;
	BoundingBox elementBoundingBox = element.get_bounding_box();
	TIntersectionType iType;
	TTetrahedron tetrahedron;

	START_TIMER("find_elements_1D");
	((BIHTree *)bihTree_)->find_bounding_box(elementBoundingBox, searchedElements_);
	END_TIMER("find_elements_1D");

        total_measure = 0.0;
	pressure_ = 0.0;
	START_TIMER("compute_pressure_1D");
	ADD_CALLS(searchedElements_.size());
	for (std::vector<unsigned int>::iterator it = searchedElements_.begin(); it!=searchedElements_.end(); it++)
	{
		int idx = *it;
		ElementFullIter ele = mesh_->element( idx );
		if (ele->dim() == 3) {
			createTetrahedron(ele, tetrahedron);
			GetIntersection(element, tetrahedron, iType, measure);
			if (iType == line) {
                                pressure_ += pressures_[ idx ] * measure;
                                total_measure += measure;
			}
		} else {
			//xprintf(Err, "Dimension of element must be 3!\n");
		}
	}
        pressure_ /= total_measure;
	END_TIMER("compute_pressure_1D");
}



template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::createTetrahedron(Element *ele, TTetrahedron &te) {
	ASSERT(( ele->dim() == 3 ), "Dimension of element must be 3!\n");

	te.SetPoints(TPoint(ele->node[0]->point()(0), ele->node[0]->point()(1), ele->node[0]->point()(2)),
				TPoint(ele->node[1]->point()(0), ele->node[1]->point()(1), ele->node[1]->point()(2)),
				TPoint(ele->node[2]->point()(0), ele->node[2]->point()(1), ele->node[2]->point()(2)),
				TPoint(ele->node[3]->point()(0), ele->node[3]->point()(1), ele->node[3]->point()(2)) );
}



template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::createTriangle(Element *ele, TTriangle &tr) {
	ASSERT(( ele->dim() == 2 ), "Dimension of element must be 2!\n");

	tr.SetPoints(TPoint(ele->node[0]->point()(0), ele->node[0]->point()(1), ele->node[0]->point()(2)),
				 TPoint(ele->node[1]->point()(0), ele->node[1]->point()(1), ele->node[1]->point()(2)),
				 TPoint(ele->node[2]->point()(0), ele->node[2]->point()(1), ele->node[2]->point()(2)) );
}



template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::createAbscissa(Element *ele, TAbscissa &ab) {
	ASSERT(( ele->dim() == 1 ), "Dimension of element must be 1!\n");

	ab.SetPoints(TPoint(ele->node[0]->point()(0), ele->node[0]->point()(1), ele->node[0]->point()(2)),
			 	 TPoint(ele->node[1]->point()(0), ele->node[1]->point()(1), ele->node[1]->point()(2)) );
}


template <int spacedim, class Value>
FieldResult FieldInterpolatedP0<spacedim, Value>::value(const Point<spacedim> &p, ElementAccessor<spacedim> &elm, typename Value::return_type &value) {

    switch (elm.dim()) {
        case 1: {
            TAbscissa abscissa;
            createAbscissa(elm.element(), abscissa);
            calculate_abscissa_pressure(abscissa);
            break;
        }
        case 2: {
            TTriangle triangle;
            createTriangle(elm.element(), triangle);
            calculate_triangle_pressure(triangle);
            break;
        }
        default:
            xprintf(Err, "Dimension of element must be 1 or 2!\n");
            break;
    }

	Value val(value);
	val(0,0)=pressure_;
	return result_other;
}





#endif
