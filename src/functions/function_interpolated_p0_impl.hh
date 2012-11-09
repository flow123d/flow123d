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

#include "functions/functions_all.hh"

#ifndef FUNCTION_INTERPOLATED_P0_IMPL_HH_
#define FUNCTION_INTERPOLATED_P0_IMPL_HH_

#include "function_interpolated_p0.hh"
#include "system/system.hh"
#include "mesh/msh_gmshreader.h"
#include "new_mesh/bih_tree.hh"
#include "new_mesh/ngh/include/intersection.h"
#include "new_mesh/ngh/include/point.h"
#include "system/tokenizer.hh"


template <int dim>
FunctionInterpolatedP0<dim>::FunctionInterpolatedP0(const unsigned int n_components, const double init_time)
: FunctionBase<dim>( n_components, init_time)
{

}



template <int dim>
Input::Type::Record &FunctionInterpolatedP0<dim>::get_input_type() {
    using namespace  Input::Type;

    static Record rec("FunctionInterpolatedP0", "Function given by P0 data on another mesh. Currently defined only on boundary.");

    if (! rec.is_finished()) {
        rec.derive_from(FunctionBase<dim>::get_input_type());
        // TODO: use mesh record here, but make it possibly of the same simplicity (the file name only)
        rec.declare_key("mesh", FileName::input(),Default::obligatory(),
                "File with the mesh from which we interpolate. (currently only GMSH supported)");
        // TODO: allow interpolation from VTK files (contains also mesh), and our own format of raw data, that includes:
        // mesh, dof_handler, and dof values
        rec.declare_key("raw_data", FileName::input(), Default::obligatory(),
                "File with raw output from flow calculation. Currently we can interpolate only pressure.");
        rec.finish();
    }
    return rec;
}



template <int dim>
void FunctionInterpolatedP0<dim>::init_from_input(Input::Record rec) {
    set_source_of_interpolation(
            rec.val<FilePath>("mesh"),
            rec.val<FilePath>("raw_data")
            );
}




template <int dim>
void FunctionInterpolatedP0<dim>::set_element(Element *element){
	element_ = element;

	switch (element_->dim()) {
		case 1: {
			TAbscissa abscissa;
			createAbscissa(element, abscissa);
			calculate_abscissa_pressure(abscissa);
			break;
		}
		case 2: {
			TTriangle triangle;
			createTriangle(element, triangle);
			calculate_triangle_pressure(triangle);
			break;
		}
		default:
			xprintf(Err, "Dimension of element must be 1 or 2!\n");
			break;
	}
}



template <int dim>
void FunctionInterpolatedP0<dim>::set_source_of_interpolation(const FilePath & mesh_file,
		const FilePath & raw_output) {

	// read mesh, create tree
    {
       mesh_ = new Mesh();
	   GmshMeshReader reader;
	   reader.read( mesh_file, mesh_);
	   // no call to mesh->setup_topology, we need only elements, no connectivity
    }
	bihTree_ = new BIHTree(mesh_);

	// read pressures
	FILE* raw_output_file = xfopen( string(raw_output).c_str(), "rt");
	read_pressures(raw_output_file);
	xfclose(raw_output_file);

	//calculate_interpolation();
}



template <int dim>
void FunctionInterpolatedP0<dim>::read_pressures(FILE* raw_output) {
	xprintf(Msg, " - FunctionInterpolatedP0->read_pressures(FILE* raw_output)\n");

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

/*

template <int dim>
void FunctionInterpolatedP0<dim>::read_element_data_from_gmsh(istream &in, const  string &field_name) {
    xprintf(Msg, "FunctionInterpolatedP0->read_element_data_from_gmsh.\n");

    string line;
    while (1) {
        skip_to(in, "$ElementData");
        std::getline( in, line);
        unsigned int n_str = lexical_cast<unsigned int> (line);
        if (n_str == 0) continue;
        std::getline( in, line);
        if (line != field_name) continue;


    skip_to(in, filed_name)
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

*/

template <int dim>
void FunctionInterpolatedP0<dim>::calculate_triangle_pressure(TTriangle &element) {
	double total_measure, measure;
	BoundingBox elementBoundingBox = element.get_bounding_box();
	
	TIntersectionType iType;
	TTetrahedron tetrahedron;

	START_TIMER("find_elements_2D");
	((BIHTree *)bihTree_)->find_elements(elementBoundingBox, searchedElements_);
	END_TIMER("find_elements_2D");

	total_measure = 0.0;
	pressure_ = 0.0;

        START_TIMER("compute_pressure_2D");
    ADD_CALLS( searchedElements_.size());
	for (std::vector<int>::iterator it = searchedElements_.begin(); it!=searchedElements_.end(); it++)
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



template <int dim>
void FunctionInterpolatedP0<dim>::calculate_abscissa_pressure(TAbscissa &element) {
        double total_measure, measure;
	BoundingBox elementBoundingBox = element.get_bounding_box();
	TIntersectionType iType;
	TTetrahedron tetrahedron;

	START_TIMER("find_elements_1D");
	((BIHTree *)bihTree_)->find_elements(elementBoundingBox, searchedElements_);
	END_TIMER("find_elements_1D");

        total_measure = 0.0;
	pressure_ = 0.0;
	START_TIMER("compute_pressure_1D");
	ADD_CALLS(searchedElements_.size());
	for (std::vector<int>::iterator it = searchedElements_.begin(); it!=searchedElements_.end(); it++)
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



template <int dim>
void FunctionInterpolatedP0<dim>::createTetrahedron(Element *ele, TTetrahedron &te) {
	ASSERT(( ele->dim() == 3 ), "Dimension of element must be 3!\n");

	te.SetPoints(TPoint(ele->node[0]->point()(0), ele->node[0]->point()(1), ele->node[0]->point()(2)),
				TPoint(ele->node[1]->point()(0), ele->node[1]->point()(1), ele->node[1]->point()(2)),
				TPoint(ele->node[2]->point()(0), ele->node[2]->point()(1), ele->node[2]->point()(2)),
				TPoint(ele->node[3]->point()(0), ele->node[3]->point()(1), ele->node[3]->point()(2)) );
}



template <int dim>
void FunctionInterpolatedP0<dim>::createTriangle(Element *ele, TTriangle &tr) {
	ASSERT(( ele->dim() == 2 ), "Dimension of element must be 2!\n");

	tr.SetPoints(TPoint(ele->node[0]->point()(0), ele->node[0]->point()(1), ele->node[0]->point()(2)),
				 TPoint(ele->node[1]->point()(0), ele->node[1]->point()(1), ele->node[1]->point()(2)),
				 TPoint(ele->node[2]->point()(0), ele->node[2]->point()(1), ele->node[2]->point()(2)) );
}



template <int dim>
void FunctionInterpolatedP0<dim>::createAbscissa(Element *ele, TAbscissa &ab) {
	ASSERT(( ele->dim() == 1 ), "Dimension of element must be 1!\n");

	ab.SetPoints(TPoint(ele->node[0]->point()(0), ele->node[0]->point()(1), ele->node[0]->point()(2)),
			 	 TPoint(ele->node[1]->point()(0), ele->node[1]->point()(1), ele->node[1]->point()(2)) );
}



template <int dim>
double FunctionInterpolatedP0<dim>::value(const Point &p, const unsigned int component) const
{
	return pressure_;
}



template <int dim>
void FunctionInterpolatedP0<dim>::vector_value(const Point &p, std::vector<double> &value) const
{
	xprintf(Msg, " - Method vector_value is not used and implemented in class FunctionInterpolatedP0\n");
}



template <int dim>
void FunctionInterpolatedP0<dim>::value_list(const std::vector<Point>  &point_list,
					  std::vector<double>         &value_list,
					  const unsigned int  component) const
{
	xprintf(Msg, " - Method value_list is not used and implemented in class FunctionInterpolatedP0\n");
}



template <int dim>
void FunctionInterpolatedP0<dim>::vector_value_list (const std::vector<Point> &point_list,
                            std::vector< std::vector<double> > &value_list) const
{
	xprintf(Msg, " - Method vector_value_list is not used and implemented in class FunctionInterpolatedP0\n");
}


#endif
