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
#include "mesh/bih_tree.hh"
#include "mesh/ngh/include/intersection.h"
#include "mesh/ngh/include/point.h"
#include "system/sys_profiler.hh"
//#include "boost/lexical_cast.hpp"
//#include "system/tokenizer.hh"
//#include "system/xio.h"


namespace it = Input::Type;
/*template <int spacedim, class Value>
it::Record FieldInterpolatedP0<spacedim, Value>::input_type
    = it::Record("FieldInterpolatedP0", "Field given by P0 data on another mesh. Currently defined only on boundary.")
	.derive_from(FieldBase<spacedim, Value>::input_type)
	// TODO: use mesh record here, but make it possibly of the same simplicity (the file name only)
	.declare_key("mesh", it::FileName::input(),it::Default::obligatory(),
			"File with the mesh from which we interpolate. (currently only GMSH supported)")
	// TODO: allow interpolation from VTK files (contains also mesh), and our own format of raw data, that includes:
	// mesh, dof_handler, and dof values
	.declare_key("raw_data", it::FileName::input(), it::Default::obligatory(),
			"File with raw output from flow calculation. Currently we can interpolate only pressure."); */

template <int spacedim, class Value>
it::Record FieldInterpolatedP0<spacedim, Value>::input_type
    = FieldInterpolatedP0<spacedim, Value>::get_input_type(FieldBase<spacedim, Value>::input_type, NULL);


template <int spacedim, class Value>
Input::Type::Record FieldInterpolatedP0<spacedim, Value>::get_input_type(
        Input::Type::AbstractRecord &a_type, typename Value::ElementInputType *eit
        )
{
    it::Record type=
        it::Record("FieldInterpolatedP0", FieldBase<spacedim,Value>::template_name()+" Field constant in space.")
        .derive_from(a_type)
        .declare_key("gmsh_file", IT::FileName::input(), IT::Default::obligatory(),
                "Input file with ASCII GMSH file format.")
        .declare_key("field_name", IT::String(), IT::Default::obligatory(),
                "The values of the Field are read from the $ElementData section with field name given by this key.")
        .close();

    return type;
}


template <int spacedim, class Value>
FieldInterpolatedP0<spacedim, Value>::FieldInterpolatedP0(const unsigned int n_comp)
: FieldBase<spacedim, Value>(n_comp)
{}




template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::init_from_input(const Input::Record &rec) {

	// read mesh, create tree
    {
       mesh_ = new Mesh();
       reader_ = new GmshMeshReader( rec.val<FilePath>("gmsh_file") );
       reader_->read_mesh( mesh_);
	   // no call to mesh->setup_topology, we need only elements, no connectivity
    }
	bih_tree_ = new BIHTree(mesh_);

    // allocate data_
	unsigned int data_size = (mesh_->element.size() + mesh_->bc_elements.size()) * (this->value_.n_rows() * this->value_.n_cols());
    data_ = new double[data_size];
    std::fill(data_, data_ + data_size, 0.0);


	field_name_ = rec.val<std::string>("field_name");
}




/*template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::set_source_of_interpolation(const FilePath & mesh_file,
		const FilePath & raw_output) {

	// read mesh, create tree
    {
       mesh_ = new Mesh();
       reader_ = new GmshMeshReader(mesh_file);
	   reader_->read_mesh( mesh_);
	   // no call to mesh->setup_topology, we need only elements, no connectivity
    }
    bih_tree_ = new BIHTree(mesh_);

	// read pressures
	Tokenizer tok(mesh_file);

	FILE* raw_output_file = xfopen( string(raw_output).c_str(), "rt");
	read_pressures(raw_output_file);
	xfclose(raw_output_file);
	//read_element_data_from_gmsh(tok, "xx");
}*/


/*template <int spacedim, class Value>
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
}*/


template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::calculate_triangle_pressure(TTriangle &element) {
	double total_measure, measure;
	BoundingBox elementBoundingBox = element.get_bounding_box();
	
	TIntersectionType iType;
	TTetrahedron tetrahedron;

	START_TIMER("find_elements_2D");
	((BIHTree *)bih_tree_)->find_bounding_box(elementBoundingBox, searched_elements_);
	END_TIMER("find_elements_2D");

	total_measure = 0.0;
	pressure_ = 0.0;

    START_TIMER("compute_pressure_2D");
    ADD_CALLS( searched_elements_.size());
	for (std::vector<unsigned int>::iterator it = searched_elements_.begin(); it!=searched_elements_.end(); it++)
	{
		int idx = *it;
		ElementFullIter ele = mesh_->element( idx );
		if (ele->dim() == 3) {
				createTetrahedron(ele, tetrahedron);
				GetIntersection(element, tetrahedron, iType, measure);
				if (iType == area) {
						pressure_ += data_[ idx ] * measure;
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
	((BIHTree *)bih_tree_)->find_bounding_box(elementBoundingBox, searched_elements_);
	END_TIMER("find_elements_1D");

        total_measure = 0.0;

    /**
     * TODO:
     * nahradit  pressure_ -> value_
     *
     * for(unsigned int i=0; i < value_.n_rows(); i ++)
     *      for( ... n_cols() )
     *          value_.at(i,j) = 0.0;
     *
     * Value tmp_value;
     * Value::from_raw(tmp_value, (typename Value::element_type *)(data_+idx));
     * for(unsigned int i=0; i < value_.n_rows(); i ++)
     *      for( ... n_cols() )
     *          value_.at(i,j) += measure * tmp_value.at(i,j);
     *
     * ?? spojit caluculate_abscissa a calculate_triangle
     */
	pressure_ = 0.0;
	START_TIMER("compute_pressure_1D");
	ADD_CALLS(searched_elements_.size());
	for (std::vector<unsigned int>::iterator it = searched_elements_.begin(); it!=searched_elements_.end(); it++)
	{
		int idx = *it;
		ElementFullIter ele = mesh_->element( idx );
		if (ele->dim() == 3) {
			createTetrahedron(ele, tetrahedron);
			GetIntersection(element, tetrahedron, iType, measure);
			if (iType == line) {
				pressure_ += data_[ idx ] * measure;
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
bool FieldInterpolatedP0<spacedim, Value>::set_time(double time) {
    ASSERT(mesh_, "Null mesh pointer of elementwise field: %s, did you call init_from_input(Input::Record)?\n", field_name_.c_str());
    ASSERT(data_, "Null data pointer.\n");
    if (reader_ == NULL) return false;

    GMSH_DataHeader search_header;
    search_header.actual = false;
    search_header.field_name = field_name_;
    search_header.n_components = this->value_.n_rows() * this->value_.n_cols();
    search_header.n_entities = mesh_->element.size() + mesh_->bc_elements.size();
    search_header.time = time;

    bool boundary_domain_=false;
    reader_->read_element_data(search_header, data_, mesh_->elements_id_maps(boundary_domain_)  );

    return search_header.actual;
}

/*
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
*/


template <int spacedim, class Value>
typename Value::return_type const &FieldInterpolatedP0<spacedim, Value>::value(const Point<spacedim> &p, const ElementAccessor<spacedim> &elm)
{
	if (&elm != computed_elm_) {
		//cout << "first computing" << endl;
		computed_elm_ = &elm;

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
	}
	//else {
	//	cout << "second computing" << endl;
	//}

	this->value_(0,0) = pressure_;
    return this->r_value_;
}


template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::value_list(const std::vector< Point<spacedim> >  &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list)
{
	// not supported yet
}





#endif
