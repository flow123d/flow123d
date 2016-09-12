/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    field_interpolated_p0.impl.hh
 * @brief   
 */

#ifndef FIELD_INTERPOLATED_P0_IMPL_HH_
#define FIELD_INTERPOLATED_P0_IMPL_HH_

#include "field_interpolated_p0.hh"
#include "system/system.hh"
#include "mesh/msh_gmshreader.h"
#include "mesh/bih_tree.hh"
#include "mesh/reader_instances.hh"
#include "mesh/ngh/include/intersection.h"
#include "mesh/ngh/include/point.h"
#include "system/sys_profiler.hh"


namespace it = Input::Type;

FLOW123D_FORCE_LINK_IN_CHILD(field_interpolated)



template <int spacedim, class Value>
const Input::Type::Record & FieldInterpolatedP0<spacedim, Value>::get_input_type()
{
    return it::Record("FieldInterpolatedP0", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field constant in space.")
        .derive_from(FieldAlgorithmBase<spacedim, Value>::get_input_type())
        .declare_key("gmsh_file", IT::FileName::input(), IT::Default::obligatory(),
                "Input file with ASCII GMSH file format.")
        .declare_key("field_name", IT::String(), IT::Default::obligatory(),
                "The values of the Field are read from the ```$ElementData``` section with field name given by this key.")
        .close();
}



template <int spacedim, class Value>
const int FieldInterpolatedP0<spacedim, Value>::registrar =
		Input::register_class< FieldInterpolatedP0<spacedim, Value>, unsigned int >("FieldInterpolatedP0") +
		FieldInterpolatedP0<spacedim, Value>::get_input_type().size();


template <int spacedim, class Value>
FieldInterpolatedP0<spacedim, Value>::FieldInterpolatedP0(const unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>(n_comp)
{}



template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::init_from_input(const Input::Record &rec) {

	// read mesh, create tree
    {
       source_mesh_ = new Mesh();
       reader_file_ = FilePath( rec.val<FilePath>("gmsh_file") );
       ReaderInstances::instance()->get_reader(reader_file_)->read_mesh( source_mesh_ );
	   // no call to mesh->setup_topology, we need only elements, no connectivity
    }
	bih_tree_ = new BIHTree( source_mesh_ );

    // allocate data_
	unsigned int data_size = source_mesh_->element.size() * (this->value_.n_rows() * this->value_.n_cols());
	data_ = std::make_shared<std::vector<typename Value::element_type>>();
	data_->resize(data_size);

	field_name_ = rec.val<std::string>("field_name");
}




template <int spacedim, class Value>
bool FieldInterpolatedP0<spacedim, Value>::set_time(const TimeStep &time) {
	OLD_ASSERT(source_mesh_, "Null mesh pointer of elementwise field: %s, did you call init_from_input(Input::Record)?\n", field_name_.c_str());
    if ( reader_file_ == FilePath() ) return false;
    
    //walkaround for the steady time governor - there is no data to be read in time==infinity
    //TODO: is it possible to check this before calling set_time?
    //if (time == numeric_limits< double >::infinity()) return false;
    
    // value of last computed element must be recalculated if time is changed
    computed_elm_idx_ = numeric_limits<unsigned int>::max();

    GMSH_DataHeader search_header;
    search_header.actual = false;
    search_header.field_name = field_name_;
    search_header.n_components = this->value_.n_rows() * this->value_.n_cols();
    search_header.n_entities = source_mesh_->element.size();
    search_header.time = time.end();
    
    bool boundary_domain_ = false;
    data_ = ReaderInstances::instance()->get_reader(reader_file_)->template get_element_data<typename Value::element_type>(search_header,
    		source_mesh_->elements_id_maps(boundary_domain_), this->component_idx_);

    return search_header.actual;
}



template <int spacedim, class Value>
typename Value::return_type const &FieldInterpolatedP0<spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
	OLD_ASSERT( elm.is_elemental(), "FieldInterpolatedP0 works only for 'elemental' ElementAccessors.\n");
	if (elm.idx() != computed_elm_idx_ || elm.is_boundary() != computed_elm_boundary_) {
		computed_elm_idx_ = elm.idx();
		computed_elm_boundary_ = elm.is_boundary();

		if (elm.dim() == 3) {
			xprintf(Err, "Dimension of element in target mesh must be 0, 1 or 2! elm.idx() = %d\n", elm.idx());
		}

		double epsilon = 4* numeric_limits<double>::epsilon() * elm.element()->measure();

		// gets suspect elements
		if (elm.dim() == 0) {
			searched_elements_.clear();
			((BIHTree *)bih_tree_)->find_point(elm.element()->node[0]->point(), searched_elements_);
		} else {
			BoundingBox bb;
			elm.element()->get_bounding_box(bb);
			searched_elements_.clear();
			((BIHTree *)bih_tree_)->find_bounding_box(bb, searched_elements_);
		}

		// set zero values of value_ object
		for (unsigned int i=0; i < this->value_.n_rows(); i++) {
			for (unsigned int j=0; j < this->value_.n_cols(); j++) {
				this->value_(i,j) = 0.0;
			}
		}

		double total_measure=0.0, measure;
		TIntersectionType iType;

		START_TIMER("compute_pressure");
		ADD_CALLS(searched_elements_.size());
		for (std::vector<unsigned int>::iterator it = searched_elements_.begin(); it!=searched_elements_.end(); it++)
		{
			ElementFullIter ele = source_mesh_->element( *it );
			if (ele->dim() == 3) {
			    ngh::set_tetrahedron_from_element(tetrahedron_, ele);
				// get intersection (set measure = 0 if intersection doesn't exist)
				switch (elm.dim()) {
					case 0: {
					    ngh::set_point_from_element(point_, elm.element());
						if ( tetrahedron_.IsInner(point_) ) {
							measure = 1.0;
						} else {
							measure = 0.0;
						}
						break;
					}
					case 1: {
					    ngh::set_abscissa_from_element(abscissa_, elm.element());
						GetIntersection(abscissa_, tetrahedron_, iType, measure);
						if (iType != line) {
							measure = 0.0;
						}
						break;
					}
			        case 2: {
			        	ngh::set_triangle_from_element(triangle_, elm.element());
						GetIntersection(triangle_, tetrahedron_, iType, measure);
						if (iType != area) {
							measure = 0.0;
						}
			            break;
			        }
			    }



				//adds values to value_ object if intersection exists
				if (measure > epsilon) {
					unsigned int index = this->value_.n_rows() * this->value_.n_cols() * (*it);
			        std::vector<typename Value::element_type> &vec = *( data_.get() );
			        typename Value::return_type & ret_type_value = const_cast<typename Value::return_type &>( Value::from_raw(this->r_value_,  (typename Value::element_type *)(&vec[index])) );
					Value tmp_value = Value( ret_type_value );

					for (unsigned int i=0; i < this->value_.n_rows(); i++) {
						for (unsigned int j=0; j < this->value_.n_cols(); j++) {
							this->value_(i,j) += tmp_value(i,j) * measure;
						}
					}
					total_measure += measure;
				}
			}
		}

		// computes weighted average
		if (total_measure > epsilon) {
			for (unsigned int i=0; i < this->value_.n_rows(); i++) {
				for (unsigned int j=0; j < this->value_.n_cols(); j++) {
					this->value_(i,j) /= total_measure;
				}
			}
		} else {
			WarningOut().fmt("Processed element with idx {} is out of source mesh!\n", elm.idx());
		}
		END_TIMER("compute_pressure");

	}
    return this->r_value_;
}



template <int spacedim, class Value>
void FieldInterpolatedP0<spacedim, Value>::value_list(const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list)
{
	OLD_ASSERT( elm.is_elemental(), "FieldInterpolatedP0 works only for 'elemental' ElementAccessors.\n");
    FieldAlgorithmBase<spacedim, Value>::value_list(point_list, elm, value_list);
}





#endif /* FIELD_INTERPOLATED_P0_IMPL_HH_ */
