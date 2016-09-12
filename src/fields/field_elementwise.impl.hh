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
 * @file    field_elementwise.impl.hh
 * @brief   
 */

#ifndef FIELD_ELEMENTWISE_IMPL_HH_
#define FIELD_ELEMENTWISE_IMPL_HH_


#include "fields/field_elementwise.hh"
#include "system/file_path.hh"
#include "input/input_type.hh"
#include "mesh/msh_gmshreader.h"
#include "mesh/reader_instances.hh"

/// Implementation.

namespace IT = Input::Type;

FLOW123D_FORCE_LINK_IN_CHILD(field_elementwise)


template <int spacedim, class Value>
const Input::Type::Record & FieldElementwise<spacedim, Value>::get_input_type()
{
    return it::Record("FieldElementwise", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field constant in space.")
        .derive_from(FieldAlgorithmBase<spacedim, Value>::get_input_type())
        .declare_key("gmsh_file", IT::FileName::input(), IT::Default::obligatory(),
                "Input file with ASCII GMSH file format.")
        .declare_key("field_name", IT::String(), IT::Default::obligatory(),
                "The values of the Field are read from the ```$ElementData``` section with field name given by this key.")
        .close();
}


template <int spacedim, class Value>
const int FieldElementwise<spacedim, Value>::registrar =
		Input::register_class< FieldElementwise<spacedim, Value>, unsigned int >("FieldElementwise") +
		FieldElementwise<spacedim, Value>::get_input_type().size();



template <int spacedim, class Value>
FieldElementwise<spacedim, Value>::FieldElementwise( unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>(n_comp),
  internal_raw_data(true), mesh_(NULL)

{
    n_components_ = this->value_.n_rows() * this->value_.n_cols();
}



template <int spacedim, class Value>
FieldElementwise<spacedim, Value>::FieldElementwise(std::shared_ptr< std::vector<typename Value::element_type> > data,
		unsigned int n_components)
: FieldAlgorithmBase<spacedim, Value>(n_components),
internal_raw_data(false), mesh_(NULL)
{
	n_components_ = this->value_.n_rows() * this->value_.n_cols();
	data_ = data;
}




template <int spacedim, class Value>
void FieldElementwise<spacedim, Value>::init_from_input(const Input::Record &rec) {
	cout << string(reader_file_) << endl;
	OLD_ASSERT( internal_raw_data, "Trying to initialize internal FieldElementwise from input.");
	OLD_ASSERT( reader_file_ == FilePath(), "Multiple call of init_from_input.\n");
    reader_file_ = FilePath( rec.val<FilePath>("gmsh_file") );
    ReaderInstances::instance()->get_reader(reader_file_);

    field_name_ = rec.val<std::string>("field_name");
}



template <int spacedim, class Value>
void FieldElementwise<spacedim, Value>::set_data_row(unsigned int boundary_idx, typename Value::return_type &value) {
    Value ref(value);
    OLD_ASSERT( this->value_.n_cols() == ref.n_cols(), "Size of variable vectors do not match.\n" );
    OLD_ASSERT( mesh_, "Null mesh pointer of elementwise field: %s, did you call set_mesh()?\n", field_name_.c_str());
    OLD_ASSERT( boundary_domain_ , "Method set_data_row can be used only for boundary fields.");
    unsigned int vec_pos = boundary_idx * n_components_;
    std::vector<typename Value::element_type> &vec = *( data_.get() );
    for(unsigned int row=0; row < ref.n_rows(); row++)
        for(unsigned int col=0; col < ref.n_cols(); col++, vec_pos++)
        	vec[vec_pos] = ref(row,col);

}


template <int spacedim, class Value>
bool FieldElementwise<spacedim, Value>::set_time(const TimeStep &time) {
	OLD_ASSERT(mesh_, "Null mesh pointer of elementwise field: %s, did you call set_mesh()?\n", field_name_.c_str());
    if ( reader_file_ == FilePath() ) return false;

    //walkaround for the steady time governor - there is no data to be read in time==infinity
    //TODO: is it possible to check this before calling set_time?
    //if (time.end() == numeric_limits< double >::infinity()) return false;
    
    GMSH_DataHeader search_header;
    search_header.actual=false;
    search_header.field_name=field_name_;
    search_header.n_components=n_components_;
    search_header.n_entities=n_entities_;
    search_header.time=time.end();


    data_ = ReaderInstances::instance()->get_reader(reader_file_)-> template get_element_data<typename Value::element_type>(search_header,
    		mesh_->elements_id_maps(boundary_domain_), this->component_idx_);
    return search_header.actual;
}



template <int spacedim, class Value>
void FieldElementwise<spacedim, Value>::set_mesh(const Mesh *mesh, bool boundary_domain) {
    // set mesh only once or to same value
	OLD_ASSERT(mesh_ == nullptr || mesh_ == mesh, "Trying to change mesh of the FieldElementwise.");
    boundary_domain_ = boundary_domain;

    mesh_=mesh;
    if (boundary_domain_) {
        n_entities_=mesh_->bc_elements.size();
    } else {
        n_entities_=mesh_->n_elements();
    }

    // allocate
    if (!data_) {
    	data_ = std::make_shared<std::vector<typename Value::element_type>>();
    	data_->resize(n_entities_ * n_components_);
    }

}



/**
 * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
 */
template <int spacedim, class Value>
typename Value::return_type const & FieldElementwise<spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
        OLD_ASSERT( elm.is_elemental(), "FieldElementwise works only for 'elemental' ElementAccessors.\n");
        OLD_ASSERT( elm.is_boundary() == boundary_domain_, "Trying to get value of FieldElementwise '%s' for wrong ElementAccessor type (boundary/bulk).\n", field_name_.c_str() );

        unsigned int idx = n_components_*elm.idx();
        std::vector<typename Value::element_type> &vec = *( data_.get() );

        return Value::from_raw(this->r_value_, (typename Value::element_type *)(&vec[idx]));
}



/**
 * Returns std::vector of scalar values in several points at once.
 */
template <int spacedim, class Value>
void FieldElementwise<spacedim, Value>::value_list (const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{
	OLD_ASSERT( elm.is_elemental(), "FieldElementwise works only for 'elemental' ElementAccessors.\n");
	OLD_ASSERT( elm.is_boundary() == boundary_domain_, "Trying to get value of FieldElementwise '%s' for wrong ElementAccessor type (boundary/bulk).\n", field_name_.c_str() );
	OLD_ASSERT_EQUAL( point_list.size(), value_list.size() );
    if (boost::is_floating_point< typename Value::element_type>::value) {
        unsigned int idx = n_components_*elm.idx();
        std::vector<typename Value::element_type> &vec = *( data_.get() );

        typename Value::return_type const &ref = Value::from_raw(this->r_value_, (typename Value::element_type *)(&vec[idx]));
        for(unsigned int i=0; i< value_list.size(); i++) {
        	OLD_ASSERT( Value(value_list[i]).n_rows()==this->value_.n_rows(),
                    "value_list[%d] has wrong number of rows: %d; should match number of components: %d\n",
                    i, Value(value_list[i]).n_rows(),this->value_.n_rows());

            value_list[i] = ref;
        }
    } else {
        xprintf(UsrErr, "FieldElementwise is not implemented for discrete return types.\n");
    }
}



template <int spacedim, class Value>
FieldElementwise<spacedim, Value>::~FieldElementwise() {}



#endif /* FIELD_ELEMENTWISE_IMPL_HH_ */
