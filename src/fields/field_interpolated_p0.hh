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
 * @file    field_interpolated_p0.hh
 * @brief   
 */

#ifndef FIELD_INTERPOLATED_P0_HH_
#define FIELD_INTERPOLATED_P0_HH_

#include <string.h>                           // for memcpy
#include <boost/exception/info.hpp>           // for operator<<, error_info:...
#include <limits>                             // for numeric_limits
#include <memory>                             // for shared_ptr
#include <string>                             // for string
#include <vector>                             // for vector
#include <armadillo>
#include "field_algo_base.hh"                 // for FieldAlgorithmBase
#include "fields/field_values.hh"             // for FieldValue<>::Enum, Fie...
#include "input/accessors.hh"                 // for ExcAccessorForNullStorage
#include "input/accessors_impl.hh"            // for Record::val
#include "input/storage.hh"                   // for ExcStorageTypeMismatch
#include "input/type_record.hh"               // for Record::ExcRecordKeyNot...
#include "mesh/element_impls.hh"              // for Element::dim
#include "system/exceptions.hh"               // for ExcAssertMsg::~ExcAsser...
#include "system/file_path.hh"                // for FilePath
#include "tools/time_governor.hh"             // for TimeStep
class BIHTree;
class Mesh;
template <int spacedim> class ElementAccessor;


template <int spacedim, class Value>
class FieldInterpolatedP0: public FieldAlgorithmBase<spacedim, Value> {
public:

    typedef typename FieldAlgorithmBase<spacedim, Value>::Point Point;
    typedef FieldAlgorithmBase<spacedim, Value> FactoryBaseType;

	/**
	 * Constructor
	 */
	FieldInterpolatedP0(unsigned int n_comp=0);

	/**
	 * Declare Input type.
	 */
	static const Input::Type::Record & get_input_type();

	/**
	 * Initialization from the input interface.
	 */
	virtual void init_from_input(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data);

    /**
     * Update time and possibly update data from GMSH file.
     */
    bool set_time(const TimeStep &time) override;

    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void value_list(const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);

protected:
    /// mesh, which is interpolated
	std::shared_ptr<Mesh> source_mesh_;

	/// mesh reader file
	FilePath reader_file_;

    /// Raw buffer of n_entities rows each containing Value::size() doubles.
	std::shared_ptr< std::vector<typename Value::element_type> > data_;

	/// vector stored suspect elements in calculating the intersection
	std::vector<unsigned int> searched_elements_;

	/// field name read from input
	std::string field_name_;

	/// tree of mesh elements
	BIHTree* bih_tree_;

	/// stored index to last computed element
	unsigned int computed_elm_idx_ = numeric_limits<unsigned int>::max();

	/// stored flag if last computed element is boundary
	unsigned int computed_elm_boundary_;


    /// Default value of element if not set in mesh data file
    double default_value_;

    /// Accessor to Input::Record
    Input::Record in_rec_;
private:
    /// Registrar of class to factory
    static const int registrar;

};



#endif /* FUNCTION_INTERPOLATED_P0_HH_ */
