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
 * @file    type_abstract.cc
 * @brief
 */

#include "input_type.hh"
#include "type_repository.hh"
#include "attribute_lib.hh"

#include "system/system.hh"

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/functional/hash.hpp>

namespace Input {
namespace Type {

using namespace std;


/************************************************
 * implementation of Abstract
 */

Abstract::Abstract()
: child_data_( boost::make_shared<ChildData>( "EmptyAbstract", "" ) )
{
	close();
	finish();
}



Abstract::Abstract(const Abstract& other)
: TypeBase( other ), child_data_(other.child_data_)
{}



Abstract::Abstract(const string & type_name_in, const string & description)
: child_data_( boost::make_shared<ChildData>( type_name_in, description ) )
{}


TypeBase::TypeHash Abstract::content_hash() const
{
	TypeHash seed=0;
    boost::hash_combine(seed, "Abstract");
    boost::hash_combine(seed, type_name());
    boost::hash_combine(seed, child_data_->description_);
	for (auto &param : parameter_map_) {
		boost::hash_combine(seed, param.first );
		boost::hash_combine(seed, param.second );
	}
	boost::hash_combine(seed, this->generic_type_hash_);
    //boost::hash_combine(seed, child_data_->generic_content_hash_);
    //for( Record &key : child_data_->list_of_childs) {
    //    boost::hash_combine(seed, key.content_hash() );
    //}
    return seed;
}



Abstract & Abstract::allow_auto_conversion(const string &type_default) {
    if (child_data_->closed_) xprintf(PrgErr, "Can not specify default value for TYPE key as the Abstract '%s' is closed.\n", type_name().c_str());
    child_data_->selection_default_=Default("\""+type_default+"\""); // default record is closed; other constructor creates the zero item
    return *this;
}



const Record  & Abstract::get_descendant(const string& name) const
{
    ASSERT( child_data_->selection_of_childs->is_finished(), "Can not get descendant of unfinished AbstractType\n");
    return child_data_->list_of_childs[ child_data_->selection_of_childs->name_to_int(name) ];
}



const Record * Abstract::get_default_descendant() const {
    if ( have_default_descendant() ) {
    	int sel_val = child_data_->selection_default_.get_storage( child_data_->selection_of_childs )->get_int();
        return &( get_descendant( child_data_->selection_of_childs->int_to_name(sel_val) ) );
    }
    return NULL;
}



const Selection  & Abstract::get_type_selection() const
{
    return * child_data_->selection_of_childs;
}


unsigned int Abstract::child_size() const {
    return child_data_->list_of_childs.size();
}


int Abstract::add_child(const Record &subrec)
{
    ASSERT( child_data_->closed_, "Can not add descendant to Abstract that is not closed.\n");

    if (std::find(child_data_->list_of_childs.begin(), child_data_->list_of_childs.end(), subrec) == child_data_->list_of_childs.end()) {
        child_data_->selection_of_childs->add_value(child_data_->list_of_childs.size(), subrec.type_name());
        child_data_->list_of_childs.push_back(subrec);
    }

    return 1;
}


bool Abstract::finish(bool is_generic) {
	if (child_data_->finished_) return true;

	ASSERT(child_data_->closed_, "Finished Abstract '%s' must be closed!", this->type_name().c_str());

	child_data_->selection_of_childs->close();

	child_data_->finished_ = true;

	for (auto &child : child_data_->list_of_childs) {
		if (!is_generic && child.is_root_of_generic_subtree()) THROW( ExcGenericWithoutInstance() << EI_Object(child.type_name()) );
		child.add_parent(*this);
       	child_data_->finished_ = child_data_->finished_ && child.finish(is_generic);
	}

    // check validity of possible default value of TYPE key
    if ( have_default_descendant() ) {
		try {
			child_data_->selection_default_.check_validity(child_data_->selection_of_childs);
		} catch (ExcWrongDefault & e) {
			xprintf(PrgErr, "Default value '%s' for TYPE key do not match any descendant of Abstract '%s'.\n", child_data_->selection_default_.value().c_str(), type_name().c_str());
		}
    }

	return (child_data_->finished_);
}


Abstract &Abstract::close() {
	child_data_->closed_=true;
    return *( Input::TypeRepository<Abstract>::get_instance().add_type( *this ) );
}


bool Abstract::is_finished() const {
    return child_data_->finished_;
}


bool Abstract::is_closed() const {
	return child_data_->closed_;
}


string Abstract::type_name() const {
    return child_data_->type_name_;
}


Default &Abstract::get_selection_default() const {
	return child_data_->selection_default_;
}

bool Abstract::have_default_descendant() const {
	// obligatory value if default is not set, see @p selection_default_
    if ( !child_data_->selection_default_.is_obligatory() )  {
        if ( child_data_->selection_default_.has_value_at_declaration() ) {
        	return true;
        }
    }
	return false;
}



TypeBase::MakeInstanceReturnType Abstract::make_instance(std::vector<ParameterPair> vec) const {
	Abstract abstract = this->deep_copy();
	ParameterMap parameter_map;

	// Set close flag - add_child method required closed child_data
	abstract.child_data_->closed_ = true;
	// make instances of all descendant records and add them into instance of abstract
	for (auto &child : child_data_->list_of_childs) {
		MakeInstanceReturnType inst = child.make_instance(vec);
		abstract.add_child( static_cast<Record &>( *(inst.first) ) );
		ParameterMap other_map = inst.second;
		parameter_map.insert(other_map.begin(), other_map.end());
	}
	// Unset close flag - necessary for set parameters
	abstract.child_data_->closed_ = false;

	// Set parameters and generic type as attributes
	abstract.set_parameters_attribute(parameter_map);
	abstract.parameter_map_ = parameter_map;
	abstract.add_attribute(FlowAttributes::generic_type(), this->hash_str() );
	abstract.generic_type_hash_ = this->content_hash();

	return std::make_pair( boost::make_shared<Abstract>(abstract.close()), parameter_map );
}


Abstract Abstract::deep_copy() const {
	Abstract abstract = Abstract();
	abstract.child_data_ =  boost::make_shared<Abstract::ChildData>(*this->child_data_);
	abstract.child_data_->closed_ = false;
	abstract.child_data_->finished_ = false;
	abstract.child_data_->list_of_childs.clear();
	abstract.child_data_->selection_of_childs = boost::make_shared<Selection>(this->type_name() + "_TYPE_selection");
	abstract.attributes_ = boost::make_shared<attribute_map>(*attributes_);
	abstract.generic_type_hash_ = this->generic_type_hash_;
	abstract.parameter_map_ = this->parameter_map_;
	return abstract;
}


Abstract &Abstract::root_of_generic_subtree() {
	root_of_generic_subtree_ = true;
	return *this;
}


/*Abstract &Abstract::set_generic_content_hash(TypeHash generic_content_hash) {
	child_data_->generic_content_hash_ = generic_content_hash;
	return *this;
}*/


Abstract::ChildDataIter Abstract::begin_child_data() const {
    return child_data_->list_of_childs.begin();
}

Abstract::ChildDataIter Abstract::end_child_data() const {
    return child_data_->list_of_childs.end();
}


/************************************************
 * implementation of AdHocAbstract
 */

AdHocAbstract::AdHocAbstract(const Abstract &ancestor)
: Abstract("Derived AdHocAbstract", "This description doesn't have print out."),
  ancestor_(ancestor)
{
	//test default descendant of ancestor
	const Record * default_desc = ancestor.get_default_descendant();
	if (default_desc) {
		allow_auto_conversion( default_desc->type_name() );
	}

	this->close();

}


AdHocAbstract &AdHocAbstract::add_child(const Record &subrec)
{
	Abstract::add_child(subrec);

	return *this;
}


bool AdHocAbstract::finish(bool is_generic)
{
	if (child_data_->finished_) return true;

	const_cast<Abstract &>(ancestor_).finish(is_generic);

	//test default descendant of ancestor
	const Record * default_desc = ancestor_.get_default_descendant();
	if (default_desc) {
		allow_auto_conversion( default_desc->type_name() );
	}

	for (Abstract::ChildDataIter it = ancestor_.child_data_->list_of_childs.begin(); it != ancestor_.child_data_->list_of_childs.end(); ++it) {
	    child_data_->selection_of_childs->add_value(child_data_->list_of_childs.size(), (*it).type_name());
	    child_data_->list_of_childs.push_back(*it);
	}

	return Abstract::finish(is_generic);
}


TypeBase::TypeHash AdHocAbstract::content_hash() const {
	TypeHash seed=0;
    boost::hash_combine(seed, "AdHocAbstract");
    boost::hash_combine(seed, type_name());
    boost::hash_combine(seed, child_data_->description_);
    boost::hash_combine(seed, ancestor_.type_name());

    return seed;
}



} // closing namespace Type
} // closing namespace Input


