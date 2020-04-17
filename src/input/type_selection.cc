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
 * @file    type_selection.cc
 * @brief   
 */

#include "input/type_selection.hh"
#include "input/type_repository.hh"
#include <boost/functional/hash.hpp>

namespace Input {
namespace Type {

using std::string;

Selection::Selection()
: data_(std::make_shared<SelectionData>("EmptySelection"))
{
    close();
    finish();
}



Selection::Selection(const Selection& other)
: Scalar(other), data_(other.data_)
{ }



Selection::Selection(const string &name, const string &desc)
: data_(std::make_shared<SelectionData>(name))
{
    data_->description_=desc;
}



Selection &Selection::add_value(const int value, const std::string &key,
        const std::string &description, TypeBase::attribute_map attributes)
{
    ASSERT(!is_closed())(key)(type_name()).error("Declaration of new key in closed Selection.");

    data_->add_value(value, key, description, attributes);
    return *this;
}

Selection &Selection::add_attribute(std::string key, TypeBase::json_string value) {
    this->add_attribute_(key, value);
    return *this;
}


const Selection & Selection::close() const {
    data_->closed_=true;
    return *( Input::TypeRepository<Selection>::get_instance().add_type( *this ) );
}


FinishStatus Selection::finish(FinishStatus finish_type ) {
	ASSERT(finish_type != FinishStatus::none_).error();

	if (this->is_finished()) return data_->finish_status_;

	ASSERT(data_->closed_)(this->type_name()).error();

	data_->finish_status_ = finish_type;
	return data_->finish_status_;
}



TypeBase::TypeHash Selection::content_hash() const
{
    std::size_t seed=0;
    boost::hash_combine(seed, "Selection");
    boost::hash_combine(seed, type_name());
    boost::hash_combine(seed, data_->description_);
    for( Key &key : data_->keys_) {
        boost::hash_combine(seed, key.key_);
        boost::hash_combine(seed, key.description_);
        boost::hash_combine(seed, key.value);
    }
    return seed;
}



FinishStatus Selection::finish_status() const {
    return data_->finish_status_;
}


bool Selection::is_finished() const {
    return data_->finish_status_ != FinishStatus::none_;
}


bool Selection::is_closed() const {
	return data_->closed_;
}

string Selection::type_name() const {
   return data_->type_name_;
}


string Selection::class_name() const {
	return "Selection";
}


bool Selection::operator==(const TypeBase &other) const {
    return typeid(*this) == typeid(other) && (type_name() == static_cast<const Selection *>(&other)->type_name());
}



int Selection::name_to_int(const string &key) const {
    finished_check();
    KeyHash key_h = key_hash(key);
    SelectionData::key_to_index_const_iter it = data_->key_to_index_.find(key_h);
    if (it != data_->key_to_index_.end())
        return (data_->keys_[it->second].value);
    else
        throw ExcSelectionKeyNotFound() << EI_KeyName(key) << EI_Selection(*this);
}


string Selection::int_to_name(const int &val) const {
	finished_check();
	auto it = data_->value_to_index_.find(val);
	if (it != data_->value_to_index_.end())
		return data_->keys_[it->second].key_;
	else
		throw ExcSelectionValueNotFound() << EI_Value(val) << EI_Selection(*this);
}


Selection &Selection::copy_values(const Selection &sel)
{
	for (auto it = sel.begin(); it != sel.end(); ++it)
	{
		int value = it->value;
		while (data_->value_to_index_.find(value) != data_->value_to_index_.end()) value++;
		add_value(value, it->key_, it->description_);
	}

	return *this;
}



string Selection::key_list() const {
    ostringstream os;
    for(unsigned int i=0; i<size(); i++) os << "'" <<data_->keys_[i].key_ << "' ";
    return os.str();
}



// Implements @p TypeBase::make_instance.
TypeBase::MakeInstanceReturnType Selection::make_instance(FMT_UNUSED std::vector<ParameterPair> vec) {
	return std::make_pair( std::make_shared<Selection>(*this), ParameterMap() );
}



Selection::SelectionData::SelectionData(const string &name)
: type_name_(name), closed_(false), finish_status_(FinishStatus::none_)
{}



void Selection::SelectionData::add_value(const int value, const std::string &key,
        const std::string &description, TypeBase::attribute_map attributes) {
    KeyHash key_h = TypeBase::key_hash(key);
    ASSERT(key_to_index_.find(key_h) == key_to_index_.end())(key)(type_name_).error("Key already exists in Selection.");
    value_to_index_const_iter it = value_to_index_.find(value);
    const std::string &existing_key = keys_[it->second].key_;
    ASSERT(it == value_to_index_.end())(value)(key)(existing_key)(type_name_)
    		.error("Value of new key conflicts with value of existing key in Selection");

    unsigned int new_idx = key_to_index_.size();
    key_to_index_.insert(std::make_pair(key_h, new_idx));
    value_to_index_.insert(std::make_pair(value, new_idx));

    Key tmp_key = { new_idx, key, description, value, attributes };
    keys_.push_back(tmp_key);
}





} // closing namespace Type
} // close namespace Input
