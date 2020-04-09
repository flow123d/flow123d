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
 * @file    type_tuple.cc
 * @brief
 */

#include "input_type.hh"
#include "type_repository.hh"

#include <boost/functional/hash.hpp>


namespace Input {
namespace Type {



Tuple::Tuple()
: Record() {}


Tuple::Tuple(const Tuple & other)
: Record(other) {}


Tuple::Tuple(const string & type_name_in, const string & description)
{
	this->data_ = std::make_shared<RecordData>(type_name_in, description);
}


TypeBase::TypeHash Tuple::content_hash() const
{
	TypeHash seed=0;
    boost::hash_combine(seed, "Tuple");
    data_->content_hash(seed);
    return seed;
}


string Tuple::class_name() const {
	return "Tuple";
}


Tuple & Tuple::allow_auto_conversion(const string &)
{
	ASSERT(false)(this->type_name()).error("Call of allow_auto_conversion method is forbidden for type Tuple");
	return *this;
}


bool Tuple::operator==(const TypeBase &other) const
{ return  typeid(*this) == typeid(other) &&
                 (type_name() == static_cast<const Tuple *>(&other)->type_name() );
}


const Tuple &Tuple::close() const {
	ASSERT_GT(data_->keys.size(), 0)(this->type_name()).error("Empty Tuple!\n");
    data_->closed_=true;

    bool allow_auto_conversion = true; // if no key or only first key is obligatory, tuple is auto convertible
    for (vector<Key>::iterator it=data_->keys.begin(); it!=data_->keys.end(); it++)
		if ( it->default_.is_obligatory() ) {
			if (it != data_->keys.begin()) {
				allow_auto_conversion = false;
			}
		}
    if (allow_auto_conversion) { // Add autoconvertibility
        data_->auto_conversion_key_idx = 0;
        data_->auto_conversion_key=data_->keys.begin()->key_;
    }

    return *( Input::TypeRepository<Tuple>::get_instance().add_type( *this ) );
}


FinishStatus Tuple::finish(FinishStatus finish_type)
{
	ASSERT(finish_type != FinishStatus::none_).error();
	ASSERT(finish_type != FinishStatus::in_perform_).error();
	ASSERT(data_->finish_status_ != FinishStatus::in_perform_)(this->type_name()).error("Recursion in the IST element of type Tuple.");

	if (this->is_finished()) return data_->finish_status_;

	ASSERT(data_->closed_)(this->type_name()).error();

	data_->finish_status_ = FinishStatus::in_perform_;

    // iterates through keys
    bool obligatory_keys = true; // check order of keys (at first obligatory keys are defined, then other keys)
    for (vector<Key>::iterator it=data_->keys.begin(); it!=data_->keys.end(); it++) {
    	// Performs check order of keys and check auto-conversion
    	Default dflt = it->default_;
		if ( dflt.is_obligatory() ) {
			if ( !obligatory_keys ) {
				THROW( ExcTupleWrongKeysOrder() << EI_TupleName(this->type_name()) );
			}
		} else {
			obligatory_keys = false;
		}

		// Performs finish of keys
		if (typeid( *(it->type_.get()) ) == typeid(Instance)) {
			it->type_->finish(FinishStatus::generic_); // finish Instance object
			it->type_ = it->type_->make_instance().first;
		}
		if ((finish_type != FinishStatus::generic_) && it->type_->is_root_of_generic_subtree())
			THROW( ExcGenericWithoutInstance() << EI_Object(it->type_->type_name()) );
		it->type_->finish(finish_type);
		ASSERT(it->type_->is_finished()).error();
		if (finish_type == FinishStatus::delete_) it->type_.reset();
    }

    data_->finish_status_ = finish_type;
    return (data_->finish_status_);
}


Tuple &Tuple::derive_from(Abstract &)
{
	ASSERT(false)(this->type_name()).error("Call of derive_from method is forbidden for type Tuple");
	return *this;
}


unsigned int Tuple::obligatory_keys_count() const {
	unsigned int obligatory_keys_count=0;
	for ( KeyIter it= this->begin(); it != this->end(); ++it) {
		if ( it->default_.is_obligatory() ) ++obligatory_keys_count;
	}
	return obligatory_keys_count;
}


TypeBase::MakeInstanceReturnType Tuple::make_instance(std::vector<ParameterPair> vec)  {
	Tuple instance_tuple = this->deep_copy();
    auto parameter_map = this->set_instance_data(instance_tuple, vec);
    return std::make_pair( std::make_shared<Tuple>(instance_tuple.close()), parameter_map );
}


Tuple Tuple::deep_copy() const {
	Tuple tuple = Tuple();
	tuple.data_ =  std::make_shared<Record::RecordData>(*this->data_);
	tuple.data_->closed_ = false;
	tuple.data_->finish_status_ = FinishStatus::none_;
	tuple.generic_type_hash_ = this->generic_type_hash_;
	tuple.parameter_map_ = this->parameter_map_;
	tuple.copy_attributes(*attributes_);
	return tuple;
}


Tuple &Tuple::root_of_generic_subtree() {
	root_of_generic_subtree_ = true;
	return *this;
}


Tuple &Tuple::declare_key(const string &key, std::shared_ptr<TypeBase> type,
                        const Default &default_value, const string &description,
                        TypeBase::attribute_map key_attributes)
{
    Record::declare_key(key, type, default_value, description, key_attributes);
    return *this;
}


template <class KeyType>
Tuple &Tuple::declare_key(const string &key, const KeyType &type,
                        const Default &default_value, const string &description,
                        TypeBase::attribute_map key_attributes)
{
    Record::declare_key(key, type, default_value, description, key_attributes);
    return *this;
}



template <class KeyType>
Tuple &Tuple::declare_key(const string &key, const KeyType &type,
                        const string &description,
                        TypeBase::attribute_map key_attributes)
{
    Record::declare_key(key, type, description, key_attributes);
    return *this;
}



// explicit instantiation of template methods

#define TUPLE_DECLARE_KEY(TYPE) \
template Tuple & Tuple::declare_key<TYPE>(const string &key, const TYPE &type, const Default &default_value, const string &description, TypeBase::attribute_map key_attributes); \
template Tuple & Tuple::declare_key<TYPE>(const string &key, const TYPE &type, const string &description, TypeBase::attribute_map key_attributes)


TUPLE_DECLARE_KEY(String);
TUPLE_DECLARE_KEY(Integer);
TUPLE_DECLARE_KEY(Double);
TUPLE_DECLARE_KEY(Bool);
TUPLE_DECLARE_KEY(FileName);
TUPLE_DECLARE_KEY(Selection);
TUPLE_DECLARE_KEY(Array);
TUPLE_DECLARE_KEY(Record);
TUPLE_DECLARE_KEY(Abstract);
TUPLE_DECLARE_KEY(AdHocAbstract);
TUPLE_DECLARE_KEY(Parameter);
TUPLE_DECLARE_KEY(Instance);
TUPLE_DECLARE_KEY(Tuple);


} // closing namespace Type
} // closing namespace Input
