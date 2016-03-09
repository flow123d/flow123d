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
: Record(type_name_in, description) {}


TypeBase::TypeHash Tuple::content_hash() const
{
	TypeHash seed=0;
    boost::hash_combine(seed, "Tuple");
    data_->content_hash(seed);
    return seed;
}


Tuple & Tuple::allow_auto_conversion(const string &from_key)
{
	THROW( ExcForbiddenTupleMethod() << EI_Method("allow_auto_conversion") << EI_TupleName(this->type_name()) );
	return *this;
}


bool Tuple::operator==(const TypeBase &other) const
{ return  typeid(*this) == typeid(other) &&
                 (type_name() == static_cast<const Tuple *>(&other)->type_name() );
}


const Tuple &Tuple::close() const {
    data_->closed_=true;
    return *( Input::TypeRepository<Tuple>::get_instance().add_type( *this ) );
}


bool Tuple::finish(bool is_generic)
{
	if (data_->finished) return true;

	ASSERT(data_->closed_, "Finished Tuple '%s' must be closed!", this->type_name().c_str());

    data_->finished = true;

    // iterates through keys
    bool obligatory_keys = true; // check order of keys (at first obligatory keys are defined, then other keys)
    bool allow_auto_conversion = true;
    for (vector<Key>::iterator it=data_->keys.begin(); it!=data_->keys.end(); it++) {
    	// Performs check order of keys and check auto-conversion
    	Default dflt = it->default_;
		if ( dflt.is_obligatory() ) {
			if (it != data_->keys.begin()) {
				allow_auto_conversion = false;
			}
			if ( !obligatory_keys ) {
				THROW( ExcTupleWrongKeysOrder() << EI_TupleName(this->type_name()) );
			}
		} else {
			obligatory_keys = false;
		}

		// Performs finish of keys
		if (typeid( *(it->type_.get()) ) == typeid(Instance)) it->type_ = it->type_->make_instance().first;
		if (!is_generic && it->type_->is_root_of_generic_subtree()) THROW( ExcGenericWithoutInstance() << EI_Object(it->type_->type_name()) );
       	data_->finished = data_->finished && it->type_->finish(is_generic);
    }

    // Add autoconvertibility
    if (allow_auto_conversion) {
        data_->auto_conversion_key_idx = 0;
        data_->auto_conversion_key=data_->keys.begin()->key_;
    }

    return (data_->finished);
}


Tuple &Tuple::derive_from(Abstract &parent)
{
	THROW( ExcForbiddenTupleMethod() << EI_Method("derive_from") << EI_TupleName(this->type_name()) );
	return *this;
}


unsigned int Tuple::obligatory_keys_count() const {
	unsigned int obligatory_keys_count=0;
	for ( KeyIter it= this->begin(); it != this->end(); ++it) {
		if ( it->default_.is_obligatory() ) ++obligatory_keys_count;
	}
	return obligatory_keys_count;
}


Tuple &Tuple::declare_key(const string &key, boost::shared_ptr<TypeBase> type,
                        const Default &default_value, const string &description)
{
    Record::declare_key(key, type, default_value, description);
    return *this;
}


template <class KeyType>
Tuple &Tuple::declare_key(const string &key, const KeyType &type,
                        const Default &default_value, const string &description)
{
    Record::declare_key(key, type, default_value, description);
    return *this;
}



template <class KeyType>
Tuple &Tuple::declare_key(const string &key, const KeyType &type,
                        const string &description)
{
    Record::declare_key(key, type, description);
    return *this;
}



// explicit instantiation of template methods

#define TUPLE_DECLARE_KEY(TYPE) \
template Tuple & Tuple::declare_key<TYPE>(const string &key, const TYPE &type, const Default &default_value, const string &description); \
template Tuple & Tuple::declare_key<TYPE>(const string &key, const TYPE &type, const string &description)


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
