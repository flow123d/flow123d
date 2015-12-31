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
 * @file    type_record.cc
 * @brief   
 */

#include "input_type.hh"
#include "type_repository.hh"

#include "system/system.hh"
#include "input/reader_to_storage.hh"

#include <boost/typeof/typeof.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/functional/hash.hpp>

namespace Input {
namespace Type {

using namespace std;

/*******************************************************************
 * implementation of Default
 */

Default::Default()
: value_("OPTIONAL"), type_(no_default_optional_type), storage_(NULL)
{}

Default::Default(const std::string & value)
: value_(value), type_(default_at_declaration), storage_(NULL)
{
    boost::algorithm::trim(value_);
}



Default::Default(enum DefaultType type, const std::string & value)
: value_(value), type_(type), storage_(NULL)
{}

TypeBase::TypeHash Default::content_hash() const
{   TypeBase::TypeHash seed = 0;
    boost::hash_combine(seed, "Default");
    boost::hash_combine(seed, type_);
    boost::hash_combine(seed, value_);
    return seed;
}


bool Default::check_validity(boost::shared_ptr<TypeBase> type) const
{
	if ( storage_ ) return true;
	if ( !has_value_at_declaration() ) return false;

	try {
		istringstream is("[\n" + value_ + "\n]");
		Input::ReaderToStorage reader;
		reader.read_stream(is, Array(type), FileFormat::format_JSON);
		storage_ = reader.get_storage()->get_item(0);
		return true;
	} catch ( Input::ReaderToStorage::ExcNotJSONFormat &e ) {
		THROW( ExcWrongDefault() << EI_DefaultStr( value_ ) << EI_TypeName(type->type_name()));
	} catch ( Input::ReaderToStorage::ExcInputError &e ) {
		THROW( ExcWrongDefault() << EI_DefaultStr( value_ ) << EI_TypeName(type->type_name()));
	}
}


Input::StorageBase *Default::get_storage(boost::shared_ptr<TypeBase> type) const
{
	if ( !storage_ ) this->check_validity(type);
	return storage_;
}



/**********************************************************************************
 * implementation of Type::Record
 */


Record::Record()
: data_( boost::make_shared<RecordData> ("EmptyRecord","") )
{
	close();
    finish();
}



Record::Record(const Record & other)
: TypeBase( other ), data_(other.data_)
{}



Record::Record(const string & type_name_in, const string & description)
: data_( boost::make_shared<RecordData>(type_name_in, description) )

{}


TypeBase::TypeHash Record::content_hash() const
{
	TypeHash seed=0;
    boost::hash_combine(seed, "Record");
    boost::hash_combine(seed, type_name());
    boost::hash_combine(seed, data_->description_);
    boost::hash_combine(seed, data_->auto_conversion_key);
    for( Key &key : data_->keys) {
    	if (key.key_ != "TYPE") {
    		boost::hash_combine(seed, key.key_);
    		boost::hash_combine(seed, key.description_);
    		boost::hash_combine(seed, key.type_->content_hash() );
    		boost::hash_combine(seed, key.default_.content_hash() );
    	}
    }
    attribute_content_hash(seed);
    return seed;
}



Record &Record::allow_auto_conversion(const string &from_key) {
    ASSERT(data_->auto_conversion_key_idx == -1, "Can not use key %s for auto conversion, the key is already set.", from_key.c_str());
    data_->auto_conversion_key_idx = 0;
    data_->auto_conversion_key=from_key;

    return *this;
}



void Record::make_copy_keys(Record &origin) {

    ASSERT(origin.is_closed(), "Origin record is not closed!\n");

	std::vector<Key>::iterator it = data_->keys.begin();
	if (data_->keys.size() && it->key_ == "TYPE") it++; // skip TYPE key if exists

	int n_inserted = 0;
	for(KeyIter pit=origin.data_->keys.begin(); pit != origin.data_->keys.end(); ++pit) {
		Key tmp_key=*pit;    // make temporary copy of the key
		KeyHash key_h = key_hash(tmp_key.key_);

		tmp_key.derived = true;

		// we have to copy TYPE also since there should be place in storage for it
		// however we change its Default to name of actual Record
		if (tmp_key.key_=="TYPE")
			tmp_key.default_=Default( "\""+type_name()+"\"" );

		// check for duplicate keys, override keys of the parent record by the child record
		RecordData::key_to_index_const_iter kit = data_->key_to_index.find(key_h);
		if (kit != data_->key_to_index.end()) {
			// in actual record exists a key with same name as in parent record
			// use values form the child record
			Key *k = &(data_->keys[kit->second+n_inserted]); // indices in key_to_index are not yet updated

			tmp_key.key_ = k->key_;
			tmp_key.description_ = k->description_;
			tmp_key.type_ = k->type_;
			tmp_key.default_ = k->default_;
			tmp_key.derived = false;
			k->key_ = ""; // mark original key for deletion
		}

		data_->key_to_index[key_h] = tmp_key.key_index;

		it = data_->keys.insert(it, tmp_key)+1;
		n_inserted++;
	}
	// delete duplicate keys and update key indices
	for (unsigned int i=0; i<data_->keys.size(); i++) {
		if (data_->keys[i].key_.compare("") == 0) {
			data_->keys.erase( data_->keys.begin()+i);
			i--;
		} else {
			data_->keys[i].key_index = i;
			data_->key_to_index[key_hash( data_->keys[i].key_)] = i;
		}
	}
}



Record &Record::derive_from(Abstract &parent) {
	ASSERT( parent.is_closed(), "Parent Abstract '%s' must be closed!\n", parent.type_name().c_str());
	ASSERT( data_->keys.size() == 0 || (data_->keys.size() == 1 && data_->keys[0].key_ == "TYPE"),
			"Derived record '%s' can have defined only TYPE key!\n", this->type_name().c_str() );

	// add Abstract to vector of parents
	data_->parent_vec_.push_back( boost::make_shared<Abstract>(parent) );

	if (data_->keys.size() == 0) {
    	data_->declare_key("TYPE", boost::make_shared<String>(), Default( "\""+type_name()+"\"" ), "Sub-record Selection.");
    }

	return *this;
}



Record &Record::copy_keys(const Record &other) {
	ASSERT( other.is_closed(), "Record '%s' must be closed!\n", other.type_name().c_str());

   	Record tmp(other);
   	make_copy_keys(tmp);

   	return *this;
}


bool Record::is_finished() const {
    return data_->finished;
}


bool Record::is_closed() const {
	return data_->closed_;
}




bool Record::finish(bool is_generic)
{

	if (data_->finished) return true;

	ASSERT(data_->closed_, "Finished Record '%s' must be closed!", this->type_name().c_str());

    data_->finished = true;
    for (vector<Key>::iterator it=data_->keys.begin(); it!=data_->keys.end(); it++)
    {
    	if (it->key_ != "TYPE") {
			if (typeid( *(it->type_.get()) ) == typeid(Instance)) it->type_ = it->type_->make_instance().first;
			if (!is_generic && it->type_->is_root_of_generic_subtree()) THROW( ExcGenericWithoutInstance() << EI_Object(it->type_->type_name()) );
           	data_->finished = data_->finished && it->type_->finish(is_generic);
        }
    }

    // Check default values of autoconvertible records
    if (data_->auto_conversion_key_idx != -1 ) {
        data_->auto_conversion_key_idx=key_index(data_->auto_conversion_key);

        // check that all other obligatory keys have default values
        for(KeyIter it=data_->keys.begin(); it != data_->keys.end(); ++it) {
            if (it->default_.is_obligatory() && (int)(it->key_index) != data_->auto_conversion_key_idx)
                xprintf(PrgErr, "Finishing Record auto convertible from the key '%s', but other obligatory key: '%s' has no default value.\n",
                        data_->auto_conversion_key_iter()->key_.c_str(), it->key_.c_str());
        }
    }

    return (data_->finished);
}



const Record &Record::close() const {
    data_->closed_=true;
    const Record & rec = *( Input::TypeRepository<Record>::get_instance().add_type( *this ) );
    for (auto &parent : data_->parent_vec_) {
    	parent->add_child(rec);
    }
    data_->parent_vec_.clear();

    return rec;
}



string Record::type_name() const {
    return data_->type_name_;
}



bool Record::operator==(const TypeBase &other) const
{ return  typeid(*this) == typeid(other) &&
                 (type_name() == static_cast<const Record *>(&other)->type_name() );
}


Record::KeyIter Record::auto_conversion_key_iter() const {
    finished_check();
    return data_->auto_conversion_key_iter();
}


Record &Record::declare_type_key() {
	ASSERT(data_->keys.size() == 0, "Declaration of TYPE key must be carried as the first.");
	data_->declare_key("TYPE", boost::make_shared<String>(), Default::obligatory(),
			"Sub-record selection.");
	return *this;
}

Record &Record::has_obligatory_type_key() {
	ASSERT( ! data_->parent_vec_.size(), "Record with obligatory TYPE key can't be derived.\n");
	declare_type_key();
	return *this;
}


TypeBase::MakeInstanceReturnType Record::make_instance(std::vector<ParameterPair> vec) const {
	Record rec = this->deep_copy();
	ParameterMap parameter_map;
	// Replace keys of type Parameter
	for (std::vector<Key>::iterator key_it=rec.data_->keys.begin(); key_it!=rec.data_->keys.end(); key_it++) {
		if ( key_it->key_ != "TYPE" ) { // TYPE key isn't substituted
			MakeInstanceReturnType inst = key_it->type_->make_instance(vec);
			key_it->type_ = inst.first;
			ParameterMap other_map = inst.second;
			parameter_map.insert(other_map.begin(), other_map.end());
		}
	}
	// Set attributes
	rec.set_parameters_attribute(parameter_map);
	rec.add_attribute("generic_type", this->hash_str());

	return std::make_pair( boost::make_shared<Record>(rec.close()), parameter_map );
}


Record Record::deep_copy() const {
	Record rec = Record();
	rec.data_ =  boost::make_shared<Record::RecordData>(*this->data_);
	rec.data_->closed_ = false;
	rec.data_->finished = false;
	rec.attributes_ = boost::make_shared<attribute_map>(*attributes_);
	return rec;
}


const Record &Record::add_parent(Abstract &parent) const {
	ASSERT( parent.is_closed(), "Parent Abstract '%s' must be closed!\n", parent.type_name().c_str());

	// check if parent exists in parent_vec_ vector
	TypeHash hash = parent.content_hash();
	for (auto &parent : data_->parent_vec_) {
		if ( parent->content_hash() == hash ) {
			return *this;
		}
	}

	data_->parent_vec_.push_back( boost::make_shared<Abstract>(parent) );

	// finish inheritance
	ASSERT( data_->keys.size() > 0 && data_->keys[0].key_ == "TYPE",
				"Derived record '%s' must have defined TYPE key!\n", this->type_name().c_str() );
	data_->keys[0].default_ = Default( "\""+type_name()+"\"" );

	return *this;
}


Record &Record::root_of_generic_subtree() {
	root_of_generic_subtree_ = true;
	return *this;
}



/**********************************************************************************
 * implementation of Type::Record::RecordData
 */

Record::RecordData::RecordData(const string & type_name_in, const string & description)
:description_(description),
 type_name_(type_name_in),
 finished(false),
 closed_(false),
 derived_(false),
 auto_conversion_key_idx(-1)    // auto conversion turned off
{

}


Record::KeyIter Record::RecordData::auto_conversion_key_iter() const {
	if (auto_conversion_key_idx >= 0) return keys.begin() + auto_conversion_key_idx;
  	else return keys.end();
}



void Record::RecordData::declare_key(const string &key,
                         boost::shared_ptr<TypeBase> type,
                         const Default &default_value, const string &description)
{
    ASSERT(!closed_, "Can not add key '%s' into closed record '%s'.\n", key.c_str(), type_name_.c_str());
    // validity test of default value
    try {
    	default_value.check_validity(type);
    } catch (ExcWrongDefault & e) {
        e << EI_KeyName(key);
        throw;
    }
    if (finished) xprintf(PrgErr, "Declaration of key: %s in finished Record type: %s\n", key.c_str(), type_name_.c_str());

    if (key!="TYPE" && ! TypeBase::is_valid_identifier(key))
        xprintf(PrgErr, "Invalid key identifier %s in declaration of Record type: %s\n", key.c_str(), type_name_.c_str());

    KeyHash key_h = key_hash(key);
    key_to_index_const_iter it = key_to_index.find(key_h);
    if ( it == key_to_index.end() ) {
       key_to_index.insert( std::make_pair(key_h, keys.size()) );
       Key tmp_key = { (unsigned int)keys.size(), key, description, type, default_value, false};
       keys.push_back(tmp_key);
    } else {
       if (keys[it->second].derived) {
        Key tmp_key = { it->second, key, description, type, default_value, false};
        keys[ it->second ] = tmp_key;
       } else {
           xprintf(Err,"Re-declaration of the key: %s in Record type: %s\n", key.c_str(), type_name_.c_str() );
       }
    }

}

Record &Record::declare_key(const string &key, boost::shared_ptr<TypeBase> type,
                        const Default &default_value, const string &description)
{
    data_->declare_key(key, type, default_value, description);
    return *this;
}


template <class KeyType>
Record &Record::declare_key(const string &key, const KeyType &type,
                        const Default &default_value, const string &description)
// this accept only lvalues - we assume that these are not local variables
{
    // ASSERT MESSAGE: The type of declared keys has to be a class derived from TypeBase.
    BOOST_STATIC_ASSERT( (boost::is_base_of<TypeBase, KeyType>::value) );
	boost::shared_ptr<TypeBase> type_copy = boost::make_shared<KeyType>(type);
	return declare_key(key, type_copy, default_value, description);
}



template <class KeyType>
Record &Record::declare_key(const string &key, const KeyType &type,
                        const string &description)
{
    return declare_key(key,type, Default::optional(), description);
}



// explicit instantiation of template methods

#define RECORD_DECLARE_KEY(TYPE) \
template Record & Record::declare_key<TYPE>(const string &key, const TYPE &type, const Default &default_value, const string &description); \
template Record & Record::declare_key<TYPE>(const string &key, const TYPE &type, const string &description)


RECORD_DECLARE_KEY(String);
RECORD_DECLARE_KEY(Integer);
RECORD_DECLARE_KEY(Double);
RECORD_DECLARE_KEY(Bool);
RECORD_DECLARE_KEY(FileName);
RECORD_DECLARE_KEY(Selection);
RECORD_DECLARE_KEY(Array);
RECORD_DECLARE_KEY(Record);
RECORD_DECLARE_KEY(Abstract);
RECORD_DECLARE_KEY(AdHocAbstract);
RECORD_DECLARE_KEY(Parameter);
RECORD_DECLARE_KEY(Instance);




} // closing namespace Type
} // closing namespace Input


