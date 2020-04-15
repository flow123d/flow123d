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
#include "input/reader_to_storage.hh"
#include "input/reader_internal_base.hh"
#include "attribute_lib.hh"

#include <boost/typeof/typeof.hpp>
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


bool Default::check_validity(std::shared_ptr<TypeBase> type) const
{
	if ( storage_ ) return true;
	if ( !has_value_at_declaration() ) return false;
    
    type->finish();

	try {
		istringstream is("[\n" + value_ + "\n]");
		Input::ReaderToStorage reader;
		reader.read_stream(is, Array(type), FileFormat::format_JSON);
		storage_ = reader.get_storage()->get_item(0);
		return true;
	} catch ( Input::ReaderInternalBase::ExcNotJSONFormat &e ) {
		THROW( ExcWrongDefaultJSON() << EI_DefaultStr( value_ ) << EI_TypeName(type->type_name())
				<< make_nested_ei(e) );
	} catch ( Input::ReaderInternalBase::ExcInputError &e ) {
		THROW( ExcWrongDefault() << EI_DefaultStr( value_ ) << EI_TypeName(type->type_name())
				<< make_nested_ei(e) );
	}
}


Input::StorageBase *Default::get_storage(std::shared_ptr<TypeBase> type) const
{
	if ( !storage_ ) this->check_validity(type);
	return storage_;
}



/**********************************************************************************
 * implementation of Type::Record
 */


Record::Record()
: data_( std::make_shared<RecordData> ("EmptyRecord","") )
{
	declare_key("empty_key", Integer(), "Only for correct 'close()'.");
	close();
    finish();
}



Record::Record(const Record & other)
: TypeBase( other ), data_(other.data_)
{}



Record::Record(const string & type_name_in, const string & description)
: data_( std::make_shared<RecordData>(type_name_in, description) )
{
	data_->declare_key("TYPE", std::make_shared<String>(), Default( "\""+type_name()+"\"" ),
	        "Sub-record Selection.", TypeBase::attribute_map());
}


TypeBase::TypeHash Record::content_hash() const
{
	TypeHash seed=0;
    boost::hash_combine(seed, "Record");
    data_->content_hash(seed);
    return seed;
}



Record &Record::allow_auto_conversion(const string &from_key) {
	ASSERT(data_->auto_conversion_key_idx == -1)(from_key).error("auto conversion key is already set");
    data_->auto_conversion_key_idx = 0;
    data_->auto_conversion_key=from_key;

    return *this;
}



void Record::make_copy_keys(Record &origin) {

	ASSERT( origin.is_closed() ).error();

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
	ASSERT( parent.is_closed() )(parent.type_name()).error();
	ASSERT( data_->keys.size() > 0 && data_->keys[0].key_ == "TYPE" )(this->type_name())
			.error("Derived record must have defined TYPE key!");

	// check if parent exists in parent_vec_ vector
	TypeHash hash = parent.content_hash();
	for (auto &parent : data_->parent_vec_) {
		if ( parent->content_hash() == hash ) {
			return *this;
		}
	}

	// add Abstract to vector of parents
	data_->parent_vec_.push_back( std::make_shared<Abstract>(parent) );

	return *this;
}


Record &Record::copy_keys(const Record &other) {
	ASSERT( other.is_closed() )(other.type_name()).error();

   	Record tmp(other);
   	make_copy_keys(tmp);

   	return *this;
}


FinishStatus Record::finish_status() const {
    return data_->finish_status_;
}


bool Record::is_finished() const {
    return (data_->finish_status_ != FinishStatus::none_) && (data_->finish_status_ != FinishStatus::in_perform_);
}


bool Record::is_closed() const {
	return data_->closed_;
}




FinishStatus Record::finish(FinishStatus finish_type)
{
	ASSERT(finish_type != FinishStatus::none_).error();
	ASSERT(finish_type != FinishStatus::in_perform_).error();
	ASSERT(data_->finish_status_ != FinishStatus::in_perform_)(this->type_name())(this->type_name()).error("Recursion in the IST element of type Record.");

	if (this->is_finished()) return data_->finish_status_;

	ASSERT(data_->closed_)(this->type_name()).error();

    data_->finish_status_ = FinishStatus::in_perform_;
    for (vector<Key>::iterator it=data_->keys.begin(); it!=data_->keys.end(); it++)
    {

      	if (it->key_ != "TYPE") {
			if (typeid( *(it->type_.get()) ) == typeid(Instance)) {
				it->type_->finish(FinishStatus::generic_); // finish Instance object
				it->type_ = it->type_->make_instance().first;
			}
			if ((finish_type != FinishStatus::generic_) && it->type_->is_root_of_generic_subtree())
			    THROW( ExcGenericWithoutInstance()
			            << EI_Object(it->type_->type_name())
			            << EI_TypeName(this->type_name()));
			it->type_->finish(finish_type);
			ASSERT(it->type_->is_finished()).error();
			if (finish_type == FinishStatus::delete_) it->type_.reset();
        }

        if (finish_type == FinishStatus::regular_) {
            try {
                it->default_.check_validity(it->type_);
            } catch (ExcWrongDefaultJSON & e) {
                e << EI_KeyName(it->key_);
                throw;
            } catch (ExcWrongDefault & e) {
                e << EI_KeyName(it->key_);
                throw;
            }
        }

    }

    // Check default values of autoconvertible records
    if (data_->auto_conversion_key_idx != -1 ) {
        data_->auto_conversion_key_idx=key_index(data_->auto_conversion_key);

        // check that all other obligatory keys have default values
        for(KeyIter it=data_->keys.begin(); it != data_->keys.end(); ++it) {
        	const string &other_key = it->key_;
        	ASSERT(!it->default_.is_obligatory() || (int)(it->key_index) == data_->auto_conversion_key_idx)
        			   (data_->auto_conversion_key_iter()->key_)(other_key)
					   .error("Finishing auto convertible Record from given key, but other obligatory key has no default value.");
        }
    }

    data_->finish_status_ = finish_type;
    return (data_->finish_status_);
}



Record &Record::close() const {
	ASSERT_GT(data_->keys.size(), 0)(this->type_name()).error("Empty Record!\n");
    data_->closed_=true;
    Record & rec = *( Input::TypeRepository<Record>::get_instance().add_type( *this ) );
    for (auto &parent : data_->parent_vec_) {
    	parent->add_child(rec);
    }

    return rec;
}



string Record::type_name() const {
    return data_->type_name_;
}



string Record::class_name() const {
	return "Record";
}



bool Record::operator==(const TypeBase &other) const
{ return  typeid(*this) == typeid(other) &&
                 (type_name() == static_cast<const Record *>(&other)->type_name() );
}


Record::KeyIter Record::auto_conversion_key_iter() const {
    finished_check();
    return data_->auto_conversion_key_iter();
}


/*Record &Record::declare_type_key() {
	ASSERT(data_->keys.size() == 0).error("Declaration of TYPE key must be carried as the first.");
	data_->declare_key("TYPE", std::make_shared<String>(), Default::obligatory(),
			"Sub-record selection.", TypeBase::attribute_map());
	return *this;
}*/

/*Record &Record::has_obligatory_type_key() {
	ASSERT(! data_->parent_vec_.size()).error("Record with obligatory TYPE key can't be derived");
	declare_type_key();
	return *this;
}*/

Record &Record::add_attribute(std::string key, TypeBase::json_string value) {
    this->add_attribute_(key, value);
    return *this;
}


TypeBase::MakeInstanceReturnType Record::make_instance(std::vector<ParameterPair> vec) {
    Record instance_rec = this->deep_copy();
    auto parameter_map = this->set_instance_data(instance_rec, vec);
	return std::make_pair( std::make_shared<Record>(instance_rec.close()), parameter_map );
}


TypeBase::ParameterMap Record::set_instance_data(Record &instance_rec,  std::vector<ParameterPair> vec) {
    ParameterMap p_map;

    // Replace keys of type Parameter
	for (std::vector<Key>::iterator key_it=instance_rec.data_->keys.begin(); key_it!=instance_rec.data_->keys.end(); key_it++) {
		if ( key_it->key_ != "TYPE" ) { // TYPE key isn't substituted
			MakeInstanceReturnType inst = key_it->type_->make_instance(vec);
			key_it->type_ = inst.first;
			ParameterMap other_map = inst.second;
			p_map.insert(other_map.begin(), other_map.end());
		}
	}

    // Set attributes
    instance_rec.parameter_map_ = p_map;
    this->set_generic_attributes(p_map);

    // Set reference to generic type of the instance.
    bool param_substituted = (p_map.size() > 0);
    auto generic_hash = this->content_hash();
	auto instance_hash = instance_rec.content_hash();
	ASSERT(param_substituted == (generic_hash != instance_hash));
	if (param_substituted) {
	    // Avoid self reference for non-parametrized types.
	    instance_rec.generic_type_hash_ = generic_hash;
	}
	return p_map;
}



Record Record::deep_copy() const {
	Record rec = Record();
	rec.data_ =  std::make_shared<Record::RecordData>(*this->data_);
	rec.data_->closed_ = false;
	rec.data_->finish_status_ = FinishStatus::none_;
	rec.copy_attributes(*attributes_);
	rec.generic_type_hash_ = this->generic_type_hash_;
	rec.parameter_map_ = this->parameter_map_;
	return rec;
}


/*const Record &Record::add_parent(Abstract &parent) const {
	ASSERT( parent.is_closed() )(parent.type_name()).error();

	// check if parent exists in parent_vec_ vector
	TypeHash hash = parent.content_hash();
	for (auto &parent : data_->parent_vec_) {
		if ( parent->content_hash() == hash ) {
			return *this;
		}
	}

	data_->parent_vec_.push_back( std::make_shared<Abstract>(parent) );

	// finish inheritance
	ASSERT( data_->keys.size() > 0 && data_->keys[0].key_ == "TYPE" )(this->type_name())
			.error("Derived record must have defined TYPE key!");
	data_->keys[0].default_ = Default( "\""+type_name()+"\"" );

	return *this;
}*/


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
 finish_status_(FinishStatus::none_),
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
                         std::shared_ptr<TypeBase> type,
                         const Default &default_value,
                         const string &description,
                         TypeBase::attribute_map key_attributes)
{
	ASSERT(!closed_)(key)(this->type_name_).error();
    // validity test of default value

    ASSERT( finish_status_ == FinishStatus::none_ )(key)(type_name_).error("Declaration of key in finished Record");
    ASSERT( key=="TYPE" || TypeBase::is_valid_identifier(key) )(key)(type_name_).error("Invalid key identifier in declaration of Record");

    KeyHash key_h = key_hash(key);
    key_to_index_const_iter it = key_to_index.find(key_h);
    if ( it == key_to_index.end() ) {
       key_to_index.insert( std::make_pair(key_h, keys.size()) );
       Key tmp_key = { (unsigned int)keys.size(), key, description, type, default_value, false, key_attributes };
       keys.push_back(tmp_key);
    } else {
    	ASSERT( keys[it->second].derived )(key)(type_name_).error("Re-declaration of the key in Record");
        Key tmp_key = { it->second, key, description, type, default_value, false, {}};
        keys[ it->second ] = tmp_key;
    }

}

void Record::RecordData::content_hash(TypeBase::TypeHash &seed) const {
    boost::hash_combine(seed, type_name_);
    boost::hash_combine(seed, description_);
    boost::hash_combine(seed, auto_conversion_key);
    for( const Key &key : keys) {
    	if (key.key_ != "TYPE") {
    		boost::hash_combine(seed, key.key_);
    		boost::hash_combine(seed, key.description_);
    		boost::hash_combine(seed, key.type_->content_hash() );
    		boost::hash_combine(seed, key.default_.content_hash() );
    	}
    }
}


Record &Record::declare_key(const string &key, std::shared_ptr<TypeBase> type,
                        const Default &default_value, const string &description,
                        TypeBase::attribute_map key_attributes)
{
    data_->declare_key(key, type, default_value, description, key_attributes);
    return *this;
}


template <class KeyType>
Record &Record::declare_key(const string &key, const KeyType &type,
                        const Default &default_value, const string &description,
                        TypeBase::attribute_map key_attributes)
// this accept only lvalues - we assume that these are not local variables
{
    // ASSERT MESSAGE: The type of declared keys has to be a class derived from TypeBase.
    BOOST_STATIC_ASSERT( (boost::is_base_of<TypeBase, KeyType>::value) );
	std::shared_ptr<TypeBase> type_copy = std::make_shared<KeyType>(type);
	return declare_key(key, type_copy, default_value, description, key_attributes);
}



template <class KeyType>
Record &Record::declare_key(const string &key, const KeyType &type,
                        const string &description,
                        TypeBase::attribute_map key_attributes)
{
    return declare_key(key,type, Default::optional(), description, key_attributes);
}



// explicit instantiation of template methods

#define RECORD_DECLARE_KEY(TYPE) \
template Record & Record::declare_key<TYPE>(const string &key, const TYPE &type, const Default &default_value, const string &description, TypeBase::attribute_map key_attributes); \
template Record & Record::declare_key<TYPE>(const string &key, const TYPE &type, const string &description, TypeBase::attribute_map key_attributes)


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
RECORD_DECLARE_KEY(Tuple);




} // closing namespace Type
} // closing namespace Input


