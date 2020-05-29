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
 * @file    reader_internal_base.cc
 * @brief
 */


#include "system/asserts.hh"                           // for Assert, ASSERT
#include "system/file_path.hh"                         // for FilePath, File...
#include "system/logger.hh"                            // for operator<<

#include "input/reader_internal_base.hh"
#include "input/reader_to_storage.hh"
#include "input/input_type.hh"

namespace Input {

using namespace std;

/*******************************************************************
 * implementation of ReaderInternalBase
 */

ReaderInternalBase::ReaderInternalBase()
{}

StorageBase * ReaderInternalBase::make_storage(PathBase &p, const Type::TypeBase *type)
{
	ASSERT_PTR(type).error("Can not dispatch, NULL pointer to TypeBase.");

    // find reference node, if doesn't exist return NULL
    PathBase * ref_path = p.find_ref_node();
    if (ref_path) {
        // todo: mark passed references and check cyclic references

        // dereference and take data from there
    	StorageBase * storage = make_storage( *ref_path, type );
    	delete ref_path;
        return storage;
    }

    // dispatch types - complex types
    if (typeid(*type) == typeid(Type::Tuple)) {
        return make_sub_storage(p, static_cast<const Type::Tuple *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Record)) {
        return make_sub_storage(p, static_cast<const Type::Record *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Array)) {
        return make_sub_storage(p, static_cast<const Type::Array *>(type) );
    } else {
    	const Type::Abstract * abstract_record_type = dynamic_cast<const Type::Abstract *>(type);
    	if (abstract_record_type != NULL ) return make_sub_storage(p, abstract_record_type );
    }

    // return Null storage if there is null on the current location
	if (p.is_null_type()) {
		return new StorageNull();
	}

    // dispatch types - scalar types
    if (typeid(*type) == typeid(Type::Integer)) {
        return make_sub_storage(p, static_cast<const Type::Integer *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Double)) {
        return make_sub_storage(p, static_cast<const Type::Double *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Bool)) {
        return make_sub_storage(p, static_cast<const Type::Bool *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Selection)) {
        return make_sub_storage(p, static_cast<const Type::Selection *>(type) );
    } else {
        const Type::String * string_type = dynamic_cast<const Type::String *>(type);
        if (string_type != NULL ) return make_sub_storage(p, string_type );

        // default -> error
        THROW( Type::ExcUnknownDescendant() << Type::EI_TypeName(typeid(type).name()) );
    }

    return new StorageNull();
}

StorageBase * ReaderInternalBase::make_sub_storage(PathBase &p, const Type::Record *record)
{
	// control test, check correct tag (or TYPE key) if Record is derived from Abstract
	string record_name_from_tag = p.get_record_tag();
	if (record_name_from_tag == "include") {
		return make_include_storage(p, record);
	} else if (record_name_from_tag == "include_csv") {
		THROW( ExcForbiddenTag() << EI_Tag("include_csv")
			<< EI_Specification("can be used only with arrays.") << EI_Address(p.as_string()) );
	} else {
		if ( record_name_from_tag != "" ) {
			ASSERT(record_name_from_tag == record->type_name())(record_name_from_tag)(record->type_name()).error("Inconsistent tag of record.");
		}
		std::set<string> keys_to_process;
		bool effectively_null = p.is_effectively_null();
		if ( p.get_record_key_set(keys_to_process) || effectively_null ) {
	        std::set<string>::iterator set_it;

	        /*Type::Record::KeyIter key_it;
	        if ( record->has_key_iterator("TYPE", key_it) && record->auto_conversion_key_iter() != record->end() ) {
	            PathBase *type_path = p->clone();
	            if ( type_path.down( "TYPE" ) ) {
	                try {
	                	ASSERT( type_path.get_string_value() == record->type_name() )(type_path.get_string_value())(record->type_name())
	                		.error("Invalid value of TYPE key of record");
	                    make_storage(type_path, key_it->type_.get() )->get_int();
	                } catch(Type::Selection::ExcSelectionKeyNotFound &e) {
	                	return record_automatic_conversion(p, record);
	                }
	            }
	            else {  // automatic conversion
	            	return record_automatic_conversion(p, record);
	            }
	        }*/

	        StorageArray *storage_array = new StorageArray(record->size());
	        // check individual keys
	        for( Type::Record::KeyIter it= record->begin(); it != record->end(); ++it) {
	        	// remove processed key from keys_to_process
	        	set_it = keys_to_process.find(it->key_);
	        	if (set_it != keys_to_process.end()) {
	        		keys_to_process.erase(set_it);
	        	}

	            if ( !effectively_null && p.down(it->key_, it->key_index) ) {
	                // key on input => check & use it
	                // check for obsolete key

	                auto obsolete_it = it->attributes.find( Type::Attribute::obsolete() );
	                if ( obsolete_it != it->attributes.end()) {
	                    WarningOut() << "Usage of the obsolete key: '" << it->key_ << "'\n" << obsolete_it -> second;
	                }

	                StorageBase *storage = make_storage(p, it->type_.get());
	                if ( (typeid(*storage) == typeid(StorageNull)) && it->default_.has_value_at_declaration() ) {
	                	delete storage;
	                	storage = make_storage_from_default( it->default_.value(), it->type_ );
	                }
	                storage_array->new_item( it->key_index, storage );
	                p.up();
	            } else {
	                // key not on input
	                if (it->default_.is_obligatory() ) {
	                	this->generate_input_error(p, record, "Missing obligatory key '"+ it->key_ +"'.", false);
	                } else if (it->default_.has_value_at_declaration() ) {
	                   storage_array->new_item(it->key_index,
	                           make_storage_from_default( it->default_.value(), it->type_ ) );
	                } else { // defalut - optional or default at read time
	                    // set null
	                    storage_array->new_item(it->key_index, new StorageNull() );
	                }
	            }
	        }

	        for( set_it = keys_to_process.begin(); set_it != keys_to_process.end(); ++set_it) {
	        	WarningOut() << "Unprocessed key '" << (*set_it) << "' in " << record->class_name()
	        			<< " '" << p.as_string() << "'." << std::endl;
	        }

	        return storage_array;

	    } else { // automatic conversion
	    	return record_automatic_conversion(p, record);
	    }
	    // possibly construction of reduced record
	}

	return NULL;
}

StorageBase * ReaderInternalBase::make_sub_storage(PathBase &p, const Type::Tuple *tuple)
{
	int arr_size;
	if ( (arr_size = p.get_array_size()) != -1 ) {

		StorageArray *storage_array = new StorageArray(tuple->size());
        // check individual keys
		for ( Type::Record::KeyIter it= tuple->begin(); it != tuple->end(); ++it) {
        	if ( p.down(it->key_index) ) {
                // key on input => check & use it
                StorageBase *storage = make_storage(p, it->type_.get());
                if ( (typeid(*storage) == typeid(StorageNull)) && it->default_.has_value_at_declaration() ) {
                	delete storage;
                	storage = make_storage_from_default( it->default_.value(), it->type_ );
                }
                storage_array->new_item( it->key_index, storage );
                p.up();
        	} else {
                // key not on input
                if (it->default_.is_obligatory() ) {
                	stringstream ss;
                	ss << "Too small size of '" << p.get_node_type(ValueTypes::array_type) << "' defining Tuple with "
                			<< tuple->obligatory_keys_count() << " obligatory keys.";
                	this->generate_input_error(p, tuple, ss.str(), false);
                } else if (it->default_.has_value_at_declaration() ) {
                   storage_array->new_item(it->key_index,
                           make_storage_from_default( it->default_.value(), it->type_ ) );
                } else { // default - optional or default at read time
                    // set null
                    storage_array->new_item(it->key_index, new StorageNull() );
                }
        	}
        }

		if ( arr_size > (int)tuple->size() ) {
			WarningOut().fmt("Unprocessed keys in tuple '{}', tuple has {} keys but the input is specified by {} values.\n",
                    p.as_string().c_str(), tuple->size(), arr_size );
		}

        return storage_array;

	} else {
		return make_sub_storage(p, static_cast<const Type::Record *>(tuple) );
	}
}

StorageBase * ReaderInternalBase::make_sub_storage(PathBase &p, const Type::Abstract *abstr_rec)
{
	string record_tag = p.get_record_tag();
	if (record_tag == "") {
		if ( ! abstr_rec->get_selection_default().has_value_at_declaration() ) {
			this->generate_input_error(p, abstr_rec, "Can not determine type of the Abstract.", true);
		} else { // auto conversion
			return abstract_automatic_conversion(p, abstr_rec);
		}
	} else if (record_tag.substr(0,8) == "include:") { // include of abstract is predetermined with tag '!include:record_name'
		string record_name = record_tag.substr(8);
		try {
			return make_include_storage(p, &( abstr_rec->get_descendant(record_name) ) );
		} catch (Type::Selection::ExcSelectionKeyNotFound &exc) {
			this->generate_input_error(p, abstr_rec, "Wrong value '" + record_tag + "' of the Selection.", false);
		}
	} else if ((record_tag == "include") || (record_tag == "include_csv")) { // simple '!include' tag is forbidden
		THROW( ExcForbiddenTag() << EI_Tag(record_tag)
			<< EI_Specification("can't be used with abstract type.") << EI_Address(p.as_string()) );
	} else {
		try {
			return make_sub_storage(p, &( abstr_rec->get_descendant(record_tag) ) );
		} catch (Type::Selection::ExcSelectionKeyNotFound &exc) {
			this->generate_input_error(p, abstr_rec, "Wrong value '" + record_tag + "' of the Selection.", false);
		}
	}
	return NULL;
}

StorageBase * ReaderInternalBase::record_automatic_conversion(PathBase &p, const Type::Record *record)
{
	Type::Record::KeyIter auto_key_it = record->auto_conversion_key_iter();
	if ( auto_key_it != record->end() ) {
	    try {
			StorageArray *storage_array = new StorageArray(record->size());
			for( Type::Record::KeyIter it= record->begin(); it != record->end(); ++it) {
				if ( it == auto_key_it ) {
					// one key is initialized by input
					storage_array->new_item(it->key_index, make_storage(p, it->type_.get()) );
				} else if (it->default_.has_value_at_declaration() ) {
					// other key from default values
					storage_array->new_item(it->key_index,
							make_storage_from_default( it->default_.value(), it->type_ ) );
				 } else { // defalut - optional or default at read time
					 ASSERT(! it->default_.is_obligatory())(it->key_).error("Obligatory key in auto-convertible Record.");
					 // set null
					 storage_array->new_item(it->key_index, new StorageNull() );
				 }
			}

			return storage_array;
	    } catch (ExcInputError &e ) {
	        THROW( ExcAutomaticConversionError() << EI_RecordName(record->type_name())
	        	<< EI_InputErrorMessage(e.what()) );
	    }

	} else {
		this->generate_input_error(p, record, "The value should be '" + p.get_node_type(ValueTypes::obj_type) + "', but we found: ", true);
	}

	return NULL;
}

StorageBase * ReaderInternalBase::abstract_automatic_conversion(PathBase &p, const Type::Abstract *abstr_rec)
{
    // perform automatic conversion
    const Type::Record *default_child = abstr_rec->get_default_descendant();
    if (! default_child)
    	this->generate_input_error(p, abstr_rec, "Auto conversion of Abstract not allowed.\n", false);
    return make_sub_storage(p, default_child );
}

StorageBase * ReaderInternalBase::make_array_storage(PathBase &p, const Type::Array *array, int arr_size)
{
	ASSERT(p.is_array_type()).error();

    if ( array->match_size( arr_size ) ) {
      // copy the array and check type of values
      StorageArray *storage_array = new StorageArray(arr_size);
      for( int idx=0; idx < arr_size; idx++)  {
          p.down(idx);
          const Type::TypeBase &sub_type = array->get_sub_type();
          storage_array->new_item(idx, make_storage(p, &sub_type) );
          p.up();
      }
      return storage_array;

    } else {
    	stringstream ss;
    	ss << "Do not fit the size " << arr_size << " of the Array.";
    	this->generate_input_error(p, array, ss.str(), false);
    }

	return NULL; // suppress warning for non-void function
}

StorageBase * ReaderInternalBase::make_storage_from_default(const string &dflt_str, std::shared_ptr<Type::TypeBase> type) {
    try {
    	// default strings must be valid JSON
    	Type::Default dflt(dflt_str);
    	return dflt.get_storage(type);

    } catch (Input::Type::ExcWrongDefault & e) {
        // message to distinguish exceptions thrown during Default value check at declaration
    	e << Type::EI_Desc("Wrong default value while reading an input stream:\n");
        e << Type::EI_KeyName("UNKNOWN KEY");
        throw;
    } catch (Input::Type::ExcWrongDefaultJSON & e) {
        e << Type::EI_KeyName("UNKNOWN KEY");
        throw;
    }

    return NULL;
}

StorageBase * ReaderInternalBase::make_include_storage(PathBase &p, const Type::TypeBase *type)
{
    std::string included_path;
    if ( p.is_record_type() ) {
        // include is set as record with tag and file key
        if ( p.down("file") ) {
        	included_path = get_included_file(p);
            p.up();
        } else {
        	this->generate_input_error(p, type, "Missing key 'file' defines including input file.", false);
        }
    } else {
    	// include is set only with name of file (similarly as auto conversion)
    	// this case may occur only for YAML input
    	included_path = get_included_file(p);
    }

    FilePath fpath(included_path, FilePath::FileType::input_file);
	ReaderToStorage include_reader(fpath, *(const_cast<Type::TypeBase *>(type)) );
	return include_reader.get_storage();
}

bool ReaderInternalBase::read_bool_value(PathBase &p, const Type::TypeBase *type)
{
	bool value;
	try {
		value = p.get_bool_value();
	}
	catch (ExcInputError & e) {
		complete_input_error(e, p, ValueTypes::bool_type);
		e << EI_InputType(type->desc());
		throw;
	}

	return value;
}

std::int64_t ReaderInternalBase::read_int_value(PathBase &p, const Type::TypeBase *type)
{
	std::int64_t value;
	try {
		value = p.get_int_value();
	}
	catch (ExcInputError & e) {
		complete_input_error(e, p, ValueTypes::int_type);
		e << EI_InputType(type->desc());
		throw;
	}

	return value;
}

double ReaderInternalBase::read_double_value(PathBase &p, const Type::TypeBase *type)
{
    double value;
	try {
		value = p.get_double_value();
	}
	catch (ExcInputError & e) {
		complete_input_error(e, p, ValueTypes::real_type);
		e << EI_InputType(type->desc());
		throw;
	}

	return value;
}

std::string ReaderInternalBase::read_string_value(PathBase &p, const Type::TypeBase *type)
{
	string value;
	try {
		value = p.get_string_value();
	} catch (ExcInputError & e) {
		complete_input_error(e, p, ValueTypes::str_type);
		e << EI_InputType(type->desc());
		throw;
	}

    return value;
}

std::string ReaderInternalBase::get_included_file(PathBase &p)
{
	try {
		return p.get_string_value();
	}
	catch (ExcInputError & e) {
		complete_input_error(e, p, ValueTypes::str_type);
		e << EI_InputType("path to included file");
		throw;
	}
}

void ReaderInternalBase::generate_input_error(PathBase &p, const Type::TypeBase *type, std::string spec, bool add_type)
{
	if (add_type)
		THROW( ExcInputError() << EI_Specification(spec) << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) )
            << EI_ErrorAddress(p.as_string()) << EI_InputType(type->desc()) );
	else
		THROW( ExcInputError() << EI_Specification(spec) << EI_JSON_Type( "" )
            << EI_ErrorAddress(p.as_string()) << EI_InputType(type->desc()) );
}

void ReaderInternalBase::complete_input_error(ExcInputError & e, PathBase &p, ValueTypes value_type)
{
	e << EI_Specification("The value should be '" + p.get_node_type(value_type) + "', but we found: ");
    e << EI_ErrorAddress(p.as_string());
    e << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) );
}


} // namespace Input
