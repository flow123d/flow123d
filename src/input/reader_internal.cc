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
 * @file    reader_internal.cc
 * @brief
 */


#include "input/reader_internal.hh"
#include "input/reader_to_storage.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/csv_tokenizer.hh"

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
			<< EI_Specification("can be used only with arrays.") );
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
	string record_name = p.get_record_tag();
	if (record_name == "") {
		if ( ! abstr_rec->get_selection_default().has_value_at_declaration() ) {
			this->generate_input_error(p, abstr_rec, "Can not determine type of the Abstract.", true);
		} else { // auto conversion
			return abstract_automatic_conversion(p, abstr_rec);
		}
	} else if ((record_name == "include") || (record_name == "include_csv")) {
		THROW( ExcForbiddenTag() << EI_Tag(record_name)
			<< EI_Specification("can't be used with abstract type.") );
	} else {
		try {
			return make_sub_storage(p, &( abstr_rec->get_descendant(record_name) ) );
		} catch (Type::Selection::ExcSelectionKeyNotFound &exc) {
			this->generate_input_error(p, abstr_rec, "Wrong value '" + record_name + "' of the Selection.", false);
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
}

StorageBase * ReaderInternalBase::make_storage_from_default(const string &dflt_str, std::shared_ptr<Type::TypeBase> type) {
    try {
    	// default strings must be valid JSON
    	Type::Default dflt(dflt_str);
    	return dflt.get_storage(type);

    } catch (Input::Type::ExcWrongDefault & e) {
        // message to distinguish exceptions thrown during Default value check at declaration
    	e << Type::EI_Desc("Wrong default value while reading an input stream:\n");
        e << EI_KeyName("UNKNOWN KEY");
        throw;
    } catch (Input::Type::ExcWrongDefaultJSON & e) {
        e << EI_KeyName("UNKNOWN KEY");
        throw;
    }

    return NULL;
}

StorageBase * ReaderInternalBase::make_include_storage(PathBase &p, const Type::Record *record)
{
    std::string included_path;
    if ( p.is_record_type() ) {
        // include is set as record with tag and file key
        if ( p.down("file") ) {
        	included_path = get_included_file(p);
            p.up();
        } else {
        	this->generate_input_error(p, record, "Missing key 'file' defines including input file.", false);
        }
    } else {
    	// include is set only with name of file (similarly as auto conversion)
    	// this case may occur only for YAML input
    	included_path = get_included_file(p);
    }

    FilePath fpath(included_path, FilePath::FileType::input_file);
    try {
    	ReaderToStorage include_reader(fpath, *(const_cast<Type::Record *>(record)) );
        return include_reader.get_storage();
    } catch (ExcInputError &e ) {
      e << EI_File(fpath); throw;
    } catch (ExcNotJSONFormat &e) {
      e << EI_File(fpath); throw;
    }

	return NULL;
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


/*******************************************************************
 * implementation of ReaderInternal
 */

ReaderInternal::ReaderInternal()
{}

StorageBase * ReaderInternal::read_storage(PathBase &p, const Type::TypeBase *type)
{
	return this->make_storage(p, type);
}

StorageBase * ReaderInternal::make_sub_storage(PathBase &p, const Type::Array *array)
{
	int arr_size;
	if ( (arr_size = p.get_array_size()) != -1 ) {
		return this->make_array_storage(p, array, arr_size);
	} else if (p.get_record_tag() == "include_csv") {
		ReaderInternalCsvInclude reader_csv;
		return reader_csv.read_storage(p, array);
	} else {
		ReaderInternalTranspose reader_transpose;
		return reader_transpose.read_storage(p, array);
	}
}

StorageBase * ReaderInternal::make_sub_storage(PathBase &p, const Type::Selection *selection)
{
    string item_name = read_string_value(p, selection);
	try {
		int value = selection->name_to_int( item_name );
		return new StorageInt( value );
	} catch (Type::Selection::ExcSelectionKeyNotFound &exc) {
		this->generate_input_error(p, selection, "Wrong value '" + item_name + "' of the Selection.", false);
	}

    return NULL;
}

StorageBase * ReaderInternal::make_sub_storage(PathBase &p, const Type::Bool *bool_type)
{
	return new StorageBool( read_bool_value(p, bool_type) );
}

StorageBase * ReaderInternal::make_sub_storage(PathBase &p, const Type::Integer *int_type)
{
	std::int64_t value = read_int_value(p, int_type);

	if ( int_type->match(value) )
	{
		return new StorageInt( value );
	} else {
		this->generate_input_error(p, int_type, "Value out of bounds.", false);
	}
	return NULL;
}

StorageBase * ReaderInternal::make_sub_storage(PathBase &p, const Type::Double *double_type)
{
    double value = read_double_value(p, double_type);

    if (double_type->match(value)) {
        return new StorageDouble( value );
    } else {
    	this->generate_input_error(p, double_type, "Value out of bounds.", false);
    }

    return NULL;
}

StorageBase * ReaderInternal::make_sub_storage(PathBase &p, const Type::String *string_type)
{
	string value = read_string_value(p, string_type);

	if (string_type->match(value))
		return new StorageString( value );
	else
		THROW( ExcInputError() << EI_Specification("Output file can not be given by absolute path: '" + value + "'")
						<< EI_ErrorAddress(p.as_string()) << EI_JSON_Type("") << EI_InputType(string_type->desc()) );

	return NULL;
}


/*******************************************************************
 * implementation of ReaderInternalTranspose
 */

ReaderInternalTranspose::ReaderInternalTranspose()
: transpose_index_(0)
{}

StorageBase * ReaderInternalTranspose::read_storage(PathBase &p, const Type::Array *array)
{
	const Type::TypeBase &sub_type = array->get_sub_type();
	StorageBase *first_item_storage;
	try {
		first_item_storage = make_storage(p, &sub_type);
	} catch (ExcInputError &e) {
		if ( !array->match_size(1) ) {
			e << EI_Specification("The value should be '" + p.get_node_type(ValueTypes::array_type) + "', but we found: ");
		}
		e << EI_TransposeIndex(transpose_index_);
		e << EI_TransposeAddress(p.as_string());
		throw;
	}

	// automatic conversion to array with one element
	if (transpose_array_sizes_.size() == 0) {
		return make_autoconversion_array_storage(p, array, first_item_storage);
	} else {

		// check sizes of arrays stored in transpose_array_sizes_
		transpose_array_sizes_.erase( unique( transpose_array_sizes_.begin(), transpose_array_sizes_.end() ),
									  transpose_array_sizes_.end() );
		if (transpose_array_sizes_.size() == 1) {
			unsigned int sizes = transpose_array_sizes_[0]; // sizes of transposed

			// array size out of bounds
			if ( !array->match_size( sizes ) ) {
				stringstream ss;
				ss << "Result of transpose auto-conversion do not fit the size " << sizes << " of the Array.";
				this->generate_input_error(p, array, ss.str(), false);
			}

			// create storage of array
			StorageArray *storage_array = new StorageArray(sizes);
			storage_array->new_item(0, first_item_storage);
			if (sizes>1) {
				++transpose_index_;
				while (transpose_index_ < sizes) {
					try {
						storage_array->new_item(transpose_index_, make_storage(p, &sub_type));
					} catch (ExcInputError &e) {
						e << EI_TransposeIndex(transpose_index_);
						e << EI_TransposeAddress(p.as_string());
						throw;
					}
					++transpose_index_;
				}
			}

			return storage_array;
		} else {
			this->generate_input_error(p, array, "Unequal sizes of sub-arrays during transpose auto-conversion of '" + p.get_node_type(ValueTypes::array_type) + "'", false);
		}
	}

	return NULL;
}

StorageBase * ReaderInternalTranspose::make_sub_storage(PathBase &p, const Type::Array *array)
{
	int arr_size;
	if ( (arr_size = p.get_array_size()) != -1 ) {
		return this->make_array_storage(p, array, arr_size);
	} else {
		// if transposition is carried, only conversion to array with one element is allowed
		// try automatic conversion to array with one element
		const Type::TypeBase &sub_type = array->get_sub_type();
		StorageBase *one_element_storage = this->make_storage(p, &sub_type);
		return make_autoconversion_array_storage(p, array, one_element_storage);
	}
}

StorageBase * ReaderInternalTranspose::make_sub_storage(PathBase &p, const Type::Selection *selection)
{
	if ( p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, selection);
	} else {
	    string item_name = read_string_value(p, selection);
		try {
			int value = selection->name_to_int( item_name );
			return new StorageInt( value );
		} catch (Type::Selection::ExcSelectionKeyNotFound &exc) {
			this->generate_input_error(p, selection, "Wrong value '" + item_name + "' of the Selection.", false);
		}
	}

    return NULL;
}

StorageBase * ReaderInternalTranspose::make_sub_storage(PathBase &p, const Type::Bool *bool_type)
{
	if ( p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, bool_type);
	} else {
		return new StorageBool( read_bool_value(p, bool_type) );
	}
}

StorageBase * ReaderInternalTranspose::make_sub_storage(PathBase &p, const Type::Integer *int_type)
{
	if ( p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, int_type);
	}
	std::int64_t value = read_int_value(p, int_type);

	if ( int_type->match(value) )
	{
		return new StorageInt( value );
	} else {
		this->generate_input_error(p, int_type, "Value out of bounds.", false);
	}
	return NULL;
}

StorageBase * ReaderInternalTranspose::make_sub_storage(PathBase &p, const Type::Double *double_type)
{
	if ( p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, double_type);
	} else {
	    double value = read_double_value(p, double_type);

	    if (double_type->match(value)) {
	        return new StorageDouble( value );
	    } else {
	    	this->generate_input_error(p, double_type, "Value out of bounds.", false);
	    }
	}

    return NULL;
}

StorageBase * ReaderInternalTranspose::make_sub_storage(PathBase &p, const Type::String *string_type)
{
	if ( p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, string_type);
	} else {
		return new StorageString( read_string_value(p, string_type) );
	}
}

StorageBase * ReaderInternalTranspose::make_transposed_storage(PathBase &p, const Type::TypeBase *type) {
	ASSERT(p.is_array_type()).error();

	int arr_size = p.get_array_size();
	if ( arr_size == 0 ) {
		this->generate_input_error(p, type, "Empty array during transpose auto-conversion.", false);
	} else {
		if (transpose_index_ == 0) transpose_array_sizes_.push_back( arr_size );
		p.down(transpose_index_);
		StorageBase *storage = make_storage(p, type);
		p.up();
		return storage;
	}

	return NULL;
}

StorageBase * ReaderInternalTranspose::make_autoconversion_array_storage(PathBase &p, const Type::Array *array, StorageBase *item)
{
	if ( array->match_size( 1 ) ) {
		StorageArray *storage_array = new StorageArray(1);
		storage_array->new_item(0, item);

		return storage_array;
	} else {
		this->generate_input_error(p, array, "During transpose auto-conversion, the conversion to the single element array not allowed. Require type: '" + p.get_node_type(ValueTypes::array_type) + "'\nFound on input: ", true);
	}

	return NULL;
}




/*******************************************************************
 * implementation of ReaderInternalCsvInclude
 */

ReaderInternalCsvInclude::ReaderInternalCsvInclude()
{}

StorageBase * ReaderInternalCsvInclude::read_storage(PathBase &p, const Type::Array *array)
{
	if ( p.is_record_type() ) { // sub-type must be record type
		// load path to CSV file
		std::string included_file;
        if ( p.down("file") ) {
       		included_file = get_included_file(p);
            p.up();
        } else {
        	this->generate_input_error(p, array, "Missing key 'file' defines including input file.", false);
        }

        // number of head lines to skip
        unsigned int n_head_lines = 0;
        if ( p.down("n_head_lines") ) {
        	try {
        		n_head_lines = p.get_int_value();
        	}
			catch (ExcInputError & e) {
				complete_input_error(e, p, ValueTypes::int_type);
				e << EI_InputType("number of lines to skip");
				throw;
			}
        	p.up();
        }

        // open CSV file, get number of lines, skip head lines
        FilePath fp((included_file), FilePath::input_file);
        CSVTokenizer tok( fp );

        const Type::TypeBase &sub_type = array->get_sub_type(); // sub-type of array
        StorageBase *item_storage; // storage of sub-type record of included array
        csv_columns_map_.clear();
        if ( p.down("format") ) {
			try {
				csv_subtree_depth_ = p.path_.size();
				item_storage = make_storage(p, &sub_type);
			} catch (ExcMultipleDefinitionCsvColumn &e) {
				e << EI_File(tok.f_name());
				throw;
			}
            p.up();
        } else {
        	this->generate_input_error(p, array, "Missing key 'format' defines mapping column of CSV file to input subtree.", false);
        }

        // get value of maximal column index in map
        map<unsigned int, IncludeCsvData>::iterator it;
        it = csv_columns_map_.end(); --it;
        unsigned int max_column_index = it->first;

        unsigned int n_lines = tok.get_n_lines() - n_head_lines;
        tok.skip_header(n_head_lines);
        StorageArray *storage_array = new StorageArray(n_lines);
        std::set<unsigned int> unused_columns;
        for( unsigned int arr_item=0; arr_item < n_lines; ++arr_item) {
        	unsigned int i_col;
        	tok.next_line();
        	for (i_col=0; !tok.eol(); ++i_col, ++tok) {
        		it = csv_columns_map_.find(i_col);
        		if (it != csv_columns_map_.end()) {
        			switch (it->second.data_type) {
						case IncludeDataTypes::type_int: {
							int val;
							try {
								val = tok.get_int_val();
							} catch (ExcWrongCsvFormat &e) {
								e << EI_Specification("Wrong integer value");
								e << EI_ErrorAddress(p.as_string());
								throw;
							}

							const Type::Integer *int_type = static_cast<const Type::Integer *>(it->second.type);
							if ( !int_type->match(val) ) {
								THROW( ExcWrongCsvFormat() << EI_Specification("Integer value out of bounds")
										<< EI_TokenizerMsg(tok.position_msg()) << EI_ErrorAddress(p.as_string()) );
							}
							set_storage_from_csv( i_col, item_storage, new StorageInt(val) );
							break;
						}
						case IncludeDataTypes::type_double: {
							double val;
							try {
								val = tok.get_double_val();
							} catch (ExcWrongCsvFormat &e) {
								e << EI_ErrorAddress(p.as_string());
								throw;
							}

							const Type::Double *double_type = static_cast<const Type::Double *>(it->second.type);
							if ( !double_type->match(val) ) {
								THROW( ExcWrongCsvFormat() << EI_Specification("Double value out of bounds")
										<< EI_TokenizerMsg(tok.position_msg()) << EI_ErrorAddress(p.as_string()) );
							}
							set_storage_from_csv( i_col, item_storage, new StorageDouble(val) );
							break;
						}
						case IncludeDataTypes::type_bool: {
							int val;
							try {
								val = tok.get_int_val();
							} catch (ExcWrongCsvFormat &e) {
								e << EI_Specification("Wrong boolean value");
								e << EI_ErrorAddress(p.as_string());
								throw;
							}
							set_storage_from_csv( i_col, item_storage, new StorageBool(val) );
							break;
						}
						case IncludeDataTypes::type_string: {
							try {
								set_storage_from_csv( i_col, item_storage, new StorageString(tok.get_string_val()) );
							} catch (ExcWrongCsvFormat &e) {
								e << EI_Specification("Wrong string value");
								e << EI_ErrorAddress(p.as_string());
								throw;
							}
							break;
						}
						case IncludeDataTypes::type_sel: {
							const Type::Selection *selection = static_cast<const Type::Selection *>(it->second.type);
							try {
								std::string item_name = tok.get_string_val();
								int val = selection->name_to_int( item_name );
								set_storage_from_csv( i_col, item_storage, new StorageInt(val) );
							} catch (ExcWrongCsvFormat &e) {
								e << EI_Specification("Wrong selection value");
								e << EI_ErrorAddress(p.as_string());
								throw;
							} catch (Type::Selection::ExcSelectionKeyNotFound &exc) {
								THROW( ExcWrongCsvFormat() << EI_Specification("Wrong selection value")
										<< EI_TokenizerMsg(tok.position_msg()) << EI_ErrorAddress(p.as_string()) );
							}
							break;
						}
        			}
        		} else {
        			// add index of unused column
        			unused_columns.insert(i_col);
        		}
        	}
    		if ( max_column_index > (i_col-1) ) {
    			this->generate_input_error(p, array, "Count of columns in CSV file is less than expected index, defined on input.", false);
    		}
            storage_array->new_item(arr_item, item_storage->deep_copy() );
        }

        if (unused_columns.size()) { // print warning with indexes of unused columns
        	stringstream ss;
        	for (std::set<unsigned int>::iterator it=unused_columns.begin(); it!=unused_columns.end(); ++it)
        		ss << (*it) << " ";
            WarningOut().fmt("Unused columns: {}\nin imported CSV input file: {}\n", ss.str(), tok.f_name());
        }
        return storage_array;

	} else {
		this->generate_input_error(p, array, "Invalid definition of CSV include.", false);
	}
	return NULL;
}

StorageBase * ReaderInternalCsvInclude::make_sub_storage(PathBase &p, const Type::Array *array)
{
	this->generate_input_error(p, array, "Array type in CSV-included part of IST is forbidden!\n", false);
    return NULL;
}

StorageBase * ReaderInternalCsvInclude::make_sub_storage(PathBase &p, const Type::Selection *selection)
{
	const Type::Integer *format_int = new Type::Integer(0);
	std::int64_t pos = read_int_value(p, format_int);

	IncludeCsvData include_data;
	include_data.data_type = IncludeDataTypes::type_sel;
	include_data.storage_indexes = create_indexes_vector(p);
	include_data.type = selection;
	if (csv_columns_map_.find(pos)!=csv_columns_map_.end()) {
		THROW( ExcMultipleDefinitionCsvColumn() << EI_ColumnIndex(pos) << EI_ErrorAddress(p.as_string()) );
	} else {
		csv_columns_map_[pos] = include_data;
	}
	delete format_int;
	return new StorageInt( 0 );
}

StorageBase * ReaderInternalCsvInclude::make_sub_storage(PathBase &p, const Type::Bool *bool_type)
{
	const Type::Integer *format_int = new Type::Integer(0);
	std::int64_t pos = read_int_value(p, format_int);

	IncludeCsvData include_data;
	include_data.data_type = IncludeDataTypes::type_bool;
	include_data.storage_indexes = create_indexes_vector(p);
	include_data.type = bool_type;
	if (csv_columns_map_.find(pos)!=csv_columns_map_.end()) {
		THROW( ExcMultipleDefinitionCsvColumn() << EI_ColumnIndex(pos) << EI_ErrorAddress(p.as_string()) );
	} else {
		csv_columns_map_[pos] = include_data;
	}
	delete format_int;

	return new StorageBool( false );
}

StorageBase * ReaderInternalCsvInclude::make_sub_storage(PathBase &p, const Type::Integer *int_type)
{
	const Type::Integer *format_int = new Type::Integer(0);
	std::int64_t value = read_int_value(p, format_int);

	IncludeCsvData include_data;
	include_data.data_type = IncludeDataTypes::type_int;
	include_data.storage_indexes = create_indexes_vector(p);
	include_data.type = int_type;
	if (csv_columns_map_.find(value)!=csv_columns_map_.end()) {
		THROW( ExcMultipleDefinitionCsvColumn() << EI_ColumnIndex(value) << EI_ErrorAddress(p.as_string()) );
	} else {
		csv_columns_map_[value] = include_data;
	}
	delete format_int;

	return new StorageInt( 0 );
}

StorageBase * ReaderInternalCsvInclude::make_sub_storage(PathBase &p, const Type::Double *double_type)
{
	const Type::Integer *format_int = new Type::Integer(0);
	std::int64_t pos = read_int_value(p, format_int);

	IncludeCsvData include_data;
	include_data.data_type = IncludeDataTypes::type_double;
	include_data.storage_indexes = create_indexes_vector(p);
	include_data.type = double_type;
	if (csv_columns_map_.find(pos)!=csv_columns_map_.end()) {
		THROW( ExcMultipleDefinitionCsvColumn() << EI_ColumnIndex(pos) << EI_ErrorAddress(p.as_string()) );
	} else {
		csv_columns_map_[pos] = include_data;
	}
	delete format_int;

	return new StorageDouble( 0.0 );
}

StorageBase * ReaderInternalCsvInclude::make_sub_storage(PathBase &p, const Type::String *string_type)
{
	try {
		const Type::Integer *format_int = new Type::Integer(0);
		std::int64_t pos = read_int_value(p, format_int);

		IncludeCsvData include_data;
		include_data.data_type = IncludeDataTypes::type_string;
		include_data.storage_indexes = create_indexes_vector(p);
		include_data.type = string_type;
		if (csv_columns_map_.find(pos)!=csv_columns_map_.end()) {
			THROW( ExcMultipleDefinitionCsvColumn() << EI_ColumnIndex(pos) << EI_ErrorAddress(p.as_string()) );
		} else {
			csv_columns_map_[pos] = include_data;
		}
		delete format_int;

		return new StorageString("");
	} catch (ExcInputError & e) {
		// no error, string value is not forbidden in CSV include
		return new StorageString( read_string_value(p, string_type) );
	}
}

vector<unsigned int> ReaderInternalCsvInclude::create_indexes_vector(PathBase &p)
{
	vector<unsigned int> csv_storage_indexes( p.path_.size()-csv_subtree_depth_ );
	for (unsigned int i_source=csv_subtree_depth_, i_target=0; i_source<p.path_.size(); ++i_source, ++i_target ) {
		ASSERT_GE(p.path_[i_source].first, 0).error();
		csv_storage_indexes[i_target] = p.path_[i_source].first;
	}
	return csv_storage_indexes;
}

void ReaderInternalCsvInclude::set_storage_from_csv(unsigned int column_index, StorageBase * item_storage, StorageBase * new_storage)
{
	map<unsigned int, IncludeCsvData>::iterator it = csv_columns_map_.find(column_index);
	ASSERT(it!=csv_columns_map_.end()).error();

	unsigned int i;
	StorageBase *loop_storage = item_storage;
	for (i=0; i<it->second.storage_indexes.size()-1; ++i) loop_storage = loop_storage->get_item( it->second.storage_indexes[i] );
	loop_storage->set_item( it->second.storage_indexes[i], new_storage );
}

} // namespace Input
