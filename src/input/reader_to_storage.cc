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
 * @file    reader_to_storage.cc
 * @brief   
 */

#include <cstdint>
#include <limits>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "reader_to_storage.hh"
#include "input/path_json.hh"
#include "input/path_yaml.hh"
#include "input/accessors.hh"

namespace Input {
using namespace std;
using namespace internal;



/********************************************
 * Implementation of public part of ReaderToStorage
 */

ReaderToStorage::ReaderToStorage()
: storage_(nullptr),
  root_type_(nullptr),
  try_transpose_read_(false)
{}



ReaderToStorage::ReaderToStorage(const FilePath &in_file, const Type::TypeBase &root_type)
: ReaderToStorage()
{
	std::string fname = in_file;
	std::string extension = fname.substr(fname.find_last_of(".") + 1);
	FileFormat format;
	if (extension == "con") {
		format = FileFormat::format_JSON;
	} else if (extension == "yaml") {
		format = FileFormat::format_YAML;
	} else {
		THROW(ExcInputMessage() << EI_Message("Invalid extension of file " + fname + ".\nMust be 'con' or 'yaml'."));
	}

	std::ifstream in(fname.c_str());
    if (! in) {
    	THROW(ExcInputMessage() << EI_Message("Can not open main input file: '" + fname + "'.\n"));
    }

    // finish all lazy input types
    Input::Type::TypeBase::lazy_finish();

	read_stream(in, root_type, format);
}



ReaderToStorage::ReaderToStorage( const string &str, const Type::TypeBase &root_type, FileFormat format)
: ReaderToStorage()
{
	// finish all lazy input types
    Input::Type::TypeBase::lazy_finish();

	try {
		istringstream is(str);
		read_stream(is, root_type, format);
	} catch (ExcNotJSONFormat &e) {
		e << EI_File("STRING: "+str); throw;
	}
}



void ReaderToStorage::read_stream(istream &in, const Type::TypeBase &root_type, FileFormat format)
{
	OLD_ASSERT(storage_==nullptr," ");

    PathBase * root_path;
	if (format == FileFormat::format_JSON) {
		root_path = new PathJSON(in);
	} else {
		root_path = new PathYAML(in);
	}

	// guarantee to delete root_path on function return even on exception
	std::unique_ptr<PathBase> root_path_ptr(root_path);

    root_type_ = &root_type;
	try {
	    storage_ = make_storage(*root_path_ptr, root_type_);
	} catch (ExcInputError &e) {
		if (format == FileFormat::format_JSON) {
			e << EI_Format("JSON");
		} else {
			e << EI_Format("YAML");
		}
		throw;
	}

	OLD_ASSERT(  storage_ != nullptr, "Internal error in Input reader, the storage pointer is NULL after reading the stream.\n");
}






/********************************************
 * Implementation of private part of ReaderToStorage - make_storage dispatch
 */


StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::TypeBase *type)
{
	OLD_ASSERT(type != NULL, "Can not dispatch, NULL pointer to TypeBase.\n");

    // find reference node, if doesn't exist return NULL
    PathBase * ref_path = p.find_ref_node();
    if (ref_path) {
        // todo: mark passed references and check cyclic references

        // dereference and take data from there
    	StorageBase * storage = make_storage( *ref_path, type );
    	delete ref_path;
        return storage;
    }

    // return Null storage if there is null on the current location
    if (p.is_null_type())
        return new StorageNull();

    // dispatch types
    if (typeid(*type) == typeid(Type::Tuple)) {
        return make_storage(p, static_cast<const Type::Tuple *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Record)) {
        return make_storage(p, static_cast<const Type::Record *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Array)) {
        return make_storage(p, static_cast<const Type::Array *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Integer)) {
        return make_storage(p, static_cast<const Type::Integer *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Double)) {
        return make_storage(p, static_cast<const Type::Double *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Bool)) {
        return make_storage(p, static_cast<const Type::Bool *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Selection)) {
        return make_storage(p, static_cast<const Type::Selection *>(type) );
    } else {
    	const Type::Abstract * abstract_record_type = dynamic_cast<const Type::Abstract *>(type);
    	if (abstract_record_type != NULL ) return make_storage(p, abstract_record_type );

        const Type::String * string_type = dynamic_cast<const Type::String *>(type);
        if (string_type != NULL ) return make_storage(p, string_type );

        // default -> error
        OLD_ASSERT(false, "Unknown descendant of TypeBase class, name: %s\n", typeid(type).name());
    }

    return new StorageNull();
}


StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::Record *record)
{
	std::set<string> keys_to_process;
	if ( p.get_record_key_set(keys_to_process) ) {
        std::set<string>::iterator set_it;

        /*Type::Record::KeyIter key_it;
        if ( record->has_key_iterator("TYPE", key_it) && record->auto_conversion_key_iter() != record->end() ) {
            PathBase *type_path = p->clone();
            if ( type_path.down( "TYPE" ) ) {
                try {
                	if ( type_path.get_string_value() != record->type_name() ) {
                		xprintf(UsrErr, "Invalid value of TYPE key of record %s.", record->type_name().c_str());
                	}
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

            if ( p.down(it->key_) ) {
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
                    THROW( ExcInputError() << EI_Specification("Missing obligatory key '"+ it->key_ +"'.")
                            << EI_ErrorAddress(p.as_string()) << EI_InputType(record->desc()) );
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
        	xprintf(Warn, "Unprocessed key '%s' in %s '%s'.\n", (*set_it).c_str(), record->class_name().c_str(), p.as_string().c_str() );
        }

        return storage_array;

    } else { // automatic conversion
    	return record_automatic_conversion(p, record);
    }
    // possibly construction of reduced record
}


StorageBase * ReaderToStorage::record_automatic_conversion(PathBase &p, const Type::Record *record)
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
					 OLD_ASSERT( ! it->default_.is_obligatory() ,
							 "Obligatory key: '%s' in auto-convertible %s, wrong check during finish().",
							 it->key_.c_str(), record->class_name().c_str() );
					 // set null
					 storage_array->new_item(it->key_index, new StorageNull() );
				 }
			}

			return storage_array;
	    } catch (ExcInputError &e ) {
	        THROW( ExcAutomaticConversionError() << EI_RecordName(record->type_name()) << EI_InputErrorMessage(e.what()) );
	    }

	} else {
	    THROW( ExcInputError() << EI_Specification("The value should be '" + p.get_node_type(ValueTypes::obj_type) + "', but we found: ")
	            << EI_ErrorAddress(p.as_string()) << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) )
				<< EI_InputType( record->desc()) );
	}

	return NULL;
}



StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::Abstract *abstr_rec)
{
	if ( p.is_record_type() ) {

		string descendant_name = p.get_descendant_name();
		if ( descendant_name == "" ) {
			if ( ! abstr_rec->get_selection_default().has_value_at_declaration() ) {
				THROW( ExcInputError() << EI_Specification("Missing key 'TYPE' in Abstract.") << EI_ErrorAddress(p.as_string()) << EI_InputType(abstr_rec->desc()) );
			} else { // auto conversion
				return abstract_automatic_conversion(p, abstr_rec);
			}
		} else {
			try {
				return make_storage(p, &( abstr_rec->get_descendant(descendant_name) ) );
			} catch (Type::Selection::ExcSelectionKeyNotFound &exc) {
				THROW( ExcInputError() << EI_Specification("Wrong value '" + descendant_name + "' of the Selection.")
						<< EI_ErrorAddress(p.as_string()) << EI_JSON_Type( "" ) << EI_InputType(abstr_rec->get_type_selection().desc()) );
			}
		}
	} else {
		if ( ! abstr_rec->get_selection_default().has_value_at_declaration() ) {
			THROW( ExcInputError() << EI_Specification("The value should be '" + p.get_node_type(ValueTypes::obj_type) + "', but we found: ")
				<< EI_ErrorAddress(p.as_string()) << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) ) << EI_InputType(abstr_rec->desc()) );
		} else { // auto conversion
			return abstract_automatic_conversion(p, abstr_rec);
		}
	}

	return NULL;
}



StorageBase * ReaderToStorage::abstract_automatic_conversion(PathBase &p, const Type::Abstract *abstr_rec)
{
    // perform automatic conversion
    const Type::Record *default_child = abstr_rec->get_default_descendant();
    if (! default_child) THROW(ExcInputError()
    		<< EI_Specification("Auto conversion of Abstract not allowed.\n")
    		<< EI_ErrorAddress(p.as_string())
    		<< EI_InputType(abstr_rec->desc())
    		);
    return make_storage(p, default_child );
}


StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::Array *array)
{
	int arr_size;
	if ( (arr_size = p.get_array_size()) != -1 ) {
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
        	ss << arr_size;
            THROW( ExcInputError()
                    << EI_Specification("Do not fit the size " + ss.str() + " of the Array.")
                    << EI_ErrorAddress(p.as_string()) << EI_InputType(array->desc()) );
        }
    } else {
    	// if transposition is carried, only conversion to array with one element is allowed
    	if (try_transpose_read_) {
			// try automatic conversion to array with one element
    		const Type::TypeBase &sub_type = array->get_sub_type();
    		StorageBase *one_element_storage = make_storage(p, &sub_type);
    		return make_autoconversion_array_storage(p, array, one_element_storage);
        } else {
			// set variables managed transposition
			try_transpose_read_ = true;
			transpose_index_ = 0;
			transpose_array_sizes_.clear();

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
				try_transpose_read_ = false;
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
						ss << sizes;
						THROW( ExcInputError() << EI_Specification("Result of transpose auto-conversion do not fit the size " + ss.str() + " of the Array.")
								<< EI_ErrorAddress(p.as_string()) << EI_InputType(array->desc()) );
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

					try_transpose_read_ = false;
					return storage_array;
				} else {
					THROW( ExcInputError()
							<< EI_Specification("Unequal sizes of sub-arrays during transpose auto-conversion of '" + p.get_node_type(ValueTypes::array_type) + "'")
							<< EI_ErrorAddress(p.as_string()) << EI_InputType(array->desc()) );
				}
        	}
    	}

    }

    return NULL;
}



StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::Tuple *tuple)
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
                	ss << tuple->obligatory_keys_count();
                    THROW( ExcInputError()
                    		<< EI_Specification("Too small size of '" + p.get_node_type(ValueTypes::array_type) + "' defining Tuple with "
                    							+ ss.str() + " obligatory keys.")
                            << EI_ErrorAddress(p.as_string())
							<< EI_InputType(tuple->desc()) );
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
            xprintf(Warn, "Unprocessed keys in tuple '%s', tuple has %d keys but the input is specified by %d values.\n",
                    p.as_string().c_str(), tuple->size(), arr_size );
        }

        return storage_array;

	} else {
		return make_storage(p, static_cast<const Type::Record *>(tuple) );
	}
}



StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::Selection *selection)
{
	if ( try_transpose_read_ && p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, selection);
	}
    string item_name;
	try {
		item_name = p.get_string_value();
		int value = selection->name_to_int( item_name  );
		return new StorageInt( value );
	} catch (ExcInputError & e) {
		e << EI_Specification("The value should be '" + p.get_node_type(ValueTypes::str_type) + "', but we found: ");
        e << EI_ErrorAddress(p.as_string());
        e << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) );
		e << EI_InputType(selection->desc());
		throw;
	} catch (Type::Selection::ExcSelectionKeyNotFound &exc) {
		THROW( ExcInputError() << EI_Specification("Wrong value '" + item_name + "' of the Selection.")
				<< EI_ErrorAddress(p.as_string()) << EI_JSON_Type( "" ) << EI_InputType(selection->desc()) );
	}

    return NULL;
}



StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::Bool *bool_type)
{
	if ( try_transpose_read_ && p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, bool_type);
	}
	try {
		return new StorageBool( p.get_bool_value() );
	}
	catch (ExcInputError & e) {
		e << EI_Specification("The value should be '" + p.get_node_type(ValueTypes::bool_type) + "', but we found: ");
		e << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) );
		e << EI_ErrorAddress(p.as_string());
		e << EI_InputType(bool_type->desc());
		throw;
	}
    return NULL;
}



StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::Integer *int_type)
{
	if ( try_transpose_read_ && p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, int_type);
	}
	std::int64_t value;
	try {
		value = p.get_int_value();
	}
	catch (ExcInputError & e) {
		e << EI_Specification("The value should be '" + p.get_node_type(ValueTypes::int_type) + "', but we found: ");
		e << EI_ErrorAddress(p.as_string());
		e << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) );
		e << EI_InputType(int_type->desc());
		throw;
	}

	if ( int_type->match(value) )
	{
		return new StorageInt( value );
	} else {
		THROW( ExcInputError() << EI_Specification("Value out of bounds.") << EI_ErrorAddress(p.as_string())
				<< EI_InputType(int_type->desc()) );
	}

    return NULL;
}



StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::Double *double_type)
{
	if ( try_transpose_read_ && p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, double_type);
	}
    double value;

	try {
		value = p.get_double_value();
	}
	catch (ExcInputError & e) {
		e << EI_Specification("The value should be '" + p.get_node_type(ValueTypes::real_type) + "', but we found: ");
		e << EI_ErrorAddress(p.as_string());
		e << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) );
		e << EI_InputType(double_type->desc());
		throw;
	}

    if (double_type->match(value)) {
        return new StorageDouble( value );
    } else {
        THROW( ExcInputError() << EI_Specification("Value out of bounds.") << EI_ErrorAddress(p.as_string())
        		<< EI_InputType(double_type->desc()) );
    }

    return NULL;
}



StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::String *string_type)
{
	if ( try_transpose_read_ && p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, string_type);
	}
	string value;
	try {
		value = p.get_string_value();
	}
	catch (ExcInputError & e) {
		if (try_transpose_read_) {
			return this->make_transposed_storage(p, string_type);
		}
		e << EI_Specification("The value should be '" + p.get_node_type(ValueTypes::str_type) + "', but we found: ");
        e << EI_ErrorAddress(p.as_string());
        e << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) );
		e << EI_InputType(string_type->desc());
		throw;
	}

	if (string_type->match(value))
		return new StorageString( value );
	else
		THROW( ExcInputError() << EI_Specification("Output file can not be given by absolute path: '" + value + "'")
						<< EI_ErrorAddress(p.as_string()) << EI_JSON_Type("") << EI_InputType(string_type->desc()) );

	return NULL;
}



StorageBase * ReaderToStorage::make_storage_from_default(const string &dflt_str, boost::shared_ptr<Type::TypeBase> type) {
    try {
    	// default strings must be valid JSON
    	Type::Default dflt(dflt_str);
    	return dflt.get_storage(type);

    } catch (Input::Type::ExcWrongDefault & e) {
        // message to distinguish exceptions thrown during Default value check at declaration
        xprintf(Msg, "Wrong default value while reading an input stream:\n");
        e << EI_KeyName("UNKNOWN KEY");
        throw;
    }

    return NULL;
}



StorageBase * ReaderToStorage::make_transposed_storage(PathBase &p, const Type::TypeBase *type) {
	OLD_ASSERT(try_transpose_read_, "Unset flag try_transpose_read_!\n");
	OLD_ASSERT(p.is_array_type(), "Head node of path must be of type array!\n");

	int arr_size = p.get_array_size();
	if ( arr_size == 0 ) {
		THROW( ExcInputError() << EI_Specification("Empty array during transpose auto-conversion.")
			<< EI_ErrorAddress(p.as_string()) << EI_InputType(type->desc()) );
	} else {
		if (transpose_index_ == 0) transpose_array_sizes_.push_back( arr_size );
		p.down(transpose_index_);
		StorageBase *storage = make_storage(p, type);
		p.up();
		return storage;
	}

	return NULL;
}



StorageBase * ReaderToStorage::make_autoconversion_array_storage(PathBase &p, const Type::Array *array, StorageBase *item)
{
	if ( array->match_size( 1 ) ) {
		StorageArray *storage_array = new StorageArray(1);
		storage_array->new_item(0, item);

		return storage_array;
	} else {
		THROW( ExcInputError()
				<< EI_Specification("During transpose auto-conversion, the conversion to the single element array not allowed. Require type: '" + p.get_node_type(ValueTypes::array_type) + "'\nFound on input: ")
				<< EI_ErrorAddress(p.as_string()) << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) ) << EI_InputType(array->desc()) );
	}

	return NULL;
}



/********************************************88
 * Implementation
 */

template <class T>
T ReaderToStorage::get_root_interface() const
{
	OLD_ASSERT(storage_, "NULL pointer to storage !!! \n");

    Address addr(storage_, root_type_);
    // try to create an iterator just to check type
    Iterator<T>( *root_type_, addr, 0);

    auto tmp_root_type = static_cast<const typename T::InputType &>(*root_type_);
    return T( addr, tmp_root_type );
}


template ::Input::Record ReaderToStorage::get_root_interface<::Input::Record>() const;
template ::Input::Array ReaderToStorage::get_root_interface<::Input::Array>() const;
template ::Input::AbstractRecord ReaderToStorage::get_root_interface<::Input::AbstractRecord>() const;
template ::Input::Tuple ReaderToStorage::get_root_interface<::Input::Tuple>() const;
//template ReaderToStorage::get_root_interface<::Input::>()->::Input::AbstractRecord const;


} // namespace Input
