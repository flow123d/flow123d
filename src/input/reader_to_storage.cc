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
#include <boost/lexical_cast.hpp>

#include "reader_to_storage.hh"
#include "input/path_json.hh"
#include "input/path_yaml.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "system/tokenizer.hh"


namespace Input {
using namespace std;
using namespace internal;



/********************************************
 * Implementation of public part of ReaderToStorage
 */

ReaderToStorage::ReaderToStorage()
: storage_(nullptr),
  root_type_(nullptr),
  try_read_(TryRead::none),
  transpose_index_(0)
{}



ReaderToStorage::ReaderToStorage(const FilePath &in_file, Type::TypeBase &root_type)
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

	std::ifstream in;
	in_file.open_stream(in);

    // finish of root_type ensures finish of whole IST
    root_type.finish();

	read_stream(in, root_type, format);
}



ReaderToStorage::ReaderToStorage( const string &str, Type::TypeBase &root_type, FileFormat format)
: ReaderToStorage()
{
    // finish of root_type ensures finish of whole IST
    root_type.finish();

	try {
		istringstream is(str);
		read_stream(is, root_type, format);
	} catch (ExcNotJSONFormat &e) {
		e << EI_File("STRING: "+str); throw;
	}
}



StorageBase *ReaderToStorage::get_storage()
{
	return storage_;
}



void ReaderToStorage::read_stream(istream &in, const Type::TypeBase &root_type, FileFormat format)
{
	ASSERT(storage_==nullptr).error();

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

	ASSERT_PTR(storage_).error();
}






/********************************************
 * Implementation of private part of ReaderToStorage - make_storage dispatch
 */


StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::TypeBase *type)
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
        return make_storage(p, static_cast<const Type::Tuple *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Record)) {
        return make_storage(p, static_cast<const Type::Record *>(type) );
    } else
    if (typeid(*type) == typeid(Type::Array)) {
        return make_storage(p, static_cast<const Type::Array *>(type) );
    } else {
    	const Type::Abstract * abstract_record_type = dynamic_cast<const Type::Abstract *>(type);
    	if (abstract_record_type != NULL ) return make_storage(p, abstract_record_type );
    }

    // return Null storage if there is null on the current location
	if (p.is_null_type()) {
		return new StorageNull();
	}

    // dispatch types - scalar types
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
        const Type::String * string_type = dynamic_cast<const Type::String *>(type);
        if (string_type != NULL ) return make_storage(p, string_type );

        // default -> error
        THROW( Type::ExcUnknownDescendant() << Type::EI_TypeName(typeid(type).name()) );
    }

    return new StorageNull();
}


StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::Record *record)
{
	// control test, check correct tag (or TYPE key) if Record is derived from Abstract
	string record_name_from_tag = p.get_record_tag();
	if (record_name_from_tag == "include") {
		return make_include_storage(p, record);
	} else if (record_name_from_tag == "include_csv") {
		ASSERT(false).error("Include of CSV is not supported yet.\n");
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

	            if ( !effectively_null && p.down(it->key_) ) {
	                // key on input => check & use it
	                // check for obsolete key

	                auto obsolete_it = it->attributes.find( Type::Attribute::obsolete() );
	                if ( obsolete_it != it->attributes.end()) {
	                    WarningOut() << "Usage of the obsolete key: '" << it->key_ << "'\n" << obsolete_it -> second;
	                }

	                if (try_read_ == TryRead::csv_include) csv_storage_indexes_.push_back(it->key_index);
	                StorageBase *storage = make_storage(p, it->type_.get());
	                if ( (typeid(*storage) == typeid(StorageNull)) && it->default_.has_value_at_declaration() ) {
	                	delete storage;
	                	storage = make_storage_from_default( it->default_.value(), it->type_ );
	                }
	                storage_array->new_item( it->key_index, storage );
	                if (try_read_ == TryRead::csv_include) csv_storage_indexes_.pop_back();
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
					 ASSERT(! it->default_.is_obligatory())(it->key_).error("Obligatory key in auto-convertible Record.");
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
	string record_name = p.get_record_tag();
	if (record_name == "") {
		if ( ! abstr_rec->get_selection_default().has_value_at_declaration() ) {
			THROW( ExcInputError() << EI_Specification("Can not determine type of the Abstract.") << EI_ErrorAddress(p.as_string())
					<< EI_JSON_Type( p.get_node_type(p.get_node_type_index()) ) << EI_InputType(abstr_rec->desc()) );
		} else { // auto conversion
			return abstract_automatic_conversion(p, abstr_rec);
		}
	} else if ((record_name == "include") || (record_name == "include_csv")) {
		THROW( ExcForbiddenAbstractTag() << EI_Tag(record_name) );
	} else {
		try {
			return make_storage(p, &( abstr_rec->get_descendant(record_name) ) );
		} catch (Type::Selection::ExcSelectionKeyNotFound &exc) {
			THROW( ExcInputError() << EI_Specification("Wrong value '" + record_name + "' of the Selection.")
					<< EI_ErrorAddress(p.as_string()) << EI_JSON_Type( "" ) << EI_InputType(abstr_rec->get_type_selection().desc()) );
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
	if ( try_read_ == TryRead::csv_include) {
		THROW( ExcInputError() << EI_Specification("Array type in CSV-included part of IST is forbidden!\n")
						   << EI_ErrorAddress(p.as_string()) << EI_InputType(array->desc()) );
	}

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
        switch (try_read_) {
			case TryRead::transposed: {
				// if transposition is carried, only conversion to array with one element is allowed
				// try automatic conversion to array with one element
				const Type::TypeBase &sub_type = array->get_sub_type();
				StorageBase *one_element_storage = make_storage(p, &sub_type);
				return make_autoconversion_array_storage(p, array, one_element_storage);
			}
			case TryRead::csv_include: {
				// Case doesn't happen. See at the begin of this method.
				break;
			}
			default: {
				if (p.get_record_tag() == "include_csv") {
					try_read_ = TryRead::csv_include;
					return this->make_include_csv_storage(p, array);
				} else {
					// set variables managed transposition
					try_read_ = TryRead::transposed;
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
						try_read_ = TryRead::none;
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

							try_read_ = TryRead::none;
							return storage_array;
						} else {
							THROW( ExcInputError()
									<< EI_Specification("Unequal sizes of sub-arrays during transpose auto-conversion of '" + p.get_node_type(ValueTypes::array_type) + "'")
									<< EI_ErrorAddress(p.as_string()) << EI_InputType(array->desc()) );
						}
					}
				}
				break;
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
        		if (try_read_ == TryRead::csv_include) csv_storage_indexes_.push_back(it->key_index);
                StorageBase *storage = make_storage(p, it->type_.get());
                if ( (typeid(*storage) == typeid(StorageNull)) && it->default_.has_value_at_declaration() ) {
                	delete storage;
                	storage = make_storage_from_default( it->default_.value(), it->type_ );
                }
                storage_array->new_item( it->key_index, storage );
                if (try_read_ == TryRead::csv_include) csv_storage_indexes_.pop_back();
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
			WarningOut().fmt("Unprocessed keys in tuple '{}', tuple has {} keys but the input is specified by {} values.\n",
                    p.as_string().c_str(), tuple->size(), arr_size );
		}

        return storage_array;

	} else {
		return make_storage(p, static_cast<const Type::Record *>(tuple) );
	}
}



StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::Selection *selection)
{
	if ( (try_read_ == TryRead::transposed) && p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, selection);
	} else if ( try_read_ == TryRead::csv_include) {
		try {
			unsigned int pos = p.get_int_value();
			IncludeCsvData include_data;
			include_data.data_type = IncludeDataTypes::type_sel;
			include_data.storage_indexes = csv_storage_indexes_;
			include_data.type = selection;
			if (csv_columns_map_.find(pos)!=csv_columns_map_.end()) {
				ASSERT(false).error("Change to exception");
			} else {
				csv_columns_map_[pos] = include_data;
			}
			return new StorageInt( 0 );
		} catch (ExcInputError & e) {
			e << EI_Specification("The value in definition of CSV format should be '" + p.get_node_type(ValueTypes::int_type) + "', but we found: ");
			e << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) );
			e << EI_ErrorAddress(p.as_string());
			e << EI_InputType(selection->desc());
			throw;
		}
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
	if ( (try_read_ == TryRead::transposed) && p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, bool_type);
	} else if (try_read_ == TryRead::csv_include)
		try {
			unsigned int pos = p.get_int_value();
			IncludeCsvData include_data;
			include_data.data_type = IncludeDataTypes::type_bool;
			include_data.storage_indexes = csv_storage_indexes_;
			include_data.type = bool_type;
			if (csv_columns_map_.find(pos)!=csv_columns_map_.end()) {
				ASSERT(false).error("Change to exception");
			} else {
				csv_columns_map_[pos] = include_data;
			}
			return new StorageBool( false );
		} catch (ExcInputError & e) {
			e << EI_Specification("The value in definition of CSV format should be '" + p.get_node_type(ValueTypes::int_type) + "', but we found: ");
			e << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) );
			e << EI_ErrorAddress(p.as_string());
			e << EI_InputType(bool_type->desc());
			throw;
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
	if ( (try_read_ == TryRead::transposed) && p.is_array_type() ) {
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

	if (try_read_ == TryRead::csv_include) {
		IncludeCsvData include_data;
		include_data.data_type = IncludeDataTypes::type_int;
		include_data.storage_indexes = csv_storage_indexes_;
		include_data.type = int_type;
		if (csv_columns_map_.find(value)!=csv_columns_map_.end()) {
			ASSERT(false).error("Change to exception");
		} else {
			csv_columns_map_[value] = include_data;
		}
		return new StorageInt( 0 );
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
	if ( (try_read_ == TryRead::transposed) && p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, double_type);
	} else if (try_read_ == TryRead::csv_include) {
		try {
			// read index of column in CSV file
			unsigned int pos = p.get_int_value();
			IncludeCsvData include_data;
			include_data.data_type = IncludeDataTypes::type_double;
			include_data.storage_indexes = csv_storage_indexes_;
			include_data.type = double_type;
			if (csv_columns_map_.find(pos)!=csv_columns_map_.end()) {
				ASSERT(false).error("Change to exception");
			} else {
				csv_columns_map_[pos] = include_data;
			}
			return new StorageDouble( 0.0 );
		} catch (ExcInputError & e) {
			e << EI_Specification("The value in definition of CSV format should be '" + p.get_node_type(ValueTypes::int_type) + "', but we found: ");
			e << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) );
			e << EI_ErrorAddress(p.as_string());
			e << EI_InputType(double_type->desc());
			throw;
		}
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
	if ( (try_read_ == TryRead::transposed) && p.is_array_type() ) {
		// transpose auto-conversion for array type
		return this->make_transposed_storage(p, string_type);
	} else if (try_read_ == TryRead::csv_include) {
		try {
			// read index of column in CSV file
			unsigned int pos = p.get_int_value();
			IncludeCsvData include_data;
			include_data.data_type = IncludeDataTypes::type_string;
			include_data.storage_indexes = csv_storage_indexes_;
			include_data.type = string_type;
			if (csv_columns_map_.find(pos)!=csv_columns_map_.end()) {
				ASSERT(false).error("Change to exception");
			} else {
				csv_columns_map_[pos] = include_data;
			}
			return new StorageString("");
		} catch (ExcInputError & e) {
			// no error, string value is not forbidden in CSV include
		}
	}
	string value;
	try {
		value = p.get_string_value();
	} catch (ExcInputError & e) {
		if ((try_read_ == TryRead::transposed)) {
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



StorageBase * ReaderToStorage::make_storage_from_default(const string &dflt_str, std::shared_ptr<Type::TypeBase> type) {
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



StorageBase * ReaderToStorage::make_transposed_storage(PathBase &p, const Type::TypeBase *type) {
	ASSERT(try_read_ == TryRead::transposed).error();
	ASSERT(p.is_array_type()).error();

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



StorageBase * ReaderToStorage::make_include_storage(PathBase &p, const Type::Record *record)
{
    std::string included_path;
    if ( p.is_record_type() ) {
        // include is set as record with tag and file key
        if ( p.down("file") ) {
        	included_path = get_included_file(p);
            p.up();
        } else {
    	    THROW( ExcInputError() << EI_Specification("Missing key 'file' defines including input file.")
                               << EI_ErrorAddress(p.as_string()) << EI_InputType(record->desc()) );
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
    } catch (ReaderToStorage::ExcInputError &e ) {
      e << ReaderToStorage::EI_File(fpath); throw;
    } catch (ReaderToStorage::ExcNotJSONFormat &e) {
      e << ReaderToStorage::EI_File(fpath); throw;
    }

	return NULL;
}


StorageBase * ReaderToStorage::make_include_csv_storage(PathBase &p, const Type::Array *array)
{
	using namespace boost;

	if ( p.is_record_type() ) { // sub-type must be record type
		// load path to CSV file
		std::string included_file;
        if ( p.down("file") ) {
       		included_file = get_included_file(p);
            p.up();
        } else {
    	    THROW( ExcInputError() << EI_Specification("Missing key 'file' defines including input file.")
                               << EI_ErrorAddress(p.as_string()) << EI_InputType(array->desc()) );
        }

        // number of head lines to skip
        unsigned int n_head_lines = 0;
        if ( p.down("n_head_lines") ) {
        	try {
        		n_head_lines = p.get_int_value();
        	}
			catch (ExcInputError & e) {
				e << EI_Specification("The value should be '" + p.get_node_type(ValueTypes::int_type) + "', but we found: ");
				e << EI_ErrorAddress(p.as_string());
				e << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) );
				e << EI_InputType("number of lines to skip");
				throw;
			}
        	p.up();
        }

        // sub-type of array
        const Type::TypeBase &sub_type = array->get_sub_type();
        // for every leaf of input subtree holds index of columns in CSV file
        StorageBase *item_storage;
        csv_storage_indexes_.clear();
        csv_columns_map_.clear();
        if ( p.down("format") ) {
            item_storage = make_storage(p, &sub_type);
            p.up();
        } else {
    	    THROW( ExcInputError() << EI_Specification("Missing key 'format' defines mapping column of CSV file to input subtree.")
                               << EI_ErrorAddress(p.as_string()) << EI_InputType(array->desc()) );
        }
        ASSERT_EQ(csv_storage_indexes_.size(), 0).error();

        // open CSV file, get number of lines, skip head lines
        FilePath fp((included_file), FilePath::input_file);
        Tokenizer tok( fp );
        unsigned int n_lines = 0; // number of lines
        while ( !tok.eof() ) {
        	tok.next_line(false);
        	n_lines++;
        }
        if (tok.line().size()==0) n_lines--; // removes last line if it is empty
        n_lines -= n_head_lines; // subtracts number of skipped head lines
        tok.set_position( Tokenizer::Position() );
        for (unsigned int i=0; i<n_head_lines; i++) { // skip head lines
        	tok.next_line(false);
        }

        StorageArray *storage_array = new StorageArray(n_lines);
        for( unsigned int arr_item=0; arr_item < n_lines; ++arr_item) {
        	tok.next_line();
        	for (unsigned int i_col=0; !tok.eol(); ++i_col, ++tok) {
        		map<unsigned int, IncludeCsvData>::iterator it = csv_columns_map_.find(i_col);
        		if (it != csv_columns_map_.end()) {
        			switch (it->second.data_type) {
						case IncludeDataTypes::type_int: {
							int val = lexical_cast<int>(*tok);
							set_storage_from_csv( i_col, item_storage, new StorageInt(val) );
							break;
						}
						case IncludeDataTypes::type_double: {
							double val = lexical_cast<double>(*tok);
							set_storage_from_csv( i_col, item_storage, new StorageDouble(val) );
							break;
						}
						case IncludeDataTypes::type_bool: {
							int val = lexical_cast<int>(*tok);
							set_storage_from_csv( i_col, item_storage, new StorageBool(val) );
							break;
						}
						case IncludeDataTypes::type_string: {
							set_storage_from_csv( i_col, item_storage, new StorageString(*tok) );
							break;
						}
						case IncludeDataTypes::type_sel: {
							std::string item_name;
							const Type::Selection *selection = static_cast<const Type::Selection *>(it->second.type);
							{
								item_name = *tok;
								int val = selection->name_to_int( item_name );
								set_storage_from_csv( i_col, item_storage, new StorageInt(val) );
							}
							/*catch (ExcInputError & e) {
								e << EI_Specification("The value should be '" + p.get_node_type(ValueTypes::str_type) + "', but we found: ");
						        e << EI_ErrorAddress(p.as_string());
						        e << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) );
								e << EI_InputType(selection->desc());
								throw;
							} catch (Type::Selection::ExcSelectionKeyNotFound &exc) {
								THROW( ExcInputError() << EI_Specification("Wrong value '" + item_name + "' of the Selection.")
										<< EI_ErrorAddress(p.as_string()) << EI_JSON_Type( "" ) << EI_InputType(selection->desc()) );
							}*/
							// TODO complete try-catch blocks in all tokenizer readings
							break;
						}
        			}
        		} else {
        			// add to warning
        		}
        	}

            storage_array->new_item(arr_item, item_storage->deep_copy() );
        }
        try_read_ = TryRead::none;
        return storage_array;

	} else {
	    THROW( ExcInputError() << EI_Specification("Invalid definition of CSV include.")
                           << EI_ErrorAddress(p.as_string()) << EI_InputType(array->desc()) );
	}
	return NULL;
}



std::string ReaderToStorage::get_included_file(PathBase &p)
{
	try {
		return p.get_string_value();
	}
	catch (ExcInputError & e) {
		e << EI_Specification("The value should be '" + p.get_node_type(ValueTypes::str_type) + "', but we found: ");
        e << EI_ErrorAddress(p.as_string());
        e << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) );
		e << EI_InputType("path to included file");
		throw;
	}
}



void ReaderToStorage::set_storage_from_csv(unsigned int column_index, StorageBase * item_storage, StorageBase * new_storage)
{
	map<unsigned int, IncludeCsvData>::iterator it = csv_columns_map_.find(column_index);
	ASSERT(it!=csv_columns_map_.end()).error();

	unsigned int i;
	StorageBase *loop_storage = item_storage;
	for (i=0; i<it->second.storage_indexes.size()-1; ++i) loop_storage = loop_storage->get_item( it->second.storage_indexes[i] );
	loop_storage->set_item( it->second.storage_indexes[i], new_storage );
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
	ASSERT_PTR(storage_).error();

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
