/*
 * reader_to_storage.cc
 *
 *  Created on: May 7, 2012
 *      Author: jb
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
  root_type_(nullptr)
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
	read_stream(in, root_type, format);
}



ReaderToStorage::ReaderToStorage( const string &str, const Type::TypeBase &root_type, FileFormat format)
: ReaderToStorage()
{
	try {
		istringstream is(str);
		read_stream(is, root_type, format);
	} catch (ExcNotJSONFormat &e) {
		e << EI_File("STRING: "+str); throw;
	}
}



void ReaderToStorage::read_stream(istream &in, const Type::TypeBase &root_type, FileFormat format)
{
    ASSERT(storage_==nullptr," ");

    // finish all lazy input types
    Input::Type::TypeBase::lazy_finish();

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

    ASSERT(  storage_ != nullptr, "Internal error in Input reader, the storage pointer is NULL after reading the stream.\n");
}






/********************************************
 * Implementation of private part of ReaderToStorage - make_storage dispatch
 */


StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::TypeBase *type)
{
    ASSERT(type != NULL, "Can not dispatch, NULL pointer to TypeBase.\n");

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
    	const Type::AbstractRecord * abstract_record_type = dynamic_cast<const Type::AbstractRecord *>(type);
    	if (abstract_record_type != NULL ) return make_storage(p, abstract_record_type );

        const Type::String * string_type = dynamic_cast<const Type::String *>(type);
        if (string_type != NULL ) return make_storage(p, string_type );

        // default -> error
        ASSERT(false, "Unknown descendant of TypeBase class, name: %s\n", typeid(type).name());
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
                storage_array->new_item(it->key_index, make_storage(p, it->type_.get()) );
                p.up();
            } else {
                // key not on input
                if (it->default_.is_obligatory() ) {
                    THROW( ExcInputError() << EI_Specification("Missing obligatory key '"+ it->key_ +"'.")
                            << EI_ErrorAddress(p.as_string()) << EI_InputType(record->desc()) );
                } else if (it->default_.has_value_at_declaration() ) {
                   storage_array->new_item(it->key_index,
                           make_storage_from_default( it->default_.value(), it->type_.get() ) );
                } else { // defalut - optional or default at read time
                    // set null
                    storage_array->new_item(it->key_index, new StorageNull() );
                }
            }
        }

        for( set_it = keys_to_process.begin(); set_it != keys_to_process.end(); ++set_it) {
        	xprintf(Warn, "Unprocessed key '%s' in record '%s'.\n", (*set_it).c_str(), p.as_string().c_str() );
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
							make_storage_from_default( it->default_.value(), it->type_.get() ) );
				 } else { // defalut - optional or default at read time
					 ASSERT( ! it->default_.is_obligatory() ,
							 "Obligatory key: '%s' in auto-convertible record, wrong check during finish().", it->key_.c_str());
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



StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::AbstractRecord *abstr_rec)
{
	if ( p.is_record_type() ) {

		string descendant_name = p.get_descendant_name();
		if ( descendant_name == "" ) {
			if ( ! abstr_rec->get_selection_default().has_value_at_declaration() ) {
				THROW( ExcInputError() << EI_Specification("Missing key 'TYPE' in AbstractRecord.") << EI_ErrorAddress(p.as_string()) << EI_InputType(abstr_rec->desc()) );
			} else { // auto conversion
				return abstract_rec_automatic_conversion(p, abstr_rec);
			}
		} else {
			try {
				unsigned int descendant_index = (unsigned int)abstr_rec->get_type_selection().name_to_int( descendant_name );
				return make_storage(p, &( abstr_rec->get_descendant(descendant_index) ) );
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
			return abstract_rec_automatic_conversion(p, abstr_rec);
		}
	}

	return NULL;
}



StorageBase * ReaderToStorage::abstract_rec_automatic_conversion(PathBase &p, const Type::AbstractRecord *abstr_rec)
{
    // perform automatic conversion
    const Type::Record *default_child = abstr_rec->get_default_descendant();
    if (! default_child) THROW(ExcInputError()
    		<< EI_Specification("Auto conversion of AbstractRecord not allowed.\n")
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
            THROW( ExcInputError()
                    << EI_Specification("Do not fit into size limits of the Array.")
                    << EI_ErrorAddress(p.as_string()) << EI_InputType(array->desc()) );
        }
    } else {
        // try automatic conversion to array with one element
        if ( array->match_size( 1 ) ) {
            StorageArray *storage_array = new StorageArray(1);
            const Type::TypeBase &sub_type = array->get_sub_type();
            storage_array->new_item(0, make_storage(p, &sub_type) );

            return storage_array;
        } else {
            THROW( ExcInputError() << EI_Specification("Automatic conversion to array not allowed. The value should be '" + p.get_node_type(ValueTypes::array_type) + "', but we found: ")
                    << EI_ErrorAddress(p.as_string()) << EI_JSON_Type( p.get_node_type(p.get_node_type_index()) ) << EI_InputType(array->desc()) );
        }
    }

    return NULL;
}



StorageBase * ReaderToStorage::make_storage(PathBase &p, const Type::Selection *selection)
{
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
	string value;
	try {
		value = p.get_string_value();
	}
	catch (ExcInputError & e) {
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



StorageBase * ReaderToStorage::make_storage_from_default(const string &dflt_str, const Type::TypeBase *type) {
    try {
    	/*
    	// Possible simplification of this method (need default strings to be valid JSON)
    	ReaderToStorage  tmp_storage(dflt_str, *type);
    	return tmp_storage.storage_;
		*/

        // an auto-convertible AbstractRecord can be initialized form default value
    	const Type::AbstractRecord *a_record = dynamic_cast<const Type::AbstractRecord *>(type);
    	if (a_record != NULL ) {
    		ASSERT( a_record->get_selection_default().has_value_at_declaration(),
    				"Can not initialize (non-auto-convertible) AbstractRecord '%s' by default value\n", type->type_name().c_str() );
            return make_storage_from_default( dflt_str, a_record->get_default_descendant() );
        } else
        if (typeid(*type) == typeid(Type::Record) ) {
            // an auto-convertible Record can be initialized form default value
            const Type::Record *record = static_cast<const Type::Record *>(type);
            Type::Record::KeyIter auto_key_it = record->auto_conversion_key_iter();

            ASSERT( auto_key_it != record->end(), "Can not initialize (non-auto-convertible) Record '%s' by default value\n",
            		type->type_name().c_str());
			StorageArray *storage_array = new StorageArray(record->size());
			for( Type::Record::KeyIter it= record->begin(); it != record->end(); ++it) {
				if ( it == auto_key_it ) {
					// one key is initialized by the record default string
					storage_array->new_item(it->key_index, make_storage_from_default(dflt_str, it->type_.get()) );
				} else {

					ASSERT( ! it->default_.is_obligatory(),
							"Missing default value for key: '%s' in auto-convertible record, wrong check during finish().", it->key_.c_str());

					if (it->default_.has_value_at_declaration() ) {
					   storage_array->new_item(it->key_index,
							   make_storage_from_default( it->default_.value(), it->type_.get() ) );
					} else { // defalut - optional or default at read time
						// set null
						storage_array->new_item(it->key_index, new StorageNull() );
					}
				}
			}

			return storage_array;
        } else
        if (typeid(*type) == typeid(Type::Array) ) {
            const Type::Array *array = static_cast<const Type::Array *>(type);
            if ( array->match_size(1) ) {
               // try auto conversion to array
                StorageArray *storage_array = new StorageArray(1);
                const Type::TypeBase &sub_type = array->get_sub_type();
                storage_array->new_item(0, make_storage_from_default(dflt_str, &sub_type) );
                return storage_array;
            } else {
            	THROW(ExcInputMessage() << EI_Message("Can not initialize Array '" + type->type_name() + "' by default value, size 1 not allowed.\n"));
            }

        } else
        if (typeid(*type) == typeid(Type::Integer)) {
            return new StorageInt( static_cast<const Type::Integer *>(type) ->from_default(dflt_str) );
        } else
        if (typeid(*type) == typeid(Type::Double)) {
            return new StorageDouble( static_cast<const Type::Double *>(type) ->from_default(dflt_str) );
        } else
        if (typeid(*type) == typeid(Type::Bool)) {
            return new StorageBool( static_cast<const Type::Bool *>(type) ->from_default(dflt_str) );
        } else
        if (typeid(*type) == typeid(Type::Selection)) {
                return new StorageInt( static_cast<const Type::Selection *>(type) ->from_default(dflt_str) );
        } else {
            const Type::String * string_type = dynamic_cast<const Type::String *>(type);
            if (string_type != NULL ) return new StorageString( string_type->from_default(dflt_str) );

            // default error
            ASSERT(false, "Can not store default value for type: %s\n", typeid(type).name());
        }


    } catch (Input::Type::ExcWrongDefault & e) {
        // message to distinguish exceptions thrown during Default value check at declaration
        xprintf(Msg, "Wrong default value while reading an input stream:\n");
        e << EI_KeyName("UNKNOWN KEY");
        throw;
    }

    return NULL;
}



/********************************************88
 * Implementation
 */

template <class T>
T ReaderToStorage::get_root_interface() const
{
    ASSERT(storage_, "NULL pointer to storage !!! \n");

    Address addr(storage_, root_type_);
    // try to create an iterator just to check type
    Iterator<T>( *root_type_, addr, 0);

    auto tmp_root_type = static_cast<const typename T::InputType &>(*root_type_);
    return T( addr, tmp_root_type );
}


template ::Input::Record ReaderToStorage::get_root_interface<::Input::Record>() const;
template ::Input::Array ReaderToStorage::get_root_interface<::Input::Array>() const;
template ::Input::AbstractRecord ReaderToStorage::get_root_interface<::Input::AbstractRecord>() const;
//template ReaderToStorage::get_root_interface<::Input::>()->::Input::AbstractRecord const;


} // namespace Input
